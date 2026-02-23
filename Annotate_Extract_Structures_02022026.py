#!/usr/bin/env python3
"""
Pipeline:
- Read a TSV with columns: Index (UniProt accession) and Gene (gene symbol)
- Annotate UniProt entries with keywords + DOMAIN features (and InterPro)
- Download AlphaFold PDBs to a folder (via EBI FTP "latest", with fallbacks)
- Run RCSB Pairwise Structure Alignment (jFATCAT-rigid) via alignment API
  aligning all structures to the first (reference) structure.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import requests

# ----------------------------
# User-editable defaults
# ----------------------------

DEFAULT_GENES_OF_INTEREST = [
    # Your forced-in set from earlier (edit as needed)
    "OXA1L", "TRAM1", "TIMM17A", "TIMM17B", "MTCH1", "GHITM", "DERL1", "BCAP29", "MBOAT7", "PTDSS1"
]

# UniProt REST base
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# AlphaFold FTP (EBI) base; we try multiple patterns
ALPHAFOLD_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold"

# RCSB Alignment API base
RCSB_ALIGN_BASE = "https://alignment.rcsb.org"  # endpoints: /submit, /results


# ----------------------------
# Helpers
# ----------------------------

def chunked(xs: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(xs), n):
        yield xs[i : i + n]


def safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def http_get_stream(url: str, timeout: int = 60) -> requests.Response:
    r = requests.get(url, stream=True, timeout=timeout)
    r.raise_for_status()
    return r


def head_ok(url: str, timeout: int = 30) -> bool:
    try:
        r = requests.head(url, timeout=timeout, allow_redirects=True)
        return 200 <= r.status_code < 300
    except Exception:
        return False


# ----------------------------
# Step 1: Read TSV
# ----------------------------

def load_gene_uniprot_from_tsv(tsv_path: Path) -> pd.DataFrame:
    """
    Expects columns like:
      - Index: UniProt accession (e.g., P07602)
      - Gene: gene symbol (e.g., PSAP)

    Returns DF with columns: gene, uniprot_acc
    """
    df = pd.read_csv(tsv_path, sep="\t")
    required = {"Index", "Gene"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"TSV missing required columns {sorted(missing)}. Found: {list(df.columns)}")

    out = df.rename(columns={"Index": "uniprot_acc", "Gene": "gene"})[["gene", "uniprot_acc"]].copy()
    # clean whitespace
    out["gene"] = out["gene"].astype(str).str.strip()
    out["uniprot_acc"] = out["uniprot_acc"].astype(str).str.strip()
    # drop empties
    out = out[(out["gene"] != "") & (out["uniprot_acc"] != "")]
    # drop duplicates (gene+acc)
    out = out.drop_duplicates(subset=["gene", "uniprot_acc"]).reset_index(drop=True)
    return out


def filter_genes(df: pd.DataFrame, genes: Optional[List[str]]) -> pd.DataFrame:
    if not genes:
        return df
    gene_set = set(genes)
    return df[df["gene"].isin(gene_set)].reset_index(drop=True)


# ----------------------------
# Step 2: UniProt annotation (keywords + domains)
# ----------------------------

@dataclass
class UniProtAnnot:
    accession: str
    gene_primary: str
    protein_name: str
    keywords: List[str]
    domains: List[str]
    interpro: List[str]


def fetch_uniprot_annotations(accessions: List[str], sleep_s: float = 0.2) -> Dict[str, UniProtAnnot]:
    """
    Uses UniProt REST API to fetch:
      - keywords
      - feature(DOMAIN)
      - InterPro cross-refs

    Returns dict: accession -> UniProtAnnot

    We do it in chunks to stay under URL/request-size limits.
    """
    results: Dict[str, UniProtAnnot] = {}

    # UniProt fields reference: return_fields / API queries docs.
    fields = ",".join([
        "accession",
        "protein_name",
        "gene_primary",
        "keyword",
        "feature(DOMAIN)",
        "xref_interpro",
    ])

    # UniProt search supports query like: (accession:P07602 OR accession:Q9H3G5 ...)
    for batch in chunked(accessions, 50):
        query = " OR ".join([f"accession:{a}" for a in batch])
        params = {
            "query": f"({query})",
            "format": "tsv",
            "fields": fields,
            "size": 500,  # enough for 50 accessions
        }

        r = requests.get(UNIPROT_SEARCH, params=params, timeout=60)
        r.raise_for_status()

        # Parse TSV response
        # Columns can vary slightly; we rely on header names.
        lines = r.text.splitlines()
        if not lines:
            continue
        header = lines[0].split("\t")
        rows = [ln.split("\t") for ln in lines[1:] if ln.strip()]

        df = pd.DataFrame(rows, columns=header)

        # Normalize expected columns
        col_map = {
            "Entry": "accession",
            "Protein names": "protein_name",
            "Gene Names (primary)": "gene_primary",
            "Keywords": "keyword",
            "Domain [FT]": "理解",  # placeholder if UniProt changes; we detect below
        }

        # Identify columns robustly
        # UniProt typically returns:
        #  - 'Entry' or 'Entry\t...'
        #  - 'Gene Names (primary)'
        #  - 'Protein names'
        #  - 'Keywords'
        #  - 'Domain [FT]'
        #  - 'Cross-reference (InterPro)'
        # We'll just look them up by substring.
        def find_col(substr: str) -> Optional[str]:
            for c in df.columns:
                if substr.lower() in c.lower():
                    return c
            return None

        c_entry = find_col("Entry")
        c_gene = find_col("Gene Names (primary)") or find_col("Gene Names")
        c_pname = find_col("Protein names") or find_col("Protein Name")
        c_kw = find_col("Keywords")
        c_dom = find_col("Domain [FT]") or find_col("Domain")
        c_ipr = find_col("InterPro")

        for _, row in df.iterrows():
            acc = (row.get(c_entry) or "").strip()
            if not acc:
                continue

            gene_primary = (row.get(c_gene) or "").strip()
            protein_name = (row.get(c_pname) or "").strip()

            keywords_raw = (row.get(c_kw) or "").strip()
            domains_raw = (row.get(c_dom) or "").strip()
            interpro_raw = (row.get(c_ipr) or "").strip()

            keywords = [k.strip() for k in keywords_raw.split(";") if k.strip()] if keywords_raw else []
            domains = [d.strip() for d in domains_raw.split(";") if d.strip()] if domains_raw else []
            interpro = [i.strip() for i in interpro_raw.split(";") if i.strip()] if interpro_raw else []

            results[acc] = UniProtAnnot(
                accession=acc,
                gene_primary=gene_primary,
                protein_name=protein_name,
                keywords=keywords,
                domains=domains,
                interpro=interpro,
            )

        time.sleep(sleep_s)

    return results


# ----------------------------
# Step 3: AlphaFold PDB downloads (FTP "latest" + fallbacks)
# ----------------------------

def alphafold_candidate_urls(uniprot_acc: str) -> List[Tuple[str, bool]]:
    """
    Returns list of (url, is_gz).
    We try 'latest' first; then older versioned patterns if needed.
    """
    acc = uniprot_acc.strip()

    # Common published filename patterns
    # - On the website, files often appear under https://alphafold.ebi.ac.uk/files/...
    # - On FTP, the 'latest' directory provides current models.
    # We attempt both FTP and web 'files' as fallbacks.
    return [
        # FTP latest PDB (gz)
        (f"{ALPHAFOLD_FTP_BASE}/latest/AF-{acc}-F1-model_v4.pdb.gz", True),
        # FTP latest mmCIF (gz) -> (we won't use for PDB, but you could extend)
        # (f"{ALPHAFOLD_FTP_BASE}/latest/AF-{acc}-F1-model_v4.cif.gz", True),

        # Website files (sometimes easier; still official)
        (f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb", False),
        (f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb.gz", True),

        # Older versions occasionally referenced in examples (v2)
        (f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb", False),
        (f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb.gz", True),
    ]


def download_alphafold_pdb(uniprot_acc: str, out_dir: Path, overwrite: bool = False) -> Optional[Path]:
    """
    Downloads AlphaFold PDB for UniProt accession.
    Saves as: out_dir/AF-{acc}-F1.pdb
    Returns local path or None if not found.
    """
    safe_mkdir(out_dir)
    acc = uniprot_acc.strip()
    out_path = out_dir / f"AF-{acc}-F1.pdb"
    if out_path.exists() and not overwrite:
        return out_path

    candidates = alphafold_candidate_urls(acc)

    for url, is_gz in candidates:
        if not head_ok(url):
            continue

        try:
            r = http_get_stream(url, timeout=120)
            if is_gz or url.endswith(".gz"):
                # stream gunzip to file
                with gzip.GzipFile(fileobj=r.raw) as gz, open(out_path, "wb") as fh:
                    fh.write(gz.read())
            else:
                with open(out_path, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            fh.write(chunk)
            return out_path
        except Exception as e:
            print(f"[WARN] download failed {acc} from {url}: {e}")

    print(f"[WARN] Could not find AlphaFold PDB for {acc} in tried locations.")
    return None


# ----------------------------
# Step 4: RCSB Alignment API (pairwise, jFATCAT)
# ----------------------------

def rcsb_csm_id_for_uniprot(acc: str) -> str:
    """
    RCSB CSM ID format for AlphaFold models:
      AlphaFold ID: AF-{ACC}-F1
      RCSB CSM ID:  AF_AF{ACC}F1   (no hyphens)
    Examples in RCSB docs: AF-Q9Y6F1-F1 -> AF_AFQ9Y6F1F1
    """
    acc = acc.strip()
    return f"AF_AF{acc}F1"


def submit_rcsb_alignment_job(
    uniprot_accessions: List[str],
    method: str = "jFATCAT-rigid",
    chain_id: str = "A",
    poll_s: float = 2.0,
    timeout_s: float = 600.0,
) -> Dict:
    """
    Submits a pairwise alignment job aligning all structures to the first (reference).
    Uses RCSB CSM IDs for AlphaFold (AF_AF{ACC}F1) and chain 'A' (typical for AF models).

    Saves and returns the final JSON result from /results.

    IMPORTANT:
    The Alignment API runs asynchronously: /submit returns a ticket; /results returns status until COMPLETE.
    """
    if len(uniprot_accessions) < 2:
        raise ValueError("Need at least 2 accessions for pairwise alignment.")

    # Build structure list
    structures = []
    for acc in uniprot_accessions:
        structures.append({
            # Use RCSB CSM ID for AlphaFold models:
            "structure_id": rcsb_csm_id_for_uniprot(acc),
            "chain_id": chain_id,
        })

    # Payload schema:
    # The public guide states /submit takes JSON query with 'options' and 'context'.
    # Exact field names may evolve; if your server returns 400, inspect error and adjust keys.
    query_obj = {
        "query": {
            "context": {
                "mode": "pairwise",
                "structures": structures
            },
            "options": {
                "method": method
            }
        }
    }

    submit_url = f"{RCSB_ALIGN_BASE}/submit"
    results_url = f"{RCSB_ALIGN_BASE}/results"

    # POST as form/multipart-ish. This generally matches their guidance about multipart/form-data.
    resp = requests.post(submit_url, data={"query": json.dumps(query_obj)}, timeout=60)
    resp.raise_for_status()

    ticket = resp.text.strip().strip('"')
    if not re.fullmatch(r"[0-9a-fA-F-]{36}", ticket):
        # Some deployments return JSON; handle it
        try:
            ticket = resp.json().get("ticket") or resp.json().get("id") or ticket
        except Exception:
            pass

    print(f"[INFO] Alignment ticket: {ticket}")

    # Poll /results until COMPLETE/ERROR
    t0 = time.time()
    while True:
        if time.time() - t0 > timeout_s:
            raise TimeoutError(f"Alignment job timed out after {timeout_s}s (ticket={ticket})")

        r = requests.get(results_url, params={"ticket": ticket}, timeout=60)
        r.raise_for_status()
        data = r.json()

        status = (data.get("status") or "").upper()
        if status == "COMPLETE":
            return data
        if status == "ERROR":
            raise RuntimeError(f"Alignment API returned ERROR: {json.dumps(data, indent=2)[:2000]}")

        time.sleep(poll_s)


# ----------------------------
# Main
# ----------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Annotate UniProt + download AlphaFold PDBs + RCSB jFATCAT align.")
    ap.add_argument("--tsv", required=True, type=Path, help="Input TSV (expects columns Index=UniProt acc, Gene=gene symbol).")
    ap.add_argument("--genes", nargs="*", default=None, help="Optional gene symbols to filter to. If omitted, uses all rows.")
    ap.add_argument("--default-genes-of-interest", action="store_true",
                    help="If set and --genes not provided, filter to built-in DEFAULT_GENES_OF_INTEREST.")
    ap.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    ap.add_argument("--overwrite-pdb", action="store_true", help="Overwrite existing downloaded PDBs.")
    ap.add_argument("--uniprot-cache", default=None, type=Path, help="Optional path to write/read UniProt annotation cache JSON.")
    ap.add_argument("--align-method", default="jFATCAT-rigid", help="Alignment method (e.g. jFATCAT-rigid, jFATCAT-flexible).")
    ap.add_argument("--chain-id", default="A", help="Chain ID for alignment (AlphaFold models usually 'A').")
    args = ap.parse_args()

    outdir: Path = args.outdir
    safe_mkdir(outdir)

    # Load TSV
    df = load_gene_uniprot_from_tsv(args.tsv)

    # Choose genes
    genes = args.genes
    if genes is None and args.default_genes_of_interest:
        genes = DEFAULT_GENES_OF_INTEREST
    df = filter_genes(df, genes)

    if df.empty:
        raise SystemExit("No proteins left after filtering. Check gene symbols / TSV content.")

    # Extract accessions
    accessions = df["uniprot_acc"].tolist()
    print(f"[INFO] Selected {len(accessions)} proteins.")

    # UniProt annotation (with optional cache)
    annot: Dict[str, UniProtAnnot] = {}
    cache_path: Optional[Path] = args.uniprot_cache

    if cache_path and cache_path.exists():
        print(f"[INFO] Loading UniProt cache: {cache_path}")
        raw = json.loads(cache_path.read_text())
        for acc, a in raw.items():
            annot[acc] = UniProtAnnot(
                accession=acc,
                gene_primary=a.get("gene_primary", ""),
                protein_name=a.get("protein_name", ""),
                keywords=a.get("keywords", []) or [],
                domains=a.get("domains", []) or [],
                interpro=a.get("interpro", []) or [],
            )

    missing = [a for a in accessions if a not in annot]
    if missing:
        print(f"[INFO] Fetching UniProt annotations for {len(missing)} accessions...")
        fetched = fetch_uniprot_annotations(missing)
        annot.update(fetched)

    if cache_path:
        cache_payload = {
            acc: {
                "gene_primary": a.gene_primary,
                "protein_name": a.protein_name,
                "keywords": a.keywords,
                "domains": a.domains,
                "interpro": a.interpro,
            }
            for acc, a in annot.items()
        }
        cache_path.write_text(json.dumps(cache_payload, indent=2))
        print(f"[INFO] Wrote UniProt cache: {cache_path}")

    # Write annotation table
    rows = []
    for _, r in df.iterrows():
        acc = r["uniprot_acc"]
        a = annot.get(acc)
        rows.append({
            "gene": r["gene"],
            "uniprot_acc": acc,
            "uniprot_gene_primary": a.gene_primary if a else "",
            "protein_name": a.protein_name if a else "",
            "keywords": "; ".join(a.keywords) if a else "",
            "domains": "; ".join(a.domains) if a else "",
            "interpro": "; ".join(a.interpro) if a else "",
        })
    anno_df = pd.DataFrame(rows)
    anno_tsv = outdir / "uniprot_annotation.tsv"
    anno_df.to_csv(anno_tsv, sep="\t", index=False)
    print(f"[INFO] Wrote annotations: {anno_tsv}")

    # Download AlphaFold PDBs
    pdb_dir = outdir / "alphafold_pdbs"
    safe_mkdir(pdb_dir)

    downloaded = []
    for acc in accessions:
        p = download_alphafold_pdb(acc, pdb_dir, overwrite=args.overwrite_pdb)
        if p:
            downloaded.append(p)

    print(f"[INFO] Downloaded {len(downloaded)}/{len(accessions)} AlphaFold PDBs to: {pdb_dir}")

    # Run alignment (requires >=2 structures)
    if len(accessions) >= 2:
        print(f"[INFO] Submitting RCSB alignment job ({args.align_method}) for {len(accessions)} proteins...")
        result = submit_rcsb_alignment_job(
            accessions,
            method=args.align_method,
            chain_id=args.chain_id,
        )

        out_json = outdir / "rcsb_alignment_results.json"
        out_json.write_text(json.dumps(result, indent=2))
        print(f"[INFO] Wrote alignment results: {out_json}")
    else:
        print("[WARN] Skipping alignment (need >= 2 proteins).")


if __name__ == "__main__":
    main()


python align_membrane_set.py \
  --tsv /mnt/data/ratio_protein_MD.tsv \
  --default-genes-of-interest \
  --outdir ./membrane_align_out \
  --uniprot-cache ./membrane_align_out/uniprot_cache.json \
  --align-method jFATCAT-rigid