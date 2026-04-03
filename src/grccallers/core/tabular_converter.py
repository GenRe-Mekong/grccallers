#!/usr/bin/env python3
"""
tabular_converter
=================
Convert a legacy tab-delimited genotype file (produced by AmpRecon's
write_genotype_files.py) into a multi-sample, bgzipped + tabix-indexed VCF
that can be consumed directly by CallerBase / K13Caller.

The genotype TSV columns are:
    ID  Amplicon  Pos  Chr  Loc  Gen  Depth  Filt

``Gen`` and ``Depth`` can be comma-separated for mixed/het calls, e.g.
``Gen=T,A  Depth=111,16``.

Usage from Python
-----------------
>>> from grccallers.core.tabular_converter import genotype_to_vcf
>>> vcf_path = genotype_to_vcf("genotype_file.tsv", "ref.fasta", "output.vcf.gz")

The function returns the path to the bgzipped, tabix-indexed VCF.
"""
from __future__ import annotations

import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="[TabConverter][%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def genotype_to_vcf(
    genotype_path: str,
    fasta_path: str,
    output_path: str,
    chrom_filter: Optional[List[str]] = None,
) -> str:
    """
    Convert a legacy genotype TSV to a multi-sample bgzipped VCF.

    Args:
        genotype_path: Path to the input tab-delimited genotype file.
        fasta_path:    Path to an indexed reference FASTA (samtools faidx).
        output_path:   Desired path for the output VCF (.vcf.gz).
                       A tabix index (.vcf.gz.tbi) is created alongside.
        chrom_filter:  If given, only rows whose ``Chr`` is in this list are
                       converted.  Contigs not present in the reference FASTA
                       are always skipped regardless of this filter.

    Returns:
        The *output_path* string (for chaining convenience).
    """
    fasta = pysam.FastaFile(fasta_path)
    fasta_contigs: Set[str] = set(fasta.references)

    # ---- 1. Load and validate the genotype file --------------------------
    df = _load_genotype(genotype_path)

    # ---- 2. Apply filters ------------------------------------------------
    if chrom_filter is not None:
        df = df[df["Chr"].isin(chrom_filter)]

    # Keep only rows whose chrom exists in the reference FASTA
    df = df[df["Chr"].isin(fasta_contigs)]
    if df.empty:
        logging.error(
            "No rows remain after filtering to reference contigs. "
            "Check that the genotype file chromosomes match the FASTA."
        )
        sys.exit(1)

    # ---- 3. Discover samples, positions, and build site index ------------
    samples: List[str] = list(dict.fromkeys(df["ID"]))
    logging.info(f"Converting {len(samples)} sample(s) across {df['Chr'].nunique()} contig(s)")

    # Group by (Chr, Loc) → list of rows
    site_rows: Dict[Tuple[str, int], Dict[str, _GenotypeRow]] = defaultdict(dict)
    for row in df.itertuples(index=False):
        chrom = str(row.Chr)
        loc   = int(row.Loc)
        site_rows[(chrom, loc)][str(row.ID)] = _GenotypeRow(
            gen=str(row.Gen),
            depth=str(row.Depth),
            filt=str(row.Filt).upper(),
        )

    # Sort sites by (chrom order in FASTA, position)
    contig_order = {c: i for i, c in enumerate(fasta.references)}
    sorted_sites = sorted(site_rows.keys(), key=lambda k: (contig_order.get(k[0], 999), k[1]))

    # ---- 4. Write VCF ----------------------------------------------------
    _write_vcf(
        output_path=output_path,
        fasta=fasta,
        samples=samples,
        sorted_sites=sorted_sites,
        site_rows=site_rows,
    )
    fasta.close()

    # ---- 5. bgzip + tabix ------------------------------------------------
    final_path = _compress_and_index(output_path)
    logging.info(f"Converted VCF written to {final_path}")
    return final_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

class _GenotypeRow:
    __slots__ = ("gen", "depth", "filt")

    def __init__(self, gen: str, depth: str, filt: str):
        self.gen   = gen
        self.depth = depth
        self.filt  = filt


def _load_genotype(path: str) -> pd.DataFrame:
    """Load and validate the genotype TSV."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"ID", "Amplicon", "Pos", "Chr", "Loc", "Gen", "Depth", "Filt"}
    missing = required - set(df.columns)
    if missing:
        logging.error(f"Genotype file missing required columns: {missing}")
        sys.exit(1)
    logging.info(f"Loaded {len(df)} rows from {path}")
    return df


def _write_vcf(
    output_path: str,
    fasta: pysam.FastaFile,
    samples: List[str],
    sorted_sites: List[Tuple[str, int]],
    site_rows: Dict[Tuple[str, int], Dict[str, _GenotypeRow]],
) -> None:
    """Write an uncompressed VCF from the parsed genotype data."""
    # Build the header
    hdr = pysam.VariantHeader()

    # Add contigs that appear in our data
    contigs_seen: Set[str] = {chrom for chrom, _ in sorted_sites}
    for contig in fasta.references:
        if contig in contigs_seen:
            hdr.add_line(
                f'##contig=<ID={contig},length={fasta.get_reference_length(contig)}>'
            )

    # FORMAT fields matching what CallerBase._allele_call() expects
    hdr.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    hdr.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">')
    hdr.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt alleles">')
    hdr.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
    hdr.add_line('##FILTER=<ID=LowQual,Description="Low quality call">')

    for sample in samples:
        hdr.add_sample(sample)

    # Strip .gz if present — we write uncompressed first, then bgzip
    raw_path = output_path.removesuffix(".gz")
    vcf_out = pysam.VariantFile(raw_path, "w", header=hdr)

    n_written = 0
    for chrom, loc in sorted_sites:
        ref_base = fasta.fetch(chrom, loc - 1, loc).upper()
        per_sample = site_rows[(chrom, loc)]

        # Collect all non-ref alleles across every sample at this site
        alt_alleles = _collect_alt_alleles(ref_base, per_sample)

        # Skip reference-only sites (no ALT alleles) — pysam requires ≥2 alleles,
        # and CallerBase only scans variant positions anyway.  We still need these
        # sites for DP coverage though, so emit them with a <NON_REF> placeholder.
        if not alt_alleles:
            alt_alleles = ["<NON_REF>"]

        # Build the VCF record
        rec = vcf_out.new_record()
        rec.contig = chrom
        rec.pos = loc  # pysam uses 1-based POS
        rec.stop = loc  # 0-based exclusive end; ensures correct INFO/END for symbolic alleles
        rec.alleles = (ref_base, *alt_alleles) if alt_alleles else (ref_base,)
        rec.filter.add("PASS")

        allele_idx = {a: i for i, a in enumerate(rec.alleles)}

        for sample in samples:
            grow = per_sample.get(sample)
            if grow is None:
                # Sample has no data at this position → missing
                rec.samples[sample]["GT"] = None
                continue

            alleles_str = [a.strip() for a in grow.gen.split(",")]
            depths_str  = [d.strip() for d in grow.depth.split(",")]

            if len(alleles_str) != len(depths_str):
                logging.warning(
                    f"[{sample}] Pos {loc}: Gen/Depth length mismatch "
                    f"({grow.gen!r} vs {grow.depth!r}), treating as missing"
                )
                rec.samples[sample]["GT"] = None
                continue

            try:
                depths = [int(d) for d in depths_str]
            except ValueError:
                logging.debug(
                    f"[{sample}] Pos {loc}: non-integer depth {grow.depth!r}, treating as missing"
                )
                rec.samples[sample]["GT"] = None
                continue

            # Build AD array in the order of rec.alleles (REF, ALT1, ALT2, ...)
            ad = [0] * len(rec.alleles)
            gt_indices = []
            for allele, dp in zip(alleles_str, depths):
                idx = allele_idx.get(allele.upper())
                if idx is not None:
                    ad[idx] += dp
                    if idx not in gt_indices:
                        gt_indices.append(idx)
                else:
                    # Allele not in this site's ALT list — shouldn't happen after _collect_alt_alleles
                    logging.debug(f"[{sample}] Pos {loc}: allele {allele!r} not in record alleles, skipping")

            total_dp = sum(ad)
            rec.samples[sample]["DP"] = total_dp
            rec.samples[sample]["AD"] = tuple(ad)

            if not gt_indices:
                gt_indices = [0]  # default to REF
            rec.samples[sample]["GT"] = tuple(gt_indices)
            rec.samples[sample].phased = False

        vcf_out.write(rec)
        n_written += 1

    vcf_out.close()
    logging.info(f"Wrote {n_written} variant sites to {raw_path}")


def _collect_alt_alleles(ref_base: str, per_sample: Dict[str, _GenotypeRow]) -> List[str]:
    """Gather all unique non-REF alleles across samples at a single site."""
    alts: List[str] = []
    seen: Set[str] = set()
    for grow in per_sample.values():
        for a in grow.gen.split(","):
            a = a.strip().upper()
            if a != ref_base and a not in seen and a != "-":
                alts.append(a)
                seen.add(a)
    return sorted(alts)


def _compress_and_index(output_path: str) -> str:
    """bgzip and tabix-index the VCF. Returns the .vcf.gz path."""
    raw_path = output_path.removesuffix(".gz")
    gz_path = raw_path + ".gz" if not output_path.endswith(".gz") else output_path

    pysam.tabix_compress(raw_path, gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)

    # Clean up the uncompressed file
    Path(raw_path).unlink(missing_ok=True)
    return gz_path
