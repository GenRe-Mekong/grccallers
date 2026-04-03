#!/usr/bin/env python3
"""CLI handler for the K13 propeller-domain mutation caller."""
import argparse
import logging
import sys
from pathlib import Path
from typing import List

from ..callers.k13 import K13Caller, K13MutationWriter
from ..core.SeqCaller import FastaWriter, QCThresholds, SeqReader
from ..core.tabular_converter import genotype_to_vcf


class _CustomFormatter(argparse.RawDescriptionHelpFormatter,
                       argparse.ArgumentDefaultsHelpFormatter):
    pass


DESCRIPTION = """\
Kelch13 Mutation Caller (Variant-First)
=======================================
Identifies non-synonymous amino-acid changes in the K13 propeller domain
from a multi-sample VCF.

Call-status values in the long-format output
--------------------------------------------
  PASS           Position passes all QC thresholds
  STOP_CODON     Variant produces a stop codon — excluded from results but
                 preserved in the long-format TSV for audit purposes.

Only mutations with status PASS appear in the compact output.
All statuses are written to the long-format output for auditability.

Output files
------------
  <prefix>_compact.tsv  one row per sample, passing mutations comma-joined
  <prefix>_long.tsv     one row per (sample, mutation, nucleotide change),
                        includes all call_status values
"""

EPILOG = (
    "Examples:\n"
    "  # Single multi-sample VCF\n"
    "  grccallers k13 \\\n"
    "      --vcf       <vcf_file> \\\n"
    "      --fasta     <fasta_file> \\\n"
    "      --gff       <gff_file> \\\n"
    "      --who-list  <who_list_file> \\\n"
    "      --output-prefix <output_prefix>\n"
    "\n"
    "  # Multiple single-sample VCFs (shell glob)\n"
    "  grccallers k13 \\\n"
    "      --vcf samples/*.vcf.gz \\\n"
    "      --fasta <fasta_file> --gff <gff_file> --who-list <who_list_file>\n"
    "\n"
    "  # List file (one VCF path per line)\n"
    "  grccallers k13 \\\n"
    "      --vcf-list vcf_paths.txt \\\n"
    "      --fasta <fasta_file> --gff <gff_file> --who-list <who_list_file>\n"
    "\n"
    "  # Legacy genotype file (tab-delimited, from AmpRecon)\n"
    "  grccallers k13 \\\n"
    "      --genotype genotype_file.tsv \\\n"
    "      --fasta <fasta_file> --gff <gff_file> --who-list <who_list_file>\n"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Register all K13-specific arguments on *parser*."""
    req = parser.add_argument_group(
        "required inputs",
        "Provide exactly one of --vcf, --vcf-list, or --genotype."
    )
    vcf_grp = req.add_mutually_exclusive_group(required=True)
    vcf_grp.add_argument(
        "--vcf", nargs="+", metavar="PATH",
        help=(
            "One or more bgzipped, tabix-indexed VCF/BCF files.  Each file may "
            "be a single-sample or multi-sample VCF; all samples across all "
            "files are processed and combined into a single output.  "
            "Shell globs are expanded by the shell before reaching this flag.  "
            "Mutually exclusive with --vcf-list and --genotype."
        ),
    )
    vcf_grp.add_argument(
        "--vcf-list", metavar="PATH",
        help=(
            "Plain-text file with one VCF path per line.  Blank lines and "
            "lines starting with '#' are ignored.  "
            "Mutually exclusive with --vcf and --genotype."
        ),
    )
    vcf_grp.add_argument(
        "--genotype", metavar="PATH",
        help=(
            "Legacy tab-delimited genotype file (from AmpRecon's "
            "write_genotype_files.py).  The file is converted to a "
            "multi-sample VCF on-the-fly before calling.  The converted "
            "VCF is written next to the output files as "
            "<output-prefix>_converted.vcf.gz for inspection.  "
            "Mutually exclusive with --vcf and --vcf-list."
        ),
    )
    req.add_argument(
        "--fasta", required=True, metavar="PATH",
        help=(
            "Reference genome FASTA indexed with 'samtools faidx'.  Only the "
            "K13 chromosome (Pf3D7_13_v3) is accessed."
        ),
    )
    req.add_argument(
        "--gff", required=True, metavar="PATH",
        help=(
            "bgzipped, tabix-indexed GFF3 annotation file.  Used to locate "
            "K13 CDS exon coordinates (gene ID PF3D7_1343700.1)."
        ),
    )
    req.add_argument(
        "--who-list", required=True, metavar="PATH",
        help=(
            "Plain-text file of WHO-validated K13 partial-resistance mutations, "
            "one per line in standard amino-acid notation (e.g. C580Y).  Lines "
            "starting with '#' and blank lines are ignored."
        ),
    )

    out = parser.add_argument_group("output")
    out.add_argument(
        "--output-prefix", default="k13_output", metavar="PREFIX",
        help=(
            "Filename prefix for the output TSV files.  For example, "
            "'k13_results' produces k13_results_compact.tsv "
            "and k13_results_long.tsv."
        ),
    )
    out.add_argument(
        "--write-fasta", action="store_true", default=False,
        help=(
            "Also write a multi-FASTA file of per-sample called sequences.  "
            "Output is written to <output-prefix>.fasta."
        ),
    )

    std = parser.add_argument_group(
        "standard QC thresholds",
        "Applied to all WHO-confirmed mutations and to the initial screen of novel mutations.",
    )
    std.add_argument(
        "--min-dp", type=int, default=10, metavar="N",
        help="Minimum total read depth (DP FORMAT field) at the position.",
    )
    std.add_argument(
        "--min-ad", type=int, default=5, metavar="N",
        help="Minimum allele depth (AD) required for an allele to be considered.",
    )
    std.add_argument(
        "--min-ratio", type=float, default=0.10, metavar="F",
        help="Minimum allele frequency (AD / DP) for an allele to pass.",
    )

    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Verbosity level.  Use DEBUG to trace individual codon/sample decisions.",
    )


def _resolve_vcf_paths(args: argparse.Namespace) -> List[Path]:
    """Return an ordered, validated list of VCF paths from --vcf, --vcf-list, or --genotype."""
    args._converted_vcf = None

    if args.genotype:
        geno_path = Path(args.genotype)
        if not geno_path.exists():
            logging.error(f"--genotype: file not found: {geno_path}")
            sys.exit(1)
        converted = f"{args.output_prefix}_converted.vcf.gz"
        logging.info(f"Converting genotype file {geno_path} → {converted}")
        genotype_to_vcf(
            genotype_path=str(geno_path),
            fasta_path=args.fasta,
            output_path=converted,
        )
        args._converted_vcf = converted
        return [Path(converted)]

    if args.vcf_list:
        list_path = Path(args.vcf_list)
        if not list_path.exists():
            logging.error(f"--vcf-list: file not found: {list_path}")
            sys.exit(1)
        paths = []
        for lineno, raw in enumerate(list_path.read_text().splitlines(), 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            p = Path(line)
            if not p.exists():
                logging.error(f"--vcf-list line {lineno}: file not found: {p}")
                sys.exit(1)
            paths.append(p)
        if not paths:
            logging.error(f"--vcf-list: no VCF paths found in {list_path}")
            sys.exit(1)
        return paths
    else:
        paths = []
        for raw in args.vcf:
            p = Path(raw)
            if not p.exists():
                logging.error(f"--vcf: file not found: {p}")
                sys.exit(1)
            paths.append(p)
        return paths


def run(args: argparse.Namespace) -> None:
    """Execute the K13 caller from a pre-parsed Namespace."""
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    for flag, path in [
        ("--fasta",    args.fasta),
        ("--gff",      args.gff),
        ("--who-list", args.who_list),
    ]:
        if not Path(path).exists():
            logging.error(f"{flag}: file not found: {path}")
            sys.exit(1)

    vcf_paths = _resolve_vcf_paths(args)
    logging.info(f"{len(vcf_paths)} VCF file(s) to process.")

    qc = QCThresholds(
        min_dp=args.min_dp,
        min_ad=args.min_ad,
        min_ratio=args.min_ratio,
    )

    seq_reader = SeqReader(fasta_path=args.fasta, gff_path=args.gff)

    merged_results = {}
    merged_seqs = {}

    for vcf_path in vcf_paths:
        logging.info(f"Calling: {vcf_path}")
        caller = K13Caller(
            vcf_path=str(vcf_path),
            who_path=args.who_list,
            qc=qc,
            seq_reader=seq_reader,
        )
        merged_results.update(caller.call())
        if args.write_fasta:
            merged_seqs.update(caller.call_seq())

    if hasattr(args, '_converted_vcf') and args._converted_vcf is not None:
        logging.info(f"Converted VCF kept at: {args._converted_vcf}")

    logging.info(f"Calling complete. {len(merged_results)} sample(s) processed.")

    writer = K13MutationWriter(merged_results)
    writer.write_compact(f"{args.output_prefix}_compact.tsv")
    writer.write_long(f"{args.output_prefix}_long.tsv")

    if args.write_fasta:
        fasta_writer = FastaWriter(merged_seqs)
        fasta_path = f"{args.output_prefix}.fasta"
        fasta_writer.write_fasta(fasta_path)
        logging.info(f"FASTA written to {fasta_path}")

    logging.info("Done.")


def _build_parser(prog: str = "grccallers k13") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=_CustomFormatter,
    )
    add_args(parser)
    return parser


def parse_args(argv=None) -> argparse.Namespace:
    return _build_parser().parse_args(argv)


def main(argv=None) -> None:
    run(parse_args(argv))


if __name__ == "__main__":
    main()
