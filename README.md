# grccallers

GRC mutation callers — variant-first callers for malaria drug-resistance genes.

Currently implements the **Kelch13 (K13) propeller-domain mutation caller**, which identifies non-synonymous amino-acid changes from VCF/BCF files or legacy AmpRecon genotype files.

## Features

- Accepts single- or multi-sample VCFs (bgzipped + tabix-indexed), a list of VCF paths, or legacy AmpRecon genotype files.
- Converts legacy genotype files to VCF on-the-fly.
- Annotates WHO-validated partial-resistance mutations.
- Outputs compact (one row per sample) and long-format (one row per mutation) TSV files.
- Optionally writes per-sample called sequences as multi-FASTA.

## Project Structure

```
src/grccallers/
├── callers/
│   └── k13.py            # K13 propeller-domain caller logic
├── cli/
│   ├── main.py            # CLI entry point (grccallers <tool>)
│   └── k13.py             # K13-specific CLI arguments & runner
├── core/
│   ├── SeqCaller.py       # Base caller framework, codon translation, VCF reading
│   └── tabular_converter.py  # Legacy genotype → VCF converter
rsc/
├── Pf3D7.fasta            # Reference genome (+ .fai index)
├── Pf3D7.gff.gz           # Gene annotations (+ .tbi index)
└── who_mutations.txt       # WHO-validated K13 mutations list
```

## Requirements

- Python ≥ 3.9
- [pysam](https://pysam.readthedocs.io/)
- [pandas](https://pandas.pydata.org/)
- A reference FASTA indexed with `samtools faidx`

## Installation

```bash
pip install .
```

This installs the `grccallers` command-line tool.

## Usage

### From a VCF file

```bash
grccallers k13 \
    --vcf       <vcf_file> \
    --fasta     rsc/Pf3D7.fasta \
    --gff       rsc/Pf3D7.gff.gz \
    --who-list  rsc/who_mutations.txt \
    --output-prefix k13_results
```

### From multiple VCFs (shell glob)

```bash
grccallers k13 \
    --vcf samples/*.vcf.gz \
    --fasta rsc/Pf3D7.fasta \
    --gff rsc/Pf3D7.gff.gz \
    --who-list rsc/who_mutations.txt
```

### From a VCF list file

```bash
grccallers k13 \
    --vcf-list vcf_paths.txt \
    --fasta rsc/Pf3D7.fasta \
    --gff rsc/Pf3D7.gff.gz \
    --who-list rsc/who_mutations.txt
```

### From a legacy AmpRecon genotype file

```bash
grccallers k13 \
    --genotype genotype_file.tsv \
    --fasta rsc/Pf3D7.fasta \
    --gff rsc/Pf3D7.gff.gz \
    --who-list rsc/who_mutations.txt
```

### QC threshold options

| Flag | Default | Description |
|------|---------|-------------|
| `--min-dp` | 10 | Minimum total read depth (DP) |
| `--min-ad` | 5 | Minimum allele depth (AD) |
| `--min-ratio` | 0.10 | Minimum allele frequency (AD/DP) |

### Other options

| Flag | Description |
|------|-------------|
| `--write-fasta` | Also output per-sample called sequences as FASTA |
| `--output-prefix` | Filename prefix for outputs (default: `k13_output`) |
| `--log-level` | Verbosity: DEBUG, INFO, WARNING, ERROR |

## Output

- `<prefix>_compact.tsv` — one row per sample, passing mutations comma-joined.
- `<prefix>_long.tsv` — one row per (sample, mutation, nucleotide change), includes all call statuses.
- `<prefix>.fasta` — (optional) per-sample called sequences.
