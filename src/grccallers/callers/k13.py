#!/usr/bin/env python3
"""
K13Caller.py
==================
Tool to do the calling of non-synonymous mutations and interpretation of the allele calls at the codon level.
This it build on top of the SeqReader and CdsReader to handle the reference sequence and coordinate mapping, 
and provides a CallerBase class that defines the general framework for variant-first calling of non-synonymous mutations.
"""
from __future__ import annotations

import logging
import sys
from collections import defaultdict
from dataclasses import dataclass, field, fields
from typing import Dict, List, Optional, Set, Tuple
from ..core.SeqCaller import CallerBase, NonSynMutation, SeqReader, CdsReader, QCThresholds, MutationWriter, FastaWriter
import pandas as pd

# ============ Setting logging ==============================================
logging.basicConfig(
    level=logging.INFO,
    format="[K13Caller][%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# ============ K13 Gene Constants ===========================================

K13_CDS_ID = "PF3D7_1343700.1"
K13_CHROM  = "Pf3D7_13_v3"
K13_START  = 1724768
K13_STOP   = 1725559
K13_STRAND = "-"

# ============ Data Class ===========================================

@dataclass
class WHOCov:
    perc_cov:    float
    missing_pos: Dict[int, List[int]] = field(default_factory=dict) # aa_pos → list of genomic positions with missing calls

@dataclass
class KelchMutation(NonSynMutation):
    is_who_validated: bool
    call_status: str

    @classmethod
    def from_nonsyn(cls, mut: NonSynMutation, is_who_validated: bool, call_status: str) -> "KelchMutation":
        return cls(
            **{f.name: getattr(mut, f.name) for f in fields(NonSynMutation)},
            is_who_validated=is_who_validated,
            call_status=call_status
        )

@dataclass
class KelchResult:
    sample_id: str
    mutations: List[KelchMutation]
    who_cov:   WHOCov


# =========== WHO-validated mutation list handling ===========================================

class WHOList:
    """
    Class to hold the WHO-validated K13 mutations loaded from a file, and provide helper functions to query and map them.

    This init by loading a plain text file where each line is a mutation in the format of "C580Y".
    Lines starting with "#" are treated as comments and ignored.
    """
    def __init__(self, mutations: frozenset):
        self._mutations = mutations

    @classmethod
    def from_file(cls, path: str) -> "WHOList":
        """
        Load WHO-validated mutations from a plain text file where each line is a mutation in the format "C580Y".
        Lines starting with "#" are treated as comments and ignored.
        """
        mutations = set()
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                mutations.add(line)
        logging.info(f"Loaded {len(mutations)} WHO-validated mutations from {path}")
        return cls(frozenset(mutations))

    def is_confirmed(self, mutation: str) -> bool:
        """
        Check if a given mutation (e.g. "C580Y") is in the WHO-validated mutation list.
        """
        return mutation in self._mutations
    
    def get_aa_mutation_map(self) -> Set[int]:
        """
        return a map of amino-acid position → set of confirmed alt AAs (e.g. {580: "Y"}) based on the WHO-validated mutation list.
        """
        aa_map: Set[int] = set()
        for mut in self._mutations:
            # Extract the amino acid position from the mutation string (e.g., "C580Y" → 580)
            try:
                pos_str = ''.join(filter(str.isdigit, mut))
                aa_pos = int(pos_str)
                aa_map.add(aa_pos)
            # TODO: Validate the mutation format more robustly early on, and handle parsing errors gracefully.
            except ValueError:
                logging.error(f"Invalid mutation format: {mut}")
                sys.exit(255)
        return aa_map
    
    def __len__(self) -> int:
        return len(self._mutations)

# ========== K13 Caller Implementation ===========================================

class K13Caller(CallerBase):

    def __init__(
        self,
        vcf_path: str,
        who_path: str,
        qc: QCThresholds,
        seq_reader: Optional[SeqReader] = None,
        cds_reader: Optional[CdsReader] = None,
    ):
        if cds_reader is None and seq_reader is not None:
            cds_reader = seq_reader.build_cds_reader(chrom=K13_CHROM, id=K13_CDS_ID)

        super().__init__(vcf_path=vcf_path, qc=qc, seq_reader=seq_reader, cds_reader=cds_reader)

        self.who_list = WHOList.from_file(who_path)

    def call(self) -> Dict[str, KelchResult]:
        """
        Perform the calling of non-synonymous mutations in the K13 gene for all samples in the VCF, and return a dict mapping sample ID → list of CodonMutation objects for that sample.

        The calling process applies the following criteria:
        1. Call non-synonymous mutation on Kelch13 gene region.
        2. Flag STOP codons (alt AA = "_") as "STOP_CODON" call status.
        3. Flag WHO validated mutations (at codon level).
        4. Calculate the coverage at WHO-associated positions 
        5. Return all called mutation to be used in the parser to write the output files.


        Returns:
            A dict mapping sample ID → list of CodonMutation objects for that sample.
        """
        samples      = self.get_sample_ids()
        logging.info(f"Calling non-synonymous mutations for {len(samples)} sample(s)...")
        aa_mutations = self.call_non_synonymous(qc=self.qc)
        who_covs     = self._who_list_coverage_all(samples)

        results: Dict[str, KelchResult] = {}

        for sample in samples:
            who_cov = who_covs[sample]
            called_muts = []
            for mut in aa_mutations[sample]:
                if mut.alt_aa == "_":
                    call_status = "STOP_CODON"
                else:
                    call_status = "PASS"

                is_who_validated = self.who_list.is_confirmed(mut.base_label)
                mut = KelchMutation.from_nonsyn(mut, is_who_validated=is_who_validated, call_status=call_status)
                called_muts.append(mut)
            
            result = KelchResult(
                sample_id = sample,
                mutations = called_muts,
                who_cov   = who_cov,
            )
            results[sample] = result

        return results
    

    def _who_list_coverage_all(self, samples: List[str]) -> Dict[str, WHOCov]:
        """
        Compute WHO-position coverage for all samples in a single pass over the
        VCF.  Each codon window is fetched exactly once (M fetches total) rather
        than once per sample per codon (N×M fetches).

        Args:
            samples: ordered list of sample IDs to evaluate.

        Returns:
            A dict mapping sample ID → WHOCov.
        """
        # --- precompute sample-invariant data once ---
        who_aa_pos: Set[int] = self.who_list.get_aa_mutation_map()
        total_positions = len(who_aa_pos) * 3

        # aa_pos → RefCodon  (codon coords are fixed, no need to re-query per sample)
        codon_map = {
            aa_pos: self._cds_reader.query_aa_num(aa_pos)  # type: ignore
            for aa_pos in who_aa_pos
        }

        # accumulators keyed by sample
        covered:     Dict[str, int]                    = {s: 0 for s in samples}
        missing_pos: Dict[str, Dict[int, List[int]]]   = {s: defaultdict(list) for s in samples}

        # --- single pass: one VCF fetch per codon window ---
        for aa_pos, codon in codon_map.items():
            coord_set: Set[int] = set(codon.codon_coords)
            records = list(self._vcf.fetch(codon.chrom, min(coord_set) - 1, max(coord_set)))  # type: ignore

            for sample in samples:
                uncovered = set(coord_set)  # reset per sample
                for rec in records:
                    dp = rec.samples[sample].get("DP")
                    if dp is not None and dp > self.qc.min_dp and rec.pos in uncovered:
                        covered[sample] += 1
                        uncovered.discard(rec.pos)
                if uncovered:
                    missing_pos[sample][aa_pos].extend(sorted(uncovered))

        return {
            sample: WHOCov(
                perc_cov    = covered[sample] / total_positions if total_positions > 0 else 0.0,
                missing_pos = dict(missing_pos[sample]),
            )
            for sample in samples
        }

# ========== Output Writer Implementation ===========================================

class K13MutationWriter(MutationWriter):
    def __init__(self, results: Dict[str, KelchResult], ignore_blank_samples: bool = False):
        self.results              = results
        self.ignore_blank_samples = ignore_blank_samples

    def write_compact(self, path: str) -> None:
        """
        Write two column of sample_id and mutations in a compact format where multiple mutations for a sample are comma-joined.
        For samples with no passing mutations, report "-" unless ignore_blank_samples is True.
        The call_status values considered as "passing mutations" for the compact output are "PASS_WHO", "PASS"

        Example:
        | sample_id | mutations       |
        |-----------|-----------------|
        | sample1   | C580Y           |
        | sample2   | A578S, R539T*   |
        | sample3   | -               |  (if ignore_blank_samples=False)

        Args:
            path: the file path to write the compact TSV output.
        """
        rows = []
        for _, k13_result in sorted(self.results.items()):
            who_missing_format = ",".join(sorted(map(str, k13_result.who_cov.missing_pos.keys()))) if k13_result.who_cov.missing_pos else "-"
            pass_muts = [m for m in k13_result.mutations]
            if self.ignore_blank_samples and not pass_muts:
                continue
            if pass_muts:
                # Sort mutations by heterozygosity status (non-het first) and then by amino acid position, and join them with commas. Het mutations are marked with an asterisk (*).
                sorted_pass_muts = [m.het_label for m in pass_muts if not m.is_het] + [m.het_label for m in pass_muts if m.is_het]
                rows.append({"sample_id": k13_result.sample_id, "mutations": ", ".join(sorted_pass_muts), "who_cov_perc": k13_result.who_cov.perc_cov, "missing_validated_muts": who_missing_format})
            else:
            # No passing mutations → report "-" (unless ignoring blank samples)
                if k13_result.who_cov.perc_cov < 1.0:
                    rows.append({"sample_id": k13_result.sample_id, "mutations": "-", "who_cov_perc": k13_result.who_cov.perc_cov, "missing_validated_muts": who_missing_format})
                else:
                    rows.append({"sample_id": k13_result.sample_id, "mutations": "WT", "who_cov_perc": k13_result.who_cov.perc_cov, "missing_validated_muts": who_missing_format})

        pd.DataFrame(rows, columns=["sample_id", "mutations", "who_cov_perc", "missing_validated_muts"]).to_csv(path, sep="\t", index=False)
        logging.info(f"Compact table written to {path}  ({len(rows)} samples)")


    def write_long(self, path: str) -> None:
        """
        Write a long-format TSV with one row per sample, mutation, and nucleotide change).
        
        The columns include:
        - sample_id: the sample ID
        - aa_position: the amino acid position of the mutation
        - ref_aa: the reference amino acid
        - alt_aa: the alternate amino acid
        - mutation: the mutation label in the format "RefAA{position}AltAA"
        - mutation_with_het: the mutation label with het status (e.g. "C580Y*" for het)
        - is_het: whether the mutation is heterozygous
        - is_who_validated: whether the mutation is in the WHO-validated list
        - call_status: the call status of the mutation (e.g. PASS_WHO, PASS, FAIL_STRICT_QC, STOP_CODON, WT)
        - chrom: the chromosome of the nucleotide change
        - nt_pos: the genomic position of the nucleotide change
        - ref_nt: the reference nucleotide of the change
        - alt_nt: the alternate nucleotide of the change 
        - dp: the total depth at the nucleotide position
        - ad: the allele depth string for all alleles at the nucleotide position
        - who_cov_perc: the coverage percentage at WHO-associated positions for this sample
        - missing_validated_muts: a comma-joined string of amino acid positions with missing calls that are associated with WHO-validated mutations, or "-" if none are missing.

        Example:
        | sample_id | aa_position | ref_aa | alt_aa | mutation | mutation_with_het | is_het | is_who_validated | call_status | chrom       | nt_pos  | ref_nt | alt_nt | dp  | ad    | who_cov_perc | missing_validated_muts |
        |-----------|-------------|--------|--------|----------|-------------------|--------|------------------|------------|--------------|---------|--------|--------|-----|-------|--------------|-----------------------|
        | sample1   | 580         | C      | Y      | C580Y    | C580Y*            | True   | True             | PASS       | Pf3D7_13_v3  | 1724816 | G      | A      | 100 | 10    | 0.95         | 527,537,538           |
        | sample2   | 578         | A      | S      | A578S    | A578S             | False  | False            | PASS       | Pf3D7_13_v3  | 1724816 | A      | G      | 80  | 12    | 1.0          | -                     |
        | sample3   | 539         | R      | _      | R539_    | R539_*            | True   | False            | STOP_CODON | Pf3D7_13_v3  | 1724817 | C      | A      | 90  | 15    | 0.85         | 441,446,449           |
        
        """
        cols = [
            "sample_id", "aa_position", "ref_aa", "alt_aa",
            "mutation", "mutation_with_het", "is_het", "is_who_validated", "call_status",
            "chrom", "nt_pos", "ref_nt",
            "alt_nt", "dp", "ad",
            "who_cov_perc", "missing_validated_muts",
        ]
        if not self.results:
            logging.warning("No mutations to write (long format).")
            pd.DataFrame(columns=cols).to_csv(path, sep="\t", index=False)
            return
        

        rows = []
        for _, k13_result in sorted(self.results.items()):
            pass_muts = [m for m in k13_result.mutations]
            for mut in [m for m in pass_muts]:  
                allele_change = mut.allele_info
                who_missing_format = ",".join(sorted(map(str, k13_result.who_cov.missing_pos.keys()))) if k13_result.who_cov.missing_pos else "-"
                base = {
                    "sample_id":         k13_result.sample_id,
                    "aa_position":       mut.aa_number,
                    "ref_aa":            mut.ref_aa,
                    "alt_aa":            mut.alt_aa,
                    "mutation":          mut.base_label,
                    "mutation_with_het": mut.het_label,
                    "is_het":            mut.is_het,
                    "is_who_validated":  mut.is_who_validated,
                    "call_status":       mut.call_status,
                    "chrom":             allele_change.chrom,
                    "nt_pos":            allele_change.pos,
                    "ref_nt":            allele_change.REF,
                    "alt_nt":            allele_change.alleles[0],
                    "dp":                allele_change.dp,
                    "ad":                allele_change.allele_depths.get(allele_change.alleles[0]),
                    "who_cov_perc":      k13_result.who_cov.perc_cov,
                    "missing_validated_muts": who_missing_format
                }
                rows.append(base)

        pd.DataFrame(rows, columns=cols).to_csv(path, sep="\t", index=False)
        logging.info(f"Long table written to {path}  ({len(rows)} rows)")


# # ---------------------------------------------------------------------------
# # Calling from tab-delimited genotype file 
# # ---------------------------------------------------------------------------

# class K13TabularCaller(K13CallerBase):
#     """
#     Kelch13 mutation caller
#     =======================
#     The K13TabularCaller calls a non-synonymous mutation in the K13 gene through a certain set of QC criteria while prioritising WHO-validated mutations list.

#     The caller takes as input a tab-delimited file generate from write_genotype_files.py with in the AmpRecon repository.
#     The file must contains at least the following specifics columns:

#     | ID       | Amplicon         | Pos      | Chr          | Loc      | Gen     | Depth   | Filt  |
#     |----------|------------------|----------|--------------|----------|---------|---------|-------|
#     | sample1  | K13_resistance_1 | 1        | Pf3D7_13_v3  | 1724817  | A       | 100     | PASS  |
#     | sample1  | K13_resistance_1 | 2        | Pf3D7_13_v3  | 1724816  | G       | 90      | PASS  |
#     | sample2  | K13_resistance_1 | 1        | Pf3D7_13_v3  | 1724817  | A,G     | 80      | PASS  |
#     | sample2  | K13_resistance_1 | 2        | Pf3D7_13_v3  | 1724816  | A       | 70      | PASS  |

#     The current implementation of select criteria for calling a non-synonymous mutation in the K13 gene is as follows:
#     1. Drop any stop codons (alt AA = "_").
#     2. Accept any non-synonymous mutation that is in the WHO-validated list.
#     3. If not in the WHO list, apply stricter QC thresholds at the allele level.

#     By default, the QC thresholds are set as follows:
#     - Standard QC: min_dp=10, min_ad=5, min_ratio=0.05
#     - Strict QC: min_dp=20, min_ad=10, min_ratio=0.1
    
#     Attributes:
#         genotype_path: the file path to the input tab-delimited genotype file.
#         fasta_path: the file path to the reference FASTA file for K13 gene.
#         qc: the QCThresholds object containing the standard and strict QC thresholds.
#         who_list: the WHOList object containing the WHO-validated K13 mutations.
#     """

#     def __init__(self, genotype_path: str, fasta_path: str, qc: QCThresholds,
#                  who_list: WHOList, **kwargs):
#         super().__init__(fasta_path, qc, who_list, **kwargs)
#         self.genotype_path         = genotype_path
#         self._samples, self._index = self._load_genotypes(genotype_path)

#     def _load_genotypes(self, path: str):
#         """
#         Load the input tab-delimited genotype file into memory and build an index for quick allele queries.
        
#         Args:
#             path: The file path to the input tab-delimited genotype file.
#         Returns:
#             A tuple containing:
#             - A list of sample IDs in the order they appear in the input file.
#             - A dict mapping (sample ID, genomic position) → (genotype string, depth string, filter string).
#         """
#         df = pd.read_csv(path, sep="\t", dtype=str)
#         required = {"ID", "Amplicon", "Pos", "Chr", "Loc", "Gen", "Depth", "Filt"}
#         missing  = required - set(df.columns)
#         if missing:
#             logging.error(f"Genotypes file missing required columns: {missing}")
#             sys.exit(255)

#         samples: List[str] = list(dict.fromkeys(df["ID"].tolist()))

#         index: Dict[tuple, tuple] = {}
#         k13_df = df[df["Chr"] == self.chrom]
#         for row in k13_df.itertuples(index=False):
#             key         = (str(row.ID), int(row.Loc))
#             index[key]  = (str(row.Gen), str(row.Depth), str(row.Filt).upper())

#         logging.info(
#             f"Loaded {len(df)} genotype rows; "
#             f"{len(samples)} unique sample(s); "
#             f"{len(k13_df)} row(s) on {self.chrom}."
#         )
#         return samples, index

#     def _get_sample_ids(self) -> List[str]:
#         """
#         Return a stable-ordered list of sample IDs from the input genotype file.

#         Returns:
#             A list of sample IDs in the order they appear in the input file.
#         """
#         return self._samples

#     def _get_variant_positions(self) -> Set[int]:
#         """
#         Scan all gene-region positions and get only those with alt variants.

#         Returns:
#             A set of genomic positions (1-based) within the gene region that have at least one non-REF allele in any sample.
#         """
#         gene_locs: Set[int] = {
#             loc for (_, loc) in self._index
#             if self.gene_stop <= loc <= self.gene_start
#         }
#         variant_locs: Set[int] = set()
#         for loc in gene_locs:
#             ref_base = self._fasta.fetch(self.chrom, loc - 1, loc).upper()
#             for sample in self._samples:
#                 row = self._index.get((sample, loc))
#                 if row is None:
#                     continue
#                 gen_str, _, filt = row
#                 if filt != "PASS" or gen_str.strip() == "-":
#                     continue
#                 if any(a.strip() != ref_base for a in gen_str.split(",")):
#                     variant_locs.add(loc)
#                     break
#         logging.info(
#             f"Tabular: {len(gene_locs)} position(s) in gene region → "
#             f"{len(variant_locs)} with ≥1 non-REF allele."
#         )
#         return variant_locs
        
#     # Mark for review: This function filtering allele under each QC threshold 
#     # This function is called again with more stringent QC threshold when mutation is not within the WHO mutation list.
#     def _allele_call(self, sample: str, loc: int, ref_base: str, strict: bool = False) -> AlleleCall:
#         """
#         Get an AlleleCall for a given sample at a specific genomic position by applying QC criteria to the genotype info from the input file.
#         The QC criteria are applied in the following order:
#         1. Filter out positions with total depth below the min_dp threshold.
#         2. Filter alleles by min_ad and min_ratio thresholds.
#         3. Determine if the call is heterozygous.
#         Args:
#             sample: the sample ID to evaluate
#             loc: the genomic position (1-based) to evaluate
#             ref_base: the reference base at the given genomic position
#             strict: whether to apply strict QC thresholds (True) or standard QC thresholds (False)
#         Returns:
#             An AlleleCall object representing the allele call for the given sample and position.
#         """
#         # Set QC thresholds based on whether strict mode is enabled
#         min_dp    = self.qc.strict_min_dp    if strict else self.qc.min_dp
#         min_ad    = self.qc.strict_min_ad    if strict else self.qc.min_ad
#         min_ratio = self.qc.strict_min_ratio if strict else self.qc.min_ratio
    
#         # get the genotype info for this sample and specific position
#         row = self._index.get((sample, loc))
#         if row is None:
#             return AlleleCall(alleles=[ref_base], is_missing=True)

#         gen_str, depth_str, filt = row

#         # Expect non-PASS filter to filtered out beforehand.
#         if gen_str.strip() == "-":
#             return AlleleCall(alleles=[ref_base], is_missing=True)


#         alleles    = [a.strip() for a in gen_str.split(",")]
#         raw_depths = [d.strip() for d in depth_str.split(",")]

#         if len(alleles) != len(raw_depths):
#             logging.error(
#                 f"[{sample}] Position {loc}: Gen and Depth column lengths differ "
#                 f"({gen_str!r} vs {depth_str!r})"
#             )
#             sys.exit(255)

#         # Convert raw depth strings to integers
#         try:
#             ads = [int(d) for d in raw_depths]
#         except ValueError:
#             logging.error(f"[{sample}] Position {loc}: Non-integer depth values in '{depth_str}'")
#             sys.exit(255)

#         # QC steps -- any failure results in REF call with depth info.

#         # Total Depth
#         total_dp = sum(ads)
#         # Build allele → depth mapping 
#         all_depths = {a: d for a, d in zip(alleles, ads)}

#         # 1. Check total depth against min_dp threshold
#         if total_dp < min_dp:
#             logging.debug(f"[{sample}] Position {loc}: total DP {total_dp} below threshold {min_dp}")
#             return AlleleCall(
#                 alleles=[ref_base], 
#                 is_missing=False,
#                 is_het=False,
#                 dp=total_dp,
#                 allele_depths=all_depths
#                 )

#         # 2. Filter alleles by min_ad and min_ratio thresholds
#         passing_alleles = [
#             a for a, ad in zip(alleles, ads)
#             if ad >= min_ad and (ad / total_dp) >= min_ratio
#         ]

#         # 2.1 If no alleles pass, return REF with depth info
#         if not passing_alleles:
#             return AlleleCall(
#                 alleles=[ref_base], 
#                 is_missing=False, 
#                 is_het=False,
#                 dp=total_dp, 
#                 allele_depths=all_depths
#             )
        
#         # 3. Determine if the call is heterozygous 
#         # 3.1 If heterozygous, return all passing alleles with is_het flag.
#         is_het: bool = len(passing_alleles) > 1
#         if is_het:
#             return AlleleCall(
#                 # alleles=[ref_base] + alt_alleles, 
#                 alleles=passing_alleles,
#                 is_missing=False, 
#                 is_het=True,
#                 dp=total_dp, 
#                 allele_depths=all_depths
#             )
#         # 3.2 If homozygous ALT, return the ALT allele(s) with is_het=False.
#         return AlleleCall(
#             alleles=passing_alleles, 
#             is_missing=False, 
#             is_het=False,
#             dp=total_dp, 
#             allele_depths=all_depths
#         )
