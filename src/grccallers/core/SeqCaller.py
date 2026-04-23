#!/usr/bin/env python3
"""
SeqCaller
==================
Tool to do the calling at alleles and codon level 
This require a refernce FASTA and a GFF annotation file to build the mapping between genomic positions and codon positions, 
and a bgzipped+tabix-indexed VCF/BCF file to query the allele information at each position.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
import itertools
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Set, Tuple

import pysam
import pandas as pd
from pysam import GTFProxy

# =======Constants and string helper=========================

# Standard codon table
# https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "_", "TAG": "_",
    "TGT": "C", "TGC": "C", "TGA": "_", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", 
}

IUPAC_CODE: Dict[frozenset, str] = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AC"): "M",
    frozenset("AG"): "R",
    frozenset("AT"): "W",
    frozenset("CG"): "S",
    frozenset("CT"): "Y",
    frozenset("GT"): "K",
    frozenset("ACG"): "V",
    frozenset("ACT"): "H",
    frozenset("AGT"): "D",
    frozenset("CGT"): "B",
    frozenset("ACGT"): "N",
}

_RC_TRANS = str.maketrans("ATCGNatcgn", "TAGCNtagcn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TRANS)[::-1]


# ======= Data Structs ======================================

@dataclass
class QCThresholds:
    """
    Object to hold the QC thresholds setting for accepting alleles at codon positions:

    min_dp: minimum total depth at the position to consider any allele calls
    min_ad: minimum allelic depth for an allele to be considered as valid
    min_ratio: minimum depth ratio (allelic depth / total depth) for an allele to be considered present

    """
    min_dp:   int   = 10
    min_ad:   int   = 5
    min_ratio: float = 0.10

@dataclass
class NonSynMutation:
    """
    Struct to hold a single non-synonymous amino acid change,
    with the allele-level evidence from the variant call.
    
    aa_number:    the amino acid position number on gene (1-based) 
    ref_aa:       the reference amino acid at this position
    alt_aa:       the alternative amino acid at this position
    is_het:       whether the mutation is called as heterozygous
    allele_info:  the AlleleCall evidence for the alt allele that produced this change
    """
    aa_number:    int
    ref_aa:       str
    alt_aa:       str
    is_het:       bool
    allele_info:  AlleleCall

    # Print amino-acid change in the format "C580Y" or "C580Y*" 
    @property
    def het_label(self) -> str:
        suffix = "*" if self.is_het else ""
        return f"{self.ref_aa}{self.aa_number}{self.alt_aa}{suffix}"

    @property
    def base_label(self) -> str:
        return f"{self.ref_aa}{self.aa_number}{self.alt_aa}"

@dataclass
class SeqCall:
    """
    Struct to hold the consensus sequence call for one sample at one codon position.
        chrom:  chromosome name
        start:  start position (1-based, inclusive)
        stop:   stop position (1-based, inclusive)
        seq: str
    """
    sampleid: str
    chrom: str
    start: int
    stop: int
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def as_fasta_header(self) -> str:
        return f">{self.sampleid}:{self.chrom}:{self.start}-{self.stop}"
    
    @property
    def as_fasta(self) -> str:
        return f"{self.as_fasta_header}\n{self.seq}"

@dataclass
class AlleleCall:
    """
    Struct to hold the single allele call at one genomic position under a certain QC thresholds.
        alleles:       list of alleles that pass QC at this position (e.g. ["A"] or ["A", "T"])
        is_missing:    whether the call is missing (e.g. due to low depth)
        is_het:        whether the call is heterozygous (i.e. contains multiple alleles)
        dp:            total depth at this position (if available)
        allele_depths: dict mapping allele → depth (if available) (e.g. {"A": 10, "T": 5})
    """
    chrom:         str            
    pos:           int           
    REF:           str            
    alleles:       List[str]      = field(default_factory=list)
    is_missing:    bool           = False
    is_het:        bool           = False
    dp:            Optional[int]  = None
    allele_depths: Dict[str, int] = field(default_factory=dict)  # allele → depth



@dataclass
class Feature:
    """
    Struct to hold information about a feature from gff including sequence.
        geneid: gene identifier
        chrom:  chromosome name
        start:  start position (1-based, inclusive)
        stop:   stop position (1-based, inclusive)
        strand: "+" or "-"
        phase:  phase of the feature (0, 1, 2)
        seq:    the sequence of the feature
    """
    geneid: str
    chrom: str
    source: str
    type: str
    start: int
    stop: int
    strand: str  
    phase: Optional[int]
    seq: str
    
    @classmethod
    def from_gff(cls, data: GTFProxy, seq: str) -> "Feature": 
        """Creates a Feature instance from a GFF row."""
        return cls(
            geneid=cls._get_ID(data.attributes) or "-",
            chrom=data.contig,
            source=data.source,
            type=data.feature,
            start=int(data.start) + 1, # Convert PySam 0-based to 1-based start
            stop=int(data.end),
            strand=data.strand,
            phase=int(data.frame) if data.frame is not None and data.frame != '.' else None,
            seq=seq
        )
    
    @staticmethod
    def _get_ID(attributes: str) -> Optional[str]:
        """Extracts the gene ID from the GFF attributes field."""
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("ID="):
                return attr[3:]
        return None

@dataclass
class RefCodon: 
    """
    Struct to hold the reference codon information for a given genomic position.
    """
    chrom: str
    aa_number: int
    ref_aa: str
    codon_seq: str
    codon_coords: List[int]
    strand: str


# ========= Sequence Handling ===============================

class SeqReader():
    def __init__(self, fasta_path: str, gff_path: Optional[str]):
        self._fasta = pysam.FastaFile(fasta_path)
        self._gff   = None if gff_path is None else pysam.TabixFile(gff_path)

    def query_fasta(self, chrom: str, start: int, end: int) -> str:
        """Query reference FASTA for a given region (1-based, inclusive)."""
        # PySam fasta.fetch uses 0-based start, 0-based exclusive end
        return self._fasta.fetch(chrom, start - 1, end).upper()
    
    def query_gff(self, chrom: str, start: int, end: int) -> List[Feature]:
        """Query features overlapping a region (1-based, inclusive)."""
        results = []
        if self._gff is None:
            return []
        try:
            for gff in self._gff.fetch(reference=chrom, start=start - 1, end=end, parser=pysam.asGFF3()):
                feat_seq = self.query_fasta(chrom, int(gff.start) + 1, int(gff.end))
                results.append(Feature.from_gff(gff, feat_seq))
        except ValueError:
            pass
        return results
    
    def fetch_geneid(self, chrom: Optional[str], gene_id: str) -> Optional[Feature]:
        """Get the start/stop coordinates for a given gene ID."""
        if self._gff is None:
            return None
            
        if chrom is None:
            for contig in self._gff.contigs:
                 gene_feat = self.fetch_geneid(contig, gene_id)
                 if gene_feat is not None:
                     return gene_feat
            return None
            
        for gff in self._gff.fetch(reference=chrom, parser=pysam.asGFF3()):
            attr_id = Feature._get_ID(gff.attributes)
            if gff.feature == "gene" and attr_id == gene_id:
                feat_seq = self.query_fasta(gff.contig, int(gff.start) + 1, int(gff.end))
                return Feature.from_gff(gff, feat_seq)
        return None
    
    def fetch_cds(self, chrom: Optional[str], id_str: str) -> Optional[List[Feature]]:
        """
        Get all CDS exons for a given gene or mRNA ID in a single GFF contig scan.

        Accepts either:
          - a gene ID (e.g. "PF3D7_1343700"): mRNA records whose Parent matches
            id_str are collected first; CDS records whose Parent is one of those
            mRNA IDs are then accepted.
          - an mRNA ID (e.g. "PF3D7_1343700.1"): CDS records whose Parent matches
            id_str directly are accepted.

        Both cases are resolved in one pass over the contig, so there is no
        secondary lookup.  When chrom is None the method iterates all contigs and
        returns on the first match.

        Returns a list of Feature objects sorted by genomic start position, or
        None if no matching CDS records are found.
        """
        if self._gff is None:
            return None

        if chrom is None:
            for contig in self._gff.contigs:
                result = self.fetch_cds(contig, id_str)
                if result is not None:
                    return result
            return None

        # Seed with id_str so it acts as a direct mRNA match as well as a gene match.
        transcript_ids: Set[str] = {id_str}
        cds_feats: List[Feature] = []

        for gff in self._gff.fetch(reference=chrom, parser=pysam.asGFF3()):
            attrs = {k: v for k, v in (a.split("=", 1) for a in gff.attributes.split(";") if "=" in a)}

            if gff.feature == "mRNA":
                # Collect mRNA children of the given gene ID.
                if attrs.get("Parent") == id_str:
                    t_id = attrs.get("ID")
                    if t_id:
                        transcript_ids.add(t_id)

            elif gff.feature == "CDS":
                parent   = attrs.get("Parent", "")
                attr_id  = attrs.get("ID", "")
                # Accept if CDS Parent is a known transcript, or if the CDS ID
                # prefix (before first ":") matches a known transcript ID.
                if parent in transcript_ids or attr_id.split(":")[0] in transcript_ids:
                    feat_seq = self.query_fasta(gff.contig, int(gff.start) + 1, int(gff.end))
                    cds_feats.append(Feature.from_gff(gff, feat_seq))

        if cds_feats:
            return sorted(cds_feats, key=lambda f: f.start)
        return None

    def build_cds_reader(
        self,
        chrom: Optional[str] = None,
        id: Optional[str]    = None,
        start: Optional[int] = None,
        end:   Optional[int] = None,
        strand: Optional[str] = None,
        phase:  int           = 0,
    ) -> "CdsReader":
        """
        Build a CdsReader using one of two strategies:

        GFF-lookup mode (id is given):
            Resolves CDS exons from the GFF by gene or mRNA ID via fetch_cds().
            start, end, strand, and phase are not required and are ignored.
            Requires gff_path to have been supplied at SeqReader construction.

        Region mode (strand is given):
            Builds a synthetic single-exon CDS directly from the FASTA region
            [chrom : start–end] using the provided strand and phase offset.
            No GFF is required; start, end, and strand must all be supplied.

        Exactly one of id or strand must be provided.
        """
        if id is not None and strand is not None:
            raise ValueError("Provide either id (GFF lookup) or strand (region mode), not both.")
        if id is None and strand is None:
            raise ValueError("Provide either id (GFF lookup) or strand (region mode).")

        if id is not None:
            if self._gff is None:
                raise ValueError("gff_path is required when using id-based CDS lookup.")
            cds_feats = self.fetch_cds(chrom, id)
            if not cds_feats:
                raise ValueError(f"No CDS features found for id={id!r} on chrom={chrom!r}.")
            return CdsReader(cds_feats)

        # Region mode
        if chrom is None or start is None or end is None:
            raise ValueError("chrom, start, and end are required in region/strand mode.")
        seq = self.query_fasta(chrom, start, end)
        feat = Feature(
            geneid=f"{chrom}:{start}-{end}",
            chrom=chrom,
            source="user",
            type="CDS",
            start=start,
            stop=end,
            strand=strand, # type: ignore
            phase=phase,
            seq=seq,
        )
        return CdsReader([feat])

# ======= CDS Handling ===============================

class CdsReader:
    """
    Class for reading and stitching multi-exonic CDS features.
    Maintains reference-oriented genomic coordinates for VCF mapping.
    """
    def __init__(self, features: List[Feature]):
        if not features:
            raise ValueError("Feature list cannot be empty.")
            
        # 1. Sort features physically (5' to 3' on the reference genome)
        self.features = sorted(features, key=lambda f: f.start)
        self.strand = self.features[0].strand
        
        # 2. Extract base transcript ID (handles "name:exon:1" format)
        self.transcript_id = self.features[0].geneid.split(":")[0] 
        
        # 3. Build spliced sequence and exact coordinate map
        self.spliced_seq, self.genomic_coords = self._build_spliced_model()
        
        # 4. Determine initial reading frame and trim seq/coords
        self._apply_phase_offset()
        
        # 5. Translate to protein
        self.aa = self._translate()

    def _build_spliced_model(self) -> Tuple[str, List[int]]:
        """Stitches exon sequences and maps each base to a genomic coordinate."""
        seq = ""
        coords = []
        for f in self.features:
            seq += f.seq
            coords.extend(range(f.start, f.stop + 1))
        return seq, coords

    def _apply_phase_offset(self):
        """Applies the phase offset from the first transcribed exon."""
        first_exon = self.features[0] if self.strand == "+" else self.features[-1]
        phase = int(first_exon.phase) if first_exon.phase is not None else 0
        
        # To avoid slice with negative zero index
        if phase == 0:
            self.coding_seq = self.spliced_seq
            self.coding_coords = self.genomic_coords
        else:
            if self.strand == "+":
                self.coding_seq = self.spliced_seq[phase:]
                self.coding_coords = self.genomic_coords[phase:]
            else:
                self.coding_seq = self.spliced_seq[:-phase]
                self.coding_coords = self.genomic_coords[:-phase]

    def _translate(self) -> str:
        """Translates the phase-adjusted spliced sequence."""
        seq_to_translate = self.coding_seq if self.strand == "+" else reverse_complement(self.coding_seq)
        aa_seq = ""
        
        for i in range(0, len(seq_to_translate), 3):
            codon = seq_to_translate[i:i+3]
            if len(codon) != 3:
                logging.debug(f"Incomplete codon at end of {self.transcript_id}")
                break
            aa_seq += CODON_TABLE.get(codon.upper(), "?")
            
        return aa_seq

    def genomic_to_aa(self, genomic_pos: int) -> RefCodon:
        """
        Maps a genomic coordinate to its Amino Acid number and codon position.
        Returns: RefCodon dataclass instance with:
            - aa_number: the amino acid position on the gene (1-based)
            - ref_aa: the reference amino acid at this position
            - codon_seq: the reference codon sequence
            - codon_coords: the list of three genomic coordinates for this codon
            - strand: the strand of the gene ("+" or "-")
        """
        if genomic_pos not in self.coding_coords:
            raise ValueError(f"Position {genomic_pos} is not in the coding sequence (may be intronic).")
            
        idx = self.coding_coords.index(genomic_pos)
        
        if self.strand == "+":
            aa_num = (idx // 3) + 1
            start_idx = (aa_num - 1) * 3
            codon_seq = self.coding_seq[start_idx : start_idx + 3]
            codon_coords = self.coding_coords[start_idx : start_idx + 3]
        else:
            transcript_idx = (len(self.coding_coords) - 1) - idx
            aa_num = (transcript_idx // 3) + 1
            end_idx = len(self.coding_seq) - ((aa_num - 1) * 3)
            codon_seq = self.coding_seq[end_idx - 3 : end_idx]
            codon_coords = self.coding_coords[end_idx - 3 : end_idx]
        ref_aa = self.aa[aa_num - 1]

        return RefCodon(
            chrom=self.features[0].chrom,
            aa_number=aa_num,
            ref_aa=ref_aa,
            codon_seq=codon_seq,
            codon_coords=codon_coords,
            strand=self.strand
        )
    
    def query_aa_num(self, aa_num: int) -> RefCodon:
        """Get the reference AA and the 3 genomic coordinates for its codon."""
        if self.strand == "+":
            start_idx = (aa_num - 1) * 3
            codon_seq = self.coding_seq[start_idx : start_idx + 3]
            codon_coords = self.coding_coords[start_idx : start_idx + 3]
            return RefCodon(
                chrom=self.features[0].chrom,
                aa_number=aa_num,
                ref_aa=CODON_TABLE.get(codon_seq.upper(), "?"),
                codon_seq=codon_seq.upper(),
                codon_coords=codon_coords,
                strand=self.strand
            )
        else:
            end_idx = len(self.coding_seq) - ((aa_num - 1) * 3)
            codon_seq = self.coding_seq[end_idx - 3 : end_idx]
            codon_coords = self.coding_coords[end_idx - 3 : end_idx]
            return RefCodon(
                chrom=self.features[0].chrom,
                aa_number=aa_num,
                ref_aa=CODON_TABLE.get(reverse_complement(codon_seq).upper(), "?"),
                codon_seq=codon_seq.upper(),
                codon_coords=codon_coords,
                strand=self.strand
            )
    
# ======== Base Caller ===============================

class CallerBase():
    """
    Abstract base for variant-first K13 mutation callers.

    Subclasses implement three methods:
      _get_sample_ids()          → stable-ordered list of sample IDs
      _get_variant_positions()   → set of 1-based genomic positions carrying a
                                   non-REF allele in ≥1 sample within the gene
      _allele_call()             → AlleleCall for one (sample, position)

    Shared base provides:
      - codon-grouping from variant positions (_scan_codon_mutation)
      - codon evaluation via _get_nonsyn_mutation
      - public call_non_synonymous() and call_seq() entry points
    """

    def __init__(
        self,
        vcf_path: str,
        qc: QCThresholds,
        seq_reader: Optional[SeqReader] = None,
        cds_reader: Optional["CdsReader"] = None,
    ):
        """
        Args:
            vcf_path:   path to the bgzipped, tabix-indexed VCF/BCF.
            qc:         QC thresholds for allele filtering.
            seq_reader: optional SeqReader (FASTA ± GFF) — enables FASTA-level
                        access within calling methods.
            cds_reader: optional CdsReader — required for call_non_synonymous();
                        build with SeqReader.build_cds_reader() beforehand.
        """
        self.qc          = qc
        self._seq        = seq_reader
        self._cds_reader = cds_reader

        self._vcf          = pysam.VariantFile(vcf_path)
        self._active_chrom: Optional[str]               = None
        self._gene_records: Dict[int, pysam.VariantRecord] = {}

    # ------ Region helpers -----------------------------------------

    def _resolve_region(
        self,
        chrom: Optional[str],
        start: Optional[int],
        end:   Optional[int],
    ) -> Tuple[str, int, int]:
        """
        Resolve a (chrom, start, end) triple for a call method.

        Priority:
          1. Explicit values passed by the caller.
          2. CDS span from cds_reader (when available).

        Any mix of explicit and defaulted values is supported — e.g. you can
        override only start while letting end default to the CDS boundary.

        Raises ValueError when a value cannot be resolved.
        """
        if self._cds_reader is not None:
            cds_chrom = self._cds_reader.features[0].chrom
            cds_start = self._cds_reader.features[0].start
            cds_end   = self._cds_reader.features[-1].stop
            chrom = chrom if chrom is not None else cds_chrom
            start = start if start is not None else cds_start
            end   = end   if end   is not None else cds_end
        if chrom is None or start is None or end is None:
            raise ValueError(
                "chrom, start, and end must be provided explicitly when no "
                "cds_reader is set."
            )
        return chrom, start, end

    def _load_records(self, chrom: str, start: int, end: int) -> None:
        """Pre-load VCF records for (chrom, start, end) into _gene_records."""
        self._active_chrom = chrom
        self._gene_records = {
            rec.pos: rec
            for rec in self._vcf.fetch(chrom, start - 1, end)
        }

    # ------ Abstract method ----------------------------------------

    # @abstractmethod
    def get_sample_ids(self) -> List[str]:
        """Return a stable-ordered list of sample IDs in VCF records."""
        return list(self._vcf.header.samples)

    # @abstractmethod
    def get_variant_positions(self ,chrom: Optional[str] = None, start: Optional[int] = None, end: Optional[int] = None ) -> Set[int]:
        """
        Return the set of 1-based genomic positions within
        [start, end].
        """
        if chrom is not None and start is not None and end is not None:
            self._load_records(chrom, start, end)
        return { pos for pos, rec in self._gene_records.items()
                 if rec.alts and rec.alts[0] not in (None, '<NON_REF>', '.') 
               }

    # @abstractmethod
    def _allele_call(self, sample: str, pos: int, qc: Optional[QCThresholds] = None) -> AlleleCall:
        """
        Return an AlleleCall for *sample* at genomic position *loc*.
        *ref_base* is the reference nucleotide at that position (forward strand).
        """
        # 1. Get the VCF record for this position, if it exists. If not, return a missing call with the reference base as the allele.
        rec = self._gene_records.get(pos)
        if rec is None:
            return AlleleCall(
                chrom=self._active_chrom, # type: ignore
                pos=pos,
                REF="-",       # type: ignore
                alleles=[None], # type: ignore
                is_missing=True)

        # Set QC thresholds based on whether QC is applied or not.
        min_dp    =  qc.min_dp if qc else self.qc.min_dp
        min_ad    =  qc.min_ad if qc else self.qc.min_ad
        min_ratio =  qc.min_ratio if qc else self.qc.min_ratio

        # Get the total depth (DP) and allele depth (AD) for this sample at this position from the VCF record.
        dp = rec.samples[sample].get("DP")
        ad = rec.samples[sample].get("AD")

        # 2. Filter this position by the total depth
        if dp is None or dp < min_dp or ad is None:
            return AlleleCall(
                chrom=self._active_chrom, # type: ignore
                pos=pos, 
                REF=rec.ref,       # type: ignore
                alleles=[rec.ref], # type: ignore
                is_missing=True)

        # 3. Filter alleles by min_ad and min_ratio thresholds, and determine if the call is heterozygous.

        ## Create a dict mapping alleles to its allelic depth
        all_alleles = rec.alleles  # This is a tuple of all alleles (REF, ALT1, ALT2, ...)
        norm_ad     = [x if x is not None else 0 for x in ad] # Normalize AD to be 0 for any None values (e.g. missing AD for some alleles)
        all_depths  = {a: d for a, d in zip(all_alleles, norm_ad)}  # type: ignore

        ## Filter alleles by min_ad and min_ratio thresholds
        passing_idx = [
            i for i, cnt in enumerate(norm_ad)
            if cnt >= min_ad and (cnt / dp) >= min_ratio
        ]

        # If no alleles pass the QC thresholds, return a call with the reference allele and is_missing=True.
        if not passing_idx:
            return AlleleCall(
                    chrom=self._active_chrom, # type: ignore
                    pos=pos,
                    REF=rec.ref,        # type: ignore
                    alleles=[rec.ref],  # type: ignore
                    is_missing=True, 
                    is_het=False,
                    dp=dp, 
                    allele_depths=all_depths
                )

        # Determine if the call is heterozygous (i.e. contains multiple alleles passing QC).
        is_het = len(passing_idx) > 1
        passing_alleles = [all_alleles[i] for i in passing_idx] # type: ignore

        if is_het:
            return AlleleCall(
                chrom=self._active_chrom, # type: ignore
                pos=pos,
                REF=rec.ref,              # type: ignore
                alleles=passing_alleles,
                is_missing=False,
                is_het=True,
                dp=dp,
                allele_depths=all_depths,
            )

        # Return a call with the passing allele(s) and is_het status.
        return AlleleCall(
            chrom=self._active_chrom, # type: ignore
            pos=pos,
            REF=rec.ref,                 # type: ignore
            alleles=passing_alleles,
            is_missing=False,
            is_het=False,
            dp=dp,
            allele_depths=all_depths,
        )

    def _scan_codon_mutation(self) -> List[int]:
        """
        Get *variant* positions for each codon (amino-acid position).

        Returns:
            A sorted list of amino acid positions with variants.
        """
        if self._cds_reader is None:
            raise RuntimeError(
                "Codon-based operations require a CdsReader. "
                "Build one with SeqReader.build_cds_reader() and pass it as "
                "cds_reader= at CallerBase construction."
            )
        aa_numbers = set()
        for pos in self.get_variant_positions():
            aa_num = self._cds_reader.genomic_to_aa(pos).aa_number
            aa_numbers.add(aa_num)

        return sorted(aa_numbers)

    # ------Codon translation-----------------------------

    @staticmethod
    def _get_nonsyn_mutation(
        sample: str,
        codon_call: List[AlleleCall],
        ref_aa: str,
        ref_bases: str,
        strand: str,
    ) -> Dict[str, List[AlleleCall]]:
        """
        Calling amino acid changes for each allele mutations. 
        Throw a warning when multiple allele combinations can lead to more than one possible alt AAs.

        Args:
            codon_call: list of AlleleCall objects for each codon position
            ref_aa: the reference amino acid for this codon
            ref_bases: the reference nucleotide sequence for this codon
            strand: the strand of the gene ("+" or "-")

        Returns:
            A dict mapping alt amino acid → list of AlleleCall objects that produce it.
        """
        alt_aas_check = set()
        # Try every possible allele combination in the codon to see what AA it produces, and check if there are multiple alt AAs.
        for b0, b1, b2 in itertools.product(
            codon_call[0].alleles, codon_call[1].alleles, codon_call[2].alleles,
        ):
            # Handle missing calls by using reference base at that position, which is equivalent to assuming the missing call is WT at this position.
            b0 = b0 if b0 is not None else ref_bases[0]
            b1 = b1 if b1 is not None else ref_bases[1]
            b2 = b2 if b2 is not None else ref_bases[2] 

            codon = b0 + b1 + b2
            
            if len(codon) != 3:
                logging.error(f"Unexpected allele combo '{codon}' in codon translation")
                sys.exit(255)

            if strand == "-":
                codon = reverse_complement(codon)

            aa = CODON_TABLE.get(codon, "?")

            if aa != ref_aa:
                alt_aas_check.add(aa)
            
        # Throw a warning when multiple allele combinations can lead to more than one possible alt AAs, and call AA by each allele mutation.
        if len(alt_aas_check) > 1:
            logging.warning(f"[{sample}] Alleles combination in codon can produce more than one possible alt AAs ({alt_aas_check})")
            logging.warning(f"[{sample}] Calling AA by single alleles mutation..")

        alt_aas: Dict[str, List[AlleleCall]] = {}
        for pos, call in enumerate(codon_call): # Iterate through each codon position and its AlleleCall
            for alt_nt in call.alleles: # Iterate through each allele at this position that pass QC
                if alt_nt is None:
                    continue # Skip missing calls, which are treated as reference allele in the codon construction
                # Construct the codon sequence by replacing the reference base with the alt_nt at the current position, while keeping the other positions as reference.
                codon_list = list(ref_bases[i] if i != pos else alt_nt for i in range(3))
                codon = "".join(codon_list)

                if strand == "-":
                    codon = reverse_complement(codon)

                aa = CODON_TABLE.get(codon, "?")

                if aa != ref_aa:
                    # Create a new AlleleCall with only the alt allele that contributes to this AA change
                    alt_aas.setdefault(aa, []).append(replace(call, alleles=[alt_nt]))
        return alt_aas

    def call_seq(
        self,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end:   Optional[int] = None,
        qc:    Optional[QCThresholds] = None,
        iupac: bool = True,
        rc:    bool = False,
    ) -> Dict[str, SeqCall]:
        """
        Build a per-sample nucleotide sequence for a genomic region.

        Args:
            chrom:  contig name; defaults to CDS chrom when cds_reader is set.
            start:  1-based inclusive start; defaults to CDS start.
            end:    1-based inclusive end; defaults to CDS end.
            qc:     QC thresholds; falls back to the caller's default QC.
            iupac:  encode heterozygous positions with IUPAC ambiguity codes.
            rc:     reverse-complement the output sequence.
        """
        chrom, start, end = self._resolve_region(chrom, start, end)
        self._load_records(chrom, start, end)
        qc = self.qc if qc is None else qc

        samples = self.get_sample_ids()
        all_seqs: Dict[str, SeqCall] = {}

        for sample in samples:
            seq = ""
            for gpos in range(start, end + 1):
                call = self._allele_call(sample, gpos, qc=qc)

                if call.is_missing:
                    seq += "X"
                elif call.is_het:
                    seq += IUPAC_CODE.get(frozenset(call.alleles), "N") if iupac else "N"
                else:
                    seq += call.alleles[0] if call.alleles else "X"

            seq = reverse_complement(seq) if rc else seq

            all_seqs[sample] = SeqCall(
                sampleid=sample,
                chrom=chrom,
                start=start,
                stop=end,
                seq=seq
            )

        return all_seqs

    def call_non_synonymous(
        self,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end:   Optional[int] = None,
        qc:    Optional[QCThresholds] = None,
    ) -> Dict[str, List[NonSynMutation]]:
        """
        Call non-synonymous mutations across all samples.

        For each codon with at least one variant position, evaluate every sample
        for amino-acid changes under the given (or default) QC thresholds.

        Requires a CdsReader (pass cds_reader= at construction time via
        SeqReader.build_cds_reader()).

        Returns:
            A dict mapping sample ID → list of NonSynMutation objects.
            Samples with no detected mutations have an empty list.
        """
        if self._cds_reader is None:
            raise RuntimeError(
                "Codon-based operations require a CdsReader. "
                "Build one with SeqReader.build_cds_reader() and pass it as "
                "cds_reader= at CallerBase construction."
            )

        chrom, start, end = self._resolve_region(chrom, start, end)
        self._load_records(chrom, start, end)

        samples    = self.get_sample_ids()
        # Get the list of amino acid positions that have at least one variant position in their codon.
        codon_list = self._scan_codon_mutation()   # [aa_number]
        all_results: Dict[str, List[NonSynMutation]] = {s: [] for s in samples}
        qc = qc if qc is not None else self.qc

        for aa_number in codon_list:
            ref_codon: RefCodon = self._cds_reader.query_aa_num(aa_number)  

            logging.debug(
                f"Codon {aa_number} ({ref_codon.ref_aa}) {ref_codon.codon_coords[0]}-{ref_codon.codon_coords[2]}"
            )

            for sample in samples:
                codon_call: List[AlleleCall] = []
                for gpos in ref_codon.codon_coords:
                    # Evaluate the three codon positions under standard QC and translate to get alt AAs.
                    al = self._allele_call(sample, gpos, qc)
                    codon_call.append(al)
                
                # Get the alt AAs for this sample at this codon, and the information on the corresponding allele calls.
                alt_aas: Dict[str, List[AlleleCall]] = self._get_nonsyn_mutation(
                    sample=sample,
                    codon_call=codon_call,
                    ref_aa=ref_codon.ref_aa,
                    ref_bases=ref_codon.codon_seq,
                    strand=ref_codon.strand,
                )

                mutations: List[NonSynMutation] = []
                # If no amino-acid change at standard QC → skip the codon
                if not alt_aas:
                    continue

                # If there are amino-acid changes, push NonSynMutation into the list 
                for alt_aa, call_list in alt_aas.items():
                    for alleles_info in call_list:
                        mutations.append(NonSynMutation(
                            aa_number=ref_codon.aa_number,
                            ref_aa=ref_codon.ref_aa,
                            alt_aa=alt_aa,
                            is_het=alleles_info.is_het,
                            allele_info=alleles_info,
                        ))
                
                if mutations:
                    all_results[sample].extend(mutations)

        return all_results

class MutationWriter():
    """
    Abstract base class for writing mutation results to output files.
    """

    def __init__(self, results: Dict[str, List[NonSynMutation]], ignore_blank_samples: bool = False):
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
        for sample, muts in sorted(self.results.items()):
            pass_muts = [m for m in muts]
            if self.ignore_blank_samples and not pass_muts:
                continue
            if pass_muts:
                    # Sort mutations by heterozygosity status (non-het first) and then by amino acid position, and join them with commas. Het mutations are marked with an asterisk (*).
                    sorted_pass_muts = [m.het_label for m in pass_muts if not m.is_het] + [m.het_label for m in pass_muts if m.is_het]
                    rows.append({"sample_id": sample, "mutations": ", ".join(sorted_pass_muts)})
            else:
                # No passing mutations → report "-" (unless ignoring blank samples)
                rows.append({"sample_id": sample, "mutations": "-"})

        pd.DataFrame(rows, columns=["sample_id", "mutations"]).to_csv(path, sep="\t", index=False)
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
        - codon_pos: the position within the codon (1, 2, or 3) of the nucleotide change
        - ref_nt: the reference nucleotide of the change
        - alt_nt: the alternate nucleotide of the change 
        - dp: the total depth at the nucleotide position
        - ad: the allele depth string for all alleles at the nucleotide position

        Example:
        | sample_id | aa_position | ref_aa | alt_aa | mutation | mutation_with_het | is_het |  chrom       | nt_pos  | codon_pos | ref_nt | alt_nt | dp  | ad    |
        |-----------|-------------|--------|--------|----------|-------------------|--------|--------------|---------|-----------|--------|--------|-----|-------|
        | sample1   | 580         | C      | Y      | C580Y    | C580Y*            | True   | Pf3D7_13_v3  | 1724816 | 2         | G      | A      | 100 | 10    |
        | sample2   | 578         | A      | S      | A578S    | A578S             | False  | Pf3D7_13_v3  | 1724816 | 1         | A      | G      | 80  | 12    |     
        
        """
        cols = [
            "sample_id", "aa_position", "ref_aa", "alt_aa",
            "mutation", "mutation_with_het", "is_het",
            "chrom", "nt_pos", "ref_nt", 
            "alt_nt", "dp", "ad",
        ]
        if not self.results:
            logging.warning("No mutations to write (long format).")
            pd.DataFrame(columns=cols).to_csv(path, sep="\t", index=False)
            return

        rows = []
        for sample, muts in sorted(self.results.items()):
            for mut in [m for m in muts]:  
                allele_change = mut.allele_info
                base = {
                    "sample_id":         sample,
                    "aa_position":       mut.aa_number,
                    "ref_aa":            mut.ref_aa,
                    "alt_aa":            mut.alt_aa,
                    "mutation":          mut.base_label,
                    "mutation_with_het": mut.het_label,
                    "is_het":            mut.is_het,
                    "chrom":             allele_change.chrom,
                    "nt_pos":            allele_change.pos,
                    "ref_nt":            allele_change.REF,
                    "alt_nt":            allele_change.alleles[0],
                    "dp":                allele_change.dp,
                    "ad":                allele_change.allele_depths.get(allele_change.alleles[0]),
                }
                rows.append(base)

        pd.DataFrame(rows, columns=cols).to_csv(path, sep="\t", index=False)
        logging.info(f"Long table written to {path}  ({len(rows)} rows)")


class FastaWriter:
    """
    Class for writing per-sample FASTA sequences to output files.
    """

    def __init__(self, seqs: Dict[str, SeqCall]):
        self.seqs = seqs

    def write_fasta(self, path: str) -> None:
        """
        Write per-sample sequences to a multi-FASTA file.

        Args:
            path: The file path to write the multi-FASTA file.
        """
        with open(path, "w") as f:
            for sample, seq_call in self.seqs.items():
                f.write(f"{seq_call.as_fasta}\n")


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

    
# # ---------------------------------------------------------------------------
# # Output Writer
# # ---------------------------------------------------------------------------


# # ---------------------------------------------------------------------------
# # CLI
# # ---------------------------------------------------------------------------

# class CustomFormatter(argparse.RawDescriptionHelpFormatter,
#                       argparse.ArgumentDefaultsHelpFormatter):
#     pass


# def parse_args(argv=None) -> argparse.Namespace:
#     parser = argparse.ArgumentParser(
#         prog="k13_caller.py",
#         description="""
# Kelch13 Mutation Caller (Variant-First)
# =======================================
# Identifies non-synonymous amino-acid changes in the K13 propeller domain
# from either a multi-sample VCF or a pre-called tabular genotype TSV.

# Call-status values in the long-format output
# --------------------------------------------
#   PASS_WHO       WHO-confirmed mutation accepted at standard QC thresholds.
#   PASS           Novel mutation that survived the stricter QC re-evaluation.
#   FAIL_STRICT_QC Novel mutation visible at standard QC but below strict QC.
#   STOP_CODON     Variant produces a stop codon — excluded from results but
#                  preserved in the long-format TSV for audit purposes.

# Only mutations with status PASS_WHO or PASS appear in the compact output.
# All statuses are written to the long-format output for auditibility.

# Output files
# ------------
#   <prefix>_compact.tsv  one row per sample, passing mutations comma-joined
#   <prefix>_long.tsv     one row per (sample, mutation, nucleotide change),
#                         includes all call_status values
#         """,
#         epilog=(
#             "Examples:\n"
#             "  # VCF input\n"
#             "  python k13_caller.py \\\n"
#             "      --vcf       <vcf_file> \\\n"
#             "      --fasta     <fasta_file> \\\n"
#             "      --who-list  <who_list_file> \\\n"
#             "      --output-prefix <output_prefix>\n"
#             "\n"
#             "  # Pre-called tabular genotype file\n"
#             "  python k13_caller.py \\\n"
#             "      --genotype-file  <genotype_file> \\\n"
#             "      --fasta         <fasta_file> \\\n"
#             "      --who-list      <who_list_file> \\\n"
#             "      --output-prefix <output_prefix>\n"
#         ),
#         formatter_class=CustomFormatter,
#     )

#     req = parser.add_argument_group(
#         "required inputs",
#         "Provide exactly one of --vcf or --genotype-file."
#     )
#     req.add_argument(
#         "--vcf", default=None, metavar="PATH",
#         help=(
#             "bgzipped, tabix-indexed multi-sample VCF/BCF covering the K13 "
#             "amplicon region.  All samples in the VCF header are processed. "
#             "Mutually exclusive with --genotype-file."
#         ),
#     )
#     req.add_argument(
#         "--genotype-file", default=None, metavar="PATH",
#         help=(
#             "Tab-delimited pre-called genotype TSV with columns: "
#             "ID / Amplicon / Pos / Chr / Loc / Gen / Depth / Filt. "
#             "Gen and Depth must be comma-separated per-allele values (e.g. "
#             "Gen=T,A  Depth=120,15 for a heterozygous call). "
#             "Filt must be 'PASS' or '-' for a missing call. "
#             "Mutually exclusive with --vcf."
#         ),
#     )
#     req.add_argument(
#         "--fasta", required=True, metavar="PATH",
#         help=(
#             "Reference genome FASTA indexed with 'samtools faidx'.  Only the "
#             "K13 chromosome (default Pf3D7_13_v3) is accessed."
#         ),
#     )
#     req.add_argument(
#         "--who-list", required=True, metavar="PATH",
#         help=(
#             "Plain-text file of WHO-validated K13 partial-resistance mutations, "
#             "one per line in standard amino-acid notation (e.g. C580Y).  Lines "
#             "starting with '#' and blank lines are ignored.  Mutations on this "
#             "list are accepted under standard QC; all others require strict QC."
#         ),
#     )

#     out = parser.add_argument_group("output")
#     out.add_argument(
#         "--output-prefix", default="k13_output", metavar="PREFIX",
#         help=(
#             "Filename prefix for the output TSV files.  For example, "
#             "'k13_results' produces k13_results_compact.tsv "
#             "and k13_results_long.tsv."
#         ),
#     )

#     std = parser.add_argument_group(
#         "standard QC thresholds",
#         "Applied to all WHO-confirmed mutations and to the initial screen of novel mutations.",
#     )
#     std.add_argument(
#         "--min-dp", type=int, default=10, metavar="N",
#         help="Minimum total read depth (DP FORMAT field) at the position.",
#     )
#     std.add_argument(
#         "--min-ad", type=int, default=5, metavar="N",
#         help="Minimum allele depth (AD) required for an allele to be considered.",
#     )
#     std.add_argument(
#         "--min-ratio", type=float, default=0.10, metavar="F",
#         help="Minimum allele frequency (AD / DP) for an allele to pass.",
#     )

#     strict = parser.add_argument_group(
#         "strict QC thresholds",
#         (
#             "Applied as a second-pass filter to non-WHO mutations.  A novel "
#             "mutation that passes standard QC but fails strict QC is labelled "
#             "FAIL_STRICT_QC in the long-format output and excluded from "
#             "compact/wide outputs."
#         ),
#     )
#     strict.add_argument(
#         "--strict-min-dp", type=int, default=20, metavar="N",
#         help="Minimum total depth for strict QC re-evaluation.",
#     )
#     strict.add_argument(
#         "--strict-min-ad", type=int, default=10, metavar="N",
#         help="Minimum allele depth for strict QC re-evaluation.",
#     )
#     strict.add_argument(
#         "--strict-min-ratio", type=float, default=0.20, metavar="F",
#         help="Minimum allele frequency for strict QC re-evaluation.",
#     )

#     gene = parser.add_argument_group(
#         "K13 gene coordinates",
#         (
#             "1-based genomic coordinates of the K13 propeller domain in the "
#             "reference assembly.  gene_start is the HIGH (5') boundary and "
#             "gene_stop is the LOW (3') boundary; the gene is on the reverse strand."
#         ),
#     )
#     gene.add_argument(
#         "--chrom", default="Pf3D7_13_v3",
#         help="Chromosome/contig name in the FASTA and VCF headers.",
#     )
#     gene.add_argument(
#         "--gene-start", type=int, default=1726997,
#         help="High (5' / upstream) boundary of the propeller domain, 1-based.",
#     )
#     gene.add_argument(
#         "--gene-stop", type=int, default=1724817,
#         help="Low (3' / downstream) boundary of the propeller domain, 1-based.",
#     )

#     parser.add_argument(
#         "--log-level", default="INFO",
#         choices=["DEBUG", "INFO", "WARNING", "ERROR"],
#         help="Verbosity level.  Use DEBUG to trace individual codon/sample decisions.",
#     )

#     return parser.parse_args(argv)


# def main(argv=None) -> None:
#     args = parse_args(argv)

#     logging.basicConfig(
#         level=getattr(logging, args.log_level),
#         format="%(asctime)s [%(levelname)s] %(message)s",
#         datefmt="%H:%M:%S",
#     )

#     if args.vcf and args.genotype_file:
#         logging.error("Provide either --vcf or --genotype-file, not both.")
#         sys.exit(1)
#     if not args.vcf and not args.genotype_file:
#         logging.error("One of --vcf or --genotype-file is required.")
#         sys.exit(1)

#     for flag, path in [("--fasta", args.fasta), ("--who-list", args.who_list)]:
#         if not Path(path).exists():
#             logging.error(f"{flag}: file not found: {path}")
#             sys.exit(1)

#     input_path = args.vcf or args.genotype_file
#     input_flag = "--vcf" if args.vcf else "--genotype-file"
#     if not Path(input_path).exists():
#         logging.error(f"{input_flag}: file not found: {input_path}")
#         sys.exit(1)

#     qc = QCThresholds(
#         min_dp=args.min_dp, min_ad=args.min_ad, min_ratio=args.min_ratio,
#         strict_min_dp=args.strict_min_dp, strict_min_ad=args.strict_min_ad,
#         strict_min_ratio=args.strict_min_ratio,
#     )
#     who_list = WHOList.from_file(args.who_list)

#     shared_kwargs = dict(
#         fasta_path=args.fasta,
#         qc=qc,
#         who_list=who_list,
#         chrom=args.chrom,
#         gene_start=args.gene_start,
#         gene_stop=args.gene_stop,
#     )

#     if args.genotype_file:
#         logging.info("Input mode: tabular genotype file (variant-first)")
#         caller = K13TabularCaller(genotype_path=args.genotype_file, **shared_kwargs)
#     else:
#         logging.info("Input mode: VCF (variant-first)")
#         caller = K13VcfCaller(vcf_path=args.vcf, **shared_kwargs)

#     logging.info("Starting K13 mutation calling …")
#     results = caller.call()
#     logging.info(f"Calling complete. {len(results)} sample(s) processed.")

#     writer = OutputWriter(results)
#     writer.write_compact(f"{args.output_prefix}_compact.tsv")
#     writer.write_long(f"{args.output_prefix}_long.tsv")

#     logging.info("Done.")


# if __name__ == "__main__":
#     main()
