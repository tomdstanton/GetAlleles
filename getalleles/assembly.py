"""
Copyright 2024 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/GetAlleles

This file is part of GetAlleles. GetAlleles is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. GetAlleles is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with GetAlleles.
If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import annotations

import os
from pathlib import Path
from subprocess import Popen, PIPE, DEVNULL
from typing import Generator, TextIO, BinaryIO
from re import compile

from getalleles.alignment import Alignment
from getalleles.utils import check_file, parse_fasta, translate, load_ttable, reverse_complement
from getalleles.log import log, warning, quit_with_error

# Constants -----------------------------------------------------------------------------------------------------------
_ASSEMBLY_FASTA_REGEX = compile(r'\.(fasta|fa|fna|ffn)(\.(gz|xz|bz2))?$')
_IUPAC_AA = set('ACDEFGHIKLMNPQRSTVWY')  # IUPAC standard protein alphabet
_IUPAC_UNAMBIGUOUS_DNA = set('GATC')  # IUPAC standard unambiguous nucleotide alphabet


# Classes -------------------------------------------------------------------------------------------------------------
class Assembly:
    def __init__(self, path: Path | None = None, name: str | None = None, contigs: dict[str: Contig] | None = None):
        self.path = path or Path()
        self.name = name or path.name.strip('.gz').rsplit('.', 1)[0] if path else ''
        self.contigs = contigs or {}
        self.alleles = []

    def __str__(self):
        return self.name

    def __len__(self):
        return sum(len(i) for i in self.contigs.values())

    def format(self, format_spec, **kwargs):
        if format_spec == 'fna':
            return ''.join(i.format(format_spec) for i in self.contigs.values())
        if format_spec in {'ffn', 'faa', 'tsv'}:
            return ''.join(i.format(format_spec, **kwargs) for i in self.alleles)

    def map(self, refs: References, threads: int, extra_args: str = '', verbose: bool = False
            ) -> Generator[Alignment, None, None]:
        fmt, cmd = ('ffn', "minimap2 -c ") if refs.mol == 'DNA' else ('faa', f"miniprot -T {refs.T} ")
        cmd += f"{extra_args} " if extra_args else ''
        cmd += f'-t {threads} {self.path} -'
        log(f"{cmd=}", verbose=verbose)
        with Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=DEVNULL, universal_newlines=True) as p:  # Init process
            for paf_line in p.communicate(refs.format(fmt))[0].splitlines():
                yield Alignment.from_paf_line(paf_line)  # Yield the alignment as an Alignment object

    def write(self, tsv: TextIO | None = None, ffn: Path | TextIO | None = None, faa: Path | TextIO | None = None, **kwargs):
        if self.alleles:
            for fh, fmt in ((tsv, 'tsv'), (ffn, 'ffn'), (faa, 'faa')):
                if fh:
                    (fh / f'{self}_alleles.{fmt}').write_text(self.format(fmt, **kwargs)) if isinstance(fh, Path) else \
                        fh.write(self.format(fmt, **kwargs))
        else:
            return warning(f"No alleles found in {self}")


class Contig:
    """
    This class describes a contig in an assembly: the name, length, and sequence.
    """

    def __init__(self, id_: str | None = None, desc: str | None = None, seq: str | None = None):
        self.id = id_ or ''
        self.desc = desc or ''
        self.seq = seq
        # self.neighbours_L = neighbours_L or {}  # For assembly graphs
        # self.neighbours_R = neighbours_R or {}  # For assembly graphs

    def __str__(self):
        return self.id

    def __len__(self):
        return len(self.seq)

    def format(self, format_spec):
        if format_spec == 'fna':
            return f">{self.id} {self.desc}\n{self.seq}\n"


class Allele:
    def __init__(self, ref: Reference, assembly: Assembly, contig: Contig, start: int, end: int, strand: int,
                 perc_id: float, perc_cov: float, cigar: str, dna_hash: callable, aa_hash: callable, **kwargs):
        self.ref = ref
        self.assembly = assembly
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.perc_id = perc_id  # Percent identity of the allele
        self.perc_cov = perc_cov  # Percent coverage of the allele
        self.cigar = cigar
        self.dna_seq = self.contig.seq[self.start:self.end] if self.strand == "+" else reverse_complement(
            self.contig.seq[self.start:self.end])
        if len(self.dna_seq) % 3 != 0:
            warning(f"DNA sequence for {self} is not a multiple of 3")
        self.aa_seq = translate(self.dna_seq, **kwargs)
        if len(self.aa_seq) < 2:
            warning(f"Protein sequence for {self} is less than 2 residues")
        if len(self.dna_seq) < len(self.ref.dna_seq) and any((self.start == 0, self.end == len(self.contig))):
            self.problems = 'Partial'
        elif len(self.aa_seq) < len(self.ref.aa_seq):
            self.problems = 'Truncated'
        else:
            self.problems = ''
        self.dna_hash = dna_hash(self.dna_seq.encode(), usedforsecurity=False).hexdigest()
        self.aa_hash = aa_hash(self.aa_seq.encode(), usedforsecurity=False).hexdigest()

    def __str__(self):
        return f"{self.ref.id} {self.assembly}|{self.contig}:{self.start}-{self.end}{self.strand}"

    def __len__(self):
        return self.end - self.start

    def format(self, format_spec, alt_header: bool = False):
        if format_spec == 'tsv':
            return (f"{self.ref}\t{self.assembly}\t{self.contig}\t{self.start}\t{self.end}\t{self.strand}\t"
                    f"{self.perc_id:.2f}\t{self.perc_cov:.2f}\t{self.cigar}\t{self.problems}\t{len(self.dna_seq)}\t"
                    f"{len(self.ref.dna_seq)}\t{len(self.aa_seq)}\t{len(self.ref.aa_seq)}\t{self.dna_hash}\t"
                    f"{self.aa_hash}\n")
        elif format_spec == 'ffn':
            return f">{self}\n{self.dna_seq}\n" if alt_header else f">{self.dna_hash}\n{self.dna_seq}\n"
        elif format_spec == 'faa':
            return f">{self}\n{self.aa_seq}\n" if alt_header else f">{self.aa_hash}\n{self.aa_seq}\n"


class References:
    def __init__(self, handle: BinaryIO, table: int, verbose: bool = False):
        self.name = handle.name
        self.T = table
        self.table = load_ttable(table)
        self.references = {n: Reference(n, d, s, table=self.table) for n, d, s in parse_fasta(handle.read(), verbose)}
        if not self.references:
            quit_with_error(f"Could not load references from {self.name}")
        if len(mol := {i.mol for i in self.references.values()}) > 1:
            quit_with_error("All references must be either DNA or Amino Acid")
        else:
            self.mol = mol.pop()

    def __str__(self):
        return self.name

    def __iter__(self):
        return iter(self.references.values())

    def __len__(self):
        return len(self.references)

    def __getitem__(self, item):
        return self.references[item]

    def format(self, format_spec):
        return ''.join(i.format(format_spec) for i in self)


class Reference:
    def __init__(self, id_: str | None = None, desc: str | None = None, seq: str | None = None, **kwargs):
        self.id = id_ or ''
        self.desc = desc or ''
        self.mol = 'DNA' if set(seq.upper()).issubset(_IUPAC_UNAMBIGUOUS_DNA) else 'AA'
        if self.mol == 'DNA':
            self.dna_seq = seq
            if len(self.dna_seq) % 3 != 0:
                warning(f"DNA sequence for {self.id} is not a multiple of 3")
            self.aa_seq = translate(self.dna_seq, **kwargs)
        else:
            self.aa_seq = seq
            self.dna_seq = ''
        if len(self.aa_seq) < 2:
            warning(f"Protein sequence for {self.id} is less than 2 residues")

    def __str__(self):
        return self.id

    def format(self, format_spec):
        if format_spec == 'ffn' and self.dna_seq:
            return f">{self.id}\n{self.dna_seq}\n"
        if format_spec == 'faa' and self.aa_seq:
            return f">{self.id}\n{self.aa_seq}\n"


# Functions -----------------------------------------------------------------------------------------------------------
def load_assembly(file: os.PathLike, verbose: bool = False) -> Assembly | None:
    """Parse an assembly file and return an Assembly object"""
    if file := check_file(file):
        if _ASSEMBLY_FASTA_REGEX.search(str(file)):
            assembly = Assembly(file, contigs={i: Contig(i, d, s) for i, d, s in parse_fasta(file, verbose)})
            log(f"Parsed {assembly} as FASTA", verbose=verbose)
            return assembly
        return warning(f"File extension must match {_ASSEMBLY_FASTA_REGEX.pattern}: {file.name}")


def write_headers(tsv: TextIO | None = None, no_header: bool = False):
    """Write the headers to a file handle."""
    if tsv and not no_header and (tsv.name == '<stdout>' or not Path(tsv.name).stat().st_size):
        tsv.write(f"Reference\tAssembly\tContig\tStart\tEnd\tStrand\tIdentity\tCoverage\tCigar\tProblems\t"
                  f"DNA_length\tRef_DNA_length\tAA_length\tRef_AA_length\tDNA_hash\tAA_hash\n")

