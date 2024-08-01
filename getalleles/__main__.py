#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

from io import TextIOWrapper
from pathlib import Path
import sys
import argparse
from typing import TextIO, Generator
import hashlib
from itertools import chain

from getalleles.version import __version__
from getalleles.log import log, warning, quit_with_error, bold
from getalleles.utils import check_python_version, check_out, check_file, check_cpus, translate, load_ttable

# Dependencies --------------------------------------------------------------------------------------------------------
try:
    import mappy as mp
except ImportError:
    quit_with_error("mappy is required for GetAlleles. Please install it with 'pip install mappy'")


# Functions ------------------------------------------------------------------------------------------------------------
def parse_args(a: list[str]):
    parser = argparse.ArgumentParser(
        usage="%(prog)s <reference> <assembly> [<assembly> ...] [options]",
        description='Extract alleles from genome assemblies', prog='getalleles', epilog=f"%(prog)s v{__version__}",
        formatter_class=argparse.RawTextHelpFormatter, add_help=False
    )
    input_parser = parser.add_argument_group(bold("Input"), "")
    input_parser.add_argument("reference", help="Reference genes in ffn(.gz) format",
                              type=check_file, metavar="reference")
    input_parser.add_argument("assembly", nargs='+', help="Assembly file(s) in fna(.gz) format",
                              type=check_file, metavar="assembly")

    allele_parser = parser.add_argument_group(bold("Allele options"), "")
    allele_parser.add_argument("-i", "--min-id", type=float, default=80, metavar='80',
                               help="Minimum identity percentage for alignment")
    allele_parser.add_argument("-c", "--min-cov", type=float, default=80, metavar='80',
                               help="Minimum coverage percentage for alignment")
    allele_parser.add_argument("--table", type=lambda i: load_ttable(int(i)),
                               default=load_ttable(11), help="Codon table to use for translation", metavar='11')
    allele_parser.add_argument("--dna-hashfunc", choices=hashlib.algorithms_available,
                               help="Algorithm for hashing the allele DNA sequence ", metavar='sha1', default='sha1')
    allele_parser.add_argument("--pro-hashfunc", choices=hashlib.algorithms_available,
                               help="Algorithm for hashing the allele protein sequence", metavar='md5', default='md5')

    output_parser = parser.add_argument_group(bold("Output options"), "")
    output_parser.add_argument("-o", "--tsv", help="Write/append tsv report to file (default: stdout)",
                               type=argparse.FileType('at'), default=sys.stdout, metavar='')
    output_parser.add_argument("--ffn",
                               default=None, const="alleles.ffn", type=check_out, nargs='?', metavar="alleles.ffn",
                               help="Output allele nucleotide sequences (single file or directory)")
    output_parser.add_argument("--faa",
                               default=None, const="alleles.faa", type=check_out, nargs='?', metavar="alleles.faa",
                               help="Output allele amino acid sequences (single file or directory)")
    output_parser.add_argument("--no-header", action='store_true', help='Suppress header line')

    other_parser = parser.add_argument_group(bold("Other options"), "")
    other_parser.add_argument('-t', '--threads', type=check_cpus, default=check_cpus(), metavar='0',
                              help="Number of indexing threads or 0 for all available")
    other_parser.add_argument("-h", "--help", action='help', help="Print help and exit")
    other_parser.add_argument("-v", "--verbose", action='store_true', help="Verbose messages")
    other_parser.add_argument("--version", action='version', version=__version__,
                              help="Print version and exit")
    if len(a) == 0:  # No arguments, print help message
        parser.print_help(sys.stderr)
    return parser.parse_args(a)


def write_headers(tsv: TextIO | None = None, no_header: bool = False):
    """Write the headers to a file handle."""
    if tsv and not no_header and (tsv.name == '<stdout>' or not Path(tsv.name).stat().st_size):
        tsv.write(f"Reference\tAssembly\tContig\tStart\tEnd\tStrand\tIdentity\tCoverage\tCigar\tProblems\t"
                  f"Copy_number\tProtein_length\tDNA_hash\tProtein_hash\n")


# Classes --------------------------------------------------------------------------------------------------------------
class AlleleError(Exception):
    pass


class Allele:
    def __init__(self, ref: Reference, assembly: str, contig: str, start: int, end: int, strand: int,
                 perc_id: float = 100, perc_cov: float = 100, cigar: str = None, problems: str = None,
                 dna_seq: str = None, dna_hash: str = None, protein_seq: str = None, protein_hash: str = None,
                 copy_n: int = 1):
        self.ref = ref
        self.assembly = assembly
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.perc_id = perc_id  # Percent identity of the allele
        self.perc_cov = perc_cov  # Percent coverage of the allele
        self.problems = problems or ''
        self.cigar = cigar or ''
        self.dna_seq = dna_seq
        self.dna_hash = dna_hash
        self.protein_seq = protein_seq
        self.protein_hash = protein_hash
        self.copy_n = copy_n

    def __str__(self):
        return f"{self.assembly}__{self.contig}__{self.ref}__{self.copy_n}"

    def __len__(self):
        return self.end - self.start

    def format(self, format_spec):
        if format_spec == 'tsv':
            return (f"{self.ref}\t{self.assembly}\t{self.contig}\t{self.start}\t{self.end}\t{self.strand}\t"
                    f"{self.perc_id:.2f}\t{self.perc_cov:.2f}\t{self.cigar}\t{self.problems}\t{self.copy_n}\t"
                    f"{len(self.protein_seq)}\t{self.dna_hash}\t{self.protein_hash}\n")
        elif format_spec == 'ffn':
            return f">{self.dna_hash}\n{self.dna_seq}\n"
        elif format_spec == 'faa':
            return f">{self.protein_hash}\n{self.protein_seq}\n"


class Reference:
    def __init__(self, name: str, seq: str, **kwargs):
        self.name = name
        self.dna_seq = seq
        if len(seq) % 3 != 0:
            warning(f"DNA sequence for {self} is not a multiple of 3")
        self.protein_seq = translate(seq, **kwargs)
        if len(self.protein_seq) < 2:
            warning(f"Protein sequence for {self} is less than 2 residues")

    def __str__(self):
        return self.name


def extract_alleles(aligner: mp.Aligner, reference: Reference, assembly_name: str, args: argparse.Namespace) -> \
        Generator[Allele, None, None]:
    copy_n = 0
    for alignment in aligner.map(reference.dna_seq):
        if ((perc_id := (alignment.mlen / alignment.blen) * 100) >= args.min_id and
                (perc_cov := (alignment.blen / len(reference.dna_seq)) * 100) >= args.min_cov):
            copy_n += 1
            allele = Allele(reference, assembly_name, alignment.ctg, alignment.r_st, alignment.r_en, alignment.strand,
                            perc_id, perc_cov, alignment.cigar_str, copy_n=copy_n)
            allele.dna_seq = aligner.seq(allele.contig, allele.start, allele.end)
            if allele.strand != 1:
                allele.dna_seq = mp.revcomp(allele.dna_seq)
            allele.protein_seq = translate(allele.dna_seq, table=args.table)
            if len(allele.dna_seq) < len(reference.dna_seq) and any((allele.start == 0, allele.end == alignment.ctg_len)):
                allele.problems = 'Partial'
            elif len(allele.protein_seq) < len(reference.protein_seq):
                allele.problems = 'Truncated'
            allele.dna_hash = args.dna_hashfunc(allele.dna_seq.encode(), usedforsecurity=False).hexdigest()
            allele.protein_hash = args.pro_hashfunc(allele.protein_seq.encode(), usedforsecurity=False).hexdigest()
            yield allele
    # if copy_n == 0:
    #     warning(f"No alignment for {reference} in {assembly_name}")


def process_assembly(assembly: Path, args: argparse.Namespace):
    assembly_name = assembly.name.strip('.gz').rsplit('.', 1)[0]
    try:
        aligner = mp.Aligner(str(assembly), n_threads=args.threads)  # Set up aligner and index assembly
    except Exception as e:
        return warning(f"Could not index {assembly_name}\n{e}")
    if alleles := list(
            chain(*[extract_alleles(aligner, reference, assembly_name, args) for reference in args.reference])):
        if args.tsv:
            args.tsv.write(''.join([i.format('tsv') for i in alleles]))
        for fh, fmt in ((args.ffn, 'ffn'), (args.faa, 'faa')):
            if fh:
                if isinstance(fh, Path):
                    (fh / f'{assembly_name}_alleles.{fmt}').write_text(''.join([i.format(fmt) for i in alleles]))
                else:
                    fh.write(''.join([i.format(fmt) for i in alleles]))


def main():
    check_python_version(3, 9)
    args = parse_args(sys.argv[1:])
    args.dna_hashfunc, args.pro_hashfunc = getattr(hashlib, args.dna_hashfunc), getattr(hashlib, args.pro_hashfunc)
    try:
        args.reference = [Reference(n, s, table=args.table) for n, s, _ in mp.fastx_read(str(args.reference))]
    except Exception as e:
        quit_with_error(f"Could not parse references\n{e}")
    write_headers(args.tsv, args.no_header)
    [process_assembly(assembly, args) for assembly in args.assembly]

    # Cleanup ----------------------------------------------------------------------------------------------------------
    for attr in vars(args):  # Close all open files in the args namespace if they aren't sys.stdout or sys.stdin
        if (x := getattr(args, attr, None)) and isinstance(x, TextIOWrapper) and x not in {sys.stdout, sys.stdin}:
            x.close()  # Close the file

    log("Done!", verbose=args.verbose)


if __name__ == '__main__':
    main()
