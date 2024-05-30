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
from pathlib import Path
from datetime import datetime
import sys
from inspect import currentframe
import argparse
from os import cpu_count
from typing import Generator, TextIO
from hashlib import md5

# Constants -----------------------------------------------------------------------------------------------------------
__VERSION__ = '0.0.1b0'
__PROG__ = 'GetAlleles'
__DESC__ = 'Extract alleles for a target gene from assemblies'


# Logging functions ---------------------------------------------------------------------------------------------------
def log(message: str = '', verbose: bool = True, rjust: int = 15):
    if verbose:
        prefix = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {currentframe().f_back.f_code.co_name:>{rjust}}]"
        sys.stderr.write(f"{prefix} {message}\n")
        # sys.stderr.flush()


def warning(message: str):
    # The currentframe() in log() will tell user the message is a warning, so no need to prefix message
    for line in message.splitlines():
        log(f"\033[1;33m{line}\033[0m")  # Bold yellow


def quit_with_error(message: str):
    for line in message.splitlines():
        log(f"\033[1;31m{line}\033[0m")  # Bold red
    sys.exit(1)


def bold(text: str):
    return f"\033[1m{text}\033[0m"


# Dependencies --------------------------------------------------------------------------------------------------------
try:
    import mappy as mp
except ImportError:
    quit_with_error("mappy is required for GetAlleles. Please install it with 'pip install mappy'")
try:
    from biotite.sequence import NucleotideSequence, ProteinSequence, align, codon
    from biotite.sequence.io import fasta
except ImportError:
    quit_with_error("Biotite is required for GetAlleles. Please install it with 'pip install biotite'")
try:
    from numpy import savetxt  # numpy should come with biotite
except ImportError:
    quit_with_error("NumPy is required for GetAlleles. Please install it with 'pip install numpy'")


# Classes -------------------------------------------------------------------------------------------------------------
class OrfError(Exception):
    pass


class Orf:
    def __init__(self, ref_id: str, source: str, seq_id: str, start: int, end: int, strand: int,
                 partial: bool = False, perc_id: float = 100, perc_cov: float = 100,
                 insertion: bool = False, deletion: bool = False):
        self.ref_id = ref_id
        self.source = source
        self.seq_id = seq_id
        self.start = start
        self.end = end
        self.strand = strand
        self.partial = partial  # Whether the allele is partial
        self.perc_id = perc_id  # Percent identity of the allele
        self.perc_cov = perc_cov  # Percent coverage of the allele
        self.insertion = insertion
        self.deletion = deletion

    def __str__(self):
        return self.seq_id

    def __repr__(self):
        return self.seq_id

    def __len__(self):
        return self.end - self.start

    def __format__(self, format_spec):
        if format_spec == 'tsv':
            return (f"{self.ref_id}\t{self.source}\t{self.seq_id}\t{self.start}\t{self.end}\t{self.strand}\t"
                    f"{self.partial}\t{self.perc_id:.2f}\t{self.perc_cov:.2f}\t{self.insertion}\t{self.deletion}")
        return str(self)


class Reference(Orf):
    """Subclass of Orf that is allowed to store a seq attribute"""

    def __init__(self, seq: NucleotideSequence, *args, **kwargs):
        self.seq = seq
        super().__init__(*args, **kwargs)


class Assembly:
    """Convenience wrapper for the index (`mp.Aligner`), `Path` and name (ID) of the assembly"""

    def __init__(self, path: Path | None = None, name: str | None = None, max_n: int = 1):
        self.path = path or Path()
        self.name = name or path.name.strip('.gz').rsplit('.', 1)[0] if path else ''
        self.aligner = mp.Aligner(str(self.path), best_n=max_n, n_threads=cpu_count())  # Set up aligner and index assembly

    def __str__(self):
        return self.name

    def seq(self, ctg: str, start: int, end: int, strand: str = 1) -> str:
        return self.aligner.seq(ctg, start, end) if strand == 1 else mp.revcomp(self.aligner.seq(ctg, start, end))

    def map(self, *args, **kwargs):
        return self.aligner.map(*args, **kwargs)


class AlleleError(Exception):
    pass


class Allele:
    """A group of ORFs with the same nucleotide sequence"""

    def __init__(self, orfs: list[Orf], representative: Orf, seq: NucleotideSequence):
        self.orfs = orfs  # List of ORFs for the allele
        self.representative = representative
        self.seq = seq
        self.truncated = False

    def __hash__(self):
        return md5(self.seq.code)

    def __str__(self):
        return self.__hash__().hexdigest()

    def __iter__(self):
        return iter(self.orfs)

    def __len__(self):
        return len(self.orfs)

    def translate(self, table: codon.CodonTable) -> ProteinSequence:
        assert isinstance(self.seq, NucleotideSequence), TypeError("Allele must have a NucleotideSequence")
        try:
            return self.seq.translate(complete=True, codon_table=table)
        except ValueError:
            warning(f"{self} is not a complete CDS")
            self.truncated = True  # Return the longest ORF if the sequence is not a complete CDS
            return max(self.seq.translate(complete=False, codon_table=table)[0], key=len)


class Ortholog:
    """A group of Alleles with the same protein sequence"""

    def __init__(self, alleles: list[Allele], representative: Allele, seq: ProteinSequence):
        self.alleles = alleles  # List of ORFs for the allele
        self.representative = representative
        self.seq = seq

    def __hash__(self):
        return md5(self.seq.code)

    def __str__(self):
        return self.__hash__().hexdigest()

    def __iter__(self):
        return iter(self.alleles)

    def __len__(self):
        return len(self.alleles)

    def __format__(self, format_spec):
        if format_spec == 'tsv':
            return ''.join([f"{format(o, format_spec)}\t{a.truncated}\t{a}\t{self}\n" for a in self for o in a])
        elif format_spec == 'ffn':
            return ''.join([f">{allele}\n{allele.seq}\n" for allele in self])
        elif format_spec == 'faa':
            return f">{self}\n{self.seq}\n"
        return str(self)


# Functions -----------------------------------------------------------------------------------------------------------
def parse_args(a: list[str]):
    parser = argparse.ArgumentParser(
        usage="%(prog)s <reference> <assembly> [<assembly> ...] [options]", epilog=f"Version {__VERSION__}",
        description=__DESC__, prog=__PROG__, formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False
    )
    input_parser = parser.add_argument_group(bold("Input"), "")
    input_parser.add_argument("reference", help="Reference gene in ffn(.gz) format",
                              type=check_file, metavar="reference")
    input_parser.add_argument("assemblies", nargs='+', help="Assemblies in fna(.gz) format",
                              type=check_file, metavar="assemblies")

    allele_parser = parser.add_argument_group(bold("Allele options"), "")
    allele_parser.add_argument("-i", "--min-id", type=float, default=90, metavar='90',
                               help="Minimum identity percentage for alignment")
    allele_parser.add_argument("-c", "--min-cov", type=float, default=90, metavar='90',
                               help="Minimum coverage percentage for alignment")
    allele_parser.add_argument("-t", "--table", type=int, default=11, help="Codon table to use for translation",
                               metavar='11')
    allele_parser.add_argument("-p", "--allow_partial", action='store_true', help="Allow partial alleles")
    allele_parser.add_argument("-n", "--max_n", type=int, default=1, help="Maximum number of alleles per assembly",
                               metavar='1')

    output_parser = parser.add_argument_group(bold("Output options"), "")
    output_parser.add_argument("-o", "--output", help="Write/append tsv report to file (default: stdout)",
                               type=argparse.FileType('at'), default=sys.stdout, metavar='')
    output_parser.add_argument("--all", metavar='alleles', const='alleles', default=None, nargs="?",
                               help='Outputs ffn, faa, aln, nwk and mat with optional prefix')
    output_parser.add_argument("--ffn", help="Output allele nucleotide sequences in FASTA format",
                               default=None, const="alleles.ffn", type=argparse.FileType('wt'), nargs='?',
                               metavar="alleles.ffn")
    output_parser.add_argument("--faa", help="Output allele protein sequences in FASTA format",
                               default=None, const="alleles.faa", type=argparse.FileType('wt'), nargs='?',
                               metavar="alleles.faa")
    output_parser.add_argument("--no-header", action='store_true', help='Suppress header line')

    msa_parser = parser.add_argument_group(bold("MSA options"), "\nAny of aln, nwk or mat will trigger MSA")
    msa_parser.add_argument("--aln", help="Write alignment to fasta file",
                            default=None, const="alleles.aln", type=argparse.FileType('wt'), nargs='?',
                            metavar="alleles.aln")
    msa_parser.add_argument("--nwk", help="Write MSA guide tree to newick file", default=None,
                            const="alleles.nwk", type=argparse.FileType('wt'), nargs='?', metavar="alleles.nwk")
    msa_parser.add_argument("--mat", help="Write MSA distance matrix to TSV file", default=None,
                            const="alleles.dist", type=argparse.FileType('wt'), nargs='?', metavar="alleles.dist")
    msa_parser.add_argument('--gap', help="Gap penalty for MSA", metavar="-10",
                            type=int, default=-10)
    msa_parser.add_argument('--terminal', help=" Apply terminal penalty for MSA", metavar="True",
                            type=bool, default=True)
    msa_parser.add_argument("--protein", help="Use protein sequence for MSA", action='store_true')

    other_parser = parser.add_argument_group(bold("Other options"), "")
    other_parser.add_argument("-h", "--help", action='help', help="Print help and exit")
    other_parser.add_argument("-v", "--verbose", action='store_true', help="Verbose messages")
    other_parser.add_argument("--version", action='version', version=__VERSION__,
                              help="Print version and exit")
    return parser.parse_args(a)


def check_file(file: str) -> Path:
    if not (file := Path(file)).exists() or not file.is_file() or not file.stat().st_size:
        quit_with_error(f"'{file.name}' not found or empty")  # Throws a warning so pipeline can continue
    return file


def write_headers(tsv: TextIO | None = None, no_header: bool = False):
    """Write the headers to a file handle."""
    if tsv:
        if tsv.name != '<stdout>' and tsv.tell() != 0:  # If file is path and not already written to
            no_header = True  # Headers already written, useful for running on HPC
        if not no_header:
            tsv.write(
                f"Reference\tAssembly\tContig\tStart\tEnd\tStrand\tPartial\tIdentity\tCoverage\tInsertion\tDeletion\t"
                f"Truncation\tAllele\tOrtholog\n"
            )


def adjust_codon_range(query_start: int, query_end: int, ref_start: int, ref_end: int, strand: int, clip=False,
                       verbose: bool = False) -> tuple[int, int, int, int]:
    """Adjust the alignment ranges to the nearest codon boundaries on the query"""
    if not clip:
        start = next(
            i for i in (0, 1, 2) if (query_start - i) % 3 == 0)  # Subtract frame to find the start of the codon
        end = next(i for i in (0, 1, 2) if (query_end + i) % 3 == 0)  # Add frame to find the end of the codon
        new_query_start = query_start - start  # Extend the start of the codon
        new_query_end = query_end + end  # Extend the end of the codon
        new_ref_start = ref_start - start if strand == 1 else ref_start - end
        new_ref_end = ref_end + end if strand == 1 else ref_end + start
    else:
        start = next(i for i in (0, 1, 2) if (query_start + i) % 3 == 0)  # Add frame to find the start of the codon
        end = next(i for i in (0, 1, 2) if (query_end - i) % 3 == 0)  # Subtract frame to find the end of the codon
        new_query_start = query_start + start  # Clip the start of the codon
        new_query_end = query_end - end  # Clip the end of the codon
        new_ref_start = ref_start + start if strand == 1 else ref_start + end
        new_ref_end = ref_end - end if strand == 1 else ref_end - start

    if query_start != new_query_start or query_end != new_query_end:
        log(f"{'Clipping' if clip else 'Extending'} query: {query_start}-{query_end} -> "
            f"{new_query_start}-{new_query_end}; Ref: {ref_start}-{ref_end} -> {new_ref_start}-{new_ref_end}",
            verbose=verbose)

    return new_query_start, new_query_end, new_ref_start, new_ref_end


def check_partial(query_length: int, query_start: int, query_end: int, ref_length: int, ref_start: int,
                  ref_end: int) -> bool:
    """Check if the query is partial or would be partial if alignment was full length"""
    # First two cases just check for overlap
    if ref_length < query_length:  # E.g. Contig length 99 and gene length 100 means 1bp is missing
        return True
    if ref_start < query_start:  # E.g. Contig position 0 and gene position 1 means 1bp is missing
        return True
    if ref_end < query_length:  # E.g. Contig position 99 and gene length 100 means 1bp is missing
        return True
    if (ref_length - ref_end) < (query_length - query_end):  # Where the alignment is partial at the end
        return True
    return False


def extract_orfs(assembly: Path, ref: Reference, min_id: float, min_cov: float, allow_partial: bool = True, max_n: int = 1,
                verbose: bool = False) -> Generator[tuple[Orf, str], None, None]:
    try:
        assembly = Assembly(assembly, max_n=max_n)  # This will trigger the build of the index
    except Exception as e:
        warning(f"Failed to parse {assembly}: {e}")  # Throws a warning so pipeline can continue
        yield None
    log(f"Aligning {ref} to {assembly}", verbose=verbose)

    n = 0
    for a in assembly.map(str(ref.seq)):
        if (perc_id := a.mlen / a.blen * 100) >= min_id and (perc_cov := len(ref) / a.blen * 100) >= min_cov:
            new_query_start, new_query_end, new_ref_start, new_ref_end = adjust_codon_range(
                a.q_st, a.q_en, a.r_st, a.r_en, a.strand, verbose=verbose)
        
            if partial := (check_partial(len(ref), new_query_start, new_query_end, a.ctg_len, new_ref_start, new_ref_end) and
                           not allow_partial):  # Skip partial alleles if not allowed
                warning(f"Partial ORF found in {assembly}")
                yield None
            n += 1
            yield (Orf(ref.seq_id, assembly.name, a.ctg, new_ref_start, new_ref_end, a.strand, partial, perc_id, perc_cov,
                        # TODO: Check if these are correct: https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
                        deletion="I" in a.cigar_str,  # Insertion in query == deletion in reference
                        insertion="D" in a.cigar_str  # Deletion in query == insertion in reference
                        ),
                    assembly.seq(a.ctg, new_ref_start, new_ref_end, a.strand))
    if n == 0:
        warning(f"No alignment found for {ref} in {assembly}")


def get_alleles(assemblies: list[Path], ref: Reference, out: TextIO, ffn: TextIO | None, faa: TextIO | None,
                table: codon.CodonTable, no_header: bool, verbose: bool, **kwargs) -> Generator[Ortholog, None, None]:
    """
    Extract alleles from assemblies for a target gene.
    This function purposefully does not use `itertools.groupby` for grouping objects by sequence as only one copy of
    each non-redundant sequence is kept in memory so this program can scale to large datasets.
    """
    alleles, orthologs = {}, {}  # Hash alleles by nucleotide sequence and orthologs by protein sequence
    for assembly in assemblies:  # Extract ORFs from each assembly
        for orf in extract_orfs(assembly, ref, verbose=verbose, **kwargs):
            if orf:
                if orf[1] not in alleles:  # orf[1] is the nucleotide sequence
                    alleles[orf[1]] = [orf[0]]  # orf[0] is the Orf object, which doesn't store the sequence
                else:  # If the nucleotide sequence is already in the dictionary, append the ORF
                    alleles[orf[1]].append(orf[0])
    if not alleles:
        quit_with_error("No alignments found")

    # Add reference to alleles
    if (seq := str(ref.seq)) in alleles:
        alleles[seq].append(ref)
    else:
        alleles[seq] = [ref]

    log(f"Found {len(alleles)} unique alleles in {len(assemblies)} assemblies", verbose=verbose)
    # If we have alleles, we can now group them by protein sequence and write tabular output
    write_headers(out, no_header)

    for seq, orfs in alleles.items():  # type: str, list[Orf]
        allele = Allele(orfs, orfs[0], NucleotideSequence(seq))
        if (protein_seq := str(allele.translate(table))) not in orthologs:
            orthologs[protein_seq] = [allele]
        else:
            orthologs[protein_seq].append(allele)

    log(f"Found {len(orthologs)} unique orthologs", verbose=verbose)
    for seq, alleles in orthologs.items():  # type: str, list[Allele]
        out.write(format(ortholog := Ortholog(alleles, alleles[0], ProteinSequence(seq)), 'tsv'))
        if faa:
            faa.write(format(ortholog, 'faa'))
        if ffn:
            ffn.write(format(ortholog, 'ffn'))
        yield ortholog


def main():
    args = parse_args(sys.argv[1:])
    try:
        name, seq, _ = next(mp.fastx_read(str(args.reference)))
        args.reference = Reference(NucleotideSequence(seq), "reference",
                                   args.reference.name.strip('.gz').rsplit('.', 1)[0], name, 0, len(seq), 1)
        args.reference.seq.translate(args.table)
    except Exception as e:
        quit_with_error(f"Failed to load reference: {e}")
    
    try:
        args.table = codon.CodonTable.load(args.table)
    except Exception as e:
        quit_with_error(f"Failed to load translation table: {e}")

    if args.all:
        args.ffn = open(f"{args.all}.ffn", 'wt')
        args.faa = open(f"{args.all}.faa", 'wt')
        args.aln = open(f"{args.all}.aln", 'wt')
        args.nwk = open(f"{args.all}.nwk", 'wt')
        args.mat = open(f"{args.all}.dist", 'wt')

    orthologs = list(get_alleles(
        args.assemblies, args.reference, args.output, args.ffn, args.faa, args.table, args.no_header,
        verbose=args.verbose, min_id=args.min_id, min_cov=args.min_cov, allow_partial=args.allow_partial, max_n=args.max_n))

    if any([args.aln, args.nwk, args.mat]):  # Perform MSA if any of these options are set
        groups = orthologs if args.protein else [allele for ortholog in orthologs for allele in ortholog.alleles]
        log(f"Performing MSA with {len(groups)} unique {'protein' if args.protein else 'nucleotide'} sequences",
            verbose=args.verbose)
        alignment, _, tree, distances = align.align_multiple(
            [i.seq for i in groups], matrix=align.SubstitutionMatrix.std_protein_matrix() if args.protein else
            align.SubstitutionMatrix.std_nucleotide_matrix(),
            terminal_penalty=args.terminal, gap_penalty=args.gap
        )  # Perform MSA
        if args.nwk:  # Write tree to file
            args.nwk.write(tree.to_newick(labels=list(map(str, groups)), include_distance=True))
        if args.aln:  # Write alignment to file
            args.aln.write(''.join([f">{x}\n{y}\n" for x, y in zip(groups, alignment.get_gapped_sequences())]))
        if args.mat:  # Write distance matrix to file
            savetxt(args.mat, distances, delimiter='\t', header='\t'.join(map(str, groups)), comments="")

    # Close all open files
    for i in [args.output, args.ffn, args.faa, args.aln, args.nwk, args.mat]:
        if i and i.name != '<stdout>':
            i.close()

    log("Done!", verbose=args.verbose)


if __name__ == '__main__':
    main()
