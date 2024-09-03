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
import sys
import argparse
import hashlib
from operator import attrgetter

from getalleles.version import __version__
from getalleles.log import log, bold
from getalleles.utils import check_python_version, check_out, check_cpus, load_ttable, check_programs
from getalleles.alignment import group_alns, cull_all
from getalleles.assembly import References, Allele, load_assembly, write_headers


# Functions ------------------------------------------------------------------------------------------------------------
def parse_args(a: list[str]):
    parser = argparse.ArgumentParser(
        usage="%(prog)s <reference> <assembly> [<assembly> ...] [options]",
        description='Extract alleles from genome assemblies', prog='getalleles', epilog=f"%(prog)s v{__version__}",
        formatter_class=argparse.RawTextHelpFormatter, add_help=False
    )
    input_parser = parser.add_argument_group(bold("Input"), "\nInput files / stdin can be compressed")
    input_parser.add_argument("reference", help="Reference genes in ffn/fna format, use - for stdin",
                              type=argparse.FileType('rb'), metavar="reference")
    input_parser.add_argument("assembly", nargs='+', help="Assembly file(s) in fna format",
                              metavar="assembly")

    alignment_parser = parser.add_argument_group(bold("Alignment options"), "")
    alignment_parser.add_argument("-i", "--min-id", type=float, default=80, metavar='80',
                                  help="Minimum identity percentage for alignment")
    alignment_parser.add_argument("-c", "--min-cov", type=float, default=80, metavar='80',
                                  help="Minimum coverage percentage for alignment")
    alignment_parser.add_argument("--best-n", type=int, default=0, metavar='0',
                                  help="Best N hits per reference or 0 to report all (that pass filters)")
    alignment_parser.add_argument("--cull", action='store_true',
                                  help="Culls overlapping references so only the best is kept (that pass filters)")
    alignment_parser.add_argument("--args", default='', metavar='',
                                  help="Extra arguments to pass to mini{map2,prot}; MUST BE WRAPPED")

    allele_parser = parser.add_argument_group(bold("Allele options"), "")
    allele_parser.add_argument("--table", type=int, default=11,
                               help="Codon table to use for translation", metavar='11')
    allele_parser.add_argument("--dna-hash", choices=hashlib.algorithms_available,
                               type=lambda i: getattr(hashlib, i), default=getattr(hashlib, 'sha1'),
                               help="Algorithm for hashing the allele DNA sequence ", metavar='sha1')
    allele_parser.add_argument("--aa-hash", choices=hashlib.algorithms_available,
                               type=lambda i: getattr(hashlib, i), default=getattr(hashlib, 'md5'),
                               help="Algorithm for hashing the allele AA sequence", metavar='md5')

    output_parser = parser.add_argument_group(bold("Output options"), "")
    output_parser.add_argument("-o", "--tsv", help="Write/append tsv report to file (default: stdout)",
                               type=argparse.FileType('at'), default=sys.stdout, metavar='')
    output_parser.add_argument("--ffn",
                               default=None, const="alleles.ffn", type=check_out, nargs='?', metavar="alleles.ffn",
                               help="Output allele DNA sequences (single file or directory)")
    output_parser.add_argument("--faa",
                               default=None, const="alleles.faa", type=check_out, nargs='?', metavar="alleles.faa",
                               help="Output allele AA sequences (single file or directory)")
    output_parser.add_argument("--alt-header", action='store_true', help='Sample-specific fasta headers')
    output_parser.add_argument("--no-header", action='store_true', help='Suppress header in TSV')

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


def main():
    check_python_version(3, 9)
    args = parse_args(sys.argv[1:])

    refs = References(args.reference, table=args.table, verbose=args.verbose)
    check_programs(['minimap2' if refs.mol == 'DNA' else 'miniprot'], verbose=args.verbose)

    write_headers(args.tsv, args.no_header)
    for assembly in args.assembly:
        if assembly := load_assembly(assembly, verbose=args.verbose):
            align, alignments = assembly.map(refs, args.threads, args.args, args.verbose), []
            # Filter alignments and group by query (the default for the `group_alns` function
            for _, alns in group_alns(filter(lambda i: i.perc_id >= args.min_id and i.perc_cov >= args.min_cov, align)):
                alignments.extend(
                    alns if not args.best_n else sorted(alns, key=attrgetter('mlen'), reverse=True)[0:args.best_n])
            for ctg, alns in group_alns(alignments, key='ctg'):  # Group alignments by contig for easy culling if needed
                ctg = assembly.contigs[ctg]  # Extract contig here, could do this in allele init but only one lookup
                assembly.alleles.extend(
                    Allele(refs[a.q], assembly, ctg, a.r_st, a.r_en, a.strand, a.perc_id, a.perc_cov, a.tags['cg'],
                           args.dna_hash, args.aa_hash, table=refs.table) for a in
                    (cull_all(alns) if args.cull else alns)
                )
            assembly.write(args.tsv, args.ffn, args.faa,
                           alt_header=args.alt_header)  # Write the results to the requested files

    # Cleanup ----------------------------------------------------------------------------------------------------------
    for attr in vars(args):  # Close all open files in the args namespace if they aren't sys.stdout or sys.stdin
        if (x := getattr(args, attr, None)) and isinstance(x, TextIOWrapper) and x not in {sys.stdout, sys.stdin}:
            x.close()  # Close the file

    log("Done!", verbose=args.verbose)


if __name__ == '__main__':
    main()
