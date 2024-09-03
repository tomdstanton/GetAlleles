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
import sys
from pathlib import Path
from typing import TextIO, Generator, BinaryIO
from zlib import decompress as gz_decompress
from gzip import open as gz_open
from bz2 import (decompress as bz2_decompress, open as bz2_open)
from lzma import (decompress as xz_decompress, open as xz_open)

from getalleles.log import log, quit_with_error, warning

# Constants ------------------------------------------------------------------------------------------------------------
_MAGIC_BYTES = {b'\x1f\x8b': 'gz', b'\x42\x5a': 'bz2', b'\xfd7zXZ\x00': 'xz'}
_OPEN = {'gz': gz_open, 'bz2': bz2_open, 'xz': xz_open}
_DECOMPRESS = {'gz': gz_decompress, 'bz2': bz2_decompress, 'xz': xz_decompress}
_MIN_N_BYTES = max(len(i) for i in _MAGIC_BYTES)  # Minimum number of bytes to read in a file to guess the compression
_COMPLEMENT = str.maketrans("ACGT", "TGCA")  # type: dict[int, int]
_TRANS_TABLES = {  # Condensed version of Biopython translation tables from stantlib.fasta
    # Can't remember where I got this from originally, make sure this gives correct translations
    1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    2: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
    3: 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    5: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
    6: 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    9: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    10: 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    11: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    12: 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    13: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
    14: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    15: 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    16: 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    21: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
    22: 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    23: 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    24: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
    25: 'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    26: 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    27: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    28: 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    29: 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    30: 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    31: 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    32: 'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    33: 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG'
}


# Functions -----------------------------------------------------------------------------------------------------------
def load_ttable(table: int) -> dict[str, str]:
    if not (t := _TRANS_TABLES.get(table)):
        quit_with_error(f'No such translation table: {table}')
    # Taken from stantlib.fasta
    return dict(zip([a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'], t))


def translate(seq: str, table: int | dict[str, str], to_stop: bool = True, stop_symbol: str = "*",
              frame: int = 0) -> str:
    if frame not in {0, 1, 2}:
        quit_with_error(f"Frame must be 0, 1 or 2, not {frame}")
    if isinstance(table, int):  # Load ttable if necessary
        table = load_ttable(table)
    protein = []  # Load ttable if necessary
    for i in range(frame, len(seq), 3):  # Iterate over seq codons (chunks of 3)
        if len(codon := seq[i:i + 3]) < 3:  # We can't translate chunks less than 3 (not codons) so break here
            break
        protein.append(aa := table.get(codon, stop_symbol))  # If lookup fails (get == None), assume stop codon
        if to_stop and aa == stop_symbol:  # Break if to_stop == True
            break
    return ''.join(protein)


def check_file(path: os.PathLike) -> Path | None:
    if not (path := Path(path)).exists():
        return warning(f'{path} does not exist')
    if not path.is_file():
        return warning(f'{path} is not a file')
    if path.stat().st_size == 0:
        return warning(f'{path} is empty')
    return path.absolute()


def check_cpus(cpus: int | None = 0, verbose: bool = False) -> int:
    cpus = os.cpu_count() if not cpus else min(cpus, os.cpu_count())
    log(f'Using {cpus} CPUs', verbose)
    return cpus


def check_out(path: str, mode: str = "at", parents: bool = True, exist_ok: bool = True) -> Path | TextIO:
    """
    Check if the user wants to create/append a file or directory.
    If it looks like/is already a file (has an extension), return the file object.
    If it looks like/is already a directory, return the directory path.
    """
    if path == '-':  # If the path is '-', return stdout
        return sys.stdout
    if (path := Path(path)).suffix:  # If the path has an extension, it's probably a file
        try:
            return path.open(mode)  # Open the file
        except Exception as e:
            quit_with_error(f'Could not open {path}: {e}')
    if not path.exists():  # Assume directory
        try:
            path.mkdir(parents=parents, exist_ok=exist_ok)  # Create the directory if it doesn't exist
        except Exception as e:
            quit_with_error(f'Could not create {path}: {e}')
    return path


def check_python_version(major: int = 3, minor: int = 9):
    if sys.version_info.major < major or sys.version_info.minor < minor:
        quit_with_error(f'Python version {major}.{minor} or greater required')


def range_overlap(range1: tuple[int, int], range2: tuple[int, int], skip_sort: bool = False) -> int:
    """
    Returns the overlap between two ranges
    :param range1: Tuple of start and end positions
    :param range2: Tuple of start and end positions
    :param skip_sort: Skip sorting each range before calculating the overlap
    :return: Integer of overlap
    """
    start1, end1 = range1 if skip_sort else sorted(range1)
    start2, end2 = range2 if skip_sort else sorted(range2)
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def decompress(data: bytes, verbose: bool = False, *args, **kwargs) -> bytes:
    """
    Decompress bytes using magic bytes at the beginning of the data
    :param data: Bytes to decompress
    :param verbose: Print log messages to stderr
    :return: Decompressed bytes
    """
    first_bytes = data[:_MIN_N_BYTES]
    for magic, compression in _MAGIC_BYTES.items():
        if first_bytes.startswith(magic):
            log(f"Assuming data is compressed with {compression}", verbose=verbose)
            try:
                return _DECOMPRESS[compression](data, *args, **kwargs)
            except Exception as e:
                return warning(f"Error decompressing with {compression}; {first_bytes=}\n{e}")
    log("Assuming data is uncompressed", verbose=verbose)
    return data


def opener(file: os.PathLike, verbose: bool = False, *args, **kwargs) -> TextIO | BinaryIO:
    """
    Opens a file with the appropriate open function based on the magic bytes at the beginning of the data
    :param file: File to open
    :param check: Check file exists and is not empty
    :param verbose: Print log messages to stderr
    :return: Decompressed bytes
    """
    with open(file, 'rb') as f:  # Open the file to read bytes
        first_bytes = f.read(_MIN_N_BYTES)  # Get the bytes necessary to guess the compression type
    for magic, compression in _MAGIC_BYTES.items():
        if first_bytes.startswith(magic):
            log(f"Assuming file is compressed with {compression}", verbose=verbose)
            try:
                return _OPEN[compression](file, *args, **kwargs)
            except Exception as e:
                return warning(f"Error opening with {compression}; {first_bytes=}\n{e}")
    log("Assuming file is uncompressed", verbose=verbose)
    return open(file, *args, **kwargs)


def parse_fasta(data: str | bytes | os.PathLike, verbose: bool = False) -> Generator[tuple[str, str, str], None, None]:
    """
    Simple multi-line fasta parser
    :param data: String, bytes or a file with the data; bytes/files will undergo automatic decompression if needed
    :param verbose: Print log messages to stderr
    :return: A generator of tuples with the name, description and sequence of each record as strings
    """
    if isinstance(data, bytes):
        data = decompress(data, verbose=verbose).decode()
    elif Path(data).is_file():
        with opener(data, mode='rt', verbose=verbose) as f:
            data = f.read()
    header, seq = '', []
    for line in data.splitlines():  # Loop over lines in chunk
        if line.startswith('>'):
            if header and seq:
                name, desc = header.split(' ', 1) if ' ' in header else (header, '')
                yield name, desc, ''.join(seq)
            header, seq = line.lstrip('>'), []
        else:
            seq.append(line)
    if header and seq:
        name, desc = header.split(' ', 1) if ' ' in header else (header, '')
        yield name, desc, ''.join(seq)
    else:
        warning("No fasta records parsed")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def check_programs(progs: list[str], verbose: bool = False):
    """Check if programs are installed and executable"""
    bins = {  # Adapted from: https://unix.stackexchange.com/a/261971/375975
        f: Path(f'{p}/{f}') for p in filter(
            os.path.isdir, os.environ["PATH"].split(os.path.pathsep)
        ) for f in os.listdir(p) if os.access(f'{p}/{f}', os.X_OK)
    }
    for program in progs:
        if program in bins:
            log(f'{program}: {bins[program]}', verbose=verbose)
        else:
            quit_with_error(f'{program} not found')
