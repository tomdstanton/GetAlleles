# GetAlleles :computer: :dna: :microbe:
**Extract alleles for a target gene from assemblies**

*by Dr. Tom Stanton (he/him)* :scientist:

### Implementation
1. This script uses minimap2 (mappy API) to align the target gene to the assemblies in 
nucleotide space. The best alignment is selected based on matching bases (`mlen`) and
is extracted from the assembly if the identity (`mlen / blen * 100`) and coverage 
(`qlen / blen * 100`).
2. The alignments are sorted by nucleotide sequence, such that only one copy of the non-redundant
sequence is stored in memory (`Allele` objects).
3. The alleles are translated and sorted by protein sequence such that only one copy of the non-redundant 
sequence is stored in memory (`Ortholog` objects).
4. The alleles/orthologs undergo an optional MSA using the `biotite.align` package.

- The current implementation has taken scalability into account, such that the script can be run
on a large number of assemblies with minimal memory usage.
- As lookups need to be performed on the alleles and orthologs to only keep non-redundant sequences in memory,
multithreading is not implemented, however multiple threads (`os.cpu_count()`) is used to index
each assembly with `mappy`.

### Setup
Currently, no installation is required as it is a standalone script. 
However, you will need to have the following libraries installed:
 - [mappy](https://github.com/lh3/minimap2/tree/master/python)
 - [biotite](https://www.biotite-python.org/)

Which can be installed with pip (recommended):
```bash
pip install mappy biotite
```
or conda:
```bash
conda install bioconda::mappy conda-forge::biotite
```

### Usage
```bash
./get_alleles.py -h
usage: GetAlleles <reference> <assembly> [<assembly> ...] [options]

Extract alleles for a target gene from assemblies

Input:

  reference             Reference gene in ffn(.gz) format
  assemblies            Assemblies in fna(.gz) format

Allele options:

  -i 90, --min-id 90    Minimum identity percentage for alignment
  -c 90, --min-cov 90   Minimum coverage percentage for alignment
  -t 11, --table 11     Codon table to use for translation
  -p, --allow_partial   Allow partial alleles

Output options:

  -o , --output         Write/append tsv report to file (default: stdout)
  --all [alleles]       Outputs ffn, faa, aln, nwk and mat with optional prefix
  --ffn [alleles.ffn]   Output allele nucleotide sequences in FASTA format
  --faa [alleles.faa]   Output allele protein sequences in FASTA format
  --no-header           Suppress header line

MSA options:

  Any of aln, nwk or mat will trigger MSA

  --aln [alleles.aln]   Write alignment to fasta file
  --nwk [alleles.nwk]   Write MSA guide tree to newick file
  --mat [alleles.dist]  Write MSA distance matrix to TSV file
  --gap -10             Gap penalty for MSA
  --terminal True       Apply terminal penalty for MSA
  --protein             Use protein sequence for MSA

Other options:

  -h, --help            Print help and exit
  -v, --verbose         Verbose messages
  --version             Print version and exit

Version 0.0.1b0
```

### TSV
By default, the program will output a tab-separated report to `<stdout>` with the following columns:
1. **Reference**: Name of the reference sequence
1. **Assembly**: Assembly name
1. **Contig**: Contig name
1. **Start**: Start position of the allele in the contig (0-based)
1. **End**: End position of the allele in the contig
1. **Strand**: Strand of the allele in the contig
1. **Partial**: True if the allele collides with the contig boundaries
1. **Identity**: Percent identity of the allele
1. **Coverage**: Percent coverage of the allele
1. **Insertion**: If the allele has insertions (`D` in cigar string)
1. **Deletion**: If the allele has deletions (`I` in cigar string)
1. **Truncation**: True if ortholog does not form a complete CDS
1. **Allele**: Allele ID (MD5 hash of the non-redundant nucleotide sequence)
1. **Ortholog**: Ortholog ID (MD5 hash of the non-redundant protein sequence)

* The `-o/--out` option will write/append the report to a file.
* The `--no-header` option will suppress the header line. If appending to a file, the header will only be written once.
* You can suppress the output by redirecting to `/dev/null` or `NUL` on Windows.

### Distance matrix
- The MSA distance matrix is a tab-separated `N+1,N` matrix with a header row of allele/ortholog IDs;
the rows from top-to-bottom are in the order of the headers from left-to-right.

### Examples

Extract an ompA gene from assemblies, perform an MSA and save all files:
```bash
./get_alleles.py ompA.ffn assemblies/*.fna --all ompA > ompA_alleles.tsv
```
OR
```bash
./get_alleles.py ompA.ffn assemblies/*.fna --all ompA -o ompA_alleles.tsv
```
To perform an MSA on the protein sequences:
```bash
./get_alleles.py ompA.ffn assemblies/*.fna --all ompA --protein
```

