# GetAlleles :computer: :dna: :microbe:
**Extract alleles from genome assemblies**

## Introduction
GetAlleles is a simple program to extract alleles of target genes from genome assemblies.
It uses alignment via [`mappy`](https://pypi.org/project/mappy/) to align reference gene
nucleotide sequences to the assembly contigs, and extracts the alleles that pass the identity
and coverage thresholds. The extracted sequences are then translated to determine whether the
protein is truncated. The nucleotide and amino acid sequences are hashed for easy allele
identification in tabular format and for submission to typing scheme databases such as
[PubMLST](https://pubmlst.org/).

## Installation
```bash
pip install git+https://github.com/tomdstanton/GetAlleles.git
```

## Usage
```
usage: getalleles <reference> <assembly> [<assembly> ...] [options]

Extract alleles from genome assemblies

Input:

  reference            Reference genes in ffn(.gz) format
  assembly             Assembly file(s) in fna(.gz) format

Allele options:

  -i 80, --min-id 80   Minimum identity percentage for alignment
  -c 80, --min-cov 80  Minimum coverage percentage for alignment
  --table 11           Codon table to use for translation
  --dna-hashfunc sha1  Algorithm for hashing the allele DNA sequence
  --pro-hashfunc md5   Algorithm for hashing the allele protein sequence

Output options:

  -o , --tsv           Write/append tsv report to file (default: stdout)
  --ffn [alleles.ffn]  Output allele nucleotide sequences (single file or directory)
  --faa [alleles.faa]  Output allele amino acid sequences (single file or directory)
  --no-header          Suppress header line

Other options:

  -t 0, --threads 0    Number of indexing threads or 0 for all available
  -h, --help           Print help and exit
  -v, --verbose        Verbose messages
  --version            Print version and exit

getalleles v0.0.1b1
```

## Output
### Tabular
By default, the program will output a BED-style TSV to `<stdout>` with the following columns:
1. **Reference**: Name of the reference sequence
1. **Assembly**: Assembly name
1. **Contig**: Contig name
1. **Start**: Start position of the allele in the contig (0-based)
1. **End**: End position of the allele in the contig
1. **Strand**: Strand of the allele in the contig
1. **Identity**: Percent identity of the allele
1. **Coverage**: Percent coverage of the allele
1. **Problems**: "Partial" if gene runs off contig edge or "Truncated" if length of translation is less than the reference.
1. **Cigar**: Alignment CIGAR string.
1. **Copy_number**: Number of copies of the _gene_ in the assembly.
1. **DNA_hash**: Hash of the DNA sequence.
1. **Protein_hash**: Hash of the protein sequence.

* This can be written/appended to a file with either the `-o/--tsv` flag or `>` and `>>` in bash.
* The `--no-header` option will suppress the header line. If appending to a file, the header will only be written once.
* You can suppress the output by redirecting to `/dev/null` or `NUL` on Windows.

### Fasta
You can output the extracted DNA and protein sequences in Fasta format with the `--ffn` 
and `--faa` flags respectively, with the defaults being `alleles.{ffn,faa}`.

These arguments take a file or directory. If the file already exists, it will be appended to;
if the argument is a directory, _one file per assembly_ will be generated.

The ID (header) of each sequence is the hash digest of the sequence.

## Examples

Extract reference genes from assemblies and save results to a file:
```bash
getalleles genes.ffn assemblies/*.fna > alleles.tsv
```
OR
```bash
getalleles genes.ffn assemblies/*.fna -o alleles.tsv
```
To do the same but output the nucleotide sequences to a fasta file.
```bash
getalleles genes.ffn assemblies/*.fna -o alleles.tsv --ffn alleles.ffn
```
OR
```bash
getalleles genes.ffn assemblies/*.fna -o alleles.tsv --ffn - > alleles.ffn
```
To output the allele protein sequences into an MSA program:
```bash
getalleles genes.ffn assemblies/*.fna -o alleles.tsv --ffn --faa - > muscle
```
To output _one sequence file per assembly_ :
```bash
getalleles genes.ffn assemblies/*.fna -o alleles.tsv --ffn dna_seqs/ --faa protein_seqs/
```

## Acknowledgements
- [Dr. Ben Vezina](https://github.com/bananabenana) for ideas and initial implementation
of multiple copies in `v0.0.1b0`.
