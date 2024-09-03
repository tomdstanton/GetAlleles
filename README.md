# GetAlleles :computer: :dna: :microbe:
**Extract alleles from genome assemblies**

## Introduction
GetAlleles is a simple program to extract alleles of target genes from genome assemblies.
It uses alignment via [`minimap2`](https://lh3.github.io/minimap2) or [`miniprot`](https://lh3.github.io/miniprot) to align reference genes to 
assembly contigs, and extracts the alleles that pass the identity
and coverage thresholds. The extracted sequences are then translated to determine whether the
protein is truncated. The nucleotide and amino acid sequences are hashed for easy allele
identification in tabular format and for submission to typing scheme databases such as
[PubMLST](https://pubmlst.org/).

## Installation
```bash
pip install git+https://github.com/tomdstanton/GetAlleles.git
```
If using nucleotide references, [`minimap2`](https://lh3.github.io/minimap2) needs to be installed in your `$PATH`.

If using amino acid references, [`miniprot`](https://lh3.github.io/miniprot) needs to be installed in your `$PATH`.

## Usage
```
usage: getalleles <reference> <assembly> [<assembly> ...] [options]

Extract alleles from genome assemblies

Input:

  Input files / stdin can be compressed

  reference            Reference genes in ffn/fna format, use - for stdin
  assembly             Assembly file(s) in fna format

Alignment options:

  -i 80, --min-id 80   Minimum identity percentage for alignment
  -c 80, --min-cov 80  Minimum coverage percentage for alignment
  --best-n 0           Best N hits per reference or 0 to report all (that pass filters)
  --cull               Culls overlapping references so only the best is kept (that pass filters)
  --args               Extra arguments to pass to mini{map2,prot}; MUST BE WRAPPED

Allele options:

  --table 11           Codon table to use for translation
  --dna-hash sha1      Algorithm for hashing the allele DNA sequence
  --aa-hash md5        Algorithm for hashing the allele AA sequence

Output options:

  -o , --tsv           Write/append tsv report to file (default: stdout)
  --ffn [alleles.ffn]  Output allele DNA sequences (single file or directory)
  --faa [alleles.faa]  Output allele AA sequences (single file or directory)
  --alt-header         Sample-specific fasta headers
  --no-header          Suppress header in TSV

Other options:

  -t 0, --threads 0    Number of indexing threads or 0 for all available
  -h, --help           Print help and exit
  -v, --verbose        Verbose messages
  --version            Print version and exit

getalleles v0.0.2b0
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
1. **DNA_length**: Length of the nucleotide sequence.
1. **Ref_DNA_length**: Length of the reference nucleotide sequence.
1. **AA_length**: Length of the protein sequence.
1. **Ref_AA_length**: Length of the reference protein sequence.
1. **DNA_hash**: Hash of the DNA sequence.
1. **AA_hash**: Hash of the protein sequence.

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
getalleles genes.fasta assemblies/*.fna > alleles.tsv
```
OR
```bash
cat genes.fasta | getalleles - assemblies/*.fna -o alleles.tsv
```
To do the same but output the nucleotide sequences to a fasta file.
```bash
getalleles genes.fasta assemblies/*.fna -o alleles.tsv --ffn alleles.ffn
```
OR
```bash
getalleles genes.fasta assemblies/*.fna -o alleles.tsv --ffn - > alleles.ffn
```
To output the allele protein sequences into an MSA program:
```bash
getalleles genes.fasta assemblies/*.fna -o alleles.tsv --ffn --faa - | muscle
```
To output _one sequence file per assembly_ :
```bash
getalleles genes.fasta assemblies/*.fna -o alleles.tsv --ffn dna_seqs/ --faa protein_seqs/
```

## Acknowledgements
- [Dr. Ben Vezina](https://github.com/bananabenana) for ideas and initial implementation
of multiple copies in `v0.0.1b0`.
