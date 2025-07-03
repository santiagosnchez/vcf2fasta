# vcf2fasta

> **this package is not actively maintained**

## Overview
`vcf2fasta` is a Python package designed to convert genomic regions and intervals into FASTA alignments using VCF files, GFF/GTF files, and a reference FASTA genome. This tool is particularly useful for bioinformaticians and researchers working with genomic data.

## Features
- Convert VCF files to FASTA format.
- Support for GFF/GTF and BED files for feature annotation.
- Options to include reference sequences and handle incomplete CDS.
- Command-line interface for easy usage.

## Installation
To install `vcf2fasta`, you can use pip:

```bash
pip install vcf2fasta
```

## Build
Build the package locally with:

```
python3 -m pip install -e .
```

## Usage
To use `vcf2fasta`, you can run the command line interface with the required arguments. Here is an example:

```bash
vcf2fasta -f genome.fas -v variants.vcf.gz -g intervals.gff -e CDS
```

### Command Line Options
- `-f`, `--fasta`: Path to the reference FASTA file (required).
- `-v`, `--vcf`: Path to the tabix-indexed VCF file (required).
- `-g`, `--gff`: Path to the GFF/GTF file (optional).
- `-b`, `--blend`: Concatenate GFF entries of the specified feature into a single alignment (default: False).
- `-i`, `--inframe`: Force the first codon of the sequence to be inframe (default: False).
- `-o`, `--out`: Name for the output directory (default: "vcf2fasta").
- `-r`, `--addref`: Include the reference sequence in the FASTA alignment (default: False).
- `-s`, `--skip`: Skip features without variants (default: False).

## Contribution
Contributions are welcome! If you would like to contribute to `vcf2fasta`, please fork the repository and submit a pull request. You can also open issues for any bugs or feature requests.

## License
This project is licensed under the MIT License. See the LICENSE file for more details.