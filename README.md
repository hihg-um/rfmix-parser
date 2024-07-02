# RFMIX-PARSER

### Usage

```
usage: rfmix-parser.py [-h] -f MSPFILE -r REGION [-m MARKERS] [-v VCF] -o OUT

Parse RFMix output

options:
  -h, --help            show this help message and exit
  -f MSPFILE, --mspFile MSPFILE
                        Path to RFMix .msp.tsv output file
  -r REGION, --region REGION
                        List with region to parse, comma-separated with CHR,SPOS,EPOS
  -m MARKERS, --markers MARKERS
                        List with markers to append to output, providing genotype and ancestry, comma-separated in format chr:pos
  -v VCF, --vcf VCF     Path to Phased VCF file, must be indexed. Must be used in conjunction with --markers
  -o OUT, --out OUT     Output file

```
