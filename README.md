# RFMIX-PARSER

### Usage

```
usage: rfmix-parser.py [-h] -f MSPFILES -r REGIONS [-s SAMPLES] -o OUT

Parse RFMix output

options:
  -h, --help            show this help message and exit
  -f MSPFILES, --mspFiles MSPFILES
                        File listing RFMix MSP output files, one per line
  -r REGIONS, --regions REGIONS
                        File listing regions to parse, tab-separated with columns CHR, SPOS, EPOS
  -s SAMPLES, --samples SAMPLES
                        File listing samples to include, one per line
  -o OUT, --out OUT     Output file
```

### Input Examples

mspFiles: `msp.txt`
```
rfmix_chr1.msp.tsv
rfmix_chr2.msp.tsv
...
rfmix_chr22.msp.tsv
```

regions: `regions.txt`
```
CHR	SPOS	EPOS
chr1	50000000	52500000
chr3	20000000	30000000
chr4	40000000	50000000
chr10	30000000	40000000
chr10	50000000	60000000
chr17	10000000	20000000
chr17	30000000	40000000
chr19	15000000	15000200
```

samples: `samples.txt`
```
samp1
samp2
samp3
...
```

### Output Example

```
ID	CHR	SPOS	EPOS	Ancestry	Proportion
samp1	chr1	50000000	52500000	EUROPE|AFRICA	1.0|1.0
samp2	chr1	50000000	52500000	AFRICA|EUROPE	1.0|1.0
samp3	chr1	50000000	52500000	EUROPE|AFRICA	1.0|1.0
samp1	chr3	20000000	30000000	AMERICA|AFRICA:EUROPE	1.0|0.55:0.45
samp2	chr3	20000000	30000000	EUROPE|EUROPE	1.0|1.0
samp3	chr3	20000000	30000000	AMERICA|AFRICA:EUROPE	1.0|0.55:0.45
samp1	chr4	40000000	50000000	EUROPE|AMERICA	1.0|1.0
samp2	chr4	40000000	50000000	EUROPE:AMERICA|EUROPE:AFRICA	0.4:0.6|0.74:0.26
samp3	chr4	40000000	50000000	AMERICA|EUROPE	1.0|1.0
```

`Ancestry` provides all reported ancestries from the rfmix msp file in the specified region from SPOS to EPOS. Haplotypes are delimited by `|` and ancestry changes within a haplotype are delimited by `:`.

`Proportion` follows the same delimiting structure and provides the proportion of each ancestry (by base length) across the specified region.