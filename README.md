# bioinfo-tools

Miscellaneous standalone tools for bioinformatic / computational biology work

## analyze_scope_cnv.R

A more user-friendly way to run the SCOPE CNV program

```
Rscript analyze_scope_cnv.R --help
```

SCOPE CNV resources: https://github.com/rujinwang/SCOPE

This program has high sensitivity and works best with bins between 1000kb and 500kb. Smaller bin sizes increase the computational time significantly and may result in errors

## kraken2_report_to_jtree.py

A way to parse the kraken2 output so it works better with other programs like the R package treeio and its function read.jtree()

```
python3 kraken2_report_to_jtree.py --help
```

Kraken2 resources: https://ccb.jhu.edu/software/kraken2/

## parse_genbank_annotations.py

A way to parse a gbff genbank annotation record and extract specific fields into a tab-delimited file

```
python3 parse_genbank_annotations.py --help
```

## filter_by_VAF.R

A tool for creating a blacklist of variants from a private cohort based on variant frequency, depth, and mapping quality, then using that blacklist to filter the variants from that private cohort

```
Rscript filter_by_VAF.R --help
```

Use the [vcf2tsv.py](https://github.com/sigven/vcf2tsv) tool to first convert VCF files to TSVs, then use my script to create blacklists and/or apply blacklists for filtering on your TSV of variants

vcf2tsv: https://github.com/sigven/vcf2tsv

This tool was inspired by this paper on filtering variants based on frequency: https://pubmed.ncbi.nlm.nih.gov/30591557/
