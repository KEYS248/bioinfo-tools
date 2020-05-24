# bioinfo-tools

Miscellaneous standalone tools for bioinformatic / computational biology work

## analyze-scope-cnv.R

A more user-friendly way to run the SCOPE CNV program

```
Rscript analyze-scope-cnv.R --help
```

SCOPE CNV resources: https://github.com/rujinwang/SCOPE

This program has high sensitivity and works best with bins between 1000kb and 500kb. Smaller bin sizes increase the computational time significantly and may result in errors

## kraken2-report-to-jtree.py

A way to parse the kraken2 output so it works better with other programs like the R package treeio and its function read.jtree()

```
python3 kraken2-report-to-jtree.py --help
```

Kraken2 resources: https://ccb.jhu.edu/software/kraken2/

## parse-genbank-annotations.py

A way to parse a gbff genbank annotation record and extract specific fields into a tab-delimited file

```
python3 parse-genbank-annotations.py --help
```