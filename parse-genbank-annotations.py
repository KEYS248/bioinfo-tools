from Bio import SeqIO
import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(
    description="""Parse a refseq genbank gbff file for specific fields and output data as tab-delimited file.
    You can find these files at NCBI FTP sites such as these: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/""")
parser.add_argument('-i', '--infile', required=True, help="input file in the gbff format")
parser.add_argument('-o', '--outfile', required=True, help="output file that will be the parsed data as a tab-delimited file")
args = parser.parse_args()

column_names = ['accession', 'refseq_version', 'genome_size', 'source', 'plasmid',
                'genes', 'cds', 'coding_genes', 'protein_cds', 'rRNAs', 'tRNAs', 'ncRNAs',
                'pseudogenes']
annotations_df = pd.DataFrame()

without_annotations = 0
for record in SeqIO.parse(args.infile, "genbank"):
    before, keyword, after = record.description.partition("plasmid ")
    plasmid = after.split(" ")[0][:-1]
    source = " ".join(before.split(" ")[:-1])
    if len(plasmid.strip()) == 0:
        plasmid = 'unnamed'
    if len(source.strip()) == 0:
        source = 'unknown'
    genus = record.description.split(" ")[0]
    species = " ".join(record.description.split(" ")[0:2])
    try:
        record_annotation = record.annotations['structured_comment']['Genome-Annotation-Data']
        genes = record_annotation['Genes (total)']
        cds = record_annotation['CDSs (total)']
        coding_genes = record_annotation['Genes (coding)']
        protein_cds = record_annotation['CDSs (with protein)']
        coding_genes = record_annotation['Genes (coding)']
        split_rRNAs = record_annotation['rRNAs'].split(" ")
        if len(split_rRNAs) > 4:
            rRNAs = int(split_rRNAs[0][:-1]) + int(split_rRNAs[1][:-1]) + int(split_rRNAs[2])
        elif len(split_rRNAs) > 2:
            rRNAs = int(split_rRNAs[0][:-1]) + int(split_rRNAs[1])
        elif len(split_rRNAs) == 2:
            rRNAs = int(split_rRNAs[0])
        else:
            rRNAs = 0
        tRNAs = record_annotation['tRNAs']
        ncRNAs = record_annotation['ncRNAs']
        pseudogenes = record_annotation['Pseudo Genes (total)']
        if len(str(pseudogenes).strip()) == 0:
            pseudogenes = 0
    except:
        genes, cds, coding_genes, protein_cds, coding_genes = 0, 0, 0, 0, 0
        protein_cds, rRNAs, tRNAs, ncRNAs, pseudogenes = 0, 0, 0, 0, 0
        without_annotations += 1
    annotations_df = annotations_df.append(pd.Series([record.name, record.id, len(record.seq), genus, species,
                                                      genes, cds, coding_genes, protein_cds, rRNAs, tRNAs, ncRNAs,
                                                      pseudogenes]), ignore_index=True)
print("{} records processed".format(len(annotations_df)))
print("\t{} of which did had some missing genome annotation data".format(without_annotations))
annotations_df.columns = column_names
annotations_df.head()
annotations_df.to_csv(args.outfile, sep="\t", index=False)#, header=False)