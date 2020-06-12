import pandas as pd
import numpy as np
import datetime
import argparse
import json
import sys

parser = argparse.ArgumentParser(
    description="""Convert a kraken2 report file to jtree format.
    Results with percent coverage less than 1.00 will be filtered out.
    The output json can then be read by programs like the R package treeio and its function read.jtree().""")
parser.add_argument('-i', '--infile', required=True, help="input file in the kraken2 report format")
parser.add_argument('-o', '--outfile', required=True, help="output file that will be in the jtree json format")
args = parser.parse_args()

kraken_df = pd.read_csv(args.infile, sep="\t")
column_names = ['percent_fragments_covered', 'fragments_covered', 'fragments_assigned', 'rank_code', 'taxid', 'sciname']
kraken_df.columns = column_names
kraken_df = kraken_df[kraken_df['percent_fragments_covered'] > 1.00]

jtree_data = []
current_newick = ''
saved_newicks = pd.DataFrame()
previous_spaces = 0
increment = 0
unique_depths = []

for row in reversed(kraken_df.index):
    increment += 1
    new_spaces = kraken_df['sciname'][row].count(' ')
    new_newick = '%s:%s{%s}' % (kraken_df['sciname'][row].strip(), str(new_spaces), str(increment))
    new_dict = {
        'edge_num': int(increment),
        'percent_fragments_covered': float(kraken_df['percent_fragments_covered'][row]),
        'fragments_covered': int(kraken_df['fragments_covered'][row]),
        'fragments_assigned': int(kraken_df['fragments_assigned'][row]),
        'rank_code': str(kraken_df['rank_code'][row]),
        'taxid': str(kraken_df['taxid'][row])
    }
    jtree_data.append(new_dict)
    if new_spaces not in unique_depths:
        unique_depths.append(new_spaces)
    if len(saved_newicks) > 0:
        if len(saved_newicks.loc[saved_newicks.iloc[:,1] == new_spaces]) > 0:
                previous_newicks = ",".join(saved_newicks.loc[saved_newicks.iloc[:, 1] == new_spaces][0])
                new_newick = '{},{}'.format(new_newick, previous_newicks)
                saved_newicks = saved_newicks[saved_newicks.iloc[:, 1] != new_spaces]
    if len(current_newick) < 1:
        current_newick = new_newick
    elif new_spaces < previous_spaces:
        current_newick = '({}){}'.format(current_newick, new_newick)
    elif new_spaces == previous_spaces:
        current_newick = '{},{}'.format(current_newick, new_newick)
    elif new_spaces > previous_spaces:
        saved_newicks = saved_newicks.append(pd.Series([current_newick, previous_spaces]), ignore_index=True)
        current_newick = new_newick
    previous_spaces = new_spaces

current_newick = '({});'.format(current_newick)
jtree_dict = {
    'tree': str(current_newick),
    'data': jtree_data,
    'metadata': {'info': 'kraken2 report to json', 
                    'data': int(len(unique_depths)), 
                    'max_depth': int(len(unique_depths))}
}
with open(args.outfile, 'w') as outfile:
    json.dump(jtree_dict, outfile)