#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np


def split_core(roary_out):
    core_genes = []
    with open(roary_out + '/gene_presence_absence.Rtab', 'r') as f:
        for line in f.readlines():
            if line.startswith('Gene\t') or len(line) < 5:
                continue
            split_line = line.split('\t', 1)
            # 1 is present, 0 is absent
            if '0' in split_line[1]:
                continue
            core_genes.append(split_line[0])
    return core_genes


def merged_core(roary_out):
    merged_data = pd.read_csv(roary_out + '/gene_presence_absence_paralogs_merged.csv',
                              index_col=0)
    print(merged_data.head())
    data_headers = 'Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC'
    data_headers = data_headers.split('","')

    for col in data_headers:
        del merged_data[col]

    merged_data.dropna(inplace=True)
    return merged_data.index


def main():
    roary_out = sys.argv[1]

    # core_genes = split_core(roary_out)
    # with open(roary_out + '/core_genes_list.txt', 'w') as f:
    #     f.write(','.join(core_genes))

    core_genes = merged_core(roary_out)
    with open(roary_out + '/core_genes_list_merged.txt', 'w') as f:
        f.write(','.join(core_genes))

main()


