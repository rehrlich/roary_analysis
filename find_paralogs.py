#!/usr/bin/env python2
"""
Author:  Rachel Ehrlich
Input is the path to a roary folder that was created using the
--dont_delete_files flag, the output directory and a nickname
for the analysis.
Outputs four .csv files
gene_presence_absence_paralogs_annotated is the original roary output
with columns added for paralog counts and lists
gene_presence_absence_paralogs_merged is formatted like the roary output
but entries for paralogs have been combined using tabs.  This is also
written to the input folder.
paralog_table maps each gene to its paralogs
summary_statistics_paralogs_merged is like the summary statistics from the
original output but has the gene counts from before splitting paralogs.
"""
import numpy as np
import sys
from itertools import izip


# Input is the folder with the roary output files
# This uses the gene_presence_absence.csv file 
# Returns the cluster names, a list of lists with the genes
# from each cluster and a list of lists with all the input data
def get_pres_abs_data(folder):
    with open(folder + "/gene_presence_absence.csv", 'rU') as f:
        lines_end_quotes = f.read().split('\n')
    lines = [x[1:-1] for x in lines_end_quotes]   
    
    # The file is a csv but the annotations contain commas...
    split_lines = [x.split('","') for x in lines if len(x) > 0]
    
    rows = split_lines[1:]
    cluster_names = [x[0] for x in rows]
    gene_data = [' '.join(row[11:]).split() for row in rows]

    return cluster_names, np.array(gene_data), np.array(split_lines)


# Input is a roary folder
# Returns the names for the clusters after splitting paralogs
# Returns two lists of lists with the genes in each cluster
# (split and unsplit)
def get_data(folder):
    cluster_names, split_clusters, all_data = get_pres_abs_data(folder)
    with open(folder + "/_inflated_unsplit_mcl_groups", 'rU') as f:
        unsplit_clusters = [x.split()
                            for x in f.read().split('\n') if len(x) > 0]
        
    return cluster_names, split_clusters, np.array(unsplit_clusters)


# Input is a list of split clusters names and a list of lists with list
# i containing the genes in cluster name i.
# Returns an inverted index that maps the gene name (prokka something) to
# the gene cluster's name
def make_split_index(cluster_names, split_clusters):
    inverted_index = dict()
    for name, genes in izip(cluster_names, split_clusters):
        for gene in genes:
            inverted_index[gene] = name
    return inverted_index


# Input: inverted index that maps the gene name (prokka something) to
# the gene cluster's name, list of lists containing gene names for the 
# unsplit clusters
# Returns a list of sets where each set is a group of paralog clusters
def map_unsplit_to_split(inverted_index, unsplit_clusters):
    paralogs = list()
    
    for cluster in unsplit_clusters:
        para_group = set()
        for gene in cluster:
            para_group.add(inverted_index[gene])
        paralogs.append(para_group)
    return paralogs


# Input is a list of sets of paralogs, and the list of all split clusters
# Assert fails if a cluster is in multiple paralog groups or if a cluster
# is missing
def check_paralogs_unique(paralogs, cluster_names):
    num_clusters = sum(len(x) for x in paralogs)
    joined_paralogs = set().union(*paralogs)
    assert(len(joined_paralogs) == num_clusters)
    assert(len(cluster_names) == num_clusters)


# Input is a list of sets of paralogs
# Outputs a list of lists where each list represents one split cluster
# The list has the number of clusters (including self) in the paralog group
# The second entry is a tab separated list of those clusters
def make_output_data(paralogs):
    out_data = dict()
    for clusters in paralogs:
        num_clusters = str(len(clusters))
        c_list = '\t'.join(clusters)
        for cluster in clusters:
            out_data[cluster] = [num_clusters, c_list]
    return out_data


# Input is a list of rows which are paralogs of each other
# Output is a single row where each column is a tab separated list of entries
# for all the input data
def merge_rows(paralog_rows):
    if len(paralog_rows) == 1:
        return paralog_rows[0]
    
    new_row = list()
    
    # This iterates over columns in the untransposed array
    for column in np.array(paralog_rows).T:
        new_entry = '\t'.join(x for x in column if len(x) > 0)
        new_row.append(new_entry)
    return new_row


# Inputs are a list of lists with all the input roary data and a 
# list of sets of paralogs
# Output is a list of lists of roary style data in which paralog
# groups have been combined.  These rows are sorted according to the 
# original file
def merge_paralogs(all_data, paralogs):
    row_num = dict()
    for index, row in enumerate(all_data[1:], start=1):
        row_num[row[0]] = index
    
    merged_output = []  
    
    for clusters in paralogs:
        paralog_rows = list()
        first_index = sys.maxint
        for cluster in clusters:
            index = row_num[cluster]
            first_index = min(first_index, index)
            paralog_rows.append(all_data[index])
        merged_output.append([first_index, merge_rows(paralog_rows)])
    
    # Sorts, removes sort index, adds heading
    merged_output.sort(key=lambda y: y[0])
    merged_output = [y[1] for y in merged_output]
    merged_output.insert(0, all_data[0])

    return merged_output


# Input: prefix - path + nickname
# file_name - suffix w/o extension, data is a list of lists
# Write data to a file with each entry being double quoted and commas and 
# newlines separating entries
def write_output(prefix, file_name, data):
    output_str = '\n'.join('"' + ('","'.join(x)) + '"' for x in data)
    with open(prefix + file_name + ".csv", 'w') as f:
        f.write(output_str)


# Inputs:
# out_data list of lists with the number of clusters (including self) 
# in the paralog group and a tab separated list of those clusters
# in_folder - roary output
# paralogs - list of sets of paralogs
# out_folder - output location
# nickname - prefix for output files
def make_output(out_data, in_folder, paralogs, out_folder, nickname):

    cluster_names, split_clusters, all_data = get_pres_abs_data(in_folder)
    
    ordered_out_data = [['num_paralogs', 'paralog_group']]
    for name in cluster_names:
        ordered_out_data.append(out_data[name])
    ordered_out_data = np.vstack(ordered_out_data)
    
    prefix = out_folder + '/' + nickname + '_'
    
    combined_output = np.hstack((all_data[:, :2], ordered_out_data,
                                 all_data[:, 2:]))
    write_output(prefix, "gene_presence_absence_paralogs_annotated",
                 combined_output)

    paralog_translation = np.hstack((all_data[:, :2], ordered_out_data))
    write_output(prefix, "paralog_table", paralog_translation)
        
    collapsed_paralogs = merge_paralogs(all_data, paralogs)
    write_output(prefix, "gene_presence_absence_paralogs_merged",
                 collapsed_paralogs)
    write_output(in_folder, "gene_presence_absence_paralogs_merged",
                 collapsed_paralogs)
    return collapsed_paralogs


# Input: a float
# Output: a string of the input as a percentage
def percent(decimal):
    return str(int(decimal * 100)) + "%"


# Inputs:
# gene_counts - a list of the number of strains that have a copy of a cluster
# (lo, hi, name) - cutoff fractions and the name of the cluster type
# total_strains - an int
# Returns a string describing the number of clusters found in some fraction
# of strains.  Ex "Core genes (99% <= strains < 100%):	137"
def get_num_clusters(gene_counts, (lo, hi, name), total_strains):
    num = sum(1 for x in gene_counts
               if x >= lo * total_strains and x < hi * total_strains)
    txt = ''.join([name, " genes (", percent(lo), " <= strains < ",
                   percent(hi), "):\t", str(num)])
    return txt


# Inputs:  a list of list with the gene +- matrix (including headings) with the
# paralogs collapsed, roary input folder
# Writes summary statistics to the roary folder
def make_summary_stats(collapsed_paralogs, in_folder, out_folder):
    gene_data = collapsed_paralogs[1:, 11:]
    total_strains = len(collapsed_paralogs[1, 11:])

    has_gene = np.vectorize(lambda x: len(x) > 0, otypes=[bool])
    gene_presence = has_gene(gene_data)
    cluster_counts = gene_presence.sum(axis=1)

    cutoffs = [(0.99, 1.01, "Core"), [0.95, 0.99, "Soft core"],
               [0.15, 0.95, "Shell"], [0.0, 0.15, "Cloud"]]

    output = [get_num_clusters(cluster_counts, cutoff, total_strains)
              for cutoff in cutoffs]
    output.append("Total genes:\t" + str(len(gene_data) - 1))
    output = '\n'.join(output).replace("< 101", "<= 100")
                  
    with open(in_folder + "/summary_statistics_paralogs_merged.csv", 'w') as f:
        f.write(output)

    with open(out_folder + "/summary_statistics_paralogs_merged.csv", 'w') as f:
        f.write(output)


# Makes gene_presence_absence.Rtab but after unsplitting paralogs
# Inputs:  a list of list with the gene +- matrix (including headings) with the
# paralogs collapsed
# Writes summary statistics to the roary folder
def make_gpa_rtab(collapsed_paralogs, in_folder):
    strain_headings = collapsed_paralogs[0, 11:]
    gene_data = collapsed_paralogs[1:, 11:]
    gene_names = np.char.replace(collapsed_paralogs[1:, 0], '\t', ',')

    has_gene = np.vectorize(lambda x: str(int(len(x) > 0)), otypes=[str])
    pres_abs = has_gene(gene_data)

    with open(in_folder + '/gene_presence_absence_paralogs_merged.Rtab', 'w') as f:
        f.write('Gene_name\t')
        f.write('\t'.join(strain_headings))
        f.write('\n')
        f.write('\n'.join('\t'.join(x) for x in
                          np.column_stack((gene_names, pres_abs))))


def main():
    in_folder = sys.argv[1]
    out_folder = sys.argv[2]
    nickname = sys.argv[3]
    
    cluster_names, split_clusters, unsplit_clusters = get_data(in_folder)
    gene_name_to_cluster_dict = make_split_index(cluster_names, split_clusters)
    paralogs = map_unsplit_to_split(gene_name_to_cluster_dict, unsplit_clusters)
    
    check_paralogs_unique(paralogs, cluster_names)
    
    cluster_to_paralog_counts = make_output_data(paralogs)

    collapsed_paralogs = np.array(make_output(cluster_to_paralog_counts,
                                              in_folder, paralogs, out_folder,
                                              nickname))

    make_gpa_rtab(collapsed_paralogs, in_folder)
    make_summary_stats(collapsed_paralogs, in_folder, out_folder)
    
    
if __name__ == "__main__":
    main()
