#!/usr/bin/env python

# Original file:  github.com/andrewjpage/Roary/blob/master/contrib/roary_plots/roary_plots.py
# Current version modified by Rachel Ehrlich

# Copyright (C) <2015> EMBL-European Bioinformatics Institute

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Neither the institution name nor the name roary_plots
# can be used to endorse or promote products derived from
# this software without prior written permission.
# For written permission, please contact <marco@ebi.ac.uk>.

# Products derived from this software may not be called roary_plots
# nor may roary_plots appear in their names without prior written
# permission of the developers. You should have received a copy
# of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.


import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import Phylo
import os

__author__ = "Marco Galardini"
__version__ = '0.1.0'


def get_options():
    # create the top-level parser
    description = "Create plots from roary outputs"
    parser = argparse.ArgumentParser(description=description,
                                     prog='roary_plots.py')

    parser.add_argument('tree', action='store',
                        help='Newick Tree file',
                        default='accessory_binary_genes.fa.newick')
    parser.add_argument('spreadsheet', action='store',
                        help='Roary gene presence/absence spreadsheet',
                        default='gene_presence_absence.csv')
    parser.add_argument('out_dir', action='store',
                        help='Output directory for figures',
                        default=os.getcwd())

    parser.add_argument('nickname', action='store',
                        help='Nickname for figures',
                        default='')

    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)

    return parser.parse_args()


def make_freq_plot(roary, outprefix):
    # Pangenome frequency plot
    plt.figure(figsize=(7, 5))

    plt.hist(roary.sum(axis=1), roary.shape[1],
             histtype="stepfilled", alpha=.7)
    plt.xlabel('Number of genomes')
    plt.ylabel('Number of genes')

    plt.savefig(outprefix + 'cluster_frequency.png')
    plt.clf()


def make_pie_chart(roary, outprefix):
    # Plot the pangenome pie chart
    plt.figure(figsize=(10, 10))
    core = roary[(roary.sum(axis=1) >= roary.shape[1] * 0.99) & (
        roary.sum(axis=1) <= roary.shape[1])].shape[0]
    softcore = roary[(roary.sum(axis=1) >= roary.shape[1] * 0.95) & (
        roary.sum(axis=1) < roary.shape[1] * 0.99)].shape[0]
    shell = roary[(roary.sum(axis=1) >= roary.shape[1] * 0.15) & (
        roary.sum(axis=1) < roary.shape[1] * 0.95)].shape[0]
    cloud = roary[roary.sum(axis=1) < roary.shape[1] * 0.15].shape[0]
    total = roary.shape[0]

    def my_autopct(pct):
        val = int(round(pct * total / 100.0))
        return '{v:d}'.format(v=val)

    a = plt.pie([core, softcore, shell, cloud],
                labels=['core\n(%d <= strains <= %d)' % (
                    roary.shape[1] * .99, roary.shape[1]),
                        'soft-core\n(%d <= strains < %d)' % (
                            roary.shape[1] * .95, roary.shape[1] * .99),
                        'shell\n(%d <= strains < %d)' % (
                            roary.shape[1] * .15, roary.shape[1] * .95),
                        'cloud\n(strains < %d)' % (roary.shape[1] * .15)],
                explode=[0.1, 0.05, 0.02, 0], radius=0.9,
                colors=[(0, 0, 1, float(x) / total) for x in
                        (core, softcore, shell, cloud)],
                autopct=my_autopct)
    plt.savefig(outprefix + 'pangenome_pie.png')
    plt.clf()


def get_roary_data(options):
    # Load roary
    roary = pd.read_table(options.spreadsheet,
                          sep=',',
                          low_memory=False)
    # Set index (group name)
    roary.set_index('Gene', inplace=True)
    # Drop the other info columns
    roary.drop(list(roary.columns[:10]), axis=1, inplace=True)
    # Transform it in a presence/absence matrix (1/0)
    roary.replace('.{2,100}', 1, regex=True, inplace=True)
    roary.replace(np.nan, 0, regex=True, inplace=True)
    return roary


def get_tree_name(options):
    # strips path
    if '/' in options.tree:
        tree_name = options.tree.split('/')[-1]
    else:
        tree_name = options.tree
    # tree_name = tree_name.rsplit('.', 1)[0]
    return tree_name


def plot_tree_heatmap(mdist, roary, roary_sorted, tree, tree_name, outprefix):

    # Plot presence/absence matrix against the tree
    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 10))

        ax1 = plt.subplot2grid((1, 40), (0, 10), colspan=30)
        a = ax1.matshow(roary_sorted.T, cmap=plt.cm.Blues,
                        vmin=0, vmax=1,
                        aspect='auto',
                        interpolation='none',
                        )

        # Creates an outline around the heatmap
        ax1.set_yticks([])
        ax1.set_xticks([])

        # Adjust colspan if strain names overlap heatmap
        ax = plt.subplot2grid((1, 40), (0, 0), colspan=7, axisbg='white')

        fig.subplots_adjust(wspace=0, hspace=0)

        ax1.set_title('Cluster possession matrix\nSorted by cluster frequency\n'
                      '(%d clusters)' % roary.shape[0])

        Phylo.draw(tree, axes=ax,
                   show_confidence=False,
                   xticks=([],), yticks=([],),
                   ylabel=('',), xlabel=('',),
                   xlim=(-0.01, mdist + 0.01),
                   axis=('off',),
                   title=('%s\n(%d strains)' % (tree_name, roary.shape[1]),),
                   do_show=False,
                   )
        plt.savefig(outprefix + 'tree_heatmap.png')
        plt.clf()


def main():
    options = get_options()
    outprefix = options.out_dir + '/' + options.nickname + '_'

    sns.set_style('white')

    tree = Phylo.read(options.tree, 'newick')

    tree_name = get_tree_name(options)

    # Max distance to create better plots
    mdist = max([tree.distance(tree.root, x) for x in tree.get_terminals()])

    roary = get_roary_data(options)

    # Sort the matrix by the sum of strains presence
    idx = roary.sum(axis=1).order(ascending=False).index
    roary_sorted = roary.ix[idx]

    make_freq_plot(roary, outprefix)

    # Sort the matrix according to tip labels in the tree
    roary_sorted = roary_sorted[[x.name for x in tree.get_terminals()]]

    plot_tree_heatmap(mdist, roary, roary_sorted, tree, tree_name,
                      outprefix)

    make_pie_chart(roary, outprefix)


if __name__ == "__main__":
    main()
