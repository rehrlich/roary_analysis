#!/usr/bin/env python3
"""
Plots of the pairwise tables make it easy to spot outliers.  This is an attempt
to use statistics to show that certain strains should be removed.
python3 pairwise_outliers.py /data/shared/homes/rachel/COPD/75_figures
"""

import sys
from glob import iglob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, deque
from statistics import median


class PairWiseComparison:
    def __init__(self, file_path):
        self.file_path = file_path
        self.pair_type, self.pair_counts = self.parse_file()
        self.stacked = self.pair_counts.stack(0)
        self.median_val, self.mad = self.calc_stats()

    def plot_counts(self, out_file):
        """
        Writes a histogram of all the data in self.pair_counts to out_file
        :param out_file:
        :return:
        """
        fig, ax = plt.subplots()
        title = "median is " + str(self.median_val) + " mad is " + str(self.mad)
        self.stacked.hist(ax=ax, bins=50)
        fig.suptitle(title)
        fig.savefig(out_file)

    def parse_file(self):
        """
        :return: the type (sim, diff, pair unique) and the upper triangular matrix
        """
        pair_type = self.file_path.split('pairwise_')[1].split('_table')[0]
        pair_counts = pd.read_table(self.file_path, index_col=0)
        return pair_type, pair_counts

    def calc_stats(self):
        median_val = self.stacked.median()
        mad = abs(self.stacked - median_val).median()
        return median_val, mad

    def find_outliers(self):
        """
        gets all scores associated with a strain.  Since the df is upper
        triangular, each score gets assigned to both the associated strains.
        returns a list of strains with low or high medians
        :return:
        """
        cluster_counts = defaultdict(list)

        for strain, counts in self.pair_counts.iteritems():
            cluster_counts[strain].extend([x for x in counts if not np.isnan(x)])

        for strain, counts in self.pair_counts.iterrows():
            cluster_counts[strain].extend([x for x in counts if not np.isnan(x)])

        thresh = self.mad * 2
        results = list()
        header = '\t'.join(['strain', 'comparison_type', 'direction',
                           'strain_median', 'pan_genome_median'])
        for strain, counts in cluster_counts.items():

            med = median(counts)
            if med < med + thresh < self.median_val:
                results.append('\t'.join([strain, self.pair_type, 'low',
                                          str(med),
                                          str(self.median_val)]))

            if med > med - thresh > self.median_val:
                results.append('\t'.join([strain, self.pair_type, 'high',
                                          str(med),
                                          str(self.median_val)]))

        return header, results


def main():
    roary_figs = sys.argv[1]
    output = deque()
    for file_name in iglob(roary_figs + '/*pairwise*i*table.txt'):
        pair_data = PairWiseComparison(file_name)
        pair_data.plot_counts(roary_figs + '/counts_hist_' +
                              pair_data.pair_type + '.pdf')
        header, results = pair_data.find_outliers()
        output.extend(results)

    output.appendleft(header)
    with open(roary_figs + '/pairwise_outliers.tsv', 'w') as f:
        f.write('\n'.join(output))


if __name__ == '__main__':
    main()