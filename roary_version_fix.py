#!/usr/bin/env python3

"""
Author:  Rache Ehrlich
Different versions of Roary format the output files differently.  Most of these
scripts were written for roary 3.2.5 but currently the cluster is running 3.5.1.
This program rewrites the gene_presence_absence.csv so that it has the number of
metadata columns found in 3.2.5.
Currently only considering 3.2.5 and 3.5.1
"""
import sys
import os
from subprocess import call
import numpy as np


class GenePresAbs:
    """
    Class invariant:  self.version is the roary version of the data in
    self.file_path
    """
    header_351 = '"Gene","Non-unique Gene name","Annotation","No. isolates",' \
                 '"No. sequences","Avg sequences per isolate","Genome Fragment",' \
                 '"Order within Fragment","Accessory Fragment",' \
                 '"Accessory Order with Fragment","QC","Min group size nuc",' \
                 '"Max group size nuc","Avg group size nuc"'
    extra_351_cols = ["Min group size nuc", "Max group size nuc",
                      "Avg group size nuc"]

    def __init__(self, roary_dir):
        self.file_path = roary_dir + '/gene_presence_absence.csv'
        assert os.path.isfile(self.file_path)

        self.original_path = self.file_path.replace('.csv', '_original.csv')
        if not os.path.isfile(self.original_path):
            call(['cp', self.file_path, self.original_path])

        self.version = self.get_version(self.original_path)

    @staticmethod
    def get_version(file_path):
        with open(file_path, 'r') as f:
            header = f.readline()
            if header.startswith(GenePresAbs.header_351):
                return '3.5.1'
            else:
                return '3.2.5'

    def get_data_325(self):
        data351 = np.genfromtxt(self.original_path, delimiter='","',
                                dtype=str)
        indices = get_column_indices(data351, col_headings=self.extra_351_cols)
        return np.delete(data351, indices, axis=1)

    def write_325(self):
        if self.version == '3.2.5':
            return
        data_325 = self.get_data_325()
        data_325 = '\n'.join('","'.join(x) for x in data_325)
        with open(self.file_path, 'w') as f:
            f.write(data_325)
            f.write('\n')
        self.version = '3.2.5'


def get_column_indices(mat, col_headings):
    positions = [np.where(mat[0, :] == heading) for heading in col_headings]
    return [x[0][0] for x in positions]


def main():
    roary_outdir = sys.argv[1]
    gpa = GenePresAbs(roary_outdir)
    if not gpa.version == '3.2.5':
        gpa.write_325()

if __name__ == "__main__":
    main()
