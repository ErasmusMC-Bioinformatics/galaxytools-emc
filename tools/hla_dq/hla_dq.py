#!/usr/bin/env python
'''
Given HLA-A and HLA-B genes annotated with BLAST IMGT/HLA database, determine type


DQA1 | DQB1 | type
--------------------
02:01| 02:02| DQ2.2
03:03| 02:02| DQ2.3
05:01| 02:01| DQ2.3
03:01| 03:02| DQ8
03:02| 03:02| DQ8
03:03| 03:02| DQ8

Annotations are of form DQ[A|B]1*xx:yy[:zz[:vv]] where xx:yy is of interest.

Example: "HLA:HLA11066 DQA1*01:05:02 768 bp"
'''

import argparse
import itertools


def to_matrix(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]


def get_list_of_associated_types(filename, column):
    with open(filename) as f:
        contents = f.readlines()

    return [':'.join(line.split('\t')[column-1].split(' ')[1].split('*')[1].split(':')[0:2]) for line in contents[1:]]


def get_associations(typesA, typesB):
    ''' Given list of genotype annotations (e.g. DQA1*02:01..) from A and B, determine possible associated serotypes '''
    ''' each combination of DQA1,DQB1,type '''
    associated_combinations = [
        ['02:01', '02:02', 'DQ2.2'],
        ['03:03', '02:02', 'DQ2.3'],
        ['05:01', '02:01', 'DQ2.5'],
        ['03:01', '03:02', 'DQ8'],
        ['03:02', '03:02', 'DQ8'],
        ['03:03', '03:02', 'DQ8']]

    return [a[2] for a in associated_combinations if a[0] in typesA and a[1] in typesB]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-A', required=True, action='append', help='BLAST hits A gene')
    parser.add_argument('-B', required=True, action='append', help='BLAST hits B gene')
    parser.add_argument('-c', '--column', default=5, type=int, help='Column number containing the BLAST annotation')
    args = parser.parse_args()

    # TODO: QC check that file A contains DQA1 annotations and B file contains DQB1 annotations?

    # find possible associated types, for all combinations of alleles
    typesA = [get_list_of_associated_types(A, args.column) for A in args.A]
    typesB = [get_list_of_associated_types(B, args.column) for B in args.B]
    associations = [get_associations(c[0], c[1]) for c in itertools.product(typesA, typesB)]
    associations = to_matrix(associations, len(args.B))

    # write output table
    headerline = '\t'+'\t'.join(['B'+str(i+1) for i in range(0, len(args.B))])+'\n'
    bcount = 0
    with open('results.tsv', 'w') as outfile:
        outfile.write(headerline)
        for line in associations:
            bcount += 1
            outfile.write('A'+str(bcount)+'\t'+'\t'.join([';'.join(sorted(set(l))) if l else '-' for l in line])+'\n')
