#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def best_blast(file):

    blast_dict = {}
    for line in read_tsv(file):
        if line[0] not in blast_dict:
            blast_dict[line[0]] = line
        else:
            if int(line[3])<int(blast_dict[line[0]][3]):
                continue
            if float(line[2])<float(blast_dict[line[0]][2]):
                continue
            blast_dict[line[0]] = line

    return blast_dict


def run_obtain_taxonomy(m6_file, taxonomy, name):

    reads_dict = {}
    blast_dict = best_blast(m6_file)
    output = open('%s.species_annotation.txt' % name, 'w')
    
    
    for line in blast_dict.values():
        if line[-1] not in reads_dict:
            reads_dict[line[-1]] = []
        reads_dict[line[-1]].append(line[0])

    for line in read_tsv(taxonomy):
        if line[0] in reads_dict:
            for i in reads_dict[line[0]]:
                output.write('%s\t%s\t%s\t%s\n' % (i, line[0], line[1], line[2]))
        else:
            print('%s does not have corresponding species classification information' % line[0])
    output.close()


def add_args(parser):

    parser.add_argument('-i', '--input', metavar='m6', type=str, required=True,
        help='Input balst comparison file.')
    parser.add_argument('-t', '--taxonomy', metavar='FILE', type=str, default= '/export2/master2/sunzy/software/species.taxonomy',
        help='Input the species classification file.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default= 'out',
        help='Output file prefix')
    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:

    obtain_taxonomy.py -- Get the classification information of reads by comparing files

attention:
    obtain_taxonomy.py -i txt.fa.m6 -n name
    obtain_taxonomy.py -i txt.fa.m6 -t /export2/master2/sunzy/software/species.taxonomy -n out
''')
    args = add_args(parser).parse_args()

    run_obtain_taxonomy(args.input, args.taxonomy, args.name)


if __name__ == "__main__":
    main()
