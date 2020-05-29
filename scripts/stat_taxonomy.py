#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


PROK = {'Archaea', 'Bacteria', 'Viruses', 'Viroids'}


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_tax(string):

    string = string.strip(';').split(';')
    euka = {'Fungi', 'Viridiplantae', 'Alveolata', 'Metazoa'}
    kind = " Other"
    classify = "Other"

    for i in string:
        if i in PROK:
            kind = "prokaryotic"
            classify = i
            break
        if i in euka:
            kind = "eukaryote"
            classify = i
            break

    return kind, classify, string[-1]


def stat_taxonomy(file, reads_number, name):

    map_number = 0
    species_dict = {}
    prok = []
    boundary_dict = {'Archaea':0, 'Bacteria':0 , 'Fungi':0, 'Viridiplantae':0,
        'Alveolata':0, 'Metazoa':0 , 'Viruses':0, 'Viroids':0, 'Other':0}

    for line in read_tsv(file, '\t'):
        taxid = line[1]
        map_number +=1

        kind, classify, species = read_tax(line[3])
        boundary_dict[classify] += 1
        if kind =="prokaryotic":
            prok.append(taxid)

        species = "%s\t%s" % (classify, line[2])
        if species not in species_dict:
            species_dict[species] = 0
        species_dict[species] += 1

    fb = open('%s.species_classify.tsv' % name, 'w')
    fb.write('#Bacteria\tArchaea\tFungi\tViridiplantae\tAlveolata\tMetazoa\tViruses\tViroids\tOther\n')
    fb.write('{0:,}\t{1:,}\t{2:,}\t{3:,}\t{4:,}\t{5:,}\t{6:,}\t{7:,}\t{8:,}\n'.format(
        boundary_dict['Bacteria'],
        boundary_dict['Archaea'],
        boundary_dict['Fungi'],
        boundary_dict['Viridiplantae'],
        boundary_dict['Alveolata'],
        boundary_dict['Metazoa'],
        boundary_dict['Viruses'],
        boundary_dict['Viroids'],
        boundary_dict['Other'])
    )

    n = 0
    cont = 0
    fs = open('%s.stat_species.tsv' % name, 'w')
    fs.write('Unknow\tUnmap\t{}\n'.format(reads_number-map_number))

    for i in sorted(species_dict.items(),key = lambda x:x[1],reverse = True):
        fs.write('{0}\t{1}\n'.format(i[0], i[1]))
        n +=1
        if n <=10 and i[0].split('\t')[0] in PROK:
            cont += 1

    fs.close()


    if len(prok)*100.0/reads_number>=20:
        LOG.info("The sample is seriously contaminated, it is recommended to use strict mode.")
    ft = open('%s.prokaryotic.taxid' % name, 'w')
    ft.write('#Prokaryote ratio:{1:.2f}\tTop 10 numbers:{1}\n'.format(len(prok)*100.0/reads_number, cont))
    for i in set(prok):
        ft.write('%s\n' % i)
    ft.close()


def add_args(parser):

    parser.add_argument('-i', '--input', metavar='m6', type=str, required=True,
        help='Input the species annotation file.')
    parser.add_argument('-rn', '--reads_number', metavar='INT', type=int, default=10000,
        help='Number of reads participating in the comparison,default=10000.')
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

    stat_taxonomy.py -- Statistically annotated species information

attention:
    stat_taxonomy.py -i txt.species_annotation.txt -n name
    stat_taxonomy.py -i txt.species_annotation.txt -rn 1000000 -n out
''')
    args = add_args(parser).parse_args()

    stat_taxonomy(args.input, args.reads_number, args.name)


if __name__ == "__main__":
    main()
