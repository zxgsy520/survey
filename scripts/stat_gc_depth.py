#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_depth(files):
    '''Read depth file'''

    for line in open(files, 'r'):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        parts = line.split('\t')

        yield [parts[0], int(parts[1]), int(parts[2])]


def read_fasta(fasta):
    '''Read fasta file'''

    seq = ''

    for line in open(fasta, 'r'):
        line = line.strip()

        if not line:
            continue
        if line.startswith('>'):
            if len(seq.split('\n'))==2:
                seq = seq.split('\n')
                yield seq[0], seq[1]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    if len(seq.split('\n'))==2:
        seq = seq.split('\n')
        yield seq[0], seq[1]


def stat_gc(string):

    string = string.upper()
    gc = string.count('G') + string.count('C')
    at = string.count('A') + string.count('T')

    if (gc+at)==0:
        gc_content = 0
    else:
        gc_content = gc*100.0/(gc+at)

    return gc_content


def run_stat_gc_depth(seq_id, seq, depth_list, output, bins, wind, error=100):

    seq_len = len(seq)

    if len(depth_list)>=seq_len-error: #Allowable depth error range is 100pb
        for i in range(0, seq_len, bins):
            if i>=seq_len:
                continue

            end = i+wind

            if end>=seq_len:
                end = seq_len

            string = seq[i:end]
            gc_content = stat_gc(string)

            if gc_content==0:
                continue
            depth = sum(depth_list[i:end])*1.0/(end-i)
            output.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(seq_id, i+1, end, gc_content, depth))

            if end>=seq_len:
                break
    else:
        LOG.debug("The length of the sequence %s is inconsistent with the length of the sequencing depth list. Use the (samtools depth -aa) parameter." % seq_id)

    return seq_id


def stat_gc_depth(depth_file, fasta, name, bins, wind, error):

    geneme_dict = {}

    for line in read_fasta(fasta):
        id = line[0].strip().split()[0].replace('>', '')
        geneme_dict[id] = line[1]

    id = ''
    depth_list = []
    output = open('%s.stat_gc_depth.tsv' % name, 'w')

    output.write('#Sequence_id\tstart\tend\tgc_content\tdepth\n')
    for line in read_depth(depth_file):
        if len(geneme_dict[line[0]])<2000 or len(geneme_dict[line[0]])<(wind/2):
            continue

        if not id:
            id = line[0]
            depth_list.append(int(line[2]))
            continue

        if id!=line[0]:
            id = run_stat_gc_depth(id, geneme_dict[id], depth_list, output, bins, wind, error)
            id = line[0]
            depth_list = []
            depth_list.append(int(line[2]))
            continue

        depth_list.append(int(line[2]))

    id = run_stat_gc_depth(id, geneme_dict[id], depth_list, output, bins, wind, error)
    output.close()


def add_gc_depth_help(parser):

    parser.add_argument('-d', '--depth', metavar='FILE', type=str, required=True,
        help='Input the coverage depth statistics file.')
    parser.add_argument('-g', '--geneme', metavar='FILE', type=str, required=True,
        help='Input the genome file (fasta).')
    parser.add_argument('-b', '--bins', metavar='INT', type=int, default=1000,
        help='Select interval size,default=1000.')
    parser.add_argument('-w', '--wind', metavar='INT', type=int, default=5000,
        help='Selected window size,default=5000.')
    parser.add_argument('-e', '--error', metavar='INT', type=int, default=1000,
        help='Allowable error range,default=1000bp.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default='out',
        help='Output file name.')

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
    stat_gc_depth.py  Statistical GC depth.

attention:
    stat_gc_depth.py -d *.depth -g geneme.fa -w 5000 -n name
    stat_gc_depth.py -d *.depth -g geneme.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args =add_gc_depth_help(parser).parse_args()

    stat_gc_depth(args.depth, args.geneme, args.name, args.bins, args.wind, args.error)


if __name__ == "__main__":

    main()
