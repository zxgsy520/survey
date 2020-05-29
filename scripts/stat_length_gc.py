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


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
       fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield seq[0], seq[1]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
       fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip("@").split()[0]
            continue
        if line.startswith('@'):
            line = line.strip("@").split()[0]
            seq = seq.split('\n')

            yield seq[0], seq[1]
            seq = ''
            seq = "%s\n" % line
        else:
            seq += "%s\n" % line

    if len(seq.split('\n'))==5:
        seq = seq.split('\n')
        yield seq[0], seq[1]


def read_depth(file):
    '''Read depth file'''

    for line in open(file, 'r'):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        parts = line.split('\t')

        yield [parts[0], int(parts[1]), int(parts[2])]


def stat_base(string):

    string = string.upper()
    g = string.count('G')
    c = string.count('C')
    a = string.count('A')
    t = string.count('T')

    return a, t, c, g


def sum_depth(file):
    
    depth_dict = {}
    
    for line in read_depth(file):
        if line[0] not in depth_dict:
            depth_dict[line[0]] = line[2]
            continue
        depth_dict[line[0]] += line[2]

    return depth_dict


def stat_base_length(genome, file, name):

     if file=="":
         depth_dict = {}
         logger.info("No depth file")
     else:
         depth_dict = sum_depth(file)

     output = open('%s.length_gc.xls' % name, 'w')
     output.write("#Id\tLength\tGC(%)\tCoverage\n")
    
     if genome.endswith(".fastq") or genome.endswith(".fq") or genome.endswith(".fastq.gz") or genome.endswith(".fq.gz"): 
         fh = read_fastq(genome)
     elif genome.endswith(".fasta") or genome.endswith(".fa") or genome.endswith(".fasta.gz") or genome.endswith(".fa.gz"):
         fh = read_fasta(genome)
     else:
         raise Exception("%r file format error" % genome)

     for seq_id, seq in fh:
         a, t, c, g = stat_base(seq)
         length = a+t+c+g
         
         if seq_id not in depth_dict:
             cover = 0
         else:
             cover = depth_dict[seq_id]*1.0/length

         output.write('{}\t{:,}\t{:.2f}\t{:.2f}\n'.format(seq_id, length, (g+c)*100.0/length, cover))
     output.close()


def add_help(parser):

    parser.add_argument('-g', '--genome', metavar='FILE', type=str, required=True,
        help='Input genome file.')
    parser.add_argument('-d', '--depth', metavar='FILE', type=str, default='',
        help='Input depth file.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default='out',
        help='Set the prefix of the output file.')

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
    splitfp.py  Split a specific format for multiple files.

attention:
    splitfp.py -i fasta.list
    splitfp.py -i fasta.list -n 100000 -o name

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    stat_base_length(args.genome, args.depth, args.name)


if __name__ == "__main__":

    main()
