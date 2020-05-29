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
        if line.startswith(">"):
            if seq!='':
                yield seq
            seq = ''

        seq += "%s\n" % line

    if seq!='':
        yield seq


def read_fastq(file):
    '''Read fastq file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''
    n = 0

    for line in fp:
        line = line.strip()

        if not line:
            continue
        n +=1
        seq += "%s\n" % line
        if n==4:
            n=0
            yield seq
            seq = ''


def read_ngs_files(files):

    r1_list = []
    r2_list = []

    for line in open(files):
        line = line.strip()

        if not line:
            continue

        line1 = line.split('\t')
        if len(line1)!=2:
            LOG.info('The %s line is not two files, check if there is a space in the middle' % line)
            continue
        if os.path.exists(line1[0]) and os.path.exists(line1[1]):
            r1_list.append(os.path.abspath(line1[0]))
            r2_list.append(os.path.abspath(line1[1]))
        else:
            LOG.info('File %s or %s may not exist' % (line1[0], line1[0]))
            continue

    return r1_list, r2_list


def read_files(files):
    '''Read file list'''

    file_list = []

    for i in open(files):
        i = i.strip()
        if not i:
            continue

        if os.path.exists(i):
            i = os.path.abspath(i)
            file_list.append(i)
        else:
            LOG.info('File %s may not exist' % i)
            continue

    return file_list


def split_fp(file_list, workdir, name, number, format='fastq'):
    '''Split the data by the number of reads'''

    n = 0
    lable = 1
    output = open('{}/{}.part_{}.{}'.format(workdir, name, lable, format), 'w')

    for i in file_list:
        LOG.info('Read %s file' % i)
        if format=='fastq':
            fh = read_fastq(i)
        else:
            fh = read_fasta(i)

        for line in fh:
            if n >= number:
                output.write(line)
                output.close()
                n = 0
                lable += 1
                output = open('{}/{}.part_{}.{}'.format(workdir, name, lable, format), 'w')

            output.write(line)
            n += 1

    output.close()


def split_data(files, workdir, name, no_ngs, number):

    if not os.path.exists(workdir):
        os.makedirs(workdir)
    else:
        LOG.info('%s file exists' % workdir)

    if no_ngs:
        file_list = read_files(files)
        name = name
     
        if file_list[0].endswith(".fastq") or file_list[0].endswith(".fq") or file_list[0].endswith(".fastq.gz") or file_list[0].endswith(".fq.gz"):
            format = 'fastq'
        else:
            format = 'fasta'

        split_fp(file_list, workdir, name, number, format)

    else:
        format = 'fastq'
        name1 = name+'.r1'
        name2 = name+'.r2'
        r1_list, r2_list = read_ngs_files(files)

        if len(r1_list)==0:
            raise Exception('Can not find the input file, please check the input file format.')

        split_fp(r1_list, workdir, name1, number, format)
        split_fp(r2_list, workdir, name2, number, format)


def split_fq_help(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Set the input file.')
    parser.add_argument('--no_ngs', action="store_true",
        help='Set the input data type.')
    parser.add_argument('-w', '--workdir', metavar='FILE', type=str, default='.',
        help='Set the output file path, default=.')
    parser.add_argument('-n', '--number', metavar='INT', type=int, default=200000,
        help='Set the number of reads after splitting the file, default=200000')
    parser.add_argument('-o', '--out', metavar='STR', type=str, default='out',
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

    args = split_fq_help(parser).parse_args()

    split_data(args.input, args.workdir, args.out, args.no_ngs, args.number)


if __name__ == "__main__":

    main()
