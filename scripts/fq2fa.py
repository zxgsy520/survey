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
            if len(seq)==5:
                yield seq[0], seq[1]
            seq = ''
            seq = "%s\n" % line
        else:
            seq += "%s\n" % line
    
    seq = seq.split('\n')
    if len(seq)==5:
        yield seq[0], seq[1]


def run_fq2fa(file, number):
    
    n = 0
    
    if number!='all':
        number=int(number) 
   
    for seq_id,seq in read_fastq(file):
        if n==number:
             break
        print('>%s\n%s' % (seq_id,seq))
        n +=1
    

def add_args(parser):

    parser.add_argument('fastq', help='Input fastq file.')
    parser.add_argument('-n', '--number', default='all',
        help='Output reads number, default=all')
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

    fq2fa.py -- Integrated database

attention:
    fq2fa.py txt.fastq >txt.fasta
''')
    args = add_args(parser).parse_args()

    run_fq2fa(args.fastq, args.number)


if __name__ == "__main__":
    main()
