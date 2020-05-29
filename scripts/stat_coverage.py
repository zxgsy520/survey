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


def read_depth(file):

    for line in open(file, 'r'):
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        line = line.split('\t')

        yield float(line[2])


def stat_depth(depth_list, file):

    depth_list = depth_list.strip().split(',')
    depth_dict = {}
    depth = 0
    n = 0

    for i in read_depth(file):
        n += 1
        depth += i
        for j in depth_list:
            j = float(j)
            if i < j:
                continue
            if j not in depth_dict:
                depth_dict[j] = 0
            depth_dict[j] += 1

    return n, depth, depth_dict


def stat_out_coverage(n, depth, depth_dict, out_file):

    out = open(out_file, 'w')
    out.write('Depth\tBase number\tCoverage ratio(%)\n')
    #http://blog.sina.com.cn/s/blog_a5d4da69010169rn.html

    for name in sorted(depth_dict):
        out.write('{:.0f}\t{:,}\t{:.2f}\n'.format(name, depth_dict[name], depth_dict[name]*100.0/n))
    out.write('\nCoverage depth(X)\t{:,.2f}\n'.format(depth*1.0/n))
    out.close()


def add_help(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Set the input file.')
    parser.add_argument('-d', '--depth', metavar='LIST', type=str, default='1,5,10,20',
        help='Set the statistics coverage depth list,default=1,5,10,20')
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
    stat_coverage.py  Count the number of bases with different coverage depths.
note:
    The format of the input file:
        The first behavior contig id, the second generation behavior base position, and the third behavior coverage depth.
        Bases covering a depth of 0 also require output.

attention:
    stat_coverage.py -i test.depth
    stat_coverage.py -i test.depth -d 1,5,10,20,50 -o name

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()
    n, depth, depth_dict = stat_depth(args.depth, args.input)

    stat_out_coverage(n, depth, depth_dict, args.out)


if __name__ == "__main__":

    main()
