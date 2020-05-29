#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

from scipy.signal import find_peaks_cwt

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def find_peak(data, size):

    r = []

    x = []
    y = []
    for k, v in sorted(data.items(), key=lambda d: d[0]):
        x.append(k)
        y.append(v)

    peakind = find_peaks_cwt(y, [size])

    for p in peakind:
        m = x.index(p)

        if m < 5:
            continue

        start = p - 2
        end = p + 2
        if end > len(y):
            end = len(y)
        r.append([x[y[start: end+1].index(max(y[start: end+1])) + start], max(y[start: end+1])])

    return r


def read_histo(file):

    r = {}

    for line in open(file):

        line = line.strip()

        if line.startswith("#Kmer_species_num:"):
            n = int(line.split(":")[1])
            continue

        if line.startswith("#Kmer_individual_num:"):
            i = int(line.split(":")[1])

        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if parts[0].startswith(">="):
            r[int(parts[0][2:])] = int(parts[1])*100.0/n
        else:
            r[int(parts[0])] = int(parts[1])*100.0/n

    return r, n, i


def add_args(parser):
    parser.add_argument("histo", help="")
    parser.add_argument("--depth", default=0, type=int, help="")

    return parser


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""


version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()

    data, n, i = read_histo(args.histo)

    if args.depth:
        main_peak = args.depth
    else:
        peaks = find_peak(data, 3)
        main_peak = sorted(peaks, key=lambda d: d[1])[-1][0]

    print("""\
kmer_num\t%s
kmer_depth\t%s
genome_size\t%s
""" % (i, main_peak, format(int(i / main_peak), ',')))


if __name__ == "__main__":
    main()
