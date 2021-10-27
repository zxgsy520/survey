#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

from scipy.signal import find_peaks

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_depth(file):
    '''Read depth file'''

    for line in open(file, 'r'):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        line = line.split('\t')

        yield float(line[3]),  float(line[4])


def stat_hist(data, bin=50):

    data = sorted(data)
    maxdata = int(max(data))
    mindata = int(min(data))
    step = int((maxdata-mindata)/bin)
    r = {}

    if step >=10:
        step = 10
    elif step < 1:
        step = 1
    else:
        pass
    for i in range(mindata, maxdata, step):
        for j in data:
            if j >= i and  j< (i+step):
                k = i+step
                if k not in r:
                    r[k] = 0
                r[k] += 1
            if j >= (i+step):
                break
            else:
                continue
    return r


def find_peak(data):

    x = list(data.keys())
    y = list(data.values())
    peaks, pros = find_peaks(y, prominence=10, width=4)

    xpeak = []
    ypeak = []
    for i in peaks:
        xpeak.append(x[i])
        ypeak.append(y[i])

    maxy = max(ypeak)
    peakx = []
    temp = 0
    for i in range(len(ypeak)):
        if ypeak[i] >= maxy/3:
            peakx.append(xpeak[i])

    return peakx


def find_peaked(file, bin=50):

    gcs = []
    depths = []
    for gc, depth in read_depth(file):
        gcs.append(gc)
        depths.append(depth)

    gcdict = stat_hist(gcs, bin)
    depthdict = stat_hist(depths, bin)
    gcpeak = find_peak(gcdict)
    depthpeak = find_peak(depthdict)

    if len(gcpeak) ==1 and len(depthpeak)==1:
        print("样本基因组GC-Depth无明显分离聚团现象，GC含量在{gc}%左右，\
样本的Depth在{depth}X左右，样本无明显污染可以进行后续分析".format(
              gc=gcpeak[0], depth=depthpeak[0])
        )
    elif len(gcpeak) == 1:
        print("样本基因组GC-Depth有明显分离聚团现象，GC含量无明显分离现象，\
GC含量在{gc}%左右，样本Depth有明显分离，样本可能为高杂合".format(
               gc=gcpeak[0])
        )
    else:
        print("样本基因组GC-Depth有明显分离聚团现象，GC含量有明显分离现象，\
样本Depth有明显分离，样本可能含有外源的污染")

    return gcpeak, depthpeak


def add_hlep_args(parser):

    parser.add_argument("gc_depth", metavar="FILE", type=str,
        help="Input the genome GC depth.")
    parser.add_argument("-b","--bin", metavar="INT", type=int, default=50,
        help='Data segmentation bin, default=50.')

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
    find_peak.py Discover the GC and deep peaks of the genome

attention:
    find_peak.py genome.fasta stat_gc_depth.tsv >gc_depth_describe.txt
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    find_peaked(args.gc_depth, args.bin)


if __name__ == "__main__":

    main()
