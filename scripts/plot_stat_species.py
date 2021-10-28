#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("Reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_species(file):

    print("#Species\tReads Number\tTaxonomy")

    n = 0
    r = {}
    species = []
    total = 0
    for line in read_tsv(file, "\t"):
        total += float(line[2])

        if n >=10:
            continue
        print("%s\t%s\t%s" % (line[1], line[2], line[0]))
        if line[1] == 'Unmap':
            continue
        n += 1
        if line[0] not in r:
            r[line[0]] = [[], []]
        species.append(line[1])
        r[line[0]][0].append(n)
        r[line[0]][1].append(float(line[2]))

    if "Bacteria" in r:
        print("#!基于污染分析分析结果显示，数据存在明显细菌的污染对后续分析可能存在一定的影响")
    else:
        print("#!基于污染分析分析结果显示，数据无明显细菌的污染可以进行后续分析")

    return r, species, total


def plot_species(r, species, total, prefix="out"):

    COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
              '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1a55FF']

    fig, ax = plt.subplots()
    n = 0
    maxy = 0
    rects = []
    for i in r:
        rects += ax.bar(r[i][0], r[i][1], color=COLORS[n], label=i)
        #rects.append(rects)
        n += 1
        if max(r[i][1]) >= maxy:
            maxy = max(r[i][1])

    plt.xticks(range(1, len(species)+1), species)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylim(0, maxy+maxy/10)

    for rect in rects:
        height = rect.get_height()
        ax.annotate('{:.2f}%'.format(height*100.0/total),
                    xy=(rect.get_x() + rect.get_width()/2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

    plt.legend(loc='upper right',
               frameon=False,
               fontsize="large",
               labelspacing=0.3,
               handlelength=0.6,
               handleheight=0.6,
               handletextpad=0.5)

    fig.tight_layout()
    plt.savefig('%s.top10_species.pdf' % prefix)
    plt.savefig('%s.top10_species.png' % prefix, dpi=700)


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input species statistics file(stat_species.tsv)')
    parser.add_argument('-p', '--prefix', metavar='FILE', type=str, default="out",
        help='Output file prefix, default=out.')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
attention:
    plot_stat_species.py stat_species.tsv -p name >name.top10_species.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()
    r, species, total = read_species(args.input)

    plot_species(r, species, total, args.prefix)
    

if __name__ == "__main__":

    main()
