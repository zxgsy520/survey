#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt


LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_data(file):

    r = {}
    tag = ""
    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")

        if line.startswith(">>"):
            if len(parts) == 2:
                tag = parts[0][2:]
                r[tag] = []
        else:
            r[tag].append(parts)

    return r


def plot_base_quality(data, ax, lr):

    x0 = [i[0] for i in data]
    x = range(len(x0))
    y1 = [float(i[1]) for i in data]  # Mean
    y2 = [float(i[2]) for i in data]  # Median
    y3 = [float(i[3]) for i in data]  # Lower Quartile
    y4 = [float(i[4]) for i in data]  # Upper Quartile
    y5 = [float(i[5]) for i in data]  # 10th Percentile
    y6 = [float(i[6]) for i in data]  # 90th Percentile

    x_num = len(x)
    ax.set_title('Quality scores across all bases in %s' % lr, fontsize=10)
    ax.set_xlabel('Position in read (bp)', fontsize=8)
    ax.set_ylim([0, 38])
    ax.set_xlim([-0.5, x_num - 0.5])

    if lr == "read1":
        ax.set_ylabel('Quality score', fontsize=8)
        ax.set_yticks([i for i in range(0, 39, 2)])
        ax.set_yticklabels([i for i in range(0, 39, 2)], fontsize=5)
    else:
        ax.set_yticks([])

    single_ticks = []
    double_ticks = []

    for i in x:
        if "-" not in x0[i]:
            single_ticks.append([i, x0[i]])
        else:
            double_ticks.append([i, x0[i]])

    ax.set_xticks([i[0] for i in single_ticks] + [i[0] for n, i in enumerate(double_ticks) if n % 3 == 1])
    ax.set_xticklabels([i[1] for i in single_ticks] + [i[1] for n, i in enumerate(double_ticks) if n % 3 == 1], fontsize=5)

    ax.axhspan(ymin=0, ymax=20, color="red", alpha=0.6, linewidth=0)
    ax.axhspan(ymin=20, ymax=28, color="orange", alpha=0.6, linewidth=0)
    ax.axhspan(ymin=28, ymax=38, color="green", alpha=0.6, linewidth=0)
    for i in range(x_num):
        if i % 2 == 0:
            alpha = 0.4
        else:
            alpha = 0.3

        ax.axvspan(xmin=i-0.5, xmax=i+0.5, alpha=alpha, color="white", linewidth=0)

    ax.vlines(x, ymin=y5, ymax=y6, linewidth=0.5)
    ax.hlines(y5, xmin=[i-0.3 for i in x], xmax=[i+0.3 for i in x], linewidth=0.5)  #
    ax.hlines(y6, xmin=[i - 0.3 for i in x], xmax=[i + 0.3 for i in x], linewidth=0.5)
    ax.hlines(y2, xmin=[i - 0.3 for i in x], xmax=[i + 0.3 for i in x], linewidth=0.5, zorder=3)
    for i in range(x_num):
        ax.add_patch(plt.Rectangle([i-0.3, y3[i]], width=0.6, height=y4[i]-y3[i], facecolor="yellow", edgecolor="black", linewidth=0.5, zorder=2))

    ax.plot(x, y1, linewidth=0.5, )

    return ax


def plot_base_content(data, ax, lr):

    x0 = [i[0] for i in data]
    x = range(len(x0))
    y1 = [float(i[1]) for i in data]  # G
    y2 = [float(i[2]) for i in data]  # A
    y3 = [float(i[3]) for i in data]  # T
    y4 = [float(i[4]) for i in data]  # C

    x_num = len(x)
    ax.set_title('Sequence content across all bases in %s' % lr, fontsize=10)
    ax.set_xlabel('Position in read (bp)', fontsize=8)
    ax.set_ylim([0, 50])
    ax.set_xlim([-0.5, x_num - 0.5])

    if lr == "read1":
        ax.set_ylabel('Percent (%)', fontsize=8)
        ax.set_yticks([i for i in range(0, 51, 5)])
        ax.set_yticklabels([i for i in range(0, 51, 5)], fontsize=5)
    else:
        ax.set_yticks([])

    single_ticks = []
    double_ticks = []

    for i in x:
        if "-" not in x0[i]:
            single_ticks.append([i, x0[i]])
        else:
            double_ticks.append([i, x0[i]])

    ax.set_xticks([i[0] for i in single_ticks] + [i[0] for n, i in enumerate(double_ticks) if n % 3 == 1])
    ax.set_xticklabels([i[1] for i in single_ticks] + [i[1] for n, i in enumerate(double_ticks) if n % 3 == 1], fontsize=5)

    for i in range(x_num):
        if i % 2 == 0:
            alpha = 0.3
        else:
            alpha = 0.1

        ax.axvspan(xmin=i-0.5, xmax=i+0.5, alpha=alpha, color="grey", linewidth=0)

    ax.plot(x, y2, linewidth=1, label="A")
    ax.plot(x, y3, linewidth=1, label="T")
    ax.plot(x, y1, linewidth=1, label="G")
    ax.plot(x, y4, linewidth=1, label="C")

    ax.legend(fontsize=8,)

    return ax


def plot_gc_distribution(data, ax, lr):

    x = [int(i[0]) for i in data]

    y0 = [float(i[1]) for i in data]  # count
    sum_y = sum(y0)
    y = [i*100.0/sum_y for i in y0]

    ax.set_title('GC distribution in %s' % lr, fontsize=10)
    ax.set_xlabel('Mean GC content (%)', fontsize=8)

    ax.set_xlim([0, 100])

    ax.grid(axis="x", linewidth=0.5, linestyle='--',)

    if lr == "read1":
        ax.set_ylabel('Percent (%)', fontsize=8)
        ax.tick_params(axis='both', which='both', labelsize=5)
    else:
        ax.set_yticks([])

    ax.set_xticks([i for i in range(0, 100, 5)])
    ax.set_xticklabels([i for i in range(0, 100, 5)], fontsize=5)

    ax.plot(x, y, linewidth=0.5, color="red")

    return ax


def add_args(parser):
    parser.add_argument("-1", "--r1", required=True, help="")
    parser.add_argument("-2", "--r2", required=True, help="")
    parser.add_argument("--name", required=True, help="")
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

    r1_data = read_data(args.r1)
    r2_data = read_data(args.r1)

    fig, ax = plt.subplots(1, 2, figsize=(8.3, 4))
    plt.subplots_adjust(top=0.9, left=0.06, right=0.98, bottom=0.15, wspace=0.02)
    plot_base_quality(r1_data["Per base sequence quality"], ax[0], lr="read1")
    plot_base_quality(r2_data["Per base sequence quality"], ax[1], lr="read2")
    plt.savefig("%s.base_quality.pdf" % args.name, dpi=700)
    plt.savefig("%s.base_quality.png" % args.name, dpi=700)

    fig, ax = plt.subplots(1, 2, figsize=(8.3, 4))
    plt.subplots_adjust(top=0.9, left=0.06, right=0.98, bottom=0.15, wspace=0.02)
    plot_base_content(r1_data["Per base sequence content"], ax[0], lr="read1")
    plot_base_content(r2_data["Per base sequence content"], ax[1], lr="read2")
    plt.savefig("%s.base_content.pdf" % args.name, dpi=700)
    plt.savefig("%s.base_content.png" % args.name, dpi=700)

    fig, ax = plt.subplots(1, 2, figsize=(8.3, 4))
    plt.subplots_adjust(top=0.9, left=0.06, right=0.98, bottom=0.15, wspace=0.02)
    plot_gc_distribution(r1_data["Per sequence GC content"], ax[0], lr="read1")
    plot_gc_distribution(r1_data["Per sequence GC content"], ax[1], lr="read2")
    plt.yticks(fontsize=5)
    plt.savefig("%s.base_gc.pdf" % args.name, dpi=700)
    plt.savefig("%s.base_gc.png" % args.name, dpi=700)


if __name__ == "__main__":
    main()

