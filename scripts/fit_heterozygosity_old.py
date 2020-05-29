#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import json
import matplotlib
matplotlib.use('Agg')

from scipy.signal import find_peaks_cwt

from matplotlib import pyplot as plt

plt.rcParams.update({
    "backend": "Agg",
    "figure.figsize": (4.2, 2.5),
    "figure.subplot.left": 0.2,
    "figure.subplot.right": 0.97,
    "figure.subplot.top": 0.92,
    "figure.subplot.bottom": 0.2,
    "font.family": "Arial",
    "lines.linewidth": 1,
    "axes.linewidth": 1,
    "xtick.major.width": 1,
    "ytick.major.width": 1,
    "axes.labelsize": 12,
})

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
        r.append([x[y[start: end + 1].index(max(y[start: end + 1])) + start], max(y[start: end + 1])])

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
            continue
        if line.startswith("#KmerSize:"):
            k = int(line.split(":")[1])
            continue
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if parts[0].startswith(">="):
            r[int(parts[0][2:])] = int(parts[1]) * 100.0 / n
        else:
            r[int(parts[0])] = int(parts[1]) * 100.0 / n

    return r, n, i, k


def m_distance(m1, m2, a):
    return sum([a[n] * (m1[n] - m2[n]) ** 2 for n in range(len(m1))]) / len(m1), m2


def fit_hetero(main_peak, data, database):
    main_peak_value = data[main_peak] * 1.0 / data[int(main_peak / 2)]
    near_database = [i for i in database if abs(i["data"][1] - main_peak) <= 1]
    if not near_database:
        raise Exception("main peak %s must > 24 && < 40" % main_peak)

    distances = sorted([m_distance([main_peak, main_peak_value], [i["data"][1], i["data"][2] * 1.0 / i["data"][4], i],
                                   [1 / main_peak ** 2, 10]) for i in near_database], key=lambda d: d[0])

    return distances[0][1], data[main_peak] * 1.0 / distances[0][1][2]["data"][2]


def read_database(file):
    with open(file) as fh:
        r = json.load(fh, encoding="utf-8")

    return r


def add_args(parser):
    parser.add_argument("histo", help="")
    parser.add_argument("--kingdom", required=True, help="")
    parser.add_argument("--depth", type=int, default=0, help="")
    parser.add_argument("--name", type=str, required=True)
    parser.add_argument("--database",metavar="FILE", type=str, default="/export/personal1/zhangxg/develop/survey/v1.0.0/scripts")

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

    data, n, i, k = read_histo(args.histo)

    database_dict = {
        "plant": "Atha.tsv",
        "animal": "Atha.tsv",
        "fungi": "Atha.tsv"
    }

    database = read_database("%s/%s" % (args.database, database_dict[args.kingdom]))

    if args.depth:
        main_peak = args.depth
    else:
        peaks = find_peak(data, 3)
        main_peak = sorted(peaks, key=lambda d: d[1])[-1][0]

    fit, alpha = fit_hetero(main_peak, data, database)
    h = float(fit[2]["data"][0]) * 100
    name = args.name.strip().split("/")[-1]

    print("""\
#sample\tkmer\tkmer_num\tkmer_depth\tgenome_size (bp)\theterozygosity
%s\t%s\t%s\t%s\t%s\t%s%%
""" % (name, k, format(i, ","), main_peak, format(int(i / main_peak), ','), h))

    x = []
    y = []

    for k in sorted(data):
        x.append(k)
        y.append(data[k])

    plt.plot(x, y, color='black', linewidth=1)

    plt.xlabel("Depth")
    plt.ylabel("Frequency (%)")
    plt.ylim([0, int(max(y[int(main_peak / 2.5):]) * 2)])
    plt.xlim([0, 150])
    plt.savefig("%s.kmer.png" % args.name, dpi=300)
    plt.savefig("%s.kmer.pdf" % args.name)

    plt.plot(x, y, color='black', label=name, linewidth=1)

    plt.plot(x, [i * alpha for i in fit[2]["y"]], "--", label="Atha %s%%" % h, linewidth=1)
    plt.xlabel("Depth")
    plt.ylabel("Frequency (%)")
    plt.ylim([0, int(max(y[int(main_peak / 2.5):]) * 2)])
    plt.xlim([0, 150])
    plt.legend(frameon=False)
    plt.savefig("%s.heterozygosity.png" % args.name, dpi=300)
    plt.savefig("%s.heterozygosity.pdf" % args.name)


if __name__ == "__main__":
    main()
