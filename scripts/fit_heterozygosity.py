#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from scipy.signal import find_peaks_cwt

LOG = logging.getLogger(__name__)

__version__ = "2.1.2"
__author__ = ("Junpeng Fan, Xingguo Zhang",)
__email__ = "jpfan@whu.edu.cn, invicoun@foxmail.com"
__all__ = []


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


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_histo(file, model="species", minn=10):

    r = {}
    n = 0
    total = 0
    total_number = 0
    total_indiv = 0

    for line in read_tsv(file):
        if line[0].startswith(">="):
            freq = int(line[0].strip(">="))
        else:
            freq = int(line[0])
        n += 1
        number = int(line[1])
        indiv = number*freq

        if n <= 255:
           total_number += number
           total_indiv += indiv
           r[freq] = [number, indiv]
        total += indiv

        if freq >=4000 and number <=minn:
           break           

    data = {}
    if model== "species":
        for i in r:
            data[i] = r[i][0]*100.0/total_number
    else:
        for i in r:
            data[i] = r[i][1]*100.0/total_indiv

    return data, total_number, total


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


def judge_main_peak(data):

    data = sorted(data, key=lambda d: d[1])
    main_peak = data[-1][0]
    if len(data)>=2:
        if data[-2][1] >= data[-1][1]/2:
            b = data[-2][0]/data[-1][0]
            if b>=1.8 and b<=2.8:
                main_peak = data[-2][0]

    return main_peak


def fit_heterozygosity(histo, model, database, kingdom, name, depth=0, length=17):

    data, n, i= read_histo(histo, model)
    database_dict = {
        "plant": "Atha.tsv",
        "animal": "Atha.tsv",
        "fungi": "Atha.tsv"
    }
    database = read_database("%s/%s" % (database, database_dict[kingdom]))

    if depth:
        main_peak = depth
    else:
        peaks = find_peak(data, 3)
        main_peak = judge_main_peak(peaks)

    fit, alpha = fit_hetero(main_peak, data, database)
    h = float(fit[2]["data"][0]) * 100
    name = name.strip().split("/")[-1]

    print("""\
#sample\tkmer\tkmer_num\tkmer_depth\tgenome_size (bp)\theterozygosity
{0}\t{1}\t{2:,}\t{3}\t{4:,.2f}\t{5:.2f}%
""".format(name, length, i, main_peak, int(i*1.0/main_peak), h))

    x = []
    y = []

    for k in sorted(data):
        x.append(k)
        y.append(data[k])

    plt.plot(x, y, color='black', linewidth=1)

    plt.xlabel("Depth")
    plt.ylabel("Frequency (%)")
    ymax = int(max(y[int(main_peak/2.5):])*2)
    if ymax <= 1:
        ymax = 1
    plt.ylim([0, ymax])
    plt.xlim([0, 150])
    plt.savefig("%s.kmer.png" % name, dpi=300)
    plt.savefig("%s.kmer.pdf" % name)

    plt.plot(x, y, color='black', label=name, linewidth=1)

    plt.plot(x, [i * alpha for i in fit[2]["y"]], "--", label="Atha {0:.2f}%".format(h), linewidth=1)
    plt.xlabel("Depth")
    plt.ylabel("Frequency (%)")
    plt.ylim([0, ymax])
    plt.xlim([0, 150])
    plt.legend(frameon=False)
    plt.savefig("%s.heterozygosity.png" % name, dpi=300)
    plt.savefig("%s.heterozygosity.pdf" % name)


def add_args(parser):

    parser.add_argument("histo",
        help="Input kmer frequency distribution.")
    parser.add_argument("-m", "--model", choices=["species", "individual"], default="species",
        help="Choose model, default=species.")
    parser.add_argument("-k", "--kingdom", choices=["plant", "animal", "fungi"], default="fungi",
        help="Choose kingdom, default=fungi.")
    parser.add_argument("-d", "--depth", metavar="INT", type=int, default=0,
        help="Input the location of the main peak of kmer, default=0")
    parser.add_argument("-n", "--name", metavar="STR", type=str, required=True,
        help="Input sample name.")
    parser.add_argument("-l", "--length", metavar="INT", type=int, default=17,
        help="Input the length of kmer, default=17")
    parser.add_argument("--database",metavar="FILE", type=str, default="/export/personal1/zhangxg/develop/survey/v1.0.0/scripts",
        help="Input the model position for heterozygosity simulation.")

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

    fit_heterozygosity(args.histo, args.model, args.database, args.kingdom, args.name, args.depth, args.length)


if __name__ == "__main__":
    main()
