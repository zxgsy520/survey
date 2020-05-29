#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
from matplotlib.pyplot import MultipleLocator

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Liu huifang",)
__email__ = "liuhuifang@grandomics.com"
__all__ = []


def get_color(species):

    sp_cl = []
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1a55FF']
    num_cl = sorted(list(set(species)))
    num_co = colors[0:len(num_cl)]

    for i in species:
        for j in range(0,len(num_cl)):
            if i == num_cl[j]:
                sp_cl.append(colors[j])
    return sp_cl, num_cl, num_co


def plot_species(args):

    sp = []
    group = []
    number = []
    f_out = open('%s.top10_species.tsv' % args.name,'w')

    f_out.write('#Species\tReads Number\tTaxonomy\n')
    with open(args.input, 'r') as f_open:
        for line in f_open:
            tem = line.strip().split('\t')
            if tem[1] == 'Unmap':
                f_out.write('%s\t%s\t%s\n' % (tem[1], tem[2], tem[0]))
            else:
                group.append(tem[0])
                sp.append(tem[1])
                number.append(float(tem[2]))

    sp_top10 = []
    group_top10 = []
    number_top10 = []
    percent_top10 = []
    for i in range(0,10):
        sp_top10.append(sp[i])
        group_top10.append(group[i])
        number_top10.append(number[i])
        per = round(100 * float(number[i])/sum(number),2)
        per = str(per) + '%'
        percent_top10.append(per)
        f_out.write('%s\t%s\t%s\n' % (sp[i], int(number[i]), group[i]))
    f_out.close()

    max_value = max(number_top10)/30
    plt.switch_backend('agg')
    fig, ax = plt.subplots()
    category_colors, x_legend, x_colors = get_color(group_top10)

    xs1 = range(len(sp_top10))
    ax.bar(xs1, number_top10, color = category_colors)
    sp_top10.insert(0, 'a')
    plt.xticks(xs1, sp_top10)

    x_major_locator = MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)
    patches = [mpatches.Patch(color=x_colors[i], label=x_legend[i]) for i in range(len(x_colors))]
    ax.legend(handles=patches)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylim(0, max(number_top10) + max(number_top10)/10)

    for i in xs1:
        ax.text(xs1[i], number_top10[i] + max_value, percent_top10[i], horizontalalignment='center', verticalalignment='center')
    fig.tight_layout()

    plt.savefig('%s.top10_species.pdf' % args.name)
    plt.savefig('%s.top10_species.png' % args.name, dpi=700)


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
        
For exmple:
        python plot_stat_species -i txt.stat_species.tsv
        
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-i', '--input',required=True,
                        help='a file,such as txt.stat_species.tsv')
    parser.add_argument('-n','--name', metavar='STR', default='out',
                        help='Output file name')

    args = parser.parse_args()

    plot_species(args)


if __name__ == "__main__":

    main()
