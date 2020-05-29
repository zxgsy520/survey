#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []

def read_depth(files):
    '''Read depth file'''

    for line in open(files, 'r'):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        line = line.split('\t')

        yield float(line[3]), float(line[4])


def beg_median(depth):

    depth = sorted(depth)
    m = int(len(depth)/2)

    return (depth[m]+depth[m-1])/2


def draw_gc_depth(files, name, minx=0, maxx=100, miny=0, maxy=500):

    gc_list = []
    depth_list = []

    for gc,depth in read_depth(files):
        gc_list.append(gc)
        depth_list.append(depth)

    maxy = beg_median(depth_list)*2

    left, width = 0.1, 0.77
    bottom, height = 0.1, 0.77
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.08]
    rect_histy = [left + width + spacing, bottom, 0.08, height]
    #Generate an image
    plt.figure(figsize=(8, 8))
#    plt.tight_layout()
    #Drawing a middle picture
    axs = plt.axes(rect_scatter)
    axs.tick_params(direction='in', top=False, right=False)
    #Set the position of the middle picture
    axx = plt.axes(rect_histx)
    axx.tick_params(direction='in', labelbottom=False)

    axy = plt.axes(rect_histy)
    axy.tick_params(direction='in', labelleft=False)

    #Draw a middle two-dimensional heat map
    #axs.hist2d(gc_list, depth_list, bins=50, range=[[minx,maxx],[miny,maxy]], cmap=plt.get_cmap('BuPu'))
    axs.scatter(gc_list, depth_list, s=20, facecolors='none', edgecolor="#107ab0", alpha=1)
    axs.set_xlabel('GC%')
    axs.set_ylabel('Sequencing Depth(X)')
    axs.set_xlim((minx, maxx))
    axs.set_ylim((miny, maxy))
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)

    axx.hist(gc_list, bins=50, range=(minx, maxx), edgecolor='white')
    axx.set_xlim(minx, maxx)
    axx.yaxis.set_visible(False)
    #axx.xaxis.set_visible(False)
    axx.spines['top'].set_visible(False)
    axx.spines['right'].set_visible(False)
    axx.spines['left'].set_visible(False)

    axy.hist(depth_list, bins=50, range=(miny, maxy), edgecolor='white', orientation='horizontal')
    axy.set_ylim(miny, maxy)
    #axy.yaxis.set_visible(False)
    axy.xaxis.set_visible(False)
    axy.spines['top'].set_visible(False)
    axy.spines['right'].set_visible(False)
    axy.spines['bottom'].set_visible(False)

    plt.savefig('%s.gc_depth.png' % name, dpi=700)
    plt.savefig("%s.gc_depth.pdf" % name)


def gc_depth_help(parser):

    parser.add_argument('-gcd', '--gc_depth', metavar='FILE', type=str, required=True,
        help='Input the coverage depth statistics file.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default='out',
        help='Output file name.')

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
    plot_depth_gc.py  Draw a gc depth map

attention:
    plot_depth_gc.py -gcd *.stat_gc_depth.tsv -n A18

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = gc_depth_help(parser).parse_args()

    draw_gc_depth(args.gc_depth, args.name)


if __name__ == "__main__":

    main()
