#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import argparse
import logging

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def json2tsv(file, name):

    """
    sample  total_reads  total_bases    clean_reads    clean_bases   Q20 rate (%)    Q30 rate (%)   GC

    :param json:
    :return:
    """
    j = json.load(open(file))

    return [name,
            j["summary"]["before_filtering"]["total_reads"],
            j["summary"]["before_filtering"]["total_bases"],
            j["summary"]["after_filtering"]["total_reads"],
            j["summary"]["after_filtering"]["total_bases"],
            j["summary"]["after_filtering"]["q20_rate"]*100,
            j["summary"]["after_filtering"]["q30_rate"]*100,
            j["summary"]["after_filtering"]["gc_content"]*100
            ]


def add_args(parser):
    parser.add_argument("fastp", help="")
    parser.add_argument("--name", required=True)

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
    name = args.name.strip().split("/")[-1]
    print("""\
#sample\ttotal reads\ttotal bases\tclean reads\tclean bases\tQ20 rate (%)\tQ30 rate (%)\tGC (%)
{}\t{:,}\t{:,}\t{:,}\t{:,}\t{:.2f}\t{:.2f}\t{:.2f}\
""".format(*json2tsv(args.fastp, name)))

if __name__ == "__main__":
    main()

