#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from survey.parser import *
from survey.ngs_qc import ngs_qc
from survey.filter_cont import filter_cont
from survey.kmer_denovo import kmer_denovo
from survey.all import run_all
from survey import __version__, __email__, __author__

def add_survey_parser(parser):

    subparsers = parser.add_subparsers(
        title='command',
        dest='commands')
    subparsers.required = True

    all_parser = subparsers.add_parser("all", help="all steps")
    all_parser = add_all_args(all_parser)
    all_parser.set_defaults(func=run_all)

    ngs_qc_parser = subparsers.add_parser('ngs_qc', help="quality control")
    ngs_qc_parser = add_ngs_qc_args(ngs_qc_parser)
    ngs_qc_parser.set_defaults(func=ngs_qc)

    filter_cont_parser = subparsers.add_parser('filter_cont', help="Sequence contamination filtering")
    filter_cont_parser = add_filter_cont_args(filter_cont_parser)
    filter_cont_parser.set_defaults(func=filter_cont)


    kmer_denovo_parser = subparsers.add_parser('kmer_denovo', help="Genomic k-mer estimation and assembly")
    kmer_denovo_parser = add_kmer_denovo_args(kmer_denovo_parser)
    kmer_denovo_parser.set_defaults(func=kmer_denovo)

    return parser   


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
URL: https://github.com/zxgsy520/survey
Survey: Eukaryotic survey analysis Pipeline

version: %s
contact: %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    parser = add_survey_parser(parser)
    args = parser.parse_args()

    args.func(args)

    return parser.parse_args()


if __name__ == "__main__":
    main()


