#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from survey.parser import *
from survey.ngs_qc import ngs_qc
from survey.filter_contamination import filter_contamination
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

    filter_contamination_parser = subparsers.add_parser('filter_contamination', help="Sequence contamination filtering")
    filter_contamination_parser = add_filter_contamination_args(filter_contamination_parser)
    filter_contamination_parser.set_defaults(func=filter_contamination)


    kmer_denovo_parser = subparsers.add_parser('kmer_denovo', help="Genomic k-mer estimation and assembly")
    kmer_denovo_parser = add_kmer_denovo_args(kmer_denovo_parser)
    kmer_denovo_parser.set_defaults(func=kmer_denovo)

    return parser   


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Nexomics Eukaryotic Genome Analysis Pipeline

version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    parser = add_survey_parser(parser)
    args = parser.parse_args()

    args.func(args)

    return parser.parse_args()


if __name__ == "__main__":
    main()


