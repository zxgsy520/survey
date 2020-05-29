#!/usr/bin/env python
# -*- coding: utf-8 -*-

from survey.config import *

__all__ = ["add_ngs_qc_args", "add_filter_contamination_args", "add_kmer_denovo_args", "add_all_args"]

def add_workflow_args(parser):
    """
    add workflow arguments to parser
    :param parser: argparse object
    :return: parser
    """

    workflow_group = parser.add_argument_group(title="Workflow arguments", )
    workflow_group.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10).")
    workflow_group.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30).")
    workflow_group.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local).")
    workflow_group.add_argument("--work_dir", metavar="DIR", type=str, default=".",
        help="Work directory (default: current directory).")
    workflow_group.add_argument("--out_dir", metavar="DIR", type=str, default=".",
        help="Output directory (default: current directory).")

    return parser


def add_ngs_qc_args(parser):

    parser.add_argument("-r1", "--read1", metavar="FILE", nargs='+', type=str, required=True, 
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-n", "--name", metavar="FILE", required=True,
        help="Sample name.")
    parser.add_argument("--trim", metavar="INT", type=int, default=5,
        help="Set trim length, default=5")
    parser.add_argument("--kingdom" , choices=["plant", "animal", "fungi"], default="fungi",
        help="Choose kingdom.")
    parser.add_argument("--thread", type=int, default=1, 
        help="threads used in genome mapping")
    parser = add_workflow_args(parser)

    return parser


def add_filter_contamination_args(parser):
    
    parser.add_argument("-r1", "--read1", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", nargs='+', type=str, required=True,
        help="Second generation sequencing of reads r1.")
    parser.add_argument("-n", "--name", metavar="FILE", type=str, default="ont",
        help="Sample name.")
    parser.add_argument("--taxid", metavar="FILE", type=str, required=True,
        help="Input the pollution ratio of the sample and the taxid of the pollutant species.")
    parser.add_argument("--kingdom" , choices=["plant", "animal", "fungi"], default="fungi",
        help="Choose kingdom.")
    parser.add_argument("-kl", "--kmer_length",  metavar="INT", type=int, default=17,
        help="Set the length of kmer.")
    parser.add_argument("-d", "--depth", metavar="INT", type=int, default=35,
        help="Set the selected kmer depth.")
    parser.add_argument("-s", "--split", choices=["split", ""], default="",
        help="Cache output, default=")
    parser.add_argument("-m", "--mode", choices=["fast", "general", "strict"], default="fast",
        help="Choose the mode of operation, default=fast.")
    parser.add_argument("-c", "--cratio", metavar="INT", type=int, default=10,
        help="Input pollution rate, valid when mode = strict, default=10.")
    parser.add_argument("-t", "--thread", type=int, default=1,
        help="Set the number of threads.")
    parser = add_workflow_args(parser)

    return parser


def add_kmer_denovo_args(parser):

    parser.add_argument("-w", "--window", metavar="INT", type=int, default=5000,
        help="Set the window size when calculating the GC depth.")
    parser.add_argument("--no_asm", action="store_true",
        help="Choose whether to assemble or not. The default is assembly.")
    parser.add_argument("-q", "--queue", metavar="STR", type=str, default="fat.q",
        help="Queue to be used for assembly tasks,default=fat.q .")
    parser = add_filter_contamination_args(parser)

    return parser


def add_all_args(parser):
    
    parser.add_argument("-kl", "--kmer_length",  metavar="INT", type=int, default=17,
        help="Set the length of kmer.")
    parser.add_argument("-d", "--depth", metavar="INT", type=int, default=35,
        help="Set the selected kmer depth.")
    parser.add_argument("-m", "--mode", choices=["fast", "general", "strict"], default="fast",
        help="Choose the mode of operation, default=fast.")
    parser.add_argument("-c", "--cratio", metavar="INT", type=int, default=10,
        help="Input pollution rate, valid when mode = strict, default=10.")
    parser.add_argument("-w", "--window", metavar="INT", type=int, default=5000,
        help="Set the window size when calculating the GC depth.")
    parser.add_argument("-s", "--split", choices=["split", ""], default="",
        help="Cache output, default=")
    parser.add_argument("--no_asm", action="store_true",
        help="Choose whether to assemble or not. The default is assembly.")
    parser.add_argument("-q", "--queue", metavar="STR", type=str, default="fat.q",
        help="Queue to be used for assembly tasks,default=fat.q .")

    parser = add_ngs_qc_args(parser)

    return parser
