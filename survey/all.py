#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os.path
import argparse
import logging
import shutil

from survey.config import *
from survey.common import check_paths, mkdir, read_tsv, read_files
from dagflow import DAG, Task, ParallelTask, do_dag
from survey.kmer_denovo import run_kmer_denovo
from survey.ngs_qc import run_ngs_qc
from survey.parser import add_all_args


LOG = logging.getLogger(__name__)

__version__ = "v1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
__all__ = []


def run_report(name, assembly, kmer_length, stat_qc, base_quality, base_content, base_gc, top10_species, top10_species_png, stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, gc_depth ,out_dir):

    if assembly.lower()=="true":
        assembly = "True"
        survey = "survey_asm"
    else:
        assembly = ""
        survey = "survey"

    shell = """
export PATH={python}:$PATH
python {script}/survey_report.py config.cfg --data {stat_qc} \\
--quality {base_quality} \\
--content {base_content} \\
--gc {base_gc} \\
--species {top10_species} \\
--species_png {top10_species_png} \\
--gse {gse_txt} \\
--scope {scope_txt} \\
--gse_png {gse_png} \\
--scope_png {scope_png} \\
--kmer {stat_heter} \\
--kmer_png {heter_png} \\
--assembly {stat_genome} \\
--gc_depth {gc_depth} \\
--out {out_dir} \\
--docx {templates}/{survey}.docx \\
--html {templates}/{survey}html
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            stat_qc=stat_qc,
            base_quality=base_quality,
            base_content=base_content,
            base_gc=base_gc,
            top10_species=top10_species,
            top10_species_png=top10_species_png,
            gse_txt=gse_txt,
            scope_txt=scope_txt,
            gse_png=gse_png,
            scope_png=scope_png,
            stat_heter=stat_heter,
            heter_png=heter_png,
            stat_genome=stat_genome,
            gc_depth=gc_depth,
            survey=survey,
            out_dir=out_dir,
            templates=TEMPLATES
            )

    config ="""
[general]
project=真菌A18基因组测序分析
id=WHWLZ-201906095A-01
name={name}
species=Unknow
strain={name}
sequencer=Illumina
author=张兴国
reviewer=尹羲农
homogeneous=True
assembly={assembly}
kmer={kmer}
pollution_description= 数据存在污染对后续的分析会有影响
depth_description= 样本的基因组碱基深度主要分布在0-120x；基因平均GC含量主要分布在20-70%。基因组GC-Depth中有明显分离的聚团现象，基因组碱基深度有明显分离，说明基因组中含有其他外源污染
""".format(name=name, kmer=kmer_length, assembly=assembly)

    out_shell = open("report.sh", 'w')
    out_config = open("config.cfg", 'w')
    out_shell.write(shell)
    out_config.write(config)
    out_shell.close()
    out_config.close()


def run_survey(r1, r2, name, trim, kingdom, mode, cratio, kmer_length, kmer_depth, thread, asm, window, job_type, queue, concurrent, refresh, work_dir, out_dir, split=""):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)
    r2 = check_paths(r2)

    clean1, clean2, taxid, stat_qc, quality, content, gc, cont_tsv, cont_png = run_ngs_qc(
        r1=r1,
        r2=r2,
        name=name,
        trim=trim,
        kingdom=kingdom,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, "01_data"),
        out_dir=os.path.join(out_dir, "01_data"))

    stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, gc_depth = run_kmer_denovo(
        r1=[clean1],
        r2=[clean2],
        taxid=taxid,
        name=name,
        mode=mode,
        cratio=cratio,
        kmer_length=kmer_length,
        kmer_depth=kmer_depth,
        kingdom=kingdom,
        asm=asm,
        window=window,
        thread=thread,
        job_type=job_type,
        queue=queue,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=out_dir,
        split=split,
        platform="illumina")

    run_report(name, asm, kmer_length, stat_qc, quality, content, gc, cont_tsv,
        cont_png, stat_heter, heter_png, scope_txt, gse_txt, scope_png,
        gse_png, stat_genome, gc_depth, out_dir)


def run_all(args):

    if args.no_asm:
        asm="false"
    else:
        asm="true"

    stat_qc, quality, content, gc, poll_png, poll_tsv ,stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome = run_survey(
        r1=args.read1,
        r2=args.read2,
        name=args.name,
        trim=args.trim,
        kingdom=args.kingdom,
        mode=args.mode,
        cratio=args.cratio,
        kmer_length=args.kmer_length,
        kmer_depth=args.depth,
        thread=args.thread,
        asm=asm,
        window=args.window,
        job_type=args.job_type,
        queue="-q %s" % args.queue,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        split=args.split
    )


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

    parser = add_all_args(parser)
    args = parser.parse_args()
    run_all(args)


if __name__ == "__main__":
    main()
