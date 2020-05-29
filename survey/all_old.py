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
from survey.bwa_mem import run_bwa_mem
#from survey.kmer_denovo import kmerfreq_task, kmer_denovo_tasks, run_gc_depth
from survey.ngs_qc import ngs_qc_tasks
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


def run_survey(r1, r2, name, trim, kingdom, kmer_length, sample_depth, thread, asm, window, job_type, queue, concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)
    r2 = check_paths(r2)

    dag = DAG("survey_qc")
    merge_task, qc_task, cont_task, result_task, clean1, clean2, quality, content, gc, stat_qc, poll_png, poll_tsv = ngs_qc_tasks(
        name=name,
        r1=r1,
        r2=r2,
        trim=trim,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    
    data_work = mkdir(os.path.join(work_dir, "choose_data"))
    freq_task1, histo1, kmer_stat, estimate1 = kmerfreq_task(
        r1=clean1,
        r2=clean2,
        name=name,
        kmer_length=kmer_length,
        thread=thread,
        job_type=job_type,
        work_dir=data_work,
        out_dir=data_work)

    dag.add_task(merge_task)
    dag.add_task(qc_task)
    qc_task.set_upstream(merge_task)
    dag.add_task(cont_task)
    dag.add_task(result_task)
    dag.add_task(freq_task1)
    freq_task1.set_upstream(qc_task)
    cont_task.set_upstream(qc_task)
    result_task.set_upstream(qc_task)
       
    do_dag(dag, concurrent, refresh)

    for line in read_tsv(kmer_stat):
        if line[0]=="kmer_depth":
             kmer_depth = int(line[1])

    if sample_depth>kmer_depth:
        LOG.debug('The amount of sequencing data may be insufficient. Sequencing depth is only %s X' % kmer_depth)
        sample_depth = kmer_depth
    proportion=sample_depth*1.0/kmer_depth
    
    dag = DAG("survey")
    choose_task, freq_task, heter_task, jellyfish_task, gse_scope_task, denovo_task, stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, genome, ngs_list = kmer_denovo_tasks(
        r1=clean1,
        r2=clean2,
        name=name,
        kmer_length=kmer_length,
        proportion=proportion,
        kingdom=kingdom,
        thread=thread,
        job_type=job_type,
        queue=queue,
        work_dir=work_dir,
        out_dir=out_dir)
    if asm=="true":
        dag.add_task(denovo_task)
    else:
        genome = "false"
        stat_genome= "false"
        ngs_list = "false"

    dag.add_task(choose_task)
    dag.add_task(freq_task)
    dag.add_task(heter_task)
    freq_task.set_upstream(choose_task)
    dag.add_task(jellyfish_task)
    jellyfish_task.set_upstream(choose_task)
    dag.add_task(gse_scope_task)
    heter_task.set_upstream(freq_task)
    gse_scope_task.set_upstream(jellyfish_task)
    do_dag(dag, concurrent, refresh)

    if ngs_list == "false":
        print("Genomics are not assembled")
        gc_depth_png = heter_png
    else:
        depth_work = mkdir(os.path.join(work_dir, "05_GC-depth"))
        depth_out = mkdir(os.path.join(out_dir, "05_GC-depth"))
        gc_depth_png = run_gc_depth(
            genome=genome,
            fastq_list=ngs_list,
            name=name,
            window=window,
            thread=thread,
            job_type=job_type,
            concurrent=concurrent,
            refresh=refresh,
            work_dir=depth_work,
            out_dir=depth_out)

    run_report(name, asm, kmer_length, stat_qc, quality, content, gc, poll_tsv, poll_png, stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, gc_depth_png ,out_dir)


    return stat_qc, quality, content, gc, poll_png, poll_tsv ,stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome


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
        kmer_length=args.kmer_length,
        sample_depth=args.depth,
        thread=args.thread,
        asm=asm,
        window=args.window,
        job_type=args.job_type,
        queue="-q %s" % args.queue,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
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
