#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

from survey.config import *
from survey.common import check_paths, mkdir, read_tsv, read_files
from dagflow import DAG, Task, ParallelTask, do_dag
from survey import __author__, __email__, __version__

LOG = logging.getLogger(__name__)


def split_ngs_task(fastq_list, name, number, data_type, job_type, work_dir, out_dir):

    task = Task(
        id="split_ngs",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{script}/splitfp.py -i {fastq_list} -w {out_dir} -o {name} -n {number} {type}
""".format(
            script=SCRIPTS,
            fastq_list=fastq_list,
            name=name,
            number=number,
            type=data_type,
            out_dir=out_dir
        )
    )
    r1_name = '%s.r1.part_*.fastq' % name
    r2_name = '%s.r2.part_*.fastq' % name

    return task, out_dir, r1_name, r2_name


def bwa_index_task(genome, name, job_type, work_dir, out_dir):

    task = Task(
        id="bwa_index",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={bwa}:$PATH
ln -sf {genome} {out_dir}/{name}.fasta
bwa index {out_dir}/{name}.fasta
""".format(
            bwa=BWA_BIN,
            genome=genome,
            name=name,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.fasta" % name)


def bwa_mem_tasks(fq_path, r1_name, r2_name, genome, name, thread, job_type, work_dir, out_dir):
     
    r1_list = read_files(fq_path, r1_name)
    r2_list = []
    number = []

    n=1
    for i in r1_list:
        r2_list.append(i.replace('.r1.part_','.r2.part_'))
        number.append(n)
        n+=1
    
    tasks = ParallelTask(
        id="bwa_mem",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={bwa}:{samtools}:$PATH

bwa mem -t {thread} {genome} {{r1}} {{r2}}| samtools view --threads {thread} -bS | samtools sort --threads {thread} -m {cpu}G -o {out_dir}/{name}.{{number}}.sort.bam
""".format(bwa=BWA_BIN,
           samtools=SAMTOOLS_BIN,
           genome=genome,
           thread=thread,
           name=name,
           cpu=int(thread)*2,
           out_dir=out_dir),
        r1=r1_list,
        r2=r2_list,
        number=number)

    return tasks, os.path.join(out_dir, "%s.*.sort.bam" % name)


def bwa_merge_bam_task(sort_bams, name, thread, job_type, work_dir, out_dir):

    task = Task(
        id="merge_bwa_bam",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={samtools}:$PATH
ls {sort_bams} >{out_dir}/bam.list
samtools merge -f -c --threads {thread} -b {out_dir}/bam.list {out_dir}/{name}.sorted.bam
samtools index {out_dir}/{name}.sorted.bam
#rm -rf {sort_bams}
""".format(
            samtools=SAMTOOLS_BIN,
            sort_bams=sort_bams,
            name=name,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.sorted.bam" % name)

    
def run_bwa_mem(fq_path, r1_name, r2_name, genome, name, thread, job_type, work_dir, out_dir):

    work_dict = {
        "index": "01_index",
        "bwa": "02_bwa",
        "merge": "03_merge"
    }
    
    index_work = mkdir(os.path.join(work_dir, work_dict["index"]))
    index_out = mkdir(os.path.join(out_dir, work_dict["index"]))
    index_task, genome = bwa_index_task(
        genome=genome,
        name=name,
        job_type=job_type,
        work_dir=index_work,
        out_dir = index_out
    )

    bwa_work = mkdir(os.path.join(work_dir, work_dict["bwa"]))
    bwa_out = mkdir(os.path.join(out_dir, work_dict["bwa"]))
    bwa_tasks, sort_bams = bwa_mem_tasks(
        fq_path=fq_path,
        r1_name=r1_name,
        r2_name=r2_name,
        genome=genome,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=bwa_work,
        out_dir=bwa_out
    )

    merge_work = mkdir(os.path.join(work_dir, work_dict["merge"]))
    merge_out = mkdir(os.path.join(out_dir, work_dict["merge"]))
    merge_task, sorted_bam = bwa_merge_bam_task(
        sort_bams=sort_bams,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=merge_work,
        out_dir=merge_out
    )

    return index_task, bwa_tasks, merge_task, sorted_bam, genome


def bwa_mem(fastq_list, genome, name, number, data_type, thread, job_type, concurrent, refresh, work_dir, out_dir):

    genome, fastq_list = check_paths([genome, fastq_list])
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    dag = DAG("split_ngs")
    split_work = mkdir(os.path.join(work_dir, "00_data"))
    split_out = mkdir(os.path.join(out_dir, "00_data"))
   
    splitfp_task, fq_path, r1_name, r2_name  = split_ngs_task(
        fastq_list=fastq_list,
        name=name,
        number=number,
        data_type=data_type,
        job_type=job_type,
        work_dir=split_work,
        out_dir=split_out
    )
    dag.add_task(splitfp_task)
    do_dag(dag, concurrent, refresh)

    dag = DAG("bwa_mem")
    index_task, bwa_tasks, merge_task, sorted_bam, genome = run_bwa_mem(
        fq_path=fq_path,
        r1_name=r1_name,
        r2_name=r2_name,
        genome=genome,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )

    dag.add_task(index_task)
    dag.add_task(*bwa_tasks)
    dag.add_task(merge_task)
    index_task.set_downstream(*bwa_tasks)
    merge_task.set_upstream(*bwa_tasks)

    do_dag(dag, concurrent, refresh)

    return sorted_bam, genome
