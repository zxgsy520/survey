#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from survey.config import *
from survey.common import check_path, check_paths, mkdir, read_tsv, read_files, get_version
from dagflow import DAG, Task, ParallelTask, do_dag


LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v1.2.0"


def split_data(r1, r2, name, number, job_type, work_dir, out_dir):

    if len(r1)!=len(r2) and len(r2)<=1:
        read = "%s.part_*.fast*" % name
        r2 = ""
    elif len(r1)==len(r2):
        read = "%s.r1.part_*.fastq" % name
    else:
        raise Exception("The input sequencing platform is abnormal.")

    dag = DAG("split_data")
    task = Task(
        id="split_data",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{script}/splitfp.py -r1 {r1} -r2 {r2} -o {name} -n {number}
#cp {name}.* {out_dir}
""".format(
            script=SCRIPTS,
            r1=r1,
            r2=r2,
            name=name,
            number=number,
            out_dir=out_dir
        )
    )

    dag.add_task(task)
    do_dag(dag, 8, 10)

    temp = read_files(work_dir, read)
    reads = []

    if len(r1)==len(r2):
        for i in temp:
            j = i.replace(".r1.part_", ".r2.part_")
            reads.append("%s %s" % (i, j))
    else:
        reads = temp

    return reads


def create_minimap_tasks(reads, genome, platform, name, thread, job_type, work_dir, out_dir, split="split"):

    number = range(len(reads))
    tmp = []

    for i in number:
        if split in "split":
            tmp.append("--split-prefix %s_%s" % (name, i))
        else:
            tmp.append("")

    tasks = ParallelTask(
        id="minimap",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={samtools}:{minimap}:$PATH
minimap2 {{tmp}} -t {thread} {x} {genome} {{query}} | samtools view --threads {thread} -bS -T | samtools sort --threads {thread} -m {memory}G -o {name}.{{number}}.sort.bam
#cp *.sort.bam {out_dir}
""".format(
            minimap=MINIMAP_BIN,
            samtools=SAMTOOLS_BIN,
            genome=genome,
            x=SEQUENCER[platform]["minimap2"],
            thread=thread,
            memory=thread*4,
            name=name,
            out_dir=out_dir
        ),
        query=reads,
        number=number,
        tmp=tmp
    )

    return tasks, os.path.join(work_dir, "%s.*.sort.bam" % name)


def merge_bam_task(bams, name, thread, job_type, work_dir, out_dir):

    task = Task(
        id="merge_bam",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={samtools}:$PATH
samtools merge -f -c --threads {thread} {name}.sorted.bam {bams}
samtools index {name}.sorted.bam
#rm {bams}
#cp {name}.sorted.bam {out_dir}
""".format(
            samtools=SAMTOOLS_BIN,
            bams=bams,
            name=name,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.sorted.bam" % name)


def run_minimap(reads, genome, platform, name, split, thread, job_type, concurrent, refresh, work_dir, out_dir):

    option = OrderedDict()
    option["minimap2"] = {
        "version": get_version(SOFTWARE_VERSION["minimap2"]),
        "option": "%s" % SEQUENCER[platform]["minimap2"]
    }

    work_dict = {
        "minimap": "01_minimap",
        "merge": "02_merge"
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))

    dag = DAG("minimap")

    minimap_tasks, bams = create_minimap_tasks(
        reads=reads,
        genome=genome,
        platform=platform,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["minimap"]),
        out_dir=out_dir,
        split=split
    )

    merge_task, bam = merge_bam_task(
        bams=bams,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["merge"]),
        out_dir=out_dir
    )

    dag.add_task(*minimap_tasks)
    dag.add_task(merge_task)
    merge_task.set_upstream(*minimap_tasks)
    do_dag(dag, concurrent, refresh)

    return bam, option


def minimap(r1, r2, genome, name, split, platform, number, thread, job_type, concurrent, refresh, work_dir, out_dir):

    genome = check_path(genome)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)

    if r2!="":
        r2 = check_paths(r2)

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    data_work = mkdir(os.path.join(work_dir, '00_data'))

    reads = split_data(
        r1=r1,
        r2=r2,
        name=name,
        number=number,
        job_type=job_type,
        work_dir=data_work,
        out_dir=out_dir
    )

    bam, option = run_minimap(
        reads=reads,
        genome=genome,
        platform=platform,
        name=name,
        split=split,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=out_dir
    )

    options["software"] = option

    with open(os.path.join(out_dir, "minimap2.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return bam


def minimap_hlep_args(parser):

    parser.add_argument("-r1", "--read1", metavar='FILE', nargs='+', type=str, required=True,
        help="Input R1 data, if it is three-generation sequencing, enter the reads sequence.")
    parser.add_argument("-r2", "--read2", metavar='FILE', nargs='+', type=str, default='',
        help="Input R2 data, not input if it is three-generation sequencing.")
    parser.add_argument("-n", "--name", metavar="STR", type=str, default="out",
        help="Input sample name.")
    parser.add_argument("-g", "--genome", metavar='FILE', type=str, required=True,
        help="Input genome file.")
    parser.add_argument("-p", "--platform", choices=["PromethION", "GridION", "RSII", "Sequel", "illumina", "mgi"], default="PromethION",
        help="Select sequencing platform, default=PromethION")
    parser.add_argument("-s", "--split", choices=["split", ""], default="",
        help="Cache output, default=")
    parser.add_argument("--number", metavar='INT', type=int, default=200000,
        help="Cache output, default=200000")
    parser.add_argument("-t", "--thread", metavar='INT', type=int, default=4,
        help="Set the running thread, default=4")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")

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
    minimap.py Compare reads to genome.

attention:
    minimap.py --read1 ngs.r1.fq --read2 ngs.r2.fq --genome genome.fa
    minimap.py --read1 tgs.read.fq --genome genome.fa

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = minimap_hlep_args(parser).parse_args()

    minimap(args.read1, args.read2, args.genome, args.name, args.split, args.platform, args.number, args.thread, args.job_type, args.concurrent, args.refresh, args.work_dir, args.out_dir)


if __name__ == "__main__":

    main()
