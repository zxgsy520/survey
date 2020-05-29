#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from collections import OrderedDict
from survey.config import *
from dagflow import DAG, Task, ParallelTask, do_dag
from survey.common import check_path, check_paths, mkdir, read_tsv, read_files, get_version
from survey.parser import add_filter_contamination_args

LOG = logging.getLogger(__name__)

__author__ = ("Xingguo Zhang",)
__email__ = "1131978210qq.com"
__version__ = "v1.0.0"


def merge_data_task(name, r1, r2, job_type, work_dir, out_dir):

    if r1[0].endswith(".gz") or r2[0].endswith(".gz"):
        suffix = ".gz"
        tools = "zcat"
    else:
        suffix = ""
        tools = "cat"

    if len(r1)<=1 and suffix == "":
        job_type = "local"
        run = """
ln -s {r1} {name}.clean.r1.fastq
ln -s {r2} {name}.clean.r2.fastq
""".format(r1=" ".join(r1), r2=" ".join(r2), name=name)
    else:
        run = """
{tools} {r1} >{name}.clean.r1.fastq
{tools} {r2} >{name}.clean.r2.fastq
""".format(tools=tools, r1=" ".join(r1), r2=" ".join(r2), name=name)

    task = Task(
        id="data_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{run}
""".format(run=run)
    )

    return task, os.path.join(work_dir, "%s.clean.r1.fastq" % name), os.path.join(work_dir, "%s.clean.r2.fastq" % name)


def create_kmerfreq_task(r1, r2, name, kmer_length, thread, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["kmerfreq"] = {
        "version": get_version(SOFTWARE_VERSION["kmerfreq"]),
        "option": "-q 33 -m 0 -k %s" % kmer_length
    }

    if kmer_length >=17:
        kmer_length=17

    task = Task(
        id="kmerfreq",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={python}:$PATH
ls {r1} >{name}.data
ls {r2} >>{name}.data
{kmerfreq} -k {kmer_length} -t {thread} -p {name} -q 33 -m 0 {name}.data > {name}.kmer.count
cp {name}.freq.stat {name}.kmerfreq.stat
python {script}/kmerfreq_stat.py {name}.freq.stat >{name}.kmer.stat
#cp {name}.kmerfreq.stat {out_dir}
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            kmerfreq=KMERFREQ,
            kmer_length=kmer_length,
            name=name,
            r1=r1,
            r2=r2,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.kmer.stat" % name), option


def stat_kmer_depth(r1, r2, name, kmer_length, thread, job_type, concurrent, refresh, work_dir, out_dir):

    dag = DAG("survey_data")
    data_task, r1, r2 = merge_data_task(
        name=name,
        r1=r1,
        r2=r2,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)

    kmerfreq_task, kmer_stat, option= create_kmerfreq_task(
        r1=r1,
        r2=r2,
        name=name,
        kmer_length=kmer_length,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)

    dag.add_task(data_task)
    dag.add_task(kmerfreq_task)
    kmerfreq_task.set_upstream(data_task)
    do_dag(dag, concurrent, refresh)

    return r1, r2, kmer_stat, option


def choose_data_task(r1, r2, name, kmer_stat, kmer_depth, job_type, work_dir, out_dir):

    for line in read_tsv(kmer_stat):
        if line[0]=="kmer_depth":
            sdepth = int(line[1])

    proportion=kmer_depth*1.0/sdepth
    if kmer_depth>sdepth:
        LOG.debug('The amount of sequencing data may be insufficient. Sequencing depth is only %s X' % sdepth)
        proportion = 1

    task = Task(
        id="choose_fastq",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
if [ {proportion} -ge 1] ; then
    ln -s {r1} {name}_choose.r1.fastq
    ln -s {r2} {name}_choose.r2.fastq
else
    {seqkit} sample -p {proportion} -2 -o {name}_choose.r1.fastq {r1}
    {seqkit} sample -p {proportion} -2 -o {name}_choose.r2.fastq {r2}
fi
#cp {name}_choose.r1.fastq {name}_choose.r2.fastq
""".format(seqkit=SEQKIT,
            r1=r1,
            r2=r2,
            name=name,
            proportion=proportion,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s_choose.r1.fastq" % name), os.path.join(work_dir, "%s_choose.r2.fastq" % name)


def choose_data(r1, r2, name, kmer_length, kmer_depth, thread, job_type, concurrent, refresh, work_dir, out_dir):

    r1, r2, kmer_stat, option = stat_kmer_depth(
        r1=r1,
        r2=r2,
        name=name,
        kmer_length=kmer_length,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=out_dir)

    dag = DAG("choose_data")
    data_task, r1, r2 = choose_data_task(
        r1=r1,
        r2=r2,
        name=name,
        kmer_stat=kmer_stat,
        kmer_depth=kmer_depth,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    dag.add_task(data_task)
    do_dag(dag, concurrent, refresh)

    return option, r1, r2


def obtain_contamination_task(taxid, name, kingdom, job_type, work_dir, out_dir, mode="general", cratio=10):

    for line in open(taxid):
        if line.startswith("#"):
             line = line.split("\t")
             prok_ratio = float(line[0].split(':')[1])
             top10 = int(line[1].split(':')[1])
        break
            

    if prok_ratio>=cratio or top10 >0:
        LOG.info("There is serious contamination of the sample, strict mode is mandatory")
        mode = "strict"

    if mode == "strict":
        run = 'blastdbcmd -db {dbase} -dbtype "nucl" -taxidlist {name}.prokaryotic.taxid -out {name}.prokaryotic.fa'.format(
            dbase=NT_TAXON["fungi"],
            name=name)
        pfa = "%s.prokaryotic.fa" % name
    else:
        run = ""
        pfa = ""

    task = Task(
        id="blastdbcmd_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
export PATH={blast}:$PATH
{script}/print.py {taxid} >{name}.prokaryotic.taxid
{run}
cat {pfa} {mbase} >{name}.ref.fa
#cp {name}.ref.fa {out_dir}
""".format(blast=BLAST_BIN,
            script=SCRIPTS,
            taxid=taxid,
            run=run,
            pfa=pfa,
            name=name,
            mbase=MC_TAXON[kingdom],
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.ref.fa" % name)


def split_data(r1, r2, name, number, job_type, concurrent, refresh, work_dir, out_dir, platform="illumina"):

    if platform in ["PromethION", "GridION" , "RSII", "Sequel"]:
        read = "%s.part_*.fast*" % name
        r2 = ""
    elif platform in ["illumina", "mgi"]:
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
    do_dag(dag, concurrent, refresh)

    temp = read_files(work_dir, read)
    reads = []

    if platform in ["illumina", "mgi"]:
        for i in temp:
            j = i.replace(".r1.part_", ".r2.part_")
            reads.append("%s %s" % (i, j))
    else:
        reads = temp

    return reads


def create_unmap_tasks(name, reference, reads, thread, job_type, work_dir, out_dir, split=""):

    option = OrderedDict()
    option["minimap2"] = {
        "version": get_version(SOFTWARE_VERSION["minimap2"]),
        "option:": "-ax sr"
    }
    option["samblaster"] = {
        "version": get_version(SOFTWARE_VERSION["samblaster"]),
        "option:": "default"
    }
    number = list(range(len(reads)))
    tmp = ""
    
    if split=="":
        tmp = ""
    else:
        tmp = "--split-prefix"

    tasks = ParallelTask(
        id="unmap_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={minimap2}:{samblaster}:$PATH
minimap2 {tmp} -t {thread} -ax sr {reference} {{reads}} |samblaster -u {name}.unmap{{number}}.fq
#cp {name}.unmap*.fq {out_dir}
""".format(
            minimap2=MINIMAP_BIN,
            samblaster=SAMBLASTER_BIN,
            tmp=tmp,
            reference=reference,
            name=name,
            thread=thread,
            out_dir=out_dir
        ),
        reads=reads,
        number=number
    )

    return tasks, os.path.join(work_dir, "%s.unmap*.fq" % name), option


def run_filter_contamination(r1, r2, name, kmer_length, kmer_depth, taxid, kingdom, thread, job_type, concurrent, refresh, work_dir, out_dir, split, mode="fast", cratio=10):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)
    r2 = check_paths(r2)
    taxid = check_path(taxid)
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    option, r1, r2 = choose_data(
        r1=r1,
        r2=r2,
        name=name,
        kmer_length=kmer_length,
        kmer_depth=kmer_depth,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=out_dir)
    options["software"].update(option)

    if mode!="fast":
        work_dict = {
            "data": "00_data",
            "ref": "01_ref",
            "ump": "02_ump"
        }
        for k, v in work_dict.items():
            mkdir(os.path.join(work_dir, v))

        reads = split_data(
            r1=r1,
            r2=r2,
            name=name,
            number=2000000,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["data"]),
            concurrent=concurrent,
            refresh=refresh,
            out_dir=out_dir,
            platform="illumina")

        dag = DAG("unmap_data")
        ref_task, ref= obtain_contamination_task(
            taxid=taxid,
            name=name,
            kingdom=kingdom,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["ref"]),
            out_dir=out_dir,
            mode=mode,
            cratio=cratio)
        dag.add_task(ref_task)

        unmap_tasks, reads, option = create_unmap_tasks(
            name=name,
            reference=ref,
            reads=reads,
            thread=thread,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["ump"]),
            out_dir=out_dir,
            split=split)
        dag.add_task(*unmap_tasks)
        ref_task.set_downstream(*unmap_tasks)
        do_dag(dag, concurrent, refresh)
        options["software"].update(option)

        reads = [reads]
    else:
        reads = [r1, r2]

    return reads, options


def filter_contamination(args):

    reads, options = run_filter_contamination(
        r1=args.read1,
        r2=args.read2,
        name=args.name,
        kmer_length=args.kmer_length,
        kmer_depth=args.depth,
        taxid=args.taxid,
        kingdom=args.kingdom,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        mode=args.mode,
        cratio=args.cratio,
        split=args.split)

    with open(os.path.join(args.out_dir, "filter_contamination.json"), "w") as fh:
        json.dump(options, fh, indent=2)


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

    parser = add_filter_contamination_args(parser)
    args = parser.parse_args()
    filter_contamination(args)


if __name__ == "__main__":
    main()
