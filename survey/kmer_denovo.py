#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from survey.config import *
from survey.common import check_path, check_paths, mkdir, read_tsv, read_files, get_version
from dagflow import DAG, Task, ParallelTask, do_dag
from survey.minimap import minimap
from survey.filter_contamination import run_filter_contamination
from survey.parser import add_kmer_denovo_args


LOG = logging.getLogger(__name__)

LOG = logging.getLogger(__name__)
__version__ = "1.1.1"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
__all__ = []


def create_kmerfreq_task(reads, name, kmer_length, thread, job_type, work_dir, out_dir):

    option = {}
    option["kmerfreq"] = {
        "version": get_version(SOFTWARE_VERSION["kmerfreq"]),
        "option": "-q 33 -m 0 -k %s" % kmer_length
    }

    if kmer_length >=17:
        kmer_length=17

    task = Task(
        id="kmerfreq_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={python}:$PATH
ls {reads} >{name}.data
{kmerfreq} -k {kmer_length} -t {thread} -p {name} -q 33 -m 0 {name}.data > {name}.kmer.count
python {script}/kmerfreq_stat.py {name}.freq.stat >{name}.kmer.stat
cp {name}.freq.stat {out_dir}/{name}.kmerfreq.stat
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            kmerfreq=KMERFREQ,
            kmer_length=kmer_length,
            name=name,
            reads=reads,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.freq.stat" % name), os.path.join(work_dir, "%s.genome_estimate" % name), option


def get_heterozygosity_task(histo, estimate, kingdom, name, job_type, work_dir, out_dir):

    task = Task(
        id="heterozygosity",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
export PATH={python}:$PATH
python {script}/fit_heterozygosity.py {histo} \
-e {estimate} --kingdom {kingdom} \
--name {name} --database {script} > {name}.heterozygosity.xls
cp {name}.heterozygosity.xls {out_dir}
cp {name}.kmer.p* {name}.heterozygosity.p* {out_dir}
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            histo=histo,
            estimate=estimate,
            name=name,
            kingdom=kingdom,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.heterozygosity.xls" % name), os.path.join(work_dir, "%s.heterozygosity.png" % name)


def kmerfreq(reads, name, kingdom, kmer_length, thread, job_type, work_dir, out_dir):

    kmerfreq_task, histo, estimate, option = create_kmerfreq_task(
        reads=reads,
        name=name,
        kmer_length=kmer_length,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)

    heter_task, stat_heter, heter_png = get_heterozygosity_task(
        histo=histo,
        estimate=estimate,
        kingdom=kingdom,
        name=name,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    heter_task.set_upstream(kmerfreq_task)

    return kmerfreq_task, heter_task, stat_heter, heter_png, option


def create_jellyfish_task(reads, name, thread, job_type, work_dir, out_dir, mode="general"):

    option = {}
    option["jellyfish"] = {
        "version": get_version(SOFTWARE_VERSION["jellyfish"]),
        "option": "-m 21 -s 1G"
    }

    if mode=="general":
        histout = "%s.histogram.txt" % name
        runh = ""
    else:
        histout = "%s.histogram_old.txt" % name
        runh = "head -n 5000 %s.histogram_old.txt >%s.histogram.txt" % (name, name)

    task = Task(
        id="jellyfish",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={jellyfish}:$PATH
jellyfish count -m 21 -s 1G -t {thread} -C {reads} -o {name}.jellyfish.jf
jellyfish histo -f {name}.jellyfish.jf -t {thread} > {histout}
{runh}
jellyfish stats -v -o {name}.stats.kmer.txt {name}.jellyfish.jf
cp {name}.histogram.txt {out_dir}
rm -rf {name}.jellyfish.jf
""".format(jellyfish=JELLYFISH_BIN,
            reads=reads,
            name=name,
            histout=histout,
            runh=runh,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.histogram.txt" % name), option


def create_gse_scope_task(histogram, name, kmer_length, job_type, work_dir, out_dir):

    task = Task(
        id="scope_gse_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
{rscript} {script}/genomescope.R {histogram} {kmer_length} 150 ./ -1 1
mv plot.png {name}.genomescope.png
mv summary.txt {name}.genomescope.txt
{rscript} {script}/findGSE.R {histogram} {kmer_length} ./
mv v1.94.est.{name}.histogram.txt.sizek{kmer_length}.curvefitted.pdf {name}.findgse.pdf
mv v1.94.est.{name}.histogram.txt.genome.size.estimated.k{kmer_length}to{kmer_length}.fitted.txt {name}.findgse.txt
{gs} -dQUIET -dNOSAFER -r300 -dBATCH -sDEVICE=pngalpha -dNOPAUSE -dNOPROMPT -sOutputFile={name}.findgse.png {name}.findgse.pdf
rm -rf round*
cp {name}.genomescope.png {name}.genomescope.txt {out_dir}
cp {name}.findgse.txt {name}.findgse.png {name}.findgse.pdf {out_dir}
""".format(rscript=RSCRIPT,
            script=SCRIPTS,
            gs=GHOSTSCRIPT,
            name=name,
            kmer_length=kmer_length,
            histogram=histogram,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.genomescope.txt" % name), os.path.join(work_dir, "%s.findgse.txt" % name), os.path.join(work_dir, "%s.genomescope.png" % name), os.path.join(work_dir, "%s.findgse.png" % name)


def gse_scope(reads, name, kmer_length, thread, job_type, work_dir, out_dir, mode="general"):

    jellyfish_task, histogram, option = create_jellyfish_task(
        reads=reads,
        name=name,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir,
        mode=mode)

    gse_scope_task, scope_txt, gse_txt, scope_png, gse_png = create_gse_scope_task(
        histogram=histogram,
        name=name,
        kmer_length=kmer_length,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    gse_scope_task.set_upstream(jellyfish_task)

    return jellyfish_task, gse_scope_task, scope_txt, gse_txt, scope_png, gse_png, option


def create_soapdenovo_task(r1, r2, name, thread, queue, job_type, work_dir, out_dir):

    option = {}
    option["soapdenovo"] = {
        "version": get_version(SOFTWARE_VERSION["soapdenovo"]),
        "option": "max_rd_len=151 avg_ins=400 reverse_seq=0 asm_flags=3 rank=1 pair_num_cutoff=3 map_len=64"
    }

    task = Task(
        id="soapdenovo_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, queue),
        script="""
export PATH={soapdenovo}:{script}:$PATH
echo -e "max_rd_len=151\n[LIB]\navg_ins=400\nreverse_seq=0\nasm_flags=3\nrank=1\npair_num_cutoff=3\nmap_len=64\nq1={r1}\nq2={r2}\n" > {name}.config
echo -e "{r1} {r2}" >{name}_ngs.list
SOAPdenovo-127mer all -s {name}.config -K 47 -p {thread} -d 2 -R 2 -o {name} >ass.log
mv {name}.scafSeq {name}.asm.fasta
stat_genome.py -s {name}.asm.fasta -r {name}.asm.tsv
cp {name}.asm.tsv {name}.asm.fasta {out_dir}
""".format(script=SCRIPTS,
            soapdenovo=SOAPDENOVO_BIN,
            r1=r1,
            r2=r2,
            name=name,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.asm.fasta" % name), os.path.join(work_dir, "%s.asm.tsv" % name), option


def stat_gc_depth_task(genome, bam, name, window, job_type, work_dir, out_dir):

    bam = check_paths(bam)

    task = Task(
        id="stat_coverage",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={samtools}:{python}:$PATH
samtools depth -aa {bam} > {name}.depth
python {script}/stat_coverage.py -i {name}.depth -d 1,5,10,20 -o {name}.coverage.xlsx
python {script}/stat_length_gc.py -d {name}.depth -g {genome} -n {name}
python {script}/stat_gc_depth.py -d {name}.depth -g {genome} -b 1000 -w 5000 -e 100 -n {name}
python {script}/draw_depth_gc.py -gcd {name}.stat_gc_depth.tsv -n {name}
#python {script}/plot_gc_depth.py -gcd {name}.stat_gc_depth.tsv -n {name}
cp {name}.coverage.xlsx {name}.length_gc.xls {name}.gc_depth.p* {out_dir}
""".format(
            samtools=SAMTOOLS_BIN,
            script=SCRIPTS,
            python=PYTHON_BIN,
            genome=genome,
            bam=bam,
            name=name,
            window=window,
            out_dir=out_dir
        )
    )

    return task , os.path.join(work_dir, "%s.gc_depth.png" % name)


def run_gc_depth(genome, r1, r2, name, platform, split, window ,thread, job_type, concurrent, refresh, work_dir, out_dir):

    genome = check_path(genome)
    r1 = check_paths(r1)
    r2 = check_paths(r2)

    sort_bam = minimap(
       r1=r1,
       r2=r2,
       genome=genome,
       name=name,
       split=split,
       platform=platform,
       number=5000000,
       thread=thread,
       job_type=job_type,
       concurrent=concurrent,
       refresh=refresh,
       work_dir=work_dir,
       out_dir=out_dir)

    sort_bam = check_paths(sort_bam)
    dag = DAG("gc_depth")

    gc_depth_task, gc_depth_png  = stat_gc_depth_task(
        genome=genome,
        bam=sort_bam,
        name=name,
        window=window,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )

    dag.add_task(gc_depth_task)
    do_dag(dag, concurrent, refresh)

    return  gc_depth_png


def run_kmer_denovo(r1, r2, taxid, name, mode, cratio, kmer_length, kmer_depth, kingdom, asm, window, thread, job_type, queue, concurrent, refresh, work_dir, out_dir, split, platform="illumina"):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)
    r2 = check_paths(r2)

    work_dict = {
        "contamination": "01_contamination",
        "gse_scope": "02_gse_scope",
        "kmerfreq": "03_Kmerfreq",
        "denovo": "04_Soapdenovo",
        "gc_depth": "05_GC-depth"
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        if k=="contamination":
            continue
        mkdir(os.path.join(out_dir, v))

    reads, options = run_filter_contamination(
        r1=r1,
        r2=r2,
        name=name,
        kmer_length=kmer_length,
        kmer_depth=kmer_depth,
        taxid=taxid,
        kingdom=kingdom,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, work_dict["contamination"]),
        out_dir=out_dir,
        mode=mode,
        cratio=cratio,
        split=split)

    dag = DAG("kmer_denovo")
    jellyfish_task, gse_scope_task, scope_txt, gse_txt, scope_png, gse_png, option = gse_scope(
        reads=" ".join(reads),
        name=name,
        kmer_length=kmer_length,
        thread=thread,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["gse_scope"]),
        out_dir=os.path.join(out_dir, work_dict["gse_scope"]),
        mode=mode)
    options["software"].update(option)
    dag.add_task(jellyfish_task)
    dag.add_task(gse_scope_task)

    kmerfreq_task, heter_task, stat_heter, heter_png, option = kmerfreq(
        reads=" ".join(reads),
        name=name,
        kingdom=kingdom,
        kmer_length=kmer_length,
        thread=thread,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["kmerfreq"]),
        out_dir=os.path.join(out_dir, work_dict["kmerfreq"]))
    options["software"].update(option)
    dag.add_task(kmerfreq_task)
    dag.add_task(heter_task)

    denovo_task, genome, stat_genome, option = create_soapdenovo_task(
        r1=" ".join(r1),
        r2=" ".join(r2),
        name=name,
        thread=thread,
        queue=queue,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["denovo"]),
        out_dir=os.path.join(out_dir, work_dict["denovo"]))
    if asm=="true":
        dag.add_task(denovo_task)
    else:
        genome = "false"
        stat_genome= "false"
    do_dag(dag, concurrent, refresh)

    if asm=="true":
        gc_depth = run_gc_depth(
            genome=genome,
            r1=" ".join(r1),
            r2=" ".join(r2),
            name=name,
            platform=platform,
            split="no_split",
            window=window,
            thread=thread,
            job_type=job_type,
            concurrent=concurrent,
            refresh=refresh,
            work_dir=os.path.join(work_dir, work_dict["gc_depth"]),
            out_dir=os.path.join(out_dir, work_dict["gc_depth"]))
    else:
        gc_depth = heter_png

    with open(os.path.join(out_dir, "kmer_denovo.json"), "w") as fh:
         json.dump(options, fh, indent=2)

    return stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, gc_depth


def kmer_denovo(args):
    if args.no_asm:
        asm="false"
    else:
        asm="true"

    run_kmer_denovo(
        r1=args.read1,
        r2=args.read2,
        name=args.name,
        taxid=args.taxid,
        mode=args.mode,
        cratio=args.cratio,
        kingdom=args.kingdom,
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

    parser = add_kmer_denovo_args(parser)
    args = parser.parse_args()
    kmer_denovo(args)


if __name__ == "__main__":
    main()
