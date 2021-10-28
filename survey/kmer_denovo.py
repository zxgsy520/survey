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
from survey.filter_cont import run_filter_cont
from survey.parser import add_kmer_denovo_args


LOG = logging.getLogger(__name__)

LOG = logging.getLogger(__name__)
__version__ = "2.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
__all__ = []

def create_jellyfish_tasks(reads, name, kingdom, kmer_length, thread, job_type, work_dir, out_dir, maximal=1000):

    option = {}
    option["jellyfish"] = {
        "version": get_version(SOFTWARE_VERSION["jellyfish"]),
        "option": "-m 21 -s 200M"
    }
    option["kmc"] = {
        "version": get_version(SOFTWARE_VERSION["kmc"]),
        "option": "off"
    }
    prefix = []
    jelly = []
    n = 1
    for i in reads:
        prefix.append(os.path.basename(i))
        jelly.append("jellyfish_%s/%s.jellyfish.jf" % (n, os.path.basename(i)))
        n += 1
    id = "jellyfish"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={jellyfish}:$PATH
jellyfish count -m {kmer_length} -s 200M -t {thread} -U {maximal} -L 0 --disk \\
-C {{reads}} -o {{prefix}}.jellyfish.jf
""".format(jellyfish=JELLYFISH_BIN,
           kmer_length=kmer_length,
           maximal=maximal,
           thread=thread),
        prefix=prefix,
        reads=reads,
    )
    
    if len(prefix)==1:
        temp = ""
        x = jelly[0]
    else:
        temp = """jellyfish merge -U {maximal} -L 0 -o {name}.jellyfish.jf \\
{hash}
""".format(hash=" ".join(jelly), name=name, maximal=maximal)
        x = "%s.jellyfish.jf" % name

    join_task = Task(
        id="merge_jellyfish",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={jellyfish}:{python}:$PATH
#jellyfish merge -U 0 -L {maximal} jellyfish*/*.jellyfish.jf -o {name}.jellyfish.jf
{temp}
jellyfish histo -f {x} -t 4 > {name}.histogram_t.txt
cat {name}.histogram_t.txt |sed 's/ /\\t/' >{name}.histogram.txt
python {script}/fit_heterozygosity.py {name}.histogram.txt --kingdom {kingdom} \\
--name {name} --length {kmer_length} --database {script} > {name}.heterozygosity.xls
cp {name}.heterozygosity.xls {out_dir}
cp {name}.kmer.p* {name}.heterozygosity.p* {out_dir}
#rm -rf jellyfish*
""".format(jellyfish=JELLYFISH_BIN,
           python=PYTHON_BIN,
           script=SCRIPTS,
           temp=temp,
           x=x,
           name=name,
           kmer_length=kmer_length,
           kingdom=kingdom,
           maximal=maximal,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)
    kmer_hist = os.path.join(work_dir, "%s.histogram.txt" % name)
    stat_heter = os.path.join(work_dir, "%s.heterozygosity.xls" % name)
    heter_png = os.path.join(work_dir, "%s.heterozygosity.png" % name)

    return tasks, join_task, option, kmer_hist, stat_heter, heter_png


def create_kmc_task(reads, name, kingdom, kmer_length, thread, job_type, work_dir, out_dir, maximal=1000):

    option = {}
    option["kmc"] = {
        "version": get_version(SOFTWARE_VERSION["kmc"]),
        "option": "-k%s -ci1 -cs%s" % (kmer_length, maximal)
    }
    option["jellyfish"] = {
        "version": get_version(SOFTWARE_VERSION["jellyfish"]),
        "option": "off"
    }
    option["survey"] = {
        "version": "v2.1.0",
        "website": "https://github.com/zxgsy520/survey"
    }

    memory = thread*4
    if memory >=20:
        memory = 20

    kmc_task = Task(
        id="kmc",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={kmc}:{python}:$PATH
ls {reads} >input.list
kmc -k{kmer_length} -t{thread} -m{memory} -ci1 -cs{maximal} -r @input.list {name} ./
kmc_tools transform {name} histogram {name}.histogram.txt -ci1 -cx{maximal}
cp {name}.histogram.txt {out_dir}
rm {name}.kmc*
""".format(kmc=KMC_BIN,
            python=PYTHON_BIN,
            reads=reads,
            name=name,
            kmer_length=kmer_length,
            maximal=maximal,
            thread=thread,
            memory=memory,
            out_dir=out_dir
        )
    )
    heter_task = Task(
        id="heterozygosity",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={python}:$PATH
python {script}/fit_heterozygosity.py {name}.histogram.txt --kingdom {kingdom} \\
--name {name} --length {kmer_length} --database {script} > {name}.heterozygosity.xls
cp {name}.heterozygosity.xls {out_dir}
cp {name}.kmer.p* {name}.heterozygosity.p* {out_dir}
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            name=name,
            kingdom=kingdom,
            kmer_length=kmer_length,
            out_dir=out_dir
        )
    )

    heter_task.set_upstream(kmc_task)
    kmer_hist = os.path.join(work_dir, "%s.histogram.txt" % name)
    stat_heter = os.path.join(work_dir, "%s.heterozygosity.xls" % name)
    heter_png = os.path.join(work_dir, "%s.heterozygosity.png" % name)

    return kmc_task, heter_task, option, kmer_hist, stat_heter, heter_png


def create_gse_scope_task(histogram, name, kmer_length, job_type, work_dir, out_dir):

    option = {}
    option["genomescope"] = {
        "version": "1.0.0",
        "option": "default"
    }
    option["findgse"] = {
        "version": "-",
        "option": "default"
    }

    scope_task = Task(
        id="genomescope",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{rscript} {script}/genomescope.R {histogram} {kmer_length} 150 ./ -1 1
mv plot.png {name}.genomescope.png
mv summary.txt {name}.genomescope.txt
cp {name}.genomescope.png {name}.genomescope.txt {out_dir}
rm -rf round*
""".format(rscript=RSCRIPT,
            script=SCRIPTS,
            name=name,
            kmer_length=kmer_length,
            histogram=histogram,
            out_dir=out_dir
        )
    )
    gse_task = Task(
        id="findgse",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{rscript} {script}/findGSE.R {histogram} {kmer_length} ./
mv v1.94.est.{name}.histogram.txt.sizek{kmer_length}.curvefitted.pdf {name}.findgse.pdf
mv v1.94.est.{name}.histogram.txt.genome.size.estimated.k{kmer_length}to{kmer_length}.fitted.txt {name}.findgse.txt
{gs} -dQUIET -dNOSAFER -r300 -dBATCH -sDEVICE=pngalpha -dNOPAUSE -dNOPROMPT -sOutputFile={name}.findgse.png {name}.findgse.pdf
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
    scope_txt = os.path.join(work_dir, "%s.genomescope.txt" % name)
    gse_txt = os.path.join(work_dir, "%s.findgse.txt" % name)
    scope_png = os.path.join(work_dir, "%s.genomescope.png" % name)
    gse_png = os.path.join(work_dir, "%s.findgse.png" % name)

    return scope_task, gse_task, option, scope_txt, gse_txt, scope_png, gse_png


def create_soapdenovo_task(r1, r2, name, thread, queue, job_type, work_dir, out_dir):

    option = {}
    option["soapdenovo"] = {
        "version": get_version(SOFTWARE_VERSION["soapdenovo"]),
        "option": "max_rd_len=151 avg_ins=400 reverse_seq=0 asm_flags=3 rank=1 pair_num_cutoff=3 map_len=64"
    }
    option["megahit"] = {
        "version": get_version(SOFTWARE_VERSION["megahit"]),
        "option:": "off"
    }

    task = Task(
        id="soapdenovo_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, queue),
        script="""
export PATH={soapdenovo}:{script}:$PATH
export LD_LIBRARY_PATH={soapdenovo}:$LD_LIBRARY_PATH
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


def create_megahit_task(read1, read2, name, thread, queue, job_type, work_dir, out_dir, memory=50):

    if queue!="":
        q = "-q %s" % queue
    else:
        q = ""
    option = {}
    option["megahit"] = {
        "version": get_version(SOFTWARE_VERSION["megahit"]),
        "option:": "--memory %s --num-cpu-threads %s" % (memory, thread)
    }
    option["soapdenovo"] = {
        "version": get_version(SOFTWARE_VERSION["soapdenovo"]),
        "option": "off"
    }
    task = Task(
        id="megahit_%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-l vf=%dG -pe smp %s %s" % (memory, thread, q),
        script="""
export PATH={megahit}:{script}:$PATH
megahit -1 {r1} -2 {r2} --memory {memory} --num-cpu-threads {thread} --out-dir 02_{name}_megahit
cp 02_{name}_megahit/final.contigs.fa {name}.asm.fasta
stat_genome.py -s {name}.asm.fasta -r {name}.asm.tsv
cp {name}.asm.tsv {name}.asm.fasta {out_dir}
rm -rf 02_{name}_megahit
""".format(megahit=MEGAHIT_BIN,
        script=SCRIPTS,
        name=name,
        r1=read1,
        r2=read2,
        thread=thread,
        memory=memory,
        out_dir=out_dir)
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
#python {script}/find_peak.py {name}.stat_gc_depth.tsv > gc_depth_describe.txt
cp {name}.coverage.xlsx {name}.length_gc.xls {name}.stat_gc_depth.tsv {name}.gc_depth.p* {out_dir}
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

    gc_depth = os.path.join(out_dir, "%s.stat_gc_depth.tsv" % name)
    gc_depth_png = os.path.join(out_dir, "%s.gc_depth.png" % name)

    return task , gc_depth, gc_depth_png


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

    gc_depth_task, gc_depth, gc_depth_png = stat_gc_depth_task(
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

    return  gc_depth, gc_depth_png


def run_kmer_denovo(r1, r2, taxid, name, mode, cratio, kmer_length, kmer_depth, kingdom,
    genome_type, asm, window, minmapq, maximal, thread, job_type, queue, concurrent,
    refresh, work_dir, out_dir, split, platform="illumina"):

    maximal = int(maximal)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)
    r2 = check_paths(r2)

    work_dict = {
        "contamination": "01_contamination",
        "kmc": "02_kmc",
        "gse_scope": "03_gse_scope",
        "denovo": "04_Soapdenovo",
        "gc_depth": "05_GC-depth"
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        if k=="contamination":
            continue
        mkdir(os.path.join(out_dir, v))

    reads, options = run_filter_cont(
        read1=r1,
        read2=r2,
        name=name,
        kmer_length=kmer_length,
        kmer_depth=kmer_depth,
        taxid=taxid,
        kingdom=kingdom,
        genome_type=genome_type,
        minmapq=minmapq,
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
    if genome_type in ["small"]:
        kmc_task, heter_task, option, kmer_hist, stat_heter, heter_png = create_kmc_task(
            reads=" ".join(reads),
            name=name,
            kmer_length=kmer_length,
            kingdom=kingdom,
            thread=thread,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["kmc"]),
            out_dir=os.path.join(out_dir, work_dict["kmc"]),
            maximal=maximal,
        )
        dag.add_task(kmc_task)
        dag.add_task(heter_task)
    else:
        jellyfish_tasks, heter_task, option, kmer_hist, stat_heter, heter_png = create_jellyfish_tasks(
            reads=reads,
            name=name,
            kmer_length=kmer_length,
            kingdom=kingdom,
            thread=thread,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["kmc"]),
            out_dir=os.path.join(out_dir, work_dict["kmc"]),
            maximal=maximal,
        )

        dag.add_task(*jellyfish_tasks)
        dag.add_task(heter_task)

    options["software"].update(option)

    scope_task, gse_task, option, scope_txt, gse_txt, scope_png, gse_png = create_gse_scope_task(
        histogram=kmer_hist,
        name=name,
        kmer_length=kmer_length,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["gse_scope"]),
        out_dir=os.path.join(out_dir, work_dict["gse_scope"])
    )
    dag.add_task(scope_task)
    dag.add_task(gse_task)
    if genome_type in ["small"]:
        scope_task.set_upstream(kmc_task)
        gse_task.set_upstream(kmc_task)
    else:
        scope_task.set_upstream(heter_task)
        gse_task.set_upstream(heter_task)
    options["software"].update(option)

    if asm=="soapdenovo":
        denovo_task, genome, stat_genome, option = create_soapdenovo_task(
            r1=" ".join(r1),
            r2=" ".join(r2),
            name=name,
            thread=thread,
            queue=queue,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["denovo"]),
            out_dir=os.path.join(out_dir, work_dict["denovo"])
        )
        dag.add_task(denovo_task)
        options["software"].update(option)
    elif asm=="megahit":
        denovo_task, genome, stat_genome, option = create_megahit_task(
            read1=" ".join(r1),
            read2=" ".join(r2),
            name=name,
            thread=thread,
            queue=queue,
            job_type=job_type,
            work_dir=os.path.join(work_dir, work_dict["denovo"]),
            out_dir=os.path.join(out_dir, work_dict["denovo"])
        )
        dag.add_task(denovo_task)
        options["software"].update(option)
    else:
        genome = "false"
        stat_genome= "false"
    do_dag(dag, concurrent, refresh)

    if asm!="":
        gc_depth_png, gc_depth = run_gc_depth(
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
        gc_depth = ""
        gc_depth_png = heter_png

    return options, stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, gc_depth_png, gc_depth


def kmer_denovo(args):

    options, stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, gc_depth_png, gc_depth = run_kmer_denovo(
        r1=args.read1,
        r2=args.read2,
        name=args.name,
        mode=args.mode,
        taxid=args.taxid,
        cratio=args.cratio,
        kingdom=args.kingdom,
        genome_type=args.genome_type,
        kmer_length=args.kmer_length,
        kmer_depth=args.depth,
        minmapq=args.minmapq,
        maximal=args.maximal,
        thread=args.thread,
        asm=args.asm,
        window=args.window,
        job_type=args.job_type,
        queue="-q %s" % args.queue,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        split=args.split
    )
    with open(os.path.join(args.out_dir, "kmer_denovo.json"), "w") as fh:
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

    parser = add_kmer_denovo_args(parser)
    args = parser.parse_args()
    kmer_denovo(args)


if __name__ == "__main__":
    main()
