#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

from survey.config import *
from survey.common import check_paths, mkdir, read_tsv, read_files
from dagflow import DAG, Task, ParallelTask, do_dag
from survey.bwa_mem import bwa_mem
from survey.parser import add_kmer_denovo_args

LOG = logging.getLogger(__name__)

LOG = logging.getLogger(__name__)
__version__ = "1.1.1"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
__all__ = []


def create_unmap_task(name, reference, r1, r2, thread, job_type, work_dir, out_dir):

    option = OrderedDict()
    option["minimap2"] = {
        "version": get_version(SOFTWARE_VERSION["minimap2"]),
        "option:": "-ax sr"
    }
    option["samblaster"] = {
        "version": get_version(SOFTWARE_VERSION["samblaster"]),
        "option:": "default"
    }

    task = Task(
        id="unmap__%s" % name,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={minimap2}:{samblaster}:$PATH
minimap2 -t {thread} -ax sr {reference} {r1} {r2} |samblaster -u {name}.unmap.fq
#cp {name}.unmap.fq {out_dir}
""".format(
            minimap2=MINIMAP_BIN,
            samblaster=SAMBLASTER_BIN,
            reference=reference,
            r1=r1,
            r2=r2,
            name=name,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.unmap.fq" % name), option


def merge_raw_data_task(name, r1, r2, tools, job_type, work_dir, out_dir):

    task = Task(
        id="merge_data",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
{tools} {r1} >{out_dir}/{name}.raw.r1.fastq
{tools} {r2} >{out_dir}/{name}.raw.r2.fastq
""".format(
            name=name,
            tools=tools,
            r1=r1,
            r2=r2,
            out_dir=out_dir
        )
    )

    r1 = os.path.join(out_dir, "%s.raw.r1.fastq" % name)
    r2 = os.path.join(out_dir, "%s.raw.r2.fastq" % name)

    return task, r1, r2


def kmerfreq_task(r1, r2, name, kmer_length, thread, job_type, work_dir, out_dir):
    """kmerfreq 统计reads的kmer频率分布"""

    if kmer_length >=17:
        kmer_length=17

    task = Task(
        id="kmerfreq",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={python}:$PATH
ls {r1} >{work_dir}/{name}.data
ls {r2} >>{work_dir}/{name}.data
{kmerfreq} -k {kmer_length} -t {thread} -p {work_dir}/{name} -q 33 -m 0 {work_dir}/{name}.data > {work_dir}/{name}.kmer.count
cp {work_dir}/{name}.freq.stat {out_dir}/{name}.kmerfreq.stat
python {script}/kmerfreq_stat.py {work_dir}/{name}.freq.stat >{work_dir}/{name}.kmer.stat
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            kmerfreq=KMERFREQ,
            kmer_length=kmer_length, #kmer数
            name=name,
            r1=r1, #文件列表
            r2=r2,
            thread=thread,
            work_dir=work_dir,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.freq.stat" % name), os.path.join(work_dir, "%s.kmer.stat" % name), os.path.join(work_dir, "%s.genome_estimate" % name)


def sample_fastq_task(r1, r2, proportion, name, job_type, work_dir):
    """fastq 取样"""

    task = Task(
        id="sample_fastq",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
{seqkit} sample -p {proportion} -2 -o {work_dir}/{name}_choose.r1.fastq {r1}
{seqkit} sample -p {proportion} -2 -o {work_dir}/{name}_choose.r2.fastq {r2}
""".format(seqkit=SEQKIT,
            r1=r1,
            r2=r2,
            name=name,
            proportion=proportion, #取样比例
            work_dir=work_dir
        )
    )

    return task, os.path.join(work_dir, "%s_choose.r1.fastq" % name), os.path.join(work_dir, "%s_choose.r2.fastq" % name), os.path.join(work_dir, "%s_choose.r*.fastq" % name)


def get_heterozygosity_task(histo, estimate, kingdom, name, job_type, work_dir, out_dir):
    """杂合度模拟"""

    task = Task(
        id="heterozygosity",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
export PATH={python}:$PATH
python {script}/fit_heterozygosity.py {histo} \
-e {estimate} --kingdom {kingdom} \
--name {out_dir}/{name} --database {script} > {out_dir}/{name}.heterozygosity.xls
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            histo=histo,
            estimate=estimate,
            name=name,
            kingdom=kingdom,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.heterozygosity.xls" % name), os.path.join(out_dir, "%s.heterozygosity.png" % name)


def get_jellyfish_task(fastq, name, depth, thread, job_type, work_dir, out_dir):
    """运行jellyfish"""

    task = Task(
        id="jellyfish",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={jellyfish}:$PATH
jellyfish count -m 21 -s 1G -t {thread} -C {fastq} -o {name}.jellyfish.jf
jellyfish histo -f {name}.jellyfish.jf -t {thread} > {name}.histogram_old.txt
head -n {depth} {name}.histogram_old.txt >{name}.histogram.txt
jellyfish stats -v -o {name}.stats.kmer.txt {name}.jellyfish.jf
cp {name}.histogram.txt {out_dir}
rm -rf {name}.jellyfish.jf
""".format(jellyfish=JELLYFISH_BIN,
            fastq=fastq,
            name=name,
            depth=depth,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, "%s.histogram.txt" % name)


def get_gse_scope_task(histogram, name, kmer_length, job_type, work_dir, out_dir):
    """运行genomescope和findGSE预估基因组大小"""

    task = Task(
        id="scope_and_gse",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
{rscript} {script}/genomescope.R {histogram} {kmer_length} 150 {work_dir} -1 1
cp {work_dir}/plot.png {out_dir}/{name}.genomescope.png
cp {work_dir}/summary.txt {out_dir}/{name}.genomescope.txt
{rscript} {script}/findGSE.R {histogram} {kmer_length} {work_dir}
mv {work_dir}/v1.94.est.{name}.histogram.txt.sizek{kmer_length}.curvefitted.pdf {work_dir}/{name}.findgse.pdf
cp {work_dir}/v1.94.est.{name}.histogram.txt.genome.size.estimated.k{kmer_length}to{kmer_length}.fitted.txt {out_dir}/{name}.findgse.txt
#convert -density 300 -quality 300 {work_dir}/{name}.findgse.pdf {work_dir}/{name}.findgse.png
{gs} -dQUIET -dNOSAFER -r300 -dBATCH -sDEVICE=pngalpha -dNOPAUSE -dNOPROMPT -sOutputFile={out_dir}/{name}.findgse.png {work_dir}/{name}.findgse.pdf
#cp {work_dir}/{name}.findgse-0.png {out_dir}/{name}.findgse.png
rm -rf {work_dir}/round*
""".format(rscript=RSCRIPT,
            script=SCRIPTS,
            gs=GHOSTSCRIPT,
            name=name,
            kmer_length=kmer_length,
            histogram=histogram,
            work_dir=work_dir,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.genomescope.txt" % name), os.path.join(out_dir, "%s.findgse.txt" % name), os.path.join(out_dir, "%s.genomescope.png" % name), os.path.join(out_dir, "%s.findgse.png" % name)


def stat_gc_depth_task(genome, bam, name, window, job_type, work_dir, out_dir):

    bam = check_paths(bam)

    task = Task(
        id="stat_coverage",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={samtools}:{python}:$PATH
samtools depth -aa {bam} > {work_dir}/{name}.depth
python {script}/stat_coverage.py -i {work_dir}/{name}.depth -d 1,5,10,20 -o {out_dir}/{name}.coverage.xlsx
python {script}/stat_length_gc.py -d {work_dir}/{name}.depth -g {genome} -n {out_dir}/{name}
python {script}/stat_gc_depth.py -d {work_dir}/{name}.depth -g {genome} -b 1000 -w 5000 -e 100 -n {work_dir}/{name}
python {script}/draw_depth_gc.py -gcd {work_dir}/{name}.stat_gc_depth.tsv -n {out_dir}/{name}
#python {script}/plot_gc_depth.py -gcd {work_dir}/{name}.stat_gc_depth.tsv -n {out_dir}/{name}
""".format(
            samtools=SAMTOOLS_BIN,
            script=SCRIPTS,
            python=PYTHON_BIN,
            genome=genome,
            bam=bam,
            name=name,
            window=window,
            work_dir=work_dir,
            out_dir=out_dir
        )
    )

    return task , os.path.join(out_dir, "%s.gc_depth.png" % name)


def soapdenovo_task(r1, r2, name, thread, queue, job_type, work_dir, out_dir):
    """soapdenove组装"""

    task = Task(
        id="soapdenovo",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, queue),
        script="""
export PATH={soapdenovo}:{script}:$PATH
echo -e "max_rd_len=151\n[LIB]\navg_ins=400\nreverse_seq=0\nasm_flags=3\nrank=1\npair_num_cutoff=3\nmap_len=64\nq1={r1}\nq2={r2}\n" > {name}.config
echo -e "{r1} {r2}" >{name}_ngs.list
SOAPdenovo-127mer all -s {name}.config -K 47 -p {thread} -d 2 -R 2 -o {name} >ass.log
cp {name}.scafSeq {out_dir}/{name}.asm.fasta
stat_genome.py -s {out_dir}/{name}.asm.fasta -r {out_dir}/{name}.asm.tsv
""".format(script=SCRIPTS,
            soapdenovo=SOAPDENOVO_BIN,
            r1=r1,
            r2=r2,
            name=name,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.asm.fasta" % name), os.path.join(out_dir, "%s.asm.tsv" % name), os.path.join(work_dir, "%s_ngs.list" % name)


def kmer_denovo_tasks(r1, r2, name, kmer_length, proportion, kingdom, thread, job_type, queue, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    clean1 = check_paths(r1)
    clean2 = check_paths(r2)
    work_dict = {
        "choose": "choose_data",
        "gse_scope": "02_gse_scope",
        "kmerfreq": "03_Kmerfreq",
        "denovo": "04_Soapdenovo",
        "gc_depth": "05_GC-depth"
    }

    choose_work = mkdir(os.path.join(work_dir, work_dict["choose"]))
    choose_task, choose_r1, choose_r2, choose_r = sample_fastq_task(
        r1=clean1,
        r2=clean2,
        proportion=proportion,
        name=name,
        job_type=job_type,
        work_dir=choose_work)

    heter_work = mkdir(os.path.join(work_dir, work_dict["kmerfreq"]))
    heter_out =  mkdir(os.path.join(out_dir, work_dict["kmerfreq"]))
    freq_task, histo, kmer_depth, estimate = kmerfreq_task(
        r1=choose_r1,
        r2=choose_r2,
        name=name,
        kmer_length=17,
        thread=thread,
        job_type=job_type,
        work_dir=heter_work,
        out_dir=heter_out
    )

    heter_task, stat_heter, heter_png = get_heterozygosity_task(
        histo=histo,
        estimate=estimate,
        kingdom=kingdom,
        name=name,
        job_type=job_type,
        work_dir=heter_work,
        out_dir=heter_out)

    scope_work = mkdir(os.path.join(work_dir, work_dict["gse_scope"]))
    scope_out =  mkdir(os.path.join(out_dir, work_dict["gse_scope"]))
    jellyfish_task, histogram = get_jellyfish_task(
        fastq=choose_r,
        name=name,
        depth=40*100,
        thread=thread,
        job_type=job_type,
        work_dir=scope_work,
        out_dir=scope_out
    )

    gse_scope_task, scope_txt, gse_txt, scope_png, gse_png = get_gse_scope_task(
        histogram=histogram,
        name=name,
        kmer_length=kmer_length,
        job_type=job_type,
        work_dir=scope_work,
        out_dir=scope_out
    )

    denovo_work = mkdir(os.path.join(work_dir, work_dict["denovo"]))
    denovo_out =  mkdir(os.path.join(out_dir, work_dict["denovo"]))
    denovo_task, genome, stat_genome, ngs_list = soapdenovo_task(
        r1=clean1,
        r2=clean2,
        name=name,
        thread=thread*2,
        queue=queue,
        job_type=job_type,
        work_dir=denovo_work,
        out_dir=denovo_out)

    return choose_task, freq_task, heter_task, jellyfish_task, gse_scope_task, denovo_task, stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome, genome, ngs_list


def run_gc_depth(genome, fastq_list, name ,window ,thread, job_type, concurrent, refresh, work_dir, out_dir):

    genome, fastq_list = check_paths([genome, fastq_list])

    sort_bam, genome = bwa_mem(
        fastq_list=fastq_list,
        genome=genome,
        name=name,
        number=5000000,
        data_type='',
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=work_dir
    )

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


def run_kmer_denovo(r1, r2, name, kingdom, kmer_length, sample_depth, thread, asm, window, job_type, queue, concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    r1 = check_paths(r1)
    r2 = check_paths(r2)

    if r1[0].endswith(".gz") or r2[0].endswith(".gz"):
        tools = "zcat"
    else:
        tools = "cat"

    dag_data = DAG("survey_data")

    data_work = mkdir(os.path.join(work_dir, "choose_data"))
    cat_data_task, clean1, clean2 = merge_raw_data_task(
        name = name,
        r1 = " ".join(r1),
        r2 = " ".join(r2),
        tools=tools,
        job_type=job_type,
        work_dir=data_work,
        out_dir=data_work
    )

    freq_task1, histo1, kmer_stat, estimate1 = kmerfreq_task(
        r1=clean1,
        r2=clean2,
        name=name,
        kmer_length=17,
        thread=thread,
        job_type=job_type,
        work_dir=data_work,
        out_dir=data_work)

    dag_data.add_task(cat_data_task)
    dag_data.add_task(freq_task1)
    freq_task1.set_upstream(cat_data_task)
    do_dag(dag_data, concurrent, refresh)

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
    return stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome


def kmer_denovo(args):

    if args.no_asm:
        asm="false"
    else:
        asm="true"

    stat_heter, heter_png, scope_txt, gse_txt, scope_png, gse_png, stat_genome = run_kmer_denovo(
        r1=args.read1,
        r2=args.read2,
        name=args.name,
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

    parser = add_kmer_denovo_args(parser)
    args = parser.parse_args()
    kmer_denovo(args)


if __name__ == "__main__":
    main()
