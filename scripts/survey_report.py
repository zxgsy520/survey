#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
[general]
project=真菌A18基因组测序分析
id=WHWLZ-201906095A-01
name=A18
species=Unknow
strain=A18
sequencer=Illumina
author=张兴国
reviewer=汤冬
homogeneous=True
assembly=True
pollution_description= 数据存在污染对后续的分析会有影响
depth_description=样本的基因组碱基深度主要分布在0-120x；基因平均GC含量主要分布在20-70%。基因组GC-Depth中有明显分离的聚团现象，基因组碱基深度有明显分离，说明基因组中含有其他外源污染
"""
import sys
import json
import shutil
import logging
import os.path
import argparse

from datetime import datetime
from jinja2 import Template
from docxtpl import DocxTemplate
try:
    from ConfigParser import RawConfigParser
except:
    from configparser import RawConfigParser

LOG = logging.getLogger(__name__)

__version__ = "2.2.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def check_path(path):

    path = os.path.abspath(path)

    if not os.path.exists(path):
        msg = "File not found '{path}'".format(**locals())
        LOG.error(msg)
        raise Exception(msg)

    return path


def check_paths(obj):
    """
    check the existence of paths
    :param obj:
    :return: abs paths
    """

    if isinstance(obj, list):
        r = []
        for path in obj:
            r.append(check_path(path))

        return r
    else:
        return check_path(obj)


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def read_picture(pictures):

    r = []

    for i in pictures:
        i = check_path(i)

        r.append(i)
        for j in ['.svg', '.pdf']:
            ni = i.replace('.png', j)

            if not os.path.exists(ni):
                continue

            r.append(ni)

    return r


def bp2mbgb(string, types="gb"):

    string = str(string).strip().replace(',', '')

    if types == "gb":
        base = 1000000000
    elif types == "mb":
        base = 1000000
    else:
        base = 1.0

    string = float(string)/base

    return "{0:,.2f}".format(string)


def read_config(cfg):
    """
    read config fron ini
    :param cfg:
    :return:
    """
    check_paths(cfg)

    r = {}
    #config = ConfigParser()
    config = RawConfigParser()
    #config.read(cfg, encoding='utf-8')
    config.read(cfg)

    for section in config.sections():
        r[section] = {}

        for option in config.options(section):
            value = config.get(section, option)

            if type(value) == type(b''):
                r[section][option] = value.decode("utf-8").strip()
            else:
                r[section][option] = value.strip()

    return r


def read_tsv(file, sep='\t'):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        line=line.replace("\t\t","\t")
        r.append(line.split(sep))

    return r


def read_table_data(file):
    
    return read_tsv(file,'\t')[0]


def read_table_species(filei, sep='\t'):
 
    r = []
    pd = ""    

    for line in open(file):
        line = line.strip()
        
        if line.startswith("#!"):
            pd = line..strip("#!")
        if not line or line.startswith("#"):
            continue
        r.append(line.split(sep))

   return r, pd


def read_table_findgse(file):

    gse = []
    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        line = line.split()

        if line[0]=="Est.":
            gse.append(line[-1])
        elif line[0]=="size_all":
            gse.append(line[-1])

    return [int(gse[0]), "-", round(float(gse[1])*100,2)]


def read_table_scope(file):

    scope = []
    title = ["Heterozygosity", "Haploid", "Repeat"]

    for line in open(file, "r"):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        line = line.split()

        if line[0]=="Heterozygosity":
            scope.append((float(line[-1].replace("%",""))+float(line[-2].replace("%","")))/2)
        elif line[1] in title:        
            scope.append(line[-2].replace(",",""))
        else:
             continue

    return [int(scope[1]), round(scope[0],2), round(int(scope[2])*100.0/int(scope[1]),2)]


def read_table_kmer(file):

    return read_tsv(file)[0]


def read_table_assembly(file):

    return read_tsv(file)


def run_report(cfg, survey_json,
               table_data, table_species, table_findgse, table_scope, table_kmer, table_assembly,
               figure_base_quality, figure_base_content, figure_base_gc, figure_species, figure_findgse, figure_genomescope, figure_kmer, figure_gc_depth,
               tpl_docx, tpl_html, out_dir
               ):

    out_dir = mkdir(out_dir)
    now = datetime.now()
    config = read_config(cfg)
    table_data = read_table_data(table_data)
    table_species, poll_desc = read_table_species(table_species, '\t')

    if config["general"]["assembly"].lower()=="true":
        if table_assembly.lower=="no":
            LOG.debug("Please enter the assembly result statistics file")

        assembly = True
        table_assembly = read_table_assembly(table_assembly)

    else:
        assembly = False
        table_assembly = [[0,0,0,0],[0],[0],[0],[0],[0],[0,0,0,0,0]]
 
    if config["general"]["homogeneous"].lower()=="true":
        if table_findgse.lower=="no" or  table_scope.lower=="no":
            LOG.debug("Please output the genomic prediction results of FindGSE and GenomeScope")

        homogeneous=True
        table_findgse = read_table_findgse(table_findgse)
        table_scope = read_table_scope(table_scope)
        findgse_genome = table_findgse[0]
        scope_genome = table_scope[0]
        max_genome = max(findgse_genome, scope_genome)
        min_genome = min(findgse_genome, scope_genome)
        average_genome = (max_genome+min_genome)*1.0/2
        software="FindGSE和GenomeScope"
        estimate_size = '%s-%s' % (bp2mbgb(min_genome, "mb"), bp2mbgb(max_genome, "mb"))
        heterozygosity = table_scope[1]
    else:
        if table_kmer.lower=="no":
            LOG.debug("Please output the genomic prediction results of Kmerfreq")
        homogeneous = False
        table_kmer = read_table_kmer(table_kmer)
        software="python的scipy包"
        estimate_size= bp2mbgb(table_kmer[4].replace(",",""), "mb")
        heterozygosity= float(table_kmer[5].replace("%",""))
        findgse_genome = 0
        scope_genome = 0
        max_genome = 0
        min_genome = 0
        average_genome = 0
    
    r = {
        "project": "",
        "id": "",
        "name": "",
        "species": "",
        "strain": "",
        "sequencer": "Illumina",
        "author": "",
        "reviewer": "",
        "year": now.year,
        "month": now.month,
        "day": now.day,
        "kmer": "",
        "raw_data": bp2mbgb(table_data[1], "gb"),
        "clean_data": bp2mbgb(table_data[3], "gb"),
        "estimate_size": estimate_size,
        "heterozygosity": '{}%'.format(heterozygosity),
        "estisoft": software,
        "homogeneous": homogeneous,
        "assembly_size": bp2mbgb(table_assembly[6][3], "mb"),
        "findgse_genome": bp2mbgb(findgse_genome, "mb"),
        "scope_genome": bp2mbgb(scope_genome, "mb"),
        "average_genome": bp2mbgb(average_genome, "mb"),
        "max_genome": bp2mbgb(max_genome, "mb"),
        "min_genome": bp2mbgb(min_genome, "mb"),
        "table_data": table_data,
        "table_species": table_species,
        "pollution_description": poll_desc,
        "table_findgse": table_findgse,
        "table_scope": table_scope,
        "table_kmer": table_kmer,
        "assembly": assembly,
        "table_assembly": table_assembly,
        "scaffold_length": bp2mbgb(table_assembly[6][3], "bp"),
        "scaffold_n50": bp2mbgb(table_assembly[0][3], "bp"),
        "scaffold_number": bp2mbgb(table_assembly[6][4], "bp"),
        "contig_length": bp2mbgb(table_assembly[6][1], "bp"),
        "contig_n50": bp2mbgb(table_assembly[0][1], "bp"),
        "contig_number": bp2mbgb(table_assembly[6][2], "bp"),
        "depth_description": "",
        "software": {},
        "database": {}
    }

    r.update(config["general"])

    with open(survey_json) as fh:
        js = json.load(fh)
        for k in js:
            r[k].update(js[k])

    tpl = DocxTemplate(tpl_docx)

    tpl.render(r)
    
    if homogeneous:
        if assembly:
            figure_name = ["base_quality.png" ,"base_content.png", "base_gc.png", "top10_species.png", "findgse.png", "genomescope.png", "gc_depth.png"]
            figure_variable = [figure_base_quality, figure_base_content, figure_base_gc, figure_species, figure_findgse, figure_genomescope, figure_gc_depth]
        else:
            figure_name = ["base_quality.png", "base_content.png", "base_gc.png", "top10_species.png", "findgse.png", "genomescope.png", "gc_depth.png"]
            figure_variable = [figure_base_quality, figure_base_content, figure_base_gc, figure_species, figure_findgse, figure_genomescope]

    else:
        if assembly:
            figure_name = ["base_quality.png", "base_content.png", "base_gc.png", "top10_species.png", "heterozygosity.png", "gc_depth.png"]
            figure_variable = [figure_base_quality, figure_base_content, figure_base_gc, figure_species, figure_kmer, figure_gc_depth]
        else:
            figure_name = ["base_quality.png", "base_content.png", "base_gc.png", "top10_species.png", "heterozygosity.png"]
            figure_variable = [figure_base_quality, figure_base_content, figure_base_gc, figure_species, figure_kmer]
    
    for i, j in zip(figure_name, figure_variable):
        tpl.replace_pic(i, j)

    tpl.save(os.path.join(out_dir, "report.docx"))

#    html_report
    for i in ["images", "static"]:
        temp = os.path.join(out_dir, i)
        if os.path.exists(temp):
            shutil.rmtree(temp)
        shutil.copytree(os.path.join(tpl_html, i), temp)

    for i in read_picture(figure_variable):
        shutil.copy(i, os.path.join(out_dir, "images/"))

    for j in ["index.html", "main.html"]:
        temp = open(os.path.join(tpl_html, j)).read()
        if type(temp) == type(b''):
            temp = temp.decode('utf-8')
        tpl = Template(temp)

        with open(os.path.join(out_dir, j), "w") as fh:
            line = tpl.render(r)
            if type(line) == type(b''):
                line = line.decode('utf-8')
            try:
                fh.write(line)
            except:
                reload(sys)
                sys.setdefaultencoding("utf-8")
                fh.write(line)

#     html_report

    return r


def report(args):
    run_report(
        cfg=args.cfg,
        survey_json=args.survey,
        table_data=args.data,
        table_species=args.species,
        table_findgse=args.gse,
        table_scope=args.scope,
        table_kmer=args.kmer,
        table_assembly=args.assembly,
        figure_base_quality=args.quality,
        figure_base_content=args.content,
        figure_base_gc=args.gc,
        figure_species=args.species_png,
        figure_findgse=args.gse_png,
        figure_genomescope=args.scope_png,
        figure_kmer=args.kmer_png,
        figure_gc_depth=args.gc_depth,
        tpl_docx=args.docx,
        tpl_html=args.html,
        out_dir=args.out
    )


def add_report_args(parser):

    parser.add_argument("cfg", help="Input configuration file.")
    parser.add_argument("-o", "--out", metavar='FILE', type=str, required=True,
        help="Output result path.")

    json_group = parser.add_argument_group(title="Json", )
    json_group.add_argument("--survey", metavar='FILE', type=str, required=True,
        help="Input the jison file for survey analysis: survey.json.")

    table_group = parser.add_argument_group(title="Tables", )
    table_group.add_argument("--data", metavar='FILE', type=str, required=True,
        help="Input QC data statistics table: QC.xls.")
    table_group.add_argument("--species", metavar='FILE', type=str, required=True,
        help="Input species distribution statistics.")
    table_group.add_argument("--gse", metavar='FILE', type=str, default="no", 
        help="Input the FindGSE genome estimate file: findgse.txt.")
    table_group.add_argument("--scope", metavar='FILE', type=str, default="no", 
        help="Input the GenomeScope estimate file: genomescope.txt.")
    table_group.add_argument("--kmer", metavar='FILE', type=str, default="no", 
        help="Input the Kmerfreq estimate file(heterozygosity.xls).")
    table_group.add_argument("--assembly", metavar='FILE', type=str, default="no", 
        help="Input genomic statistics file(asm.tsv).")

    figure_group = parser.add_argument_group(title="Figure", )
    figure_group.add_argument("--quality", metavar='FILE', type=str, required=True,
        help="Input base quality map")
    figure_group.add_argument("--content", metavar='FILE', type=str, required=True,
        help="Input base content map")
    figure_group.add_argument("--gc", metavar='FILE', type=str, required=True,
        help="Input GC profile")
    figure_group.add_argument("--species_png", metavar='FILE', type=str, required=True,
        help="Input species distribution map")
    figure_group.add_argument("--gse_png", metavar='FILE', type=str, default="no",
        help="Input the result image of FindGSE.")
    figure_group.add_argument("--scope_png", metavar='FILE', type=str, default="no",
        help="Input the result image of GenomeScope.")
    figure_group.add_argument("--kmer_png", metavar='FILE', type=str, default="no",
        help="Input the result image of Kmerfreq.")
    figure_group.add_argument("--gc_depth", metavar='FILE', type=str, default="no",
        help="Input the GC depth image.")

    template_group = parser.add_argument_group(title="Template", )
    template_group.add_argument("--docx", metavar='FILE', type=str, required=True,
        help="Input the world template file.")
    template_group.add_argument("--html", required=True, help="Input the htlm template file.")

    return parser


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

    parser = add_report_args(parser)
    args = parser.parse_args()

    report(args)


if __name__ == "__main__":
    main()
