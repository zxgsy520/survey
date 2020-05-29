#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import json
import argparse
import logging
import os.path
import shutil
from datetime import datetime

from jinja2 import Template
from docxtpl import DocxTemplate
try:
    from ConfigParser import ConfigParser
except:
    from configparser import ConfigParser

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
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


def read_config(cfg):
    """
    read config fron ini
    :param cfg:
    :return:
    """
    check_paths(cfg)

    r = {}
    config = ConfigParser()
    config.read(cfg)

    for section in config.sections():
        r[section] = {}

        for option in config.options(section):
            value = config.get(section, option).strip().decode("utf-8")
            r[section][option] = value

    return r


def read_tsv(file):

    r = []

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        r.append(line.split())

    return r


def read_table_data(file):

    return read_tsv(file)[0]


def run_report(table_data,
               figure_base_quality, figure_base_content, figure_base_gc,
               tpl_docx, out_dir
               ):
    
    out_dir = mkdir(out_dir)
    table_data = read_table_data(table_data)
    
    r = {
        "raw_data": table_data[2],
        "clean_data": table_data[4],
        "table_data": table_data
    }

    # docx report
    tpl = DocxTemplate(tpl_docx)

    tpl.render(r)

    for i, j in zip(["base_quality.png", "base_content.png", "gc_distribution.png"],
                    [figure_base_quality, figure_base_content, figure_base_gc]):
        tpl.replace_pic(i, j)

    tpl.save(os.path.join(out_dir, "report.docx"))

    return r


def ngs_qc_report(args):
    
    r = run_report(
        table_data=args.data,
        figure_base_quality=args.quality,
        figure_base_content=args.content,
        figure_base_gc=args.gc,
        tpl_docx=args.docx,
        out_dir=args.out
    )


def add_report_args(parser):

    parser.add_argument("--out", required=True, help="")

    table_group = parser.add_argument_group(title="Tables", )
    table_group.add_argument("--data", required=True, help="Input the statistical table of NGS data quality control.")

    figure_group = parser.add_argument_group(title="Figure", )
    figure_group.add_argument("--quality", required=True, help="Input base quality distribution map.")
    figure_group.add_argument("--content", required=True, help="Input base content distribution map")
    figure_group.add_argument("--gc", required=True, help="Input gc distribution map")

    template_group = parser.add_argument_group(title="Template", )
    template_group.add_argument("--docx", required=True, help="")
    
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

    ngs_qc_report(args)


if __name__ == "__main__":
    main()
