#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Junpeng Fan",)
__email__ = "jpfan@whu.edu.cn"
__all__ = []


def read_tsv(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split("\t")


def filter_hsp_with_evalue_ident(hsps, evalue, pident, outfmt):

    LOG.info("filter with evalue <= %s and pident >= %s" % (evalue, pident))

    assert "evalue" in outfmt
    assert "pident" in outfmt

    for hsp in hsps:

        if len(hsp) != len(outfmt):
            LOG.warning("%r not match %s" % (hsp, outfmt))
            continue

        _evalue = float(hsp[outfmt["evalue"]])
        _pident = float(hsp[outfmt["pident"]])

        if _evalue <= evalue and _pident >= pident:
            yield hsp


def cluster_hsp(hsps, outfmt):

    assert "qseqid" in outfmt
    assert "sseqid" in outfmt

    pqseqid = psseqid = ""

    hit = []

    for hsp in hsps:

        _qseqid = hsp[outfmt["qseqid"]]
        _sseqid = hsp[outfmt["sseqid"]]

        if _qseqid != pqseqid or _sseqid != psseqid:

            if hit:
                yield hit
                pqseqid = _qseqid
                psseqid = _sseqid
                hit = [hsp]
            else:
                pqseqid = _qseqid
                psseqid = _sseqid
                hit = [hsp]
        else:
            hit.append(hsp)

    if hit:
        yield hit


def filter_hit_with_best_evalue(hits, best_evalue, outfmt):

    LOG.info("filter with evalue <= best_evalue*%s" % best_evalue)

    assert "qseqid" in outfmt

    pqseqid = ""
    pevalue = 0

    for hit in hits:

        _qseqid = hit[0][outfmt["qseqid"]]
        _evalue = float(hit[0][outfmt["evalue"]])

        if _qseqid != pqseqid:
            yield hit
            pqseqid = _qseqid
            pevalue = _evalue
        else:
            if _evalue <= pevalue*best_evalue:
                yield hit


def filter_hit_with_best(hits, outfmt):

    LOG.info("filter with best hit")

    assert "qseqid" in outfmt

    pqseqid = ""

    for hit in hits:

        _qseqid = hit[0][outfmt["qseqid"]]

        if _qseqid != pqseqid:
            yield hit
            pqseqid = _qseqid


def _overlap(range_list):

    r = []
    p = [None, None]

    for start, end in sorted(range_list, key=lambda d: min(d)):

        if start > end:
            start, end = end, start

        if p[0] is None:
            p = [start, end]
        elif start <= p[1]:
            p[1] = end
        else:
            r.append(p)
            p = [start, end]

    if p[0] is not None:
        r.append(p)

    return r


def calculate_cov(hit, outfmt):

    qhit = [[int(i[outfmt["qstart"]]), int(i[outfmt["qend"]])] for i in hit]
    shit = [[int(i[outfmt["sstart"]]), int(i[outfmt["send"]])] for i in hit]

    qcov = sum([i[1]-i[0]+1 for i in _overlap(qhit)])*100.0/int(hit[0][outfmt["qlen"]])
    scov = sum([i[1]-i[0] + 1 for i in _overlap(shit)]) * 100.0 / int(hit[0][outfmt["slen"]])

    return qcov, scov


def filter_hit_with_cov(hits, qcov, scov, outfmt):

    for hit in hits:
        _qcov, _scov = calculate_cov(hit, outfmt)

        if _qcov >= qcov and _scov >= scov:
            yield hit


def process_blastp(file, evalue, pident, best, best_evalue, qcov, scov, outfmt, out):

    out_fmt = {}

    if "std" in outfmt:
        i = outfmt.index("std")
        outfmt = outfmt[:i] + \
                 ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"] + outfmt[i+1:]

    for i, q in enumerate(outfmt):
        out_fmt[q] = i

    outfmt = out_fmt
    data = read_tsv(file)

    data = cluster_hsp(filter_hsp_with_evalue_ident(data, evalue=evalue, pident=pident, outfmt=outfmt), outfmt)

    if best_evalue:
        data = filter_hit_with_best_evalue(data, best_evalue, outfmt=outfmt)

    if best:
        data = filter_hit_with_best(data, outfmt=outfmt)

    if qcov or scov:
        data = filter_hit_with_cov(data, qcov=qcov, scov=scov, outfmt=outfmt)

    print("#%s" % "\t".join(out))

    for hit in data:
        r = []

        for n in out:
            s = []

            for hsp in hit:
                mes = hsp[outfmt[n]]

                if mes not in s:
                    s.append(mes)
                else:
                    if n in ["pident", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]:
                        s.append(mes)

            r.append(";".join(s))

        yield r


def set_args():
    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description="""
filter for blast results

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    args.add_argument("blast", help="")
    args.add_argument("--best", action="store_true",
                      help="only get the best hit")
    args.add_argument("--best_evalue", type=float, metavar="NUM", default="1",
                      help="remain evalue <= best hit evalue * NUM (default: 1)")
    args.add_argument("--min_pident", type=float, metavar="NUM", default=30)
    args.add_argument("--min_qcov", type=float, metavar="NUM", default=30,
                      help="min coverage of query aligned (default: 30)")
    args.add_argument("--min_scov", type=float, metavar="NUM", default=30,
                      help="min coverage of subject aligned (default: 30)")
    args.add_argument("--evalue", type=float, metavar="NUM", default=0.00001,
                      help="max evalue (default: 0.00001)")
    args.add_argument("--outfmt",
                      choices=["qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "bitscore", "staxid", "stitle", "std"],
                      nargs="+",
                      default=["qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "bitscore"],
                      help="outfmt of blast result")
    args.add_argument("--out",
                      choices=["qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "sstart", "send",
                               "slen", "evalue", "bitscore", "staxid", "stitle", "std"],
                      nargs="+",
                      default=["qseqid", "sseqid"],
                      help="output field of blast result (default: qseqid sseqid)"
                      )

    return args.parse_args()


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    args = set_args()

    for i in process_blastp(args.blast, args.evalue, args.min_pident, args.best, args.best_evalue, args.min_qcov, args.min_scov, args.outfmt, args.out):
        print("\t".join(i))


if __name__ == "__main__":
    main()

