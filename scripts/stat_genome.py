#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield seq[0], seq[1]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip("@").split()[0]
            continue
        if line.startswith('@'):
            line = line.strip("@").split()[0]
            seq = seq.split('\n')

            yield seq[0], seq[1]
            seq = ''
            seq = "%s\n" % line
        else:
            seq += "%s\n" % line

    if len(seq.split('\n'))==5:
        seq = seq.split('\n')
        yield seq[0], seq[1]


def store_seq(file):

    data = {}
    
    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for seq_id, seq in fh:
        data[seq_id] = seq
    return data


def split_scaffold(id_sca,seq_sca,mingap,data_ctg,data_gap):
    first_seq,first_gap,end = 0,0,0
    check = False
    number = len(seq_sca)
    for i in range(number):
        item = seq_sca[i]
        if i == number - 1:
            if not check:
                id_ctg = id_sca + '_' + str(first_seq + 1 ) + '_' + str(i - first_seq + 1)
                data_ctg[id_ctg] = seq_sca[first_seq:]
            else:
                if item == 'N' or item == 'n':
                    last_gap = i
                else:
                    last_gap = i - 1
                if len(seq_sca[first_gap:last_gap + 1]) >= mingap:
                    id_ctg = id_sca + '_' + str(first_seq + 1) + '_' + str(first_gap - first_seq)
                    data_ctg[id_ctg] = seq_sca[first_seq:first_gap]
                    id_gap = id_sca + '_' + str(first_gap + 1) + '_' + str(last_gap + 1 - first_gap)
                    data_gap[id_gap] = seq_sca[first_gap:i]
                    if seq_sca[last_gap + 1:]:
                        id_ctg = id_sca + '_' + str(last_gap + 2) + '_' + str(i - last_gap)
                        data_ctg[id_ctg] = seq_sca[last_gap + 1:]
                else:
                    id_ctg = id_sca + '_' + str(first_seq + 1) + '_' + str(i + 1 - first_seq)
                    data_ctg[id_ctg] = seq_sca[first_seq:]

        if item == 'N' or item == 'n':
            if not check:
                first_gap = i
                check = True
        else:
            if check:
                if len(seq_sca[first_gap:i]) >= mingap:
                    #print mingap,len(seq_sca[first_gap:i])
                    #end = i
                    if first_gap - first_seq > 0:
                        id_ctg = id_sca + '_' + str(first_seq + 1) + '_' + str(first_gap - first_seq)
                        data_ctg[id_ctg] = seq_sca[first_seq:first_gap]
                    id_gap = id_sca + '_' + str(first_gap + 1) + '_' + str(i - first_gap)
                    data_gap[id_gap] = seq_sca[first_gap:i]
                    first_seq = i
                check = False


def extract_contig(data_sca,len_gap):
    data_ctg,data_gap = {},{}
    for id_sca in data_sca:
        split_scaffold(id_sca,data_sca[id_sca],len_gap,data_ctg,data_gap)
    return data_ctg,data_gap

def total(data):
    data_length = []
    num_total,num_100,num_1k,num_2k = 0,0,0,0
    seq_total,seq_100,seq_1k,seq_2k = '','','',''
    for key in data:
        seq = data[key]
        length = len(seq)
        data_length.append(length)
        if length >= 1000:
            num_100 += 1
            seq_100 += seq
        if length >= 2000:
                num_1k += 1
                seq_1k += seq
        if length >= 5000:
            num_2k += 1
            seq_2k += seq
        num_total += 1
        seq_total += seq
    return data_length,[len(seq_total),len(seq_100),len(seq_1k),len(seq_2k)],[num_total,num_100,num_1k,num_2k]

def sum_quality(total,data):
    length,longest,longest_num = 0,0,0
    data_length,data_number = [],[]
    data_check = {0:False,0.5:True,0.6:True,0.7:True,0.8:True,0.9:True}
    data_per = [0,0.5,0.6,0.7,0.8,0.9]
    number = len(data)
    data_sort = sorted(data,reverse=True)
    for i in range(number):
        if data_sort[i] >= longest:
            longest = data_sort[i]
            longest_num += 1
        length += data_sort[i]
        for j in range(1,6):
            per_now = data_per[j]
            per_fore = data_per[j -1]
            check_now = data_check[per_now]
            check_fore = data_check[per_fore]
            if length >= (total * per_now) and check_now and not check_fore:
                data_length.append(data_sort[i])
                data_number.append(i + 1)
                data_check[per_now] = False
    data_length.append(longest)
    data_number.append(longest_num)
    return data_length,data_number

def qualify(data):
    data_length,length_sum,number_sum = total(data)
    length_total = length_sum[0]
    sum_length,sum_number = sum_quality(length_total,data_length)
    result_length = sum_length + length_sum
    result_number = sum_number + number_sum
    return result_length,result_number

def result_out(file,length_sca,number_sca,length_ctg,number_ctg,length_gap,number_gap, gap):
    out = '#%13s\t%-15s\t%-15s\t%13s\t%13s' % ('StatType','ContigLength','ContigNumber','ScaffoldLength','ScaffoldNumber')
    if gap:
        out += '\t%13s\t%13s\n' % ('GapLength','GapNumber')
    else:
        out += '\n'
    type_stat = ['N50','N60','N70','N80','N90','Longest','Total','Length>=1kb','Length>=2kb','Length>=5kb']
    for i in range(10):
        out += '%13s\t%-15d\t%-15d\t%13d\t%13d\t' % (type_stat[i],length_ctg[i],number_ctg[i],length_sca[i],number_sca[i])
        if gap:
            out += '\t%13d\t%13d\n' % (length_gap[i],number_gap[i])
        else:
            out += '\n'
    #print >> sys.stdout, out
    f = open(file,'w')
    f.write(out)
    f.close()

def result_contig(file,data):
    def sort_seq(seq_one):
        seq = ''
        for i in range(len(data[id_ctg])):
            seq += data[id_ctg][i]
            if (i + 1)%100 == 0:
                seq += '\n'
        return seq

    out = ''
    for id_ctg in sorted(data):
        seq= sort_seq(data[id_ctg])
        out += '>' + id_ctg + '\n' + seq.strip('\n') + '\n'
    f = open(file,'w')
    f.write(out)
    f.close()

def main(args):

    dir_current = os.getcwd() + '/'
    name_genome = os.path.splitext(os.path.basename(args.scaffold))[0]
    if args.result:
        file_result = args.result
    else:
        file_result = dir_current + name_genome + '_result.xls'

    data_sca = store_seq(args.scaffold)
    sys.stderr.write('Storing is over!\n')
    data_ctg,data_gap = extract_contig(data_sca,args.minlen)
    sys.stderr.write('Contigs have been extracted!\n')
    length_sca,number_sca = qualify(data_sca)
    sys.stderr.write('Sacffolds have been qulified!\n')
    #print 'ctg',data_ctg.keys()
    length_ctg,number_ctg = qualify(data_ctg)
    sys.stderr.write('Contigs have been qulified!\n')
    if data_gap:
        length_gap,number_gap = qualify(data_gap)
        gap = True
    else:
        length_gap,number_gap = [], []
        gap = False
    result_out(file_result,length_sca,number_sca,length_ctg,number_ctg,length_gap,number_gap, gap)
    if args.contig:
        file_contig = dir_current + name_genome + '_contig_gap' + str(args.minlen) + '.fa'
        result_contig(file_contig,data_ctg)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Stat N50-N90 value of contigs and scaffolds.Editor:TY''')
    parser.add_argument('-s','--scaffold',metavar='fasta',required=True,help='Sequences of scaffold in *.fa/*.fasta format.')
    parser.add_argument('-m','--minlen',metavar='int',type=int,default=1,help='Define minimum length of gap.Default: 1')
    parser.add_argument('-r','--result',metavar='file',help='Define file storing results.Default:$genome_result.xls')
    parser.add_argument('-c','--contig',action='store_true',default=False,help='Define whether store contig sequences into file($genome_contig_gap*.fa).Default:False')
    args = parser.parse_args()
    #dir_current = os.getcwd() + '/'
    #name_genome = os.path.splitext(os.path.basename(args.scaffold))[0]
    #if args.result:
    #    file_result = args.result
    #else:
    #    file_result = dir_current + name_genome + '_result.xls'
    
    main(args)

