#!/usr/bin/env python
# vcf2fasta.py

import argparse
import re
import sys
import os
from string import maketrans
import multiprocessing as mp
try:
    import pysam
except:
    sys.exit("pysam is not installed. Try: pip install pysam")


def ReadGFF(file):
    gff = {}
    with open(file, "r") as f:
        for line in f:
            if re.match("^#", line):
                pass
            else:
                line = line.rstrip().split("\t")
                gname = GetGeneName(line[-1])
                if gff.get(gname):
                    gff[gname] += [line]
                else:
                    gff[gname] = [line]
    return gff

def GetGeneName(gname):
    m1 = re.search('\"(.+?)\" *;', gname)
    m2 = re.search('= *\"(.+?)\" *;', gname)
    m3 = re.search('= *(.+?) *;', gname)
    m4 = re.search('^(.+?) *;{0,1}$', gname)
    if m1:
        gname = re.sub('\"| *;*$','',m1.group())
        return gname
    elif m2:
        gname = re.sub('= *|\"| *;*$','',m2.group())
        return gname
    elif m3:
        gname = re.sub('.*= *| *;*$','',m3.group())
        return gname
    elif m4:
        gname = re.sub('.*= *| *;*$','',m4.group())
        return gname
    else:
        return gname

def UpdateAllele(vcfrec, start):
    allelelen = [ len(x) for x in vcfrec.alleles ]
    maxlen = allelelen.index(max(allelelen))
    alleles2 = tuple([ x + '-' * (allelelen[maxlen]-len(x)) for x in vcfrec.alleles ])
    return alleles2, abs(allelelen[0]-allelelen[maxlen]), allelelen[0]-1

def UpdateSeq(vcfrec, allele, addposcum, refpos, start, seq):
    pos = vcfrec.pos-1-start+addposcum
    seq = seq[:pos]+allele+seq[pos+1+refpos:]
    return seq




if __name__ == "__main__":
    seq1 = seq
    seq2 = seq
    addposcum = 0
    for rec in vcf.fetch(chrom, start, end):
        alleles,addpos,refpos = UpdateAllele(rec, start)
        seq1 = UpdateSeq(rec, alleles[0], addposcum, refpos, start, seq1)                         
        seq2 = UpdateSeq(rec, alleles[1], addposcum, refpos, start, seq2)
        addposcum += addpos
    print ">seq1\n"+seq1+"\n>seq2\n"+seq2
