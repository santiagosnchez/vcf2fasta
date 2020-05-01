#!/usr/bin/env python
# vcf2fasta.py

import argparse
import re
import sys
import os
from string import maketrans
try:
    import pysam
except:
    sys.exit("pysam is not installed. Try: pip install pysam")

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="vcf2fasta.py",
        version="0.3",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
Converts regions/intervals in the genome into FASTA alignments
provided a VCF file, a GFF file, and FASTA reference.""",
        epilog="""
All files must be indexed. So before running the code make sure 
that your reference FASTA file is indexed:

samtools faidx genome.fas

BGZIP compress and TABIX index your VCF file:

bgzip variants.vcf
tabix variants.vcf.gz

The GFF file does not need to be indexed.

examples:
python vcf2fasta.py -f genome.fas -v variants.vcf.gz -g intervals.gff -e CDS

""")
    parser.add_argument(
    '--fasta', '-f', metavar='GENOME', type=str,
    help='FASTA file with the reference genome.')
    parser.add_argument(
    '--vcf', '-v', metavar='VCF', type=str, required=True,
    help='a tabix-indexed VCF file.')
    parser.add_argument(
    '--gff', '-g', metavar='GFF', type=str,
    help='GFF file.')
    parser.add_argument(
    '--feat', '-e', metavar='FEAT', type=str,
    help='feature/annotation in the GFF file. (i.e. gene, CDS, intron)')
    parser.add_argument(
    '--dir', '-d', metavar='DIRECTORY', type=str, default='.',
    help='output directory. (default: current directory)')
    #parser.add_argument(
    #'--blend', '-b', action="store_true", default=False,
    #help='concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)')
    args = parser.parse_args()

    # read GFF file
    gff = ReadGFF(args.gff)
    # read variant file and get samples
    vcf = pysam.VariantFile(args.vcf)
    samples = vcf.fetch().next().samples.keys()
    # read genome reference file
    ref = pysam.FastaFile(args.fasta)
    # get gene keys from GFF
    genes = gff.keys()
    # create output dir
    if not os.path.exists(args.path.exists(args.dir))
        os.mkdir(args.dir)
    for gene in genes:
        genename = gene+"."+gff[gene][0][3]+"-"+gff[gene][-1][4]
        o = open(args.dir+"/"+genename+".fas", "w")
        seq1 = ''
        seq2 = ''
        for gffrec in gff[gene]:
            chrom,start,end,strand = gffrec[0],int(gffrec[3])-1,int(gffrec[4]),gffrec[6]
            s = ref.fetch(chrom, start, end).upper()
            s1 = s
            s2 = s
            addposcum = 0
            for rec in vcf.fetch(chrom, start, end):
                alleles,addpos,refpos = UpdateAllele(rec)
                s1 = UpdateSeq(rec, alleles[0], addposcum, refpos, start, s1)                         
                s2 = UpdateSeq(rec, alleles[1], addposcum, refpos, start, s2)
                addposcum += addpos
            seq1 += s1
            seq2 += s2
        if strand == "-":
            seq1 = revcomp(seq1)
            seq2 = revcomp(seq2)
        i=len(seq1) % 3
        if i == 0:
            o.write(">"+genename+"1\n"+seq1+"\n>"+genename+"2\n"+seq2+"\n")
        else:
            o.write(">"+genename+"1\n"+seq1[:-i]+"\n>"+genename+"2\n"+seq2[:-i]+"\n")
        o.close()

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

def UpdateAllele(vcfrec):
    allelelen = [ len(x) for x in vcfrec.alleles ]
    maxlen = allelelen.index(max(allelelen))
    alleles2 = tuple([ x + '-' * (allelelen[maxlen]-len(x)) for x in vcfrec.alleles ])
    return alleles2, abs(allelelen[0]-allelelen[maxlen]), allelelen[0]-1

def UpdateSeq(vcfrec, allele, addposcum, refpos, start, seq):
    pos = vcfrec.pos-1-start+addposcum
    seq = seq[:pos]+allele+seq[pos+1+refpos:]
    return seq

def revcomp(seq):
    tt = maketrans('ACGT?N','TGCA?N')
    return seq[::-1].translate(tt)

if __name__ == "__main__":
    main()
