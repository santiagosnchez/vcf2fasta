#!/usr/bin/env python
# vcf2fasta.py

import argparse, textwrap
import gzip
import sys
from re import match

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="vcf2fasta.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Converts a regions in the genome into a FASTA alignment
        based on a VCF and a GFF file.""",
        epilog=textwrap.dedent('''\
examples:
'''))
    a1 = parser.add_argument(
    '--fasta', '-f', metavar='GENOME', type=str, required=True,
    help='FASTA file with the reference genome.')
    a2 = parser.add_argument(
    '--vcf', '-v', metavar='VCF', type=str, required=True,
    help='[VCF] the vcf file.')
    a3 = parser.add_argument(
    '--gff', '-g', metavar='GFF', type=str, required=True,
    help='individual FASTA records.')
    a4 = parser.add_argument(
    '--feat', '-e', metavar='FEAT', type=str,
    help='individual FASTA records.')
    parser.add_argument(
    '--blend', '-b', action="store_true", default=False,
    help='if sequences should be printed to screen.')
    args = parser.parse_args()
    genome = readfasta(args.fasta)
    print "\n"


# functions

def readfasta(file):
    data = {}
    seqn=0
    with open(file, "r") as f:
        for line in f:
            line = line.rstrip()
            if match("^>",line):
                head = line[1:]
                data[head] = ''
                seqn += 1
                sys.stdout.write("FASTA seqs read: %s\r" % seqn),
            else:
                data[head] += line
        return data

def readgff(file):


def readvcf(file):
    

def wrapseq(seq, w):
    chunks = []
    interval = map(lambda x: x*w, range((len(seq)/w)+2))
    for i in interval:
        if i != interval[-1]:
            chunks.append(seq[i:interval[interval.index(i)+1]-1])
    return("\n".join(chunks))

def ifkeyisfound(x, key):
    try:
        x[key]
        return True
    except KeyError:
        return False
  
if __name__ == '__main__':
    main()

