#!/usr/bin/env python
# vcf2fasta.py

import argparse
import gzip
import sys
import os
from re import match, search, sub
from string import maketrans

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="vcf2fasta.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
Converts regions in the genome into FASTA alignments
provided a VCF and a GFF file.""",
        epilog="""
The VCF file can be gzipped or not. All the other files are
expected to be uncompressed.

examples:
python vcf2fasta.py -f genome.fas -v variants.vcf -g regions.gff -f CDS --blend

""")
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
    '--feat', '-e', nargs="*", metavar='FEAT', type=str, required=True,
    help='individual FASTA records.')
    parser.add_argument(
    '--blend', '-b', action="store_true", default=False,
    help='if sequences should be printed to screen.')
    args = parser.parse_args()
    # read genome
    genome = readfasta(args.fasta)
    # read gff
    gff = readgff(args.gff)
    print "Lines in GFF file: %s" % len(gff)
    # get features
    feat = args.feat
    # parse vcf file
    if search(".vcf.gz$", args.vcf):
        v = gzip.open(args.vcf, "r")
    else:
        v = open(args.vcf, "r")
    for var in v:
        if match("##", var):
            next
        elif match("#CHROM", var):
            head = var.rstrip().split("\t")
            sampstart = head.index('FORMAT')+1
            samples = head[sampstart:]
            data = {}
            notfound = []
            count = 0
        else:
            var = var.rstrip().split("\t")
            if not keyisfound(gff, var[0]):
                notfound.append(var[0])
                next
            else:
                data = vcf2fasta(var, gff, genome, head, sampstart, feat[0], data)
                count += 1
                sys.stdout.write("Reading SNP %s: %s:%s             \r" % (count,var[0],var[1])),
                sys.stdout.flush()
    # save data to disk
    genes = data.keys()
    print ""
    for fe in feat:
        if not os.path.exists(fe):
            os.makedirs(fe)
        for gene in genes:
            if args.blend:
                with open(fe+"/"+gene) as o:
                    for indiv in samples:
                        o.write(">"+indiv+"\n")
                        for region in sorted([ int(i) for i in data[gene][fe].keys() ]):
                            strand = data[gene][fe][str(region)].keys()[0]
                            seq += data[gene][fe][str(region)][strand][indiv]
                        if strand == '-':
                            seq = revcomp(seq)
                        o.write(seq+"\n")
                    sys.stdout.write("Writing: "+fe+"/"+gene+"\r"),
                    sys.stdout.flush()
            else:
                for region in sorted([ int(i) for i in data[gene][fe].keys() ]):
                    strand = data[gene][fe][str(region)].keys()[0]
                    with open(fe+"/"+gene+"."+region) as o:
                        for indiv in samples:
                            o.write(">"+indiv+"\n")
                            seq = data[gene][fe][str(region)][strand][indiv]
                            if strand == '-':
                                seq = revcomp(seq)
                            o.write(seq+"\n")
                    sys.stdout.write("Writing: "+fe+"/"+gene+"."+region+"\r"),
                    sys.stdout.flush()
    print "Done"


# functions

def vcf2fasta(var, gff, genome, head, sampstart, feat, data):
    start = [ gff[var[0]][i][3] for i in range(len(gff[var[0]])) ]
    end = [ gff[var[0]][i][4] for i in range(len(gff[var[0]])) ]
    features = [ gff[var[0]][i][2] for i in range(len(gff[var[0]])) ]
    strand = [ gff[var[0]][i][6] for i in range(len(gff[var[0]])) ]
    gname = [ gff[var[0]][i][8] for i in range(len(gff[var[0]])) ]
    for i in zip(features,start,end,strand,gname):
        if i[0] == feat:
            if int(i[1]) <= int(var[1]) <= int(i[2]):
                gname = getgnames(i[-1])
                data = startdict(data,gname)
                data[gname] = startdict(data[gname],feat)
                data[gname][feat] = startdict(data[gname],i[1])
                data[gname][feat][i[1]] = startdict(data[gname],i[3])
                for j in range(sampstart,len(head)):
                    if checkphase(var, j) == 'phased':
                        pheada = head[j]+"_a"
                        pheadb = head[j]+"_b"
                        if keyisfound(data[gname][feat][i[1]][i[3]],pheada):
                            data = updatealn(data, gname, feat, i[1], pheada, var, j, i[3])
                        elif keyisfound(data[gname][feat][i[1]][i[3]],pheadb):
                            data = updatealn(data, gname, feat, i[1], pheadb, var, j, i[3])
                        else:
                            data[gname][feat][i[1]][i[3]][pheada] = genome[var[0]][int(i[1])-1:int(i[2])]
                            data[gname][feat][i[1]][i[3]][pheadb] = genome[var[0]][int(i[1])-1:int(i[2])]
                            data = updatealn(data, gname, feat, i[1], pheada, var, j, i[3])
                            data = updatealn(data, gname, feat, i[1], pheadb, var, j, i[3])
                    else:
                        if keyisfound(data[gname][feat][i[1]][i[3]],head[j]):
                            data = updatealn(data, gname, feat, i[1], head[j], var, j, i[3])
                        else:
                            data[gname][feat][i[1]][i[3]][head[j]] = genome[var[0]][int(i[1])-1:int(i[2])]
                            data = updatealn(data, gname, feat, i[1], head[j], var, j, i[3])
    return data

def updatealn(data, gname, feat, start, samp, var, ind, strand):
    seq = data[gname][feat][start][strand][samp].upper()
    pos = int(var[1])-int(start)
    alleles = {'0':var[3], '.':'?'}
    j = 0
    if search(',', var[4]):
        alt = var[4].split(',')
        for i in alt:
            j += 1
            alleles[str(j)] = i
    else:
        alleles['1'] = var[4]
    gt = var[ind].split(":")[0]
    if checkphase(var, ind) == 'phased':
        nvar = [ alleles[i] for i in gt.split("|") ]
        if search("_a$", samp):
            data[gname][feat][start][strand][samp] = seq[:pos] + nvar[0] + seq[pos+1:]
        elif search("_b$", samp):
            data[gname][feat][start][strand][samp] = seq[:pos] + nvar[1] + seq[pos+1:]
        return data
    elif checkphase(var, ind) == 'unphased':
        nvar = [ alleles[i] for i in gt.split("/") ]
        if all([ i == nvar[0] for i in nvar ]):
            data[gname][feat][start][strand][samp] = seq[:pos] + nvar[0] + seq[pos+1:]
        else:
            data[gname][feat][start][strand][samp] = seq[:pos] + iupac(nvar) + seq[pos+1:]
        return data
    elif checkphase(var, ind) == 'haploid':
        data[gname][feat][start][strand][samp] = seq[:pos] + alleles[gt] + seq[samp][pos+1:]
        return data

def checkphase(var, ind):
    gt = var[ind].split(":")[0]
    if search('|', gt):
        return 'phased'
    elif search('/', gt):
        return 'unphased'
    else:
        return 'haploid'

def iupac(nvar):
    myiupac = {
    'AT':'W','TA':'W',
    'AC':'M','CA':'M',
    'AG':'R','GA':'R',
    'TC':'Y','CT':'Y',
    'TG':'K','GT':'K',
    'GC':'S','CG':'S',
    }
    return myiupac[''.join(nvar)]

def revcomp(seq):
    tt = maketrans('ACGT?N','TGCA?N')
    return seq[::-1].translate(tt)

def getgnames(gname):
    m1 = search('\"(.+?)\" *;{0,1}', gname)
    m2 = search('= *\"(.+?)\" *;{0,1}', gname)
    m3 = search('= *(.+?) *;{0,1}', gname)
    m4 = search('^(.+?) *;{0,1}$', gname)
    if m1:
        name = sub('\"| *;*$','',m1.group())
    elif m2:
        name = sub('= *|\"| *;*$','',m2.group())
    elif m3:
        name = sub('= *| *;*$','',m3.group())
    elif m4:
        name = sub(' *;*$','',m4.group())
    return name


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
                sys.stdout.flush()
            else:
                data[head] += line
        print ""
        return data

def readgff(file):
    gff = {}
    with open(file, "r") as g:                        
        for line in g:
            if "#" in line:
                next
            else:
                line = line.rstrip().split("\t")
                if keyisfound(gff,line[0]):
                    gff[line[0]].append(line)
                else:
                    gff[line[0]] = []
                    gff[line[0]].append(line)
    return gff

def wrapseq(seq, w):
    chunks = []
    interval = map(lambda x: x*w, range((len(seq)/w)+2))
    for i in interval:
        if i != interval[-1]:
            chunks.append(seq[i:interval[interval.index(i)+1]-1])
    return("\n".join(chunks))

def startdict(d, y):
    if keyisfound(d, y):
        return d
    else:
        d[y] = {}
        return startdict(d, y)

def keyisfound(x, key):
    try:
        x[key]
        return True
    except KeyError:
        return False
  
if __name__ == '__main__':
    main()

