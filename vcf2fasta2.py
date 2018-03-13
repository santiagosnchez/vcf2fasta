#!/usr/bin/env python
# vcf2fasta.py

import argparse
import gzip
import sys
import os
import time
from re import match, search, sub
from string import maketrans

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="vcf2fasta.py",
        version="0.1",
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
    '--vcf', '-c', metavar='VCF', type=str, required=True,
    help='[VCF] the vcf file.')
    a3 = parser.add_argument(
    '--gff', '-g', metavar='GFF', type=str, required=True,
    help='individual FASTA records.')
    a4 = parser.add_argument(
    '--feat', '-e', metavar='FEAT', type=str, required=True,
    help='individual FASTA records.')
    parser.add_argument(
    '--blend', '-b', action="store_true", default=False,
    help='if sequences should be printed to screen.')
    args = parser.parse_args()
    print " {} v{}".format(parser.prog,parser.version)
    # read genome
    genome = readfasta(args.fasta)
    # read gff
    gff = readgff(args.gff)
    sys.stdout.write(" [readgff]   lines in GFF file: {}\n".format(len(gff)))
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
            data = {}
            notfound = []
            count = 0
            t1 = time.time()
            sys.stdout.write(" [vcf2fasta] {:<10s} {:<15s} {:<10s} {:<10s}\n".format("SNP","CHROM","POS","TIME"))
        else:
            var = var.rstrip().split("\t")
            if not keyisfound(gff, var[0]):
                notfound.append(var[0])
                next
            else:
                data = vcf2fasta(var, gff, genome, head, sampstart, feat, data, args.blend)
                count += 1
                t2 = time.time()
                tnow = t2-t1
                if tnow < 60.0:
                    tnows = "{:.2f}".format(tnow)
                    sys.stdout.write(" [vcf2fasta] {:<10d} {:<15s} {:<10s} {:<6s}seconds \r".format(count,var[0],var[1],tnows)),
                    sys.stdout.flush()
                elif  60.0 < tnow < 3600.0:
                    tnows = "{:.2f}".format(tnow/60.0)
                    sys.stdout.write(" [vcf2fasta] {:<10d} {:<15s} {:<10s} {:<6s}minutes \r".format(count,var[0],var[1],tnows)),
                    sys.stdout.flush()
                elif 3600.0 < tnow < 86400.0:
                    tnows = "{:.2f}".format(tnow/3600.0)
                    sys.stdout.write(" [vcf2fasta] {:<10d} {:<15s} {:<10s} {:<6s}hours   \r".format(count,var[0],var[1],tnows)),
                    sys.stdout.flush()
                else:
                    tnows = "{:.2f}".format(tnow/86400.0)
                    sys.stdout.write(" [vcf2fasta] {:<10d} {:<15s} {:<10s} {:<6s}days    \r".format(count,var[0],var[1],tnows)),
                    sys.stdout.flush()
    if len(notfound) > 0:
        print '\n [warning] {} variants were skipped; not found in GFF'.format(len(notfound))
    else:
        print ''
    # save data to disk
    genes = data.keys()
    reg = data[genes[0]].keys()
    stra = data[genes[0]][reg[0]].keys()
    samples = sorted(data[genes[0]][reg[0]][stra[0]].keys())
    count = 0
    if not os.path.exists(feat):
        os.makedirs(feat)
    for gene in genes:
        if args.blend:
            with open(feat+"/"+gene, 'w') as o:
                for indiv in samples:
                    seq = ''
                    o.write(">"+indiv+"\n")
                    for region in sorted([ int(i) for i in data[gene].keys() ]):
                        strand = data[gene][str(region)].keys()[0]
                        seq += data[gene][str(region)][strand][indiv]
                    if strand == '-':
                        seq = revcomp(seq)
                    o.write(seq+"\n")
            count += 1
        else:
            for region in sorted([ int(i) for i in data[gene].keys() ]):
                strand = data[gene][str(region)].keys()[0]
                with open(feat+"/"+gene+"."+str(region), 'w') as o:
                    for indiv in samples:
                        o.write(">"+indiv+"\n")
                        seq = data[gene][str(region)][strand][indiv]
                        if strand == '-':
                            seq = revcomp(seq)
                        o.write(seq+"\n")
                count += 1
    print " Done writing {0} files to {1}".format(count,feat)


# functions

def vcf2fasta(var, gff, genome, head, sampstart, feat, data, blend):
    genes = gff[var[0]][feat].keys()
    if blend:
        for gname in genes:
            gmin = min([ int(k[3]) for k in gff[var[0]][feat][gname] ])
            gmax = max([ int(k[4]) for k in gff[var[0]][feat][gname] ])
            if gmin <= int(var[1]) <= gmax:
                data = startdict(data,gname)
                for k in gff[var[0]][feat][gname]:
                    data[gname] = startdict(data[gname],k[3])
                    data[gname][k[3]] = startdict(data[gname][k[3]],k[6])
                    if int(k[3]) <= int(var[1]) <= int(k[4]):
                        for j in range(sampstart,len(head)):
                            if checkphase(var, j) == 'phased':
                                pheada = head[j]+"_a"
                                pheadb = head[j]+"_b"
                                if keyisfound(data[gname][k[3]][k[6]],pheada):
                                    data = updatealn(data, gname, k[3], pheada, var, j, k[6])
                                if keyisfound(data[gname][k[3]][k[6]],pheadb):
                                    data = updatealn(data, gname, k[3], pheadb, var, j, k[6])
                                if not keyisfound(data[gname][k[3]][k[6]],pheada) and not keyisfound(data[gname][k[3]][k[6]],pheadb):
                                    data[gname][k[3]][k[6]][pheada] = genome[var[0]][int(k[3])-1:int(k[4])]
                                    data[gname][k[3]][k[6]][pheadb] = genome[var[0]][int(k[3])-1:int(k[4])]
                                    data = updatealn(data, gname, k[3], pheada, var, j, k[6])
                                    data = updatealn(data, gname, k[3], pheadb, var, j, k[6])
                            else:
                                if keyisfound(data[gname][k[3]][k[6]],head[j]):
                                    data = updatealn(data, gname, k[3], head[j], var, j, k[6])
                                else:
                                    data[gname][k[3]][k[6]][head[j]] = genome[var[0]][int(k[3])-1:int(k[4])]
                                    data = updatealn(data, gname, k[3], head[j], var, j, k[6])
                    else:
                        for j in range(sampstart,len(head)):
                            if checkphase(var, j) == 'phased':
                                pheada = head[j]+"_a"
                                pheadb = head[j]+"_b"
                                if not keyisfound(data[gname][k[3]][k[6]],pheada) and not keyisfound(data[gname][k[3]][k[6]],pheadb):
                                    data[gname][k[3]][k[6]][pheada] = genome[var[0]][int(k[3])-1:int(k[4])]
                                    data[gname][k[3]][k[6]][pheadb] = genome[var[0]][int(k[3])-1:int(k[4])]
                            else:
                                if not keyisfound(data[gname][k[3]][k[6]],head[j]):
                                    data[gname][k[3]][k[6]][head[j]] = genome[var[0]][int(k[3])-1:int(k[4])]
        return data
    else:
        for gname in genes:
            for k in gff[var[0]][feat][gname]:
                if int(k[3]) <= int(var[1]) <= int(k[4]):
                    data = startdict(data,gname)
                    data[gname] = startdict(data[gname],k[3])
                    for j in range(sampstart,len(head)):
                        if checkphase(var, j) == 'phased':
                            pheada = head[j]+"_a"
                            pheadb = head[j]+"_b"
                            if keyisfound(data[gname][k[3]][k[6]],pheada):
                                data = updatealn(data, gname, k[3], pheada, var, j, k[6])
                            if keyisfound(data[gname][k[3]][k[6]],pheadb):
                                data = updatealn(data, gname, k[3], pheada, var, j, k[6])
                            if not keyisfound(data[gname][k[3]][k[6]],pheada) and not keyisfound(data[gname][k[3]][k[6]],pheadb):
                                data[gname][k[3]][k[6]][pheada] = genome[var[0]][int(k[3])-1:int(k[4])]
                                data[gname][k[3]][k[6]][pheadb] = genome[var[0]][int(k[3])-1:int(k[4])]
                                data = updatealn(data, gname, k[3], pheada, var, j, k[6])
                                data = updatealn(data, gname, k[3], pheadb, var, j, k[6])
                        else:
                            if keyisfound(data[gname][k[3]][k[6]],head[j]):
                                data = updatealn(data, gname, k[3], head[j], var, j, k[6])
                            else:
                                data[gname][k[3]][k[6]][head[j]] = genome[var[0]][int(k[3])-1:int(k[4])]
                                data = updatealn(data, gname, k[3], head[j], var, j, k[6])



def updatealn(data, gname, start, samp, var, ind, strand):
    seq = data[gname][start][strand][samp].upper()
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
            data[gname][start][strand][samp] = seq[:pos] + nvar[0] + seq[pos+1:]
        if search("_b$", samp):
            data[gname][start][strand][samp] = seq[:pos] + nvar[1] + seq[pos+1:]
        return data
    elif checkphase(var, ind) == 'unphased':
        nvar = [ alleles[i] for i in gt.split("/") ]
        if all([ i == nvar[0] for i in nvar ]):
            data[gname][start][strand][samp] = seq[:pos] + nvar[0] + seq[pos+1:]
        else:
            data[gname][start][strand][samp] = seq[:pos] + iupac(nvar) + seq[pos+1:]
        return data
    elif checkphase(var, ind) == 'haploid':
        data[gname][start][strand][samp] = seq[:pos] + alleles[gt] + seq[samp][pos+1:]
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

#def read

def readfasta(file):
    data = {}
    c = 0
    with open(file, 'r') as f:
        lines = f.readlines()
        lines = map(lambda x: x.rstrip(), lines)
        ihead = map(lambda i: lines.index(i), filter(lambda k: ">" in k, lines))
        for i in range(len(ihead)):
            if ihead[i] != ihead[-1]:
                data[lines[ihead[i]][1:]] = ''.join(lines[ihead[i]+1:ihead[i+1]])
            else:
                data[lines[ihead[i]][1:]] = ''.join(lines[ihead[i]+1:])
            c += 1
            sys.stdout.write(" [readfasta] reading FASTA sequence: {}\r".format(c)),
            sys.stdout.flush()
    print ""
    return data

def readgff(file):
    gff = {}
    c = 0
    with open(file, "r") as g:
        lines = g.readlines()
        # get rid of comments in the GFF
        lines = filter(lambda i: match('^((?!#).)*$',i), lines)
        lines = map(lambda i: i.rstrip().split("\t"), lines)
        sys.stdout.write(" [readgff]   Reading GFF file ...")
        for line in lines:
            gff = startdict(gff,line[0])
            gname = getgnames(line[-1])
            gff[line[0]] = startdict(gff[line[0]],line[2])
            if keyisfound(gff[line[0]][line[2]],gname):
                gff[line[0]][line[2]][gname] += [line]
            else:
                gff[line[0]][line[2]][gname] = [line]
            c += 1
    print "\n [readgff]   {} entries found".format(c)
    return gff

def getgnames(gname):
    m1 = search('\"(.+?)\" *;{0,1}', gname)
    m2 = search('= *\"(.+?)\" *;{0,1}', gname)
    m3 = search('= *(.+?) *;{0,1}', gname)
    m4 = search('^(.+?) *;{0,1}$', gname)
    if m1:
        gname = sub('\"| *;*$','',m1.group())
        return gname
    elif m2:
        gname = sub('= *|\"| *;*$','',m2.group())
        return gname
    elif m3:
        gname = sub('= *| *;*$','',m3.group())
        return gname
    elif m4:
        gname = sub(' *;*$','',m4.group())
        return gname
    else:
        return gname

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

def wrapseq(seq, w):
    chunks = []
    interval = map(lambda x: x*w, range((len(seq)/w)+2))
    for i in interval:
        if i != interval[-1]:
            chunks.append(seq[i:interval[interval.index(i)+1]-1])
    return("\n".join(chunks))

if __name__ == '__main__':
    main()

