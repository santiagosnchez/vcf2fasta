#!/usr/bin/env python
# vcf2fasta.py

import argparse
import gzip
import sys
import os
import time
from re import match, search, sub
from string import maketrans
from multiprocessing import Pool

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
python vcf2fasta.py -f genome.fas -c variants.vcf -g regions.gff -e CDS -b
python vcf2fasta.py -f genome.fas -c variants.vcf -g regions.gff -e gene

""")
    parser.add_argument(
    '--fasta', '-f', metavar='GENOME', type=str,
    help='FASTA file with the reference genome.')
    parser.add_argument(
    '--vcf', '-c', metavar='VCF', type=str, required=True,
    help='[VCF] the vcf file.')
    parser.add_argument(
    '--gff', '-g', metavar='GFF', type=str,
    help='individual FASTA records.')
    parser.add_argument(
    '--feat', '-e', metavar='FEAT', type=str,
    help='individual FASTA records.')
    parser.add_argument(
    '--blend', '-b', action="store_true", default=False,
    help='concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)')
    #parser.add_argument(
    #'--wrap', '-w', const=100, nargs="?", type=int, default=False, metavar="N",
    #help='sequences will be wrapped every N characters. (default: 100)')
    args = parser.parse_args()
    print " {} v{}".format(parser.prog,parser.version)
    # read genome
    if args.fasta:
        genome = readfasta(args.fasta)
    # read gff
    if args.gff and args.feat:
        gff = readgff(args.gff, args.feat)
    # parse vcf file
    if args.fasta and args.gff and args.feat:
        variants,samples,phase = readvcf(args.vcf, gff)
        # convert to fasta
        vcf2fasta(variants, gff, genome, args.feat, args.blend, samples, phase)
    else:
        variants,samples,phase = readvcf(args.vcf, args.gff)
        vcf2fasta(variants, args.gff, args.fasta, args.feat, args.blend, samples, phase, args.vcf)
    sys.exit()

# functions

def vcf2fasta(data, gff, genome, feat, blend, samples, phase, outfile=None):
    if gff:
        chrom = gff.keys()
        gene_gff = {}
        for ch in chrom:
            genes = gff[ch].keys()
            for g in genes:
                gene_gff[g] = gff[ch][g]
        genes = data.keys()
        if phase == 'phased':
            samples = sorted(reduce(lambda x,y: x+y, map(lambda i: [i+'_a',i+'_b'], samples)))
        else:
            samples.sort()
        if not os.path.exists(feat):
            os.makedirs(feat)
        sys.stdout.write(" {:<15s}writing FASTA files ...".format("[vcf2fasta]"))
        c = 0
        for g in genes:
            if blend:
                o = open(feat+"/"+g+".fas", "w")
                concat = []
            entries = gene_gff[g]
            pos = data[g].keys()
            if len(entries) > 1:
                for e in entries:
                    tmpdata = {}
                    inregion = filter(lambda x: int(e[3]) <= x <= int(e[4]), pos)
                    if inregion:
                        for s in samples:
                            tmpdata[s] = genome[e[0]][int(e[3])-1:int(e[4])]
                        for i in inregion:
                            p = i-int(e[3])
                            for s in samples:
                                tmpdata[s] = tmpdata[s][:p] + data[g][i][s] + tmpdata[s][p+1:]
                    else:
                        if blend:
                            for s in samples:
                                tmpdata[s] = genome[e[0]][int(e[3])-1:int(e[4])]
                    if blend:
                        concat.append(tmpdata)
                    else:
                        with open(feat+"/"+g+"."+e[3]+"-"+e[4]+".fas", "w") as o:
                            for s in samples:
                                o.write(">"+s+"\n")
                                if e[6] == '-':
                                    o.write(revcomp(tmpdata[s])+"\n")
                                else:
                                    o.write(tmpdata[s]+"\n")
                if blend:
                    for s in samples:
                        o.write(">"+s+"\n")
                        seq = ''
                        for e in concat:
                            seq += e[s]
                        if entries[0][6] == '-':
                            o.write(revcomp(seq)+"\n")
                        else:
                            o.write(seq+"\n")
                    c += 1
            else:
                e = entries[0]
                tmpdata = {}
                for s in samples:
                    tmpdata[s] = genome[e[0]][int(e[3])-1:int(e[4])]
                for i in pos:
                    p = i-int(e[3])
                    for s in samples:
                        tmpdata[s] = tmpdata[s][:p] + data[g][i][s] + tmpdata[s][p+1:]
                with open(feat+"/"+g+".fas", "w") as o:
                    for s in samples:
                        o.write(">"+s+"\n")
                        if e[6] == '-':
                            o.write(revcomp(tmpdata[s])+"\n")
                        else:
                            o.write(tmpdata[s]+"\n")
                    c += 1
        sys.stdout.write("\n {:<15s}done writing {} files to {}\n".format("[vcf2fasta]",c,feat))
    else:
        if phase == 'phased':
            samples = sorted(reduce(lambda x,y: x+y, map(lambda i: [i+'_a',i+'_b'], samples)))
        else:
            samples.sort()
        out = outfile[0:outfile.index(".vcf")]+".fasta"
        sys.stdout.write(" {:<15s}writing FASTA as {}".format("[vcf2fasta]",out))
        with open(out, "w") as o:
            for ind in samples:
                o.write(">"+ind+"\n")
                o.write(wrapseq(data[ind],200)+"\n")
        sys.stdout.write("\n {:<15s}done writing FASTA\n".format("[vcf2fasta]"))

def getgene(var, gff):
    gff = gff[var[0]]
    genes = gff.keys()
    def isingene(x, p, g):
        gmin = min([ int(i[3]) for i in g ])
        gmax = max([ int(i[4]) for i in g ])
        if gmin <= p <= gmax:
            return True
        else:
            return False
    pos = [ int(i) for i in [str(var[1])] * len(genes) ]
    gff = [ gff[i] for i in genes ]
    lookup = map(isingene, genes, pos, gff)
    if not all([lookup[0] == i for i in lookup]):
        return genes[ lookup.index(True) ]

def extractvcf(data, var, head, sampstart, samples, phase, gff):
    if keyisfound(gff, var[0]):
        gene = getgene(var, gff)
        if gene:
            data = startdict(data, gene)
            pos = int(var[1])
            data[gene] = startdict(data[gene], pos)
            alleles = {'0':var[3], '.':'?'}
            j = 0
            if search(',', var[4]):
                alt = var[4].split(',')
                for i in alt:
                    j += 1
                    alleles[str(j)] = i
            else:
                alleles['1'] = var[4]
            gt = map(lambda x: var[x].split(":")[0], range(sampstart,len(head)))
            if phase == 'phased':
                for i in range(len(gt)):
                    nvar = map(lambda x: alleles[x], gt[i].split("|"))
                    data[gene][pos][samples[i]+"_a"] = nvar[0]
                    data[gene][pos][samples[i]+"_b"] = nvar[1]
                return data
            elif phase == 'unphased':
                for i in range(len(gt)):
                    nvar = map(lambda x: alleles[x], gt[i].split("/"))
                    if all([ nvar[0] == k for k in nvar ]):
                        data[gene][pos][samples[i]] = nvar[0]
                    else:
                        data[gene][pos][samples[i]] = iupac(nvar)
                return data
            elif phase == 'haploid':
                for i in range(len(gt)):
                    data[gene][pos][samples[i]] = alleles[gt[i]]
                return data
        else:
            return data
    else:
        return data

def extractvcf2(data, var, head, sampstart, phase):
    alleles = {'0':var[3], '.':'?'}
    j = 0
    if search(',', var[4]):
        alt = var[4].split(',')
        for i in alt:
            j += 1
            alleles[str(j)] = i
    else:
        alleles['1'] = var[4]
    for ind in range(sampstart,len(head)):
        gt = var[ind].split(":")[0]
        if phase == 'phased':
            nvar = map(lambda x: alleles[x], gt.split("|"))
            if data.get(head[ind]+"_a"):
                data[head[ind]+"_a"] += nvar[0]
            else:
                data[head[ind]+"_a"] = nvar[0]
            if data.get(head[ind]+"_b"):
                data[head[ind]+"_b"] += nvar[1]
            else:
                data[head[ind]+"_b"] = nvar[1]
        elif phase == 'unphased':
            nvar = map(lambda x: alleles[x], gt.split("/"))
            if data.get(head[ind]):
                if all([ nvar[0] == k for k in nvar ]):
                    data[head[ind]] += nvar[0]
                else:
                    data[head[ind]] += iupac(nvar)
            else:
                if all([ nvar[0] == k for k in nvar ]):
                    data[head[ind]] = nvar[0]
                else:
                    data[head[ind]] = iupac(nvar)
        elif phase == 'haploid':
            if data.get(head[ind]):
                data[head[ind]] += alleles[gt]
            else:
                data[head[ind]] = alleles[gt]
    return data

def readvcf(file, gff):
    data = {}
    sys.stdout.write(" {:<15s}reading VCF file ...".format("[readvcf]"))
    if search(".vcf.gz$", file):
        v = gzip.open(file, "r")
    else:
        v = open(file, "r")
    c = 0
    for var in v:
            if match("^##",var):
                next
            elif match("^#CHROM",var):
                head = var.rstrip().split("\t")
                sampstart = head.index('FORMAT')+1
                samples = head[sampstart:]
                sys.stdout.write("\n {:<15s}header has {} samples\n".format("[readvcf]",len(samples)))
            else:
                var = var.rstrip().split("\t")
                module = "[extractvcf]"
                c += 1
                if c == 1:
                    phase = checkphase(var, sampstart)
                    if phase == 'phased':
                        sys.stdout.write(" {:<15s}VCF is {}, splitting into \"a\" and \"b\" haplotypes\n".format("[readvcf]",phase))
                    elif phase == 'unphased':
                        sys.stdout.write(" {:<15s}VCF is {}, heterozygotes will be masked with IUPAC ambiguity codes\n".format("[readvcf]",phase))
                    else:
                        sys.stdout.write(" {:<15s}VCF ploidy looks {}\n".format("[readvcf]",phase))
                    sys.stdout.write(" {:<15s}{:<10s} {:<15s} {:<10s} {:<10s} {:<10s}\n".format(module,"SNP","CHROM","POS","TIME(T)","TIME(P)"))
                    t1 = time.time()
                t2 = time.time()
                if gff:
                    data = extractvcf(data, var, head, sampstart, samples, phase, gff)
                else:
                    data = extractvcf2(data, var, head, sampstart, phase)
                tpar = time.time()-t2
                tnow = time.time()-t1
                if tnow < 60.0:
                    tnows = "{:.2f}".format(tnow)
                    tpars = "{:.6f}".format(tpar)
                    sys.stdout.write(" {:<15s}{:<10d} {:<15s} {:<10s} {:<10s} {:<6s} \r".format(module,c,var[0],var[1],tnows+" s",tpars+" s")),
                    sys.stdout.flush()
                elif  60.0 < tnow < 3600.0:
                    tnows = "{:.2f}".format(tnow/60.0)
                    tpars = "{:.6f}".format(tpar)
                    sys.stdout.write(" {:<15s}{:<10d} {:<15s} {:<10s} {:<10s} {:<6s} \r".format(module,c,var[0],var[1],tnows+" m",tpars+" s")),
                    sys.stdout.flush()
                elif 3600.0 < tnow < 86400.0:
                    tnows = "{:.2f}".format(tnow/3600.0)
                    tpars = "{:.6f}".format(tpar)
                    sys.stdout.write(" {:<15s}{:<10d} {:<15s} {:<10s} {:<10s} {:<6s} \r".format(module,c,var[0],var[1],tnows+" h",tpars+" s")),
                    sys.stdout.flush()
                else:
                    tnows = "{:.2f}".format(tnow/86400.0)
                    tpars = "{:.6f}".format(tpar)
                    sys.stdout.write(" {:<15s}{:<10d} {:<15s} {:<10s} {:<10s} {:<6s} \r".format(module,c,var[0],var[1],tnows+" d",tpars+" s")),
                    sys.stdout.flush()
    v.close()
    print ""
    return data,samples,phase

def readfasta(file):
    data = {}
    c = 0
    with open(file, 'r') as f:
        lines = f.readlines()
        lines = map(lambda x: x.rstrip(), lines)
        ihead = map(lambda i: lines.index(i), filter(lambda k: ">" in k, lines))
        for i in range(len(ihead)):
            if ihead[i] != ihead[-1]:
                data[lines[ihead[i]][1:]] = ''.join(lines[ihead[i]+1:ihead[i+1]]).upper()
            else:
                data[lines[ihead[i]][1:]] = ''.join(lines[ihead[i]+1:]).upper()
            c += 1
            sys.stdout.write(" {:<15s}reading FASTA sequence: {}\r".format("[readfasta]",c)),
            sys.stdout.flush()
    print ""
    return data

def readgff(file, feat):
    gff = {}
    c = 0
    sys.stdout.write(" {:<15s}reading GFF file ...".format("[readgff]"))
    with open(file, "r") as g:
        lines = g.readlines()
        # get rid of comments in the GFF
        lines = filter(lambda i: match('^((?!#).)*$',i), lines)
        lines = map(lambda i: i.rstrip().split("\t"), lines)
        lines = filter(lambda i: i[2] == feat, lines)
        for line in lines:
            gff = startdict(gff,line[0])
            gname = getgnames(line[-1])
            if keyisfound(gff[line[0]],gname):
                gff[line[0]][gname] += [line]
            else:
                gff[line[0]][gname] = [line]
            c += 1
    print "\n {:<15s}{} {} entries found".format("[readgff]",c,feat)
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

def checkphase(var, ind):
    gt = var[ind].split(":")[0]
    if search('\|', gt):
        return 'phased'
    elif search('\/', gt):
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

def wrapseq(seq, w):
    chunks = []
    interval = map(lambda x: x*w, range((len(seq)/w)+2))
    for i in interval:
        if i != interval[-1]:
            chunks.append(seq[i:interval[interval.index(i)+1]-1])
    return("\n".join(chunks))

if __name__ == '__main__':
    main()

