#!/usr/bin/env python
# vcf2fasta

import argparse
import re
import sys
import time
import collections
import os
# test for art
try:
    import pysam
except:
    sys.exit("pysam is not installed. Try: pip install pysam")
# test for art
try:
    import art
except:
    sys.exit("art is not installed. Try: pip install art")

def main():
    print(art.text2art("vcf2fasta"))
    # parse arguments
    parser = argparse.ArgumentParser(prog="vcf2fasta.py",
        #version="0.3",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
        Converts regions/intervals in the genome into FASTA alignments
        provided a VCF file, a GFF file, and FASTA reference.\n""",
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
        \n""")
    parser.add_argument(
    '--fasta', '-f', metavar='GENOME', type=str, required=True,
    help='FASTA file with the reference genome.')
    parser.add_argument(
    '--vcf', '-v', metavar='VCF', type=str, required=True,
    help='a tabix-indexed VCF file.')
    parser.add_argument(
    '--gff', '-g', metavar='GFF', type=str, required=True,
    help='GFF file.')
    parser.add_argument(
    '--feat', '-e', metavar='FEAT', type=str, required=True,
    help='feature/annotation in the GFF file. (i.e. gene, CDS, intron)')
    parser.add_argument(
    '--blend', '-b', action="store_true", default=False,
    help='concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)')
    parser.add_argument(
    '--no-uipac', '-nu', action="store_true", default=False,
    help='selects one allele randomly if heterozygote. (default: False)')
    args = parser.parse_args()

    # read GFF file
    print('Reading VCF file [',args.gff,'] ... ', end='', sep='')
    gff = ReadGFF(args.gff)
    print('done')

    # read variant file and get samples
    print('Reading VCF file [',args.vcf,'] ... ', end='', sep='')
    vcf = pysam.VariantFile(args.vcf)
    # get a list of samples
    samples = [ x for x,y in next(vcf.fetch()).samples.items() ]
    print('done')

    # read genome reference file
    print('Reading VCF file [',args.fasta,'] ... ', end='', sep='')
    ref = pysam.FastaFile(args.fasta)
    print('done')

    # get ploidy
    ploidy = getPloidy(vcf)
    print('Ploidy is:', ploidy)

    # are genotypes phased
    phased = getPhased(vcf)
    if phased:
        print('Genotypes are phased')
    else:
        print('Genotypes are not phased')

    # output directory and print feature
    outdir = "vcf2fasta_"+args.feat
    if args.blend:
        print('Concatenating all [',args.feat,']')
    else:
        print('Writing all [',args.feat,'] separately')
    print('Setting output directory to:', outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        proceed = input(outdir + "exists. Do you want to proceed? [y|n]: ")
        if not re.match('[Yy][EEs]*', proceed):
            print('Exiting ...')
            sys.exit(parser.print_help())

    # get gene keys from GFF
    genes = gff.keys()
    print('Total number of genes found:', len(genes))

    # start counting time
    t1 = time.time()
    # start counter
    feature_counter = 0

    for gene in genes:
        feature_counter += 1
        # genename = gene+"."+gff[gene][0][3]+"-"+gff[gene][-1][4]
        sequences = getSequences(gff, gene, args.feat, args.blend, ref, vcf, ploidy, args.no_uipac, samples)
        with open(outdir + "/" + gene + ".fas", "w") as out:
            printFasta(sequences, out)
        progress = make_progress_bar(feature_counter, len(genes), t1, 70)
        print("\r", progress[0] % progress[1:], end='', flush=True)
    print('')

def getSequences(gff, gene, feat, blend, ref, vcf, ploidy, phased, no_uipac, samples):
    seqs = collections.defaultdict()
    if ploidy == 1:
        for sample in samples: seqs[sample] = ''
        if blend:
            for gffrec in gff[gene][feat]:
                #seq = ''
                #seq2 = ''
                tmpseqs = seqs.copy()
                chrom,start,end,strand = gffrec[0],int(gffrec[3])-1,int(gffrec[4]),gffrec[6]
                refseq = ref.fetch(chrom, start, end).upper()
                for sample in samples: tmpseqs[sample] = refseq
                addposcum = 0
                for rec in vcf.fetch(chrom, start, end):
                    for sample,variant in rec.samples.items():
                        alleles,addpos,refpos = UpdateAllele(variant, rec)
                        tmpseqs[sample] = UpdateSeq(rec, alleles[0], addposcum, refpos, start, tmpseqs[sample])
                        #s2 = UpdateSeq(rec, alleles[1], addposcum, refpos, start, s2)
                        addposcum += addpos
                for sample in samples: seqs[sample] = seqs[sample] + tmpseqs[sample]
            if strand == "-":
                for sample in samples: revcomp(seqs[sample])
        else:
            pass
    if ploidy == 2:
        if phased:
            for sample in samples:
                seqs[sample+"_a"] = ''
                seqs[sample+"_b"] = ''
        else:
            for sample in samples: seqs[sample] = ''
            if blend:
                for gffrec in gff[gene][feat]:
                    #seq = ''
                    #seq2 = ''
                    tmpseqs = seqs.copy()
                    chrom,start,end,strand = gffrec[0],int(gffrec[3])-1,int(gffrec[4]),gffrec[6]
                    refseq = ref.fetch(chrom, start, end).upper()
                    for sample in samples: tmpseqs[sample] = refseq
                    addposcum = 0
                    for rec in vcf.fetch(chrom, start, end):
                        for sample,variant in rec.samples.items():
                            alleles,addpos,refpos = UpdateAllele(variant, rec)
                            tmpseqs[sample] = UpdateSeq(rec, alleles[0], addposcum, refpos, start, tmpseqs[sample])
                            #s2 = UpdateSeq(rec, alleles[1], addposcum, refpos, start, s2)
                            addposcum += addpos
                    for sample in samples: seqs[sample] = seqs[sample] + tmpseqs[sample]
                if strand == "-":
                    for sample in samples: revcomp(seqs[sample])
            else:
                for gffrec in gff[gene][feat]:
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
    return seqs

def getAlleles(rec, ploidy):
    # extract all alleles for a given SNP/var. pos.
    alleles = { i[0]:i[1].alleles for i in rec.samples.items() }
    # collapse list into alleles that are segregating and are not missing
    segregating = list(set(sum([ [ x for x in alleles[i] ] for i in alleles.keys() ], [])))
    # get the length of the longest allele
    max_len = max([ len(i) for i in segregating if i is not None ])
    # make a dictionary of expanded alleles
    dict_expanded = { i:(i + '-' * (max_len - len(i))) for i in segregating if i is not None }
    # replace short alleles with expanded alleles for samples without missing data
    alleles_expanded = { i:[dict_expanded[j] for j in alleles[i]] for i in alleles.keys() if alleles[i][0] is not None }
    # add one for missing data if any, and incorporate to alleles_expanded dict
    if None in segregating:
        dict_expanded[''] = '?' * max_len
        alleles_missing = { i:[dict_expanded[''] for j in range(ploidy)] for i in alleles.keys() if alleles[i][0] is None }
        for i in alleles_missing.keys(): alleles_expanded[i] = alleles_missing[i]



def UpdateAllele(vcfrec, rec):
    allelelen = [ len(x) for x in tuple([rec.ref]) + rec.alts ]
    maxlen = allelelen.index(max(allelelen))
    if vcfrec.alleles[0]:
        alleles = tuple([ x + '-' * (allelelen[maxlen]-len(x)) for x in vcfrec.alleles ])
    else:
        alleles = tuple(['?' * allelelen[maxlen]])
    return alleles, abs(allelelen[0]-allelelen[maxlen]), allelelen[0]-1

def UpdateSeq(vcfrec, allele, addposcum, refpos, start, seq):
    pos = vcfrec.pos-1-start+addposcum
    seq = seq[:pos]+allele+seq[pos+1+refpos:]
    return seq


def printFasta(seqs, out):
    for head in seqs.keys():
        out.write(">" + head + "\n" + seqs[head] + "\n")

def getFeature(file):
    features = collections.defaultdict()
    with open(file, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            features[fields[2]] = None
    return list(features.keys())

def getGeneNames(file):
    geneNames = collections.defaultdict()
    with open(file, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            last = processGeneName(fields[8])
            if last.get('Name'):
                geneNames[last['Name']] = None
    return list(geneNames.keys())

def processGeneName(lastfield):
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        x = i.split("=")
        last[re.sub("\"| ","",x[0])] = re.sub("\"| ","",x[1])
    return last

def ReadGFF(file):
    geneNames = getGeneNames(file)
    features  = getFeature(file)
    gff = collections.defaultdict()
    for g in geneNames:
        gff[g] = collections.defaultdict()
        for f in features:
            gff[g][f] = []
    with open(file, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            last = processGeneName(fields[8])
            if last.get('Name'):
                gff[last['Name']][fields[2]].append(fields)
            else:
                gff[last['Parent']][fields[2]].append(fields)
            # if last.get('Parent'):
            #     gff[last['Parent']][fields[2]].append(fields)
            # elif last.get('Name'):
            #     gff[last['Name']][fields[2]].append(fields)
    return gff

def getPloidy(vcf):
    var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    p = sum([ len(v.get('GT')) for v in var ]) / len(var)
    return int(p)

def getPhased(vcf):
    var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    p = any([ not v.phased for v in var ])
    return not p

def revcomp(seq):
    tt = seq.maketrans('ACGT?N-','TGCA?N-')
    return seq[::-1].translate(tt)

def make_progress_bar(rec, total, t1, width):
    i = (rec/total * 100) % 100
    if i != 0:
        plus = "+" * int(i * (width/100))
        dots = "." * (width-int(i * width/100))
    else:
        plus = "+" * width
        dots = ""
        i = 100
    t2 = time.time()
    elapsed = t2-t1
    return "["+plus+dots+"] "+"%5.2f%% %7.2f s", i, elapsed

if __name__ == "__main__":
    main()
