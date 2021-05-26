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

    args = parser.parse_args()

    # read GFF file
    print('Reading GFF file [',args.gff,'] ... ', end='', sep='')
    gff = filterFeatureInGFF(ReadGFF(args.gff), args.feat)
    print('done')

    # read variant file and get samples
    print('Reading VCF file [',args.vcf,'] ... ', end='', sep='')
    vcf = pysam.VariantFile(args.vcf)
    # get a list of samples
    samples = [ x for x,y in next(vcf.fetch()).samples.items() ]
    print('done')

    # read genome reference file
    print('Reading FASTA reference file [',args.fasta,'] ... ', end='', sep='')
    ref = pysam.FastaFile(args.fasta)
    print('done')

    # get ploidy
    ploidy = getPloidy(vcf)
    print('Ploidy is:', ploidy)

    # are genotypes phased
    phased = getPhased(vcf)
    if not phased:
        print('No phased genotypes found on first variant. Treating as \"unphased\"')
    else:
        print('Phased genotypes found on first variant. Treating as \"phased\"')

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
        proceed = input(outdir + " exists. Do you want to proceed? [y|n]: ")
        if not re.match('[Yy][EEs]*', proceed):
            print('Exiting ...')
            sys.exit(parser.print_help())

    # get gene keys from GFF
    genes = list(gff.keys())
    print('Total number of genes found:', len(genes))

    # check multiple records per feat
    single = [ len(gff[i]) == 1 for i in genes ]
    if all(single):
        args.blend = True

    # start counting time
    t1 = time.time()
    # start counter
    feature_counter = 0

    for gene in genes:
        # genename = gene+"."+gff[gene][0][3]+"-"+gff[gene][-1][4]
        sequences = getSequences(gff, gene, args.feat, args.blend, ref, vcf, ploidy, phased, samples)
        if args.blend:
            with open(outdir + "/" + gene + ".fas", "w") as out:
                printFasta(sequences, out)
        else:
            for featnm in sequences.keys():
                with open(outdir + "/" + featnm + ".fas", "w") as out:
                    printFasta(sequences[featnm], out)
        feature_counter += 1
        progress = make_progress_bar(feature_counter, len(genes), t1, 70)
        print("\r", progress[0] % progress[1:], end='', flush=True)
    print('')

def getSequences(gff, gene, feat, blend, ref, vcf, ploidy, phased, samples):
    seqs = collections.defaultdict()
    if phased:
        if blend:
            for sample in samples:
                for i in range(ploidy):
                    seqs[sample + "_" + str(i)] = ''
            for gffrec in gff[gene][feat]:
                # get a new copy of seqs for every feature
                tmpseqs = collections.defaultdict()
                # extract relevant info from GFF
                chrom,start,end,strand = gffrec[0],int(gffrec[3])-1,int(gffrec[4]),gffrec[6]
                # extract sequence from reference
                refseq = ref.fetch(chrom, start, end).upper()
                # propagate reference sequence to all samples
                for sample in seqs.keys(): tmpseqs[sample] = refseq
                # initialize posiitive or negative postions to extend
                # in case of indels
                posadd = 0
                for rec in vcf.fetch(chrom, start, end):
                    # get seq position of variant
                    pos = rec.pos - start - 1 + posadd
                    # get a dict of extended alleles, including indels
                    alleles,max_len = getAlleles(rec, ploidy)
                    alleles = makePhased(alleles) # make them phased
                    ref_len = len(rec.ref)
                    for sample in alleles.keys():
                        tmpseqs[sample] = UpdateSeqPhased(alleles, sample, pos, ref_len, tmpseqs[sample])
                    # this is the cumulative number of positions to add if
                    # positions are take or are added to the sequence
                    posadd += max_len - ref_len
                for sample in seqs.keys(): seqs[sample] = seqs[sample] + tmpseqs[sample]
            # reverse complement sequence if needed
            if strand == "-":
                for sample in seqs.keys(): seqs[sample] = revcomp(seqs[sample])
        else:
            feat_ind = 0
            for gffrec in gff[gene][feat]:
                featname = gene+"_"+feat+"_"+str(feat_ind)
                seqs[featname] = collections.defaultdict()
                # extract relevant info from GFF
                chrom,start,end,strand = gffrec[0],int(gffrec[3])-1,int(gffrec[4]),gffrec[6]
                # extract sequence from reference
                refseq = ref.fetch(chrom, start, end).upper()
                for sample in samples:
                    for i in range(ploidy):
                        seqs[featname][sample + "_" + str(i)] = refseq
                # initialize posiitive or negative postions to extend
                # in case of indels
                posadd = 0
                for rec in vcf.fetch(chrom, start, end):
                    # get seq position of variant
                    pos = rec.pos - start - 1 + posadd
                    # get a dict of extended alleles, including indels
                    alleles,max_len = getAlleles(rec, ploidy)
                    alleles = makePhased(alleles) # make them phased
                    ref_len = len(rec.ref)
                    for sample in alleles.keys():
                        seqs[featname][sample] = UpdateSeqPhased(alleles, sample, pos, ref_len, seqs[featname][sample])
                    # this is the cumulative number of positions to add if
                    # positions are take or are added to the sequence
                    posadd += max_len - ref_len
                feat_ind += 1
                if strand == "-":
                    for sample in seqs.keys(): seqs[featname][sample] = revcomp(seqs[featname][sample])
    else:
        if blend:
            for sample in samples:
                seqs[sample] = ''
            for gffrec in gff[gene][feat]:
                # get a new copy of seqs for every feature
                tmpseqs = seqs.copy()
                # extract relevant info from GFF
                chrom,start,end,strand = gffrec[0],int(gffrec[3])-1,int(gffrec[4]),gffrec[6]
                # extract sequence from reference
                refseq = ref.fetch(chrom, start, end).upper()
                # propagate reference sequence to all samples
                for sample in seqs.keys(): tmpseqs[sample] = refseq
                # initialize posiitive or negative postions to extend
                # in case of indels
                posadd = 0
                for rec in vcf.fetch(chrom, start, end):
                    # get seq position of variant
                    pos = rec.pos - start - 1 + posadd
                    # get a dict of extended alleles, including indels
                    alleles,max_len = getAlleles(rec, ploidy)
                    ref_len = len(rec.ref)
                    for sample in alleles.keys():
                        tmpseqs[sample] = UpdateSeqIUPAC(alleles, sample, pos, ref_len, tmpseqs[sample])
                    # this is the cumulative number of positions to add if
                    # positions are take or are added to the sequence
                    posadd += max_len - ref_len
                for sample in seqs.keys(): seqs[sample] = seqs[sample] + tmpseqs[sample]
            # reverse complement sequence if needed
            if strand == "-":
                for sample in seqs.keys(): seqs[sample] = revcomp(seqs[sample])
        else:
            feat_ind = 0
            for gffrec in gff[gene][feat]:
                featname = gene+"_"+feat+"_"+str(feat_ind)
                seqs[featname] = collections.defaultdict()
                # extract relevant info from GFF
                chrom,start,end,strand = gffrec[0],int(gffrec[3])-1,int(gffrec[4]),gffrec[6]
                # extract sequence from reference
                refseq = ref.fetch(chrom, start, end).upper()
                for sample in samples:
                    for i in range(ploidy):
                        seqs[featname][sample + "_" + str(i)] = refseq
                # initialize posiitive or negative postions to extend
                # in case of indels
                posadd = 0
                for rec in vcf.fetch(chrom, start, end):
                    # get seq position of variant
                    pos = rec.pos - start - 1 + posadd
                    # get a dict of extended alleles, including indels
                    alleles,max_len = getAlleles(rec, ploidy)
                    ref_len = len(rec.ref)
                    for sample in alleles.keys():
                        seqs[featname][sample] = UpdateSeqIUPAC(alleles, sample, pos, ref_len, seqs[featname][sample])
                    # this is the cumulative number of positions to add if
                    # positions are take or are added to the sequence
                    posadd += max_len - ref_len
                feat_ind += 1
                if strand == "-":
                    for sample in seqs.keys(): seqs[featname][sample] = revcomp(seqs[featname][sample])
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
    return alleles_expanded, max_len

# takes a genotype dict from getAlleles and turns it into
# a phased dict with each haplotype as the sample + a zero-based index
def makePhased(alleles):
    alleles_phased = {}
    for samp in alleles.keys():
        for i in range(len(alleles[samp])):
            # append zero-based index
            alleles_phased[samp + "_" + str(i)] = alleles[samp][i]
    return alleles_phased

# function to update sequence dict based on alleles
# works on phased and haploid data
def UpdateSeqPhased(alleles, samp, pos, ref_len, seq):
    return seq[:pos] + alleles[samp] + seq[pos+ref_len:]

# same as UpdateSeqPhased but returns a collapsed genotype using IUPAC codes
def UpdateSeqIUPAC(alleles, samp, pos, ref_len, seq):
    return seq[:pos] + getIUPAC(alleles[samp]) + seq[pos+ref_len:]

# collapses list into string of IUPAC codes
### ISSUE
# produces shorter alignments
# needs fixing
###
def getIUPAC(x):
    '''
    Collapses two or more alleles into a single
    IUPAC string
    '''
    # first check if data is missing
    if x[0][0] == '?':
        return x[0]
    elif len(list(set(x))) == 1:
        return x[0]
    else:
        iupacd = ''
        # loop through positions
        for i in range(len(x[0])):
            # keep nucleotides only
            nuc = list(set([ y[i] for y in x ]))
            # if single nuc
            if len(nuc) == 1:
                iupacd += nuc[0]
            else:
                nuc = [ j for j in nuc if j != '-' ]
                if len(nuc) == 1:
                    iupacd += nuc[0]
                # if multiple nucleotides per position
                if len(nuc) == 2:
                    if   'A' in nuc and 'G' in nuc:
                        iupacd += 'R'
                    elif 'A' in nuc and 'T' in nuc:
                        iupacd += 'W'
                    elif 'A' in nuc and 'C' in nuc:
                        iupacd += 'M'
                    elif 'C' in nuc and 'T' in nuc:
                        iupacd += 'Y'
                    elif 'C' in nuc and 'G' in nuc:
                        iupacd += 'S'
                    elif 'G' in nuc and 'T' in nuc:
                        iupacd += 'K'
                    else:
                        iupacd += '?'
                elif len(nuc) == 3:
                    if 'A' in nuc and 'T' in nuc and 'C':
                        iupacd += 'H'
                    elif 'A' in nuc and 'T' in nuc and 'G':
                        iupacd += 'D'
                    elif 'G' in nuc and 'T' in nuc and 'C':
                        iupacd += 'B'
                    elif 'A' in nuc and 'G' in nuc and 'C':
                        iupacd += 'V'
                    else:
                        iupacd += '?'
                elif len(nuc) == 4:
                    if 'A' in nuc and 'G' in nuc and 'C' and 'T' in nuc:
                        iupacd += 'N'
                    else:
                        iupacd += '?'
        return iupacd

# def UpdateAllele(vcfrec, rec):
#     allelelen = [ len(x) for x in tuple([rec.ref]) + rec.alts ]
#     maxlen = allelelen.index(max(allelelen))
#     if vcfrec.alleles[0]:
#         alleles = tuple([ x + '-' * (allelelen[maxlen]-len(x)) for x in vcfrec.alleles ])
#     else:
#         alleles = tuple(['?' * allelelen[maxlen]])
#     return alleles, abs(allelelen[0]-allelelen[maxlen]), allelelen[0]-1


def printFasta(seqs, out):
    for head in seqs.keys():
        out.write(">" + head + "\n" + seqs[head] + "\n")

def getFeature(file):
    '''
    extracts a list of features from the GFF
    '''
    features = collections.defaultdict()
    with open(file, "r") as f:
        for line in f:
            if line[0] != "#":
                fields = line.rstrip().split("\t")
                features[fields[2]] = None
    return list(features.keys())

def getGeneNames(file):
    '''
    Extracts a list of all gene names in GFF
    The input is the gff file itself
    '''
    geneNames = collections.defaultdict()
    with open(file, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            last = processGeneName(fields[8])
            if last.get('Name'):
                geneNames[last['Name']] = None
            elif last.get('Parent'):
                geneNames[last['Parent']] = None
            elif last.get('ID'):
                geneNames[last['ID']] = None
    return list(geneNames.keys())

def processGeneName(lastfield):
    '''
    Makes a list of all the annotation fields in the last column [8]
    delimited by ";"
    Input is the string of the last field in GFF
    '''
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        if "=" in i:
            x = i.split("=")
            last[re.sub("\"| ","",x[0])] = re.sub("\"| ","",x[1])
    return last

def ReadGFF(file):
    '''
    returns a nested dictionary named after every feature name
    as well as every feature name [3]
    '''
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
            elif:
                gff[last['Parent']][fields[2]].append(fields)
            else:
                gff[last['ID']][fields[2]].append(fields)
            # if last.get('Parent'):
            #     gff[last['Parent']][fields[2]].append(fields)
            # elif last.get('Name'):
            #     gff[last['Name']][fields[2]].append(fields)
    return gff

def filterFeatureInGFF(gff, feat):
    '''
    Keep/filters GFF records that include specified feature
    and returns a new GFF
    '''
    filtered_gff = collections.defaultdict()
    for gene in gff.keys(): filtered_gff[gene] = collections.defaultdict()
    for gene in gff.keys():
        if len(gff[gene][feat]) != 0:
            filtered_gff[gene][feat] = gff[gene][feat]
        else:
            _ = filtered_gff.pop(gene)
    return filtered_gff

def getPloidy(vcf):
    var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    p = [ len(v.get('GT')) for v in var if v.get('GT')[0] is not None ]
    return p[0]

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
