import re
import sys
import time
import collections


# this function gets everything processed to produce a set
# of sequences in a dictionary, along with other info
def getSequences(intervals, gene, ref, vcf, ploidy, phased, samples, args):
    seqs = collections.defaultdict()
    feat_str = "_" if args.feat == "" else "_" + args.feat + "_"
    # prep sequence dictionary
    if args.blend:
        seqs[gene] = collections.defaultdict()
    else:
        feat_ind = 0
        for rec in intervals[gene]:
            featname = gene + feat_str + str(feat_ind)
            seqs[featname] = collections.defaultdict()
            feat_ind += 1

    if phased:
        for key in seqs.keys():
            for sample in samples:
                for i in range(ploidy):
                    seqs[key][sample + "_" + str(i)] = ""
            # add a dummy reference sequence if requested
            if args.addref:
                seqs[key]["REF_0"] = ""
    else:
        for key in seqs.keys():
            for sample in samples:
                seqs[key][sample] = ""
            # add a dummy reference sequence if requested
            if args.addref:
                seqs[key]["REF"] = ""
    # initial values before looping through VCF slice
    varsites = 0
    feat_ind = 0
    featname = gene
    codon_start = collections.defaultdict(list)
    for rec in intervals[gene]:
        if not args.blend:
            featname = gene + feat_str + str(feat_ind)
            feat_ind += 1
        # get a new copy of seqs for every feature
        tmpseqs = collections.defaultdict()
        # extract relevant info from GFF/BED
        if args.bed:
            chrom, start, end, strand, cs = rec[0], int(rec[1]), int(rec[2]), "+", 0
        else:
            chrom, start, end, strand, cs = (
                rec[0],
                int(rec[3]) - 1,
                int(rec[4]),
                rec[6],
                rec[7],
            )
        # add cs to codon_start
        codon_start[featname].append(cs)
        # extract sequence from reference
        refseq = ref.fetch(chrom, start, end).upper()
        # propagate reference sequence to all samples
        for sample in seqs[featname].keys():
            tmpseqs[sample] = refseq
        # initialize posiitive or negative postions to extend
        # in case of indels
        posadd = 0
        for vrec in vcf.fetch(chrom, start, end):
            # count variants within feature
            varsites += 1
            # get seq position of variant
            pos = vrec.pos - start - 1 + posadd
            # get a dict of extended alleles, including indels
            alleles, max_len = getAlleles(vrec, ploidy, phased, args.addref)
            ref_len = len(vrec.ref)
            for sample in alleles.keys():
                tmpseqs[sample] = UpdateSeq(
                    alleles, sample, pos, ref_len, tmpseqs[sample]
                )
            # this is the cumulative number of positions to add if
            # positions are take or are added to the sequence
            posadd += max_len - ref_len
        for sample in seqs[featname].keys():
            seqs[featname][sample] = seqs[featname][sample] + tmpseqs[sample]
    # reverse complement sequence if needed
    if strand == "-":
        for sample in seqs[featname].keys():
            seqs[featname][sample] = revcomp(seqs[featname][sample])
    # adjust if inframe is on
    if args.inframe:
        if args.blend and codon_start[featname][0] != ".":
            if strand == "+" and codon_start[featname][0] != "0":
                for key in seqs[featname].keys():
                    seqs[featname][key] = seqs[featname][key][
                        int(codon_start[featname][0]) :
                    ]
            elif strand == "-" and codon_start[featname][-1] != "0":
                for key in seqs[featname].keys():
                    seqs[featname][key] = seqs[featname][key][
                        int(codon_start[featname][-1]) :
                    ]
        elif not args.blend:
            for featname in seqs.keys():
                for key in seqs[featname].keys():
                    seqs[featname][key] = seqs[featname][key][
                        int(codon_start[featname][0]) :
                    ]
    return seqs, varsites


# main algorithm of vcf2fasta
# extracts alleles for each sample and expands them if indels
# deals with both phased and unphased data
def getAlleles(rec, ploidy, phased, addref):
    # extract all alleles for a given SNP/var. pos.
    alleles = {i[0]: i[1].alleles for i in rec.samples.items()}
    # collapse list into alleles that are segregating and are not missing
    segregating = list(set(sum([[x for x in alleles[i]] for i in alleles.keys()], [])))
    # add ref allele if addref is true
    if addref:
        segregating = list(set(segregating + [rec.ref]))
    # get the length of the longest allele
    max_len = max([len(i) for i in segregating if i is not None])
    # make a dictionary of expanded alleles
    dict_expanded = {
        i: (i + "-" * (max_len - len(i))) for i in segregating if i is not None
    }
    # replace short alleles with expanded alleles for samples without missing data
    alleles_expanded = {
        i: [dict_expanded[j] for j in alleles[i]]
        for i in alleles.keys()
        if alleles[i][0] is not None
    }
    # add ref allele if addref is true
    if addref:
        alleles_expanded["REF"] = [dict_expanded[rec.ref]]
    # add one for missing data if any, and incorporate to alleles_expanded dict
    if None in segregating:
        dict_expanded[""] = "?" * max_len
        alleles_missing = {
            i: [dict_expanded[""] for j in range(ploidy)]
            for i in alleles.keys()
            if alleles[i][0] is None
        }
        for i in alleles_missing.keys():
            alleles_expanded[i] = alleles_missing[i]
    if phased:
        alleles_expanded = makePhased(alleles_expanded)  # make them phased
    else:
        alleles_expanded = {
            key: getIUPAC(alleles_expanded[key]) for key in alleles_expanded.keys()
        }
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


# function to update sequence dict based on alleles or collapsed alleles
def UpdateSeq(alleles, samp, pos, ref_len, seq):
    return seq[:pos] + alleles[samp] + seq[pos + ref_len :]


# same as UpdateSeqPhased but returns a collapsed genotype using IUPAC codes
# def UpdateSeqIUPAC(alleles, samp, pos, ref_len, seq):
#     return seq[:pos] + getIUPAC(alleles[samp]) + seq[pos+ref_len:]


# collapses list into string of IUPAC codes
### ISSUE
# produces shorter alignments sometimes
# needs fixing
###
def getIUPAC(x):
    """
    Collapses two or more alleles into a single
    IUPAC string
    """
    # first check if data is missing
    if len(x) == 1 or x[0][0] == "?":
        return x[0]
    elif len(list(set(x))) == 1:
        return x[0]
    else:
        iupacd = ""
        # loop through positions
        for i in range(len(x[0])):
            # keep nucleotides only
            nuc = list(set([y[i] for y in x]))
            # if single nuc
            if len(nuc) == 1:
                iupacd += nuc[0]
            else:
                nuc = [j for j in nuc if j != "-"]
                if len(nuc) == 1:
                    iupacd += nuc[0]
                # if multiple nucleotides per position
                if len(nuc) == 2:
                    if "A" in nuc and "G" in nuc:
                        iupacd += "R"
                    elif "A" in nuc and "T" in nuc:
                        iupacd += "W"
                    elif "A" in nuc and "C" in nuc:
                        iupacd += "M"
                    elif "C" in nuc and "T" in nuc:
                        iupacd += "Y"
                    elif "C" in nuc and "G" in nuc:
                        iupacd += "S"
                    elif "G" in nuc and "T" in nuc:
                        iupacd += "K"
                    else:
                        iupacd += "?"
                elif len(nuc) == 3:
                    if "A" in nuc and "T" in nuc and "C":
                        iupacd += "H"
                    elif "A" in nuc and "T" in nuc and "G":
                        iupacd += "D"
                    elif "G" in nuc and "T" in nuc and "C":
                        iupacd += "B"
                    elif "A" in nuc and "G" in nuc and "C":
                        iupacd += "V"
                    else:
                        iupacd += "?"
                elif len(nuc) == 4:
                    if "A" in nuc and "G" in nuc and "C" and "T" in nuc:
                        iupacd += "N"
                    else:
                        iupacd += "?"
        return iupacd


def printFasta(seqs, out):
    """
    writes FASTA to file
    """
    for head in seqs.keys():
        out.write(">" + head + "\n" + seqs[head] + "\n")


def getFeature(file):
    """
    extracts a list of features from the GFF
    """
    features = collections.defaultdict()
    with open(file, "r") as f:
        for line in f:
            if line[0] != "#":
                fields = line.rstrip().split("\t")
                features[fields[2]] = None
    return list(features.keys())


def getGeneNames(file, format):
    """
    Extracts a list of all gene names in GFF/GTF
    The input is the gff/gtf file itself
    """
    geneNames = collections.defaultdict()
    with open(file, "r") as f:
        if format == "gff":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGFF(fields[8])
                    if last.get("Name"):
                        geneNames[last["Name"]] = None
                    elif last.get("Parent"):
                        geneNames[last["Parent"]] = None
                    elif last.get("ID"):
                        geneNames[last["ID"]] = None
        elif format == "gtf":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGTF(fields[8])
                    if last.get("transcript_id"):
                        geneNames[last["transcript_id"]] = None
                    elif last.get("gene_id"):
                        geneNames[last["gene_id"]] = None
    return list(geneNames.keys())


def processGeneNameGFF(lastfield):
    """
    Makes a list of all the annotation fields in the last column [8]
    delimited by ";"
    Input is the string of the last field in GFF
    """
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        if "=" in i:
            x = i.split("=")
            last[re.sub('"| ', "", x[0])] = re.sub('"| ', "", x[1])
    return last


def processGeneNameGTF(lastfield):
    """
    Makes a list of all the annotation fields in the last column [8]
    delimited by ";"
    Input is the string of the last field in GFF
    """
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        if " " in i:
            x = i.split(" ")
            last[x[0]] = re.sub('"| ', "", x[1])
    return last


def ReadBED(file):
    """
    returns a dictionary of with 0-based genomic intervals
    """
    with open(file) as f:
        lines = f.read().splitlines()
        bed = {"g" + str(i + 1): [lines[i].split("\t")] for i in range(len(lines))}
    return bed


def ReadGFF(file, parser):
    """
    returns a nested dictionary named after every feature name
    as well as every feature name [3]
    """
    # find out if GTF or GFF
    if file.split(".")[-1].lower() == "gff" or file.split(".")[-1].lower() == "gff3":
        format = "gff"
    elif file.split(".")[-1].lower() == "gtf":
        format = "gtf"
    else:
        print("Cannot figure out GFF/GTF format. File should end with .gff or .gtf")
        sys.exit(parser.print_help())
    # get gene names
    geneNames = getGeneNames(file, format)
    features = getFeature(file)
    gff = collections.defaultdict()
    for g in geneNames:
        gff[g] = collections.defaultdict()
        for f in features:
            gff[g][f] = []
    with open(file, "r") as f:
        if format == "gff":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGFF(fields[8])
                    if last.get("Name"):
                        gff[last["Name"]][fields[2]].append(fields)
                    elif last.get("Parent"):
                        gff[last["Parent"]][fields[2]].append(fields)
                    else:
                        gff[last["ID"]][fields[2]].append(fields)
                    # if last.get('Parent'):
                    #     gff[last['Parent']][fields[2]].append(fields)
                    # elif last.get('Name'):
                    #     gff[last['Name']][fields[2]].append(fields)
        elif format == "gtf":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGTF(fields[8])
                    if last.get("transcript_id"):
                        gff[last["transcript_id"]][fields[2]].append(fields)
                    elif last.get("gene_id"):
                        gff[last["gene_id"]][fields[2]].append(fields)
    return gff


def filterFeatureInGFF(gff, feat):
    """
    Keep/filters GFF records that include specified feature
    and returns a new GFF
    """
    filtered_gff = collections.defaultdict()
    for gene in gff.keys():
        if len(gff[gene][feat]) != 0:
            filtered_gff[gene] = gff[gene][feat]
        else:
            _ = filtered_gff.pop(gene)
    return filtered_gff


def getPloidy(vcf):
    var = [y for x, y in next(vcf.fetch()).samples.items()]
    p = [len(v.get("GT")) for v in var if v.get("GT")[0] is not None]
    return p[0]


def getPhased(vcf):
    var = [y for x, y in next(vcf.fetch()).samples.items()]
    p = any([not v.phased for v in var])
    return not p


def revcomp(seq):
    tt = seq.maketrans("ACGT?N-", "TGCA?N-")
    return seq[::-1].translate(tt)


def make_progress_bar(rec, total, t1, width):
    i = (rec / total * 100) % 100
    if i != 0:
        plus = "+" * int(i * (width / 100))
        dots = "." * (width - int(i * width / 100))
    else:
        plus = "+" * width
        dots = ""
        i = 100
    t2 = time.time()
    elapsed = t2 - t1
    return "[" + plus + dots + "] " + "%5.2f%% %7.2f s", i, elapsed


def get_samples(vcf):
    # Function to get a list of sample names from the VCF file using pysam
    if hasattr(vcf, "header") and hasattr(vcf.header, "samples"):
        return list(vcf.header.samples)
    return []
