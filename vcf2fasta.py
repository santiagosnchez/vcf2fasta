#!/usr/bin/env python3
# vcf2fasta

import argparse
import re
import sys
import time
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

import v2f.functions as v2f

def main():
    print(art.text2art("vcf2fasta"))
    # parse arguments
    parser = argparse.ArgumentParser(prog="vcf2fasta.py",
        #version="0.3",
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
        Converts regions/intervals in the genome into FASTA alignments
        provided a VCF file, a GFF/GTF file, and FASTA reference.\n""",
        epilog="""
        All files must be indexed. So before running the code make sure
        that your reference FASTA file is indexed:

        samtools faidx genome.fas

        BGZIP compress and TABIX index your VCF file:

        bgzip variants.vcf
        tabix variants.vcf.gz

        The GFF/GTF file does not need to be indexed.

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
    '--gff', '-g', metavar='GFF/GTF', type=str, required=False,
    help='GFF/GTF file.')
    parser.add_argument(
    '--bed', metavar='BED', type=str, required=False,
    help='BED file.')
    parser.add_argument(
    '--feat', '-e', metavar='FEAT', type=str, required=False,
    help='feature/annotation in the GFF file. (i.e. gene, CDS, intron)')
    parser.add_argument(
    '--blend', '-b', action="store_true", default=False,
    help='concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)')
    parser.add_argument(
    '--inframe', '-i', action="store_true", default=False,
    help='force the first codon of the sequence to be inframe. Useful for incomplete CDS. (default: False)')
    parser.add_argument(
    '--out', '-o', metavar='OUT', type=str, default="vcf2fasta",
    help='provide a name for the output directory (optional)')
    parser.add_argument(
    '--addref', '-r', action="store_true", default=False,
    help='include the reference sequence in the FASTA alignment (default: False)')
    parser.add_argument(
    '--skip', '-s', action="store_true", default=False,
    help='skips features without variants (default: False)')

    args = parser.parse_args()

    # resolve GFF/GTF | BED
    if args.gff and args.bed:
        print('Only --gff or --bed are allowed. Exiting ...')
        sys.exit(parser.print_help())
    elif not args.gff and not args.bed:
        print('Either --gff or --bed are required. Exiting ...')
        sys.exit(parser.print_help())
    else:
        if args.bed:
            print('Reading BED file [',args.bed,'] ... ', end='', sep='')
            intervals = v2f.ReadBED(args.bed)
        else:
            # read GFF file
            print('Reading GFF file [',args.gff,'] ... ', end='', sep='')
            intervals = v2f.filterFeatureInGFF(v2f.ReadGFF(args.gff, parser), args.feat)
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
    ploidy = v2f.getPloidy(vcf)
    print('Ploidy is:', ploidy)

    # are genotypes phased
    phased = v2f.getPhased(vcf)
    if not phased:
        print('No phased genotypes found on first variant. Treating as \"unphased\"')
    else:
        print('Phased genotypes found on first variant. Treating as \"phased\"')

    # output directory and print feature
    if args.feat:
        outdir = args.out + "_" + args.feat
    else:
        outdir = args.out
    if args.blend:
        print('Concatenating all [', args.feat, ']')
    else:
        print('Writing all [', args.feat, '] separately')
    print('Setting output directory to:', outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        proceed = input(outdir + " exists. Files will me replaced. Do you want to proceed? [y|n]: ")
        if not re.match('[Yy][EEs]*', proceed):
            print('Exiting ...')
            sys.exit(parser.print_help())


    # get gene keys from GFF
    genes = list(intervals.keys())
    print('Total number of genes found:', len(genes))

    # check multiple records per feat
    if args.gff:
        single = [ len(intervals[i]) == 1 for i in genes ]
        if all(single):
            args.blend = True
            print("Found all genes with single records. Treating as --blend")

    # start counting time
    t1 = time.time()
    # start counter
    feature_counter = 0
    # count skipped genes
    withdata = 0

    for gene in genes:
        #sequences,strand,codon_start = getSequences(gff, gene, args.feat, args.blend, ref, vcf, ploidy, phased, samples, args.addref)
        sequences, varsites = v2f.getSequences(
            intervals, 
            gene, 
            ref, 
            vcf, 
            ploidy, 
            phased, 
            samples,
            args
        )
        if args.skip and varsites != 0:
            withdata += 1
            for featname in sequences.keys():
                with open(outdir + "/" + featname + ".fas", "w") as out:
                    v2f.printFasta(sequences[featname], out)
        elif not args.skip:
            for featname in sequences.keys():
                with open(outdir + "/" + featname + ".fas", "w") as out:
                    v2f.printFasta(sequences[featname], out)
        feature_counter += 1
        progress = v2f.make_progress_bar(feature_counter, len(genes), t1, 70)
        print("\r", progress[0] % progress[1:], end='', flush=True)
    print('')
    if args.skip:
        print("Skipped", feature_counter - withdata, "genes with no variants")

if __name__ == "__main__":
    main()
