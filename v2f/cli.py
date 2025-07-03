import argparse
import sys
import os
import re
import time
import pysam

import v2f.functions as v2f_helper
import v2f.parser as v2f_parser


def main():

    parser = v2f_parser.get_parser()
    args = parser.parse_args()

    if args.gff and args.bed:
        print("Only --gff or --bed are allowed. Exiting ...")
        sys.exit(parser.print_help())
    elif not args.gff and not args.bed:
        print("Either --gff or --bed are required. Exiting ...")
        sys.exit(parser.print_help())
    else:
        if args.gff and not args.feat:
            print(
                "No feature (--feat) specified. Defaulting to 'gene'. "
                "If you want to extract all features, use --feat ''"
            )
            args.feat = "gene"
        if args.bed:
            print("Reading BED file [", args.bed, "] ... ", end="", sep="")
            intervals = v2f_helper.ReadBED(args.bed)
        else:
            # read GFF file
            print("Reading GFF file [", args.gff, "] ... ", end="", sep="")
            intervals = v2f_helper.filterFeatureInGFF(
                v2f_helper.ReadGFF(args.gff, parser), args.feat
            )
    print("done")

    print("Reading VCF file [", args.vcf, "] ... ", end="", sep="")
    vcf = pysam.VariantFile(args.vcf)
    samples = v2f_helper.get_samples(vcf)
    print("done")

    print("Reading FASTA reference file [", args.fasta, "] ... ", end="", sep="")
    ref = pysam.FastaFile(args.fasta)
    print("done")

    ploidy = v2f_helper.getPloidy(vcf)
    print("Ploidy is:", ploidy)

    phased = v2f_helper.getPhased(vcf)
    if not phased:
        print('No phased genotypes found on first variant. Treating as "unphased"')
    else:
        print('Phased genotypes found on first variant. Treating as "phased"')

    if args.feat:
        outdir = args.out + "_" + args.feat
    else:
        outdir = args.out

    if args.blend:
        print("Concatenating all [", args.feat, "]")
    else:
        print("Writing all [", "intervals" if args.feat == "" else args.feat, "] separately")
    print("Setting output directory to:", outdir)

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        proceed = input(
            outdir + " exists. Files will me replaced. Do you want to proceed? [y|n]: "
        )
        if not re.match("[Yy][EEs]*", proceed):
            print("Exiting ...")
            sys.exit(parser.print_help())

    # get gene keys from GFF
    genes = list(intervals.keys())
    print("Total number of genes found:", len(genes))

    # check multiple records per feat
    if args.gff:
        single = [len(intervals[i]) == 1 for i in genes]
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
        # sequences,strand,codon_start = getSequences(gff, gene, args.feat, args.blend, ref, vcf, ploidy, phased, samples, args.addref)
        sequences, varsites = v2f_helper.getSequences(
            intervals, gene, ref, vcf, ploidy, phased, samples, args
        )
        if args.skip and varsites != 0:
            withdata += 1
            for featname in sequences.keys():
                with open(outdir + "/" + featname + ".fas", "w") as out:
                    v2f_helper.printFasta(sequences[featname], out)
        elif not args.skip:
            for featname in sequences.keys():
                with open(outdir + "/" + featname + ".fas", "w") as out:
                    v2f_helper.printFasta(sequences[featname], out)
        feature_counter += 1
        progress = v2f_helper.make_progress_bar(feature_counter, len(genes), t1, 70)
        print("\r", progress[0] % progress[1:], end="", flush=True)
    print("")
    if args.skip:
        print("Skipped", feature_counter - withdata, "genes with no variants")


if __name__ == "__main__":
    main()
