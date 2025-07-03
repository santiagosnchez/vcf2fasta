import argparse
import art


def get_parser():
    print(art.text2art("vcf2fasta"))
    # parse arguments
    parser = argparse.ArgumentParser(
        prog="vcf2fasta.py",
        # version="0.3",
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
        \n""",
    )
    parser.add_argument(
        "--fasta",
        "-f",
        metavar="GENOME",
        type=str,
        required=True,
        help="FASTA file with the reference genome.",
    )
    parser.add_argument(
        "--vcf",
        "-v",
        metavar="VCF",
        type=str,
        required=True,
        help="a tabix-indexed VCF file.",
    )
    parser.add_argument(
        "--gff", "-g", metavar="GFF/GTF", type=str, required=False, help="GFF/GTF file."
    )
    parser.add_argument(
        "--bed", metavar="BED", type=str, required=False, help="BED file."
    )
    parser.add_argument(
        "--feat",
        "-e",
        default="",
        metavar="FEAT",
        type=str,
        required=False,
        help="feature/annotation in the GFF file. (i.e. gene, CDS, intron); default: ''.\n",
    )
    parser.add_argument(
        "--blend",
        "-b",
        action="store_true",
        default=False,
        help="concatenate GFF entries of FEAT into a single alignment. Useful for CDS. (default: False)",
    )
    parser.add_argument(
        "--inframe",
        "-i",
        action="store_true",
        default=False,
        help="force the first codon of the sequence to be inframe. Useful for incomplete CDS. (default: False)",
    )
    parser.add_argument(
        "--out",
        "-o",
        metavar="OUT",
        type=str,
        default="vcf2fasta",
        help="provide a name for the output directory (optional)",
    )
    parser.add_argument(
        "--addref",
        "-r",
        action="store_true",
        default=False,
        help="include the reference sequence in the FASTA alignment (default: False)",
    )
    parser.add_argument(
        "--skip",
        "-s",
        action="store_true",
        default=False,
        help="skips features without variants (default: False)",
    )

    return parser
