from time import localtime, strftime
import pandas as pd
import numpy as np
import argparse
from collections import namedtuple
from multiprocessing import Pool
import scipy.io as sio


annotation_categories = ["transcript", "exon", "intron", "utr5", "utr3",
    "coding", "upstream_1kb", "downstream_1kb"]

valid_chromosomes = map(str, range(1,23))

required_biomart_cols = {"Associated Gene Name":    "gene_name",
                         "Ensembl Gene ID":         "gene_id",
                         "Gene type":               "gene_type",
                         "Ensembl Transcript ID":   "transcript_id",
                         "Transcript Start (bp)":   "transcript_start",
                         "Transcript End (bp)":     "transcript_end",
                         "Chromosome Name":         "chr",
                         "Strand":                  "strand",
                         "5' UTR Start":            "utr5_start",
                         "5' UTR End":              "utr5_end",
                         "3' UTR Start":            "utr3_start",
                         "3' UTR End":              "utr3_end",
                         "Ensembl Exon ID":         "exon_id",
                         "Exon Chr Start (bp)":     "exon_start",
                         "Exon Chr End (bp)":       "exon_end",
                         "Genomic coding start":    "coding_start",
                         "Genomic coding end":      "coding_end",
                         "Exon Rank in Transcript": "exon_rank"}

required_bim_cols = {"SNP": "snp",
                     "CHR": "chr",
                     "BP":  "pos"}

snp_annotation = namedtuple("snp_annotation", annotation_categories)


def is_in_upstream_1kb(snp, pos_chr_mart_df):
    forward_upstream = ( (pos_chr_mart_df.strand == 1) &
            (snp.pos < pos_chr_mart_df.transcript_start) ).any()
    if not forward_upstream:
        reverse_upstream = ( (pos_chr_mart_df.strand == -1) &
                (snp.pos >= pos_chr_mart_df.transcript_end) ).any()
    return forward_upstream or reverse_upstream


def is_in_downstream_1kb(snp, pos_chr_mart_df):
    forward_downstream = ( (pos_chr_mart_df.strand == 1) &
            (snp.pos >= pos_chr_mart_df.transcript_end) ).any()
    if not forward_downstream:
        reverse_downstream = ( (pos_chr_mart_df.strand == -1) &
                (snp.pos < pos_chr_mart_df.transcript_start) ).any()
    return forward_downstream or reverse_downstream


def is_in_transcript(snp, pos_chr_mart_df):
    return ( (pos_chr_mart_df.transcript_start <= snp.pos) &
            (snp.pos < pos_chr_mart_df.transcript_end) ).any()


def is_in_exon(snp, pos_chr_mart_df):
    return ( (pos_chr_mart_df.exon_start <= snp.pos) &
            (snp.pos < pos_chr_mart_df.exon_end) ).any()


def is_in_utr5(snp, coding_pos_chr_mart_df):
    return ( (coding_pos_chr_mart_df.utr5_start <= snp.pos) &
            (snp.pos < coding_pos_chr_mart_df.utr5_end) ).any()


def is_in_utr3(snp, coding_pos_chr_mart_df):
    return ( (coding_pos_chr_mart_df.utr3_start <= snp.pos) &
            (snp.pos < coding_pos_chr_mart_df.utr3_end) ).any()


def is_in_coding(snp, coding_pos_chr_mart_df):
    return ( (coding_pos_chr_mart_df.coding_start <= snp.pos) &
            (snp.pos < coding_pos_chr_mart_df.coding_end) ).any()


def annotate(arg):
    chr_snp_df, chr_mart_df = arg
    annot_df = pd.DataFrame(columns=annotation_categories)
    for snp_row in chr_snp_df.itertuples():
        annot = dict.fromkeys(annotation_categories, False)
        i = ( (chr_mart_df.upstream_1kb <= snp_row.pos) &
              (snp_row.pos < chr_mart_df.downstream_1kb) )
        pos_chr_mart_df = chr_mart_df[i]
        if len(pos_chr_mart_df) != 0:
            annot["upstream_1kb"] = is_in_upstream_1kb(snp_row, pos_chr_mart_df)
            annot["downstream_1kb"] = is_in_downstream_1kb(snp_row, pos_chr_mart_df)
            annot["transcript"] = is_in_transcript(snp_row, pos_chr_mart_df)
            if annot["transcript"]:
                annot["exon"] = is_in_exon(snp_row, pos_chr_mart_df)
                annot["intron"] = not annot["exon"]
                coding_i = (pos_chr_mart_df.gene_type == "protein_coding")
                coding_pos_chr_mart_df = pos_chr_mart_df[coding_i]
                if annot["exon"] and len(pos_chr_mart_df) > 0:
                    annot["utr5"] = is_in_utr5(snp_row, coding_pos_chr_mart_df)
                    annot["utr3"] = is_in_utr3(snp_row, coding_pos_chr_mart_df)
                    annot["coding"] = is_in_coding(snp_row, coding_pos_chr_mart_df)
        annot_df.loc[snp_row.snp] = snp_annotation(**annot)
    return annot_df


def make_annotation_from_biomart(biomart_file, bim_file, out_txt, out_mat,
    test_run, test_n_snps, n_proc):
    mart_df = pd.read_table(biomart_file, usecols=list(required_biomart_cols))
    mart_df.columns = [required_biomart_cols[c] for c in mart_df.columns]
    print("%d exones in the input file" % len(mart_df))

    # change to 0-based coordinates
    mart_df.transcript_start -= 1
    mart_df.exon_start -= 1
    mart_df.coding_start -= 1

    mart_df = mart_df[mart_df.chr.isin(valid_chromosomes)]
    print("%d exones on valid chromosomes" % len(mart_df))
    mart_df["upstream_1kb"] = mart_df.transcript_start - 1000
    mart_df["downstream_1kb"] = mart_df.transcript_end + 1000

    #WARN: hardcoded CHR column
    snp_df = pd.read_table(bim_file, usecols=list(required_bim_cols),
        dtype={"CHR":str})
    snp_df.columns = [required_bim_cols[c] for c in snp_df.columns]
    print("%d snps in bim file" % len(snp_df))
    snp_df = snp_df[snp_df.chr.isin(valid_chromosomes)]
    print("%d snps on valid chromosomes" % len(snp_df))
    snp_df.drop_duplicates("snp", inplace=True)
    print("%d non duplicated snps on valid chromosomes" % len(snp_df))

    if test_run:
        # Test with random subset
        print("Taking random %d snps for testing" % test_n_snps)
        random_ind = np.random.permutation(len(snp_df))[:test_n_snps]
        snp_df = snp_df.loc[random_ind,:]

    arg_gen = ( (snp_df[snp_df.chr == c], mart_df[mart_df.chr == c])
            for c in valid_chromosomes )
    pool = Pool(processes=n_proc)
    annotation_dfs = pool.map(annotate, arg_gen)
    annot_df = pd.concat(annotation_dfs)
    print("%d SNPs were annotated" % len(annot_df))
    annot_df[annotation_categories] = annot_df[annotation_categories].astype(int)
    annot_df = annot_df.reindex(index=snp_df.snp)
    if not out_txt is None:
        annot_df.to_csv(out_txt, sep='\t', index_label="snp")
        print("%s saved" % out_txt)
    if not out_mat is None:
        mat_dict = {"annotations": annot_df.values}
        sio.savemat(out_mat, mat_dict, format="5", appendmat=False)
        print("%s saved" % out_mat)



if __name__ == "__main__":
    # Implementation notes:
    # - A file "biomart_GENCODE_basic.txt" was created using Ensembl Biomart
    #   tool. Using "Homo sapiens genes (GRCh37.p13)" dataset. Only transcripts
    #   included into GENCODE basic annotation were taken, i.e. the only fileter
    #   applied was: "GENCODE basic annotation: Only". All keys from
    #   required_biomart_cols dict presented above were taken as attributes.
    # - Definition of GENCODE basic annotation can be found here:
    #   http://grch37.ensembl.org/Help/Glossary?id=500
    # - Biomart has 1-based coordinate system, while dbSNP has 0-based
    #   coordinates. In this script everything is converted to 0-based.

    parser = argparse.ArgumentParser(description="Classify SNPs from reference "
        "template based on the biomart annotations.")
    parser.add_argument("--biomart", default="biomart_GENCODE_basic.txt.gz",
        type=str, help="File with Biomart annotations.")
    parser.add_argument("--ref", default="2558411_ref.bim",
        type=str, help="Reference template file.")
    parser.add_argument("--out-txt", default="annotations.txt", type=str,
        help="Output text file name or None.")
    parser.add_argument("--out-mat", default="annotations.mat", type=str,
        help="Output mat file name or None.")
    parser.add_argument("--test", action="store_true",
        help="Run test with randomly picked test_n_snps SNPs.")
    parser.add_argument("--test-n-snps", default=10000, type=int,
        help="Number of SNPs for testing.")
    parser.add_argument("--n-proc", default=1, type=int, help="Number of cores "
        "to use for calculation.")
    args = parser.parse_args()

    print("Started on %s" % strftime("%a, %d %b %Y %H:%M:%S", localtime()))

    make_annotation_from_biomart(args.biomart, args.ref, args.out_txt,
        args.out_mat, args.test, args.test_n_snps, args.n_proc)

    print("Finished on %s" % strftime("%a, %d %b %Y %H:%M:%S", localtime()))
