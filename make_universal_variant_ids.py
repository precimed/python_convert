"""
Input: whitespace-delimited csv file with CHR, BP, A1, A2 and ID columns.
Output: tab-separated csv file with two columns ID and UID, where ID column is taken from the input file,
    UID column is constructed as CHR:BP:AA1:AA2, where CHR and BP are from the input file, AA1 is min(A1, A2, A1_complementary, A2_complementary),
    AA2 = A2 if AA1 == A1,
    AA2 = A1 if AA1 == A2,
    AA2 = A2_complementary if AA1 == A1_complementary,
    AA2 = A1_complementary if AA1 == A2_complementary,
    min is taken based on lexicographical order.
    If either A1 or A2 contains non-ATGC char, original ID is retained, i.e. UID = ID. 

Example:
python make_universal_variant_ids.py --fname /cluster/projects/p33/users/alexeas/hrc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz \
        --chr "#CHROM" --bp POS --a1 REF --a2 ALT --id ID --out hrc.hg19.uid.txt

Resulting hrc.hg19.uid.txt file will contain two columns (no header): (1) original ID from the input file, (2) constructed univeral ID (UID)

Assume you need to map different type of variant ids between two different files:
    FILE1 has rsid
    FILE2 has chr:bp
    and both FILE1 and FILE2 contain CHR, BP, A1, A2 columns with coordinates in the same genomic build.
Then you can apply the script to each FILE1 and FILE2 separately to generate UIDs for each file. For FILE1 you will get FILE1.uid and for FILE2 you'll get FILE2.uid.
Then you can get the mapping between rsid IDs from FILE1 and chr:bp IDs from FILE2 using:
join -t$'\t' -1 2 -2 2 <(sort -k2,2 FILE1.uid) <(sort -k2,2 FILE2.uid) > FILE1_FILE2.uid
Resulting FILE1_FILE2.uid file will have three columns: (1) univeral ID (2) rsid from FILE1, (3) chr:bp from FILE2.
"""

import pandas as pd
import argparse

COMPL_DICT = {"A":"T", "T":"A", "G":"C", "C":"G"}
DEL_ACGT_TT = str.maketrans({c:"" for c in "ATGC"})
STD_COL_NAMES = ["CHR", "BP", "A1", "A2", "ID"]

def std_format(df, col_names, std_col_names):
    # standardize format:
    # - retain only relevant columns (col_names) in the specified order
    # - rename columns to standard names, std_col_names[i] should be an std col name for col_names[i]
    # - A1 and A2 to uppercase
    df = df[col_names]
    col_rename_dict = dict(zip(col_names, std_col_names))
    df.rename(columns=col_rename_dict, inplace=True)
    df["A1"] = df["A1"].str.upper()
    df["A2"] = df["A2"].str.upper()
    return df
    
def reverse_compl(seq):
    # seq is an uppercase string.
    return COMPL_DICT[seq] if len(seq) == 1 else "".join([COMPL_DICT[b] for b in seq][::-1])

def get_uid_col(df):
    # df is DataFrame with standard column names with CHR, BP, A1, A2 and ID columns in the corresponding order.
    # A1 and A2 must be capitalized.
    uid_col = []
    for chrom, bp, a1, a2, vid in df.itertuples(index=False):
        if a1.translate(DEL_ACGT_TT) == "" and a2.translate(DEL_ACGT_TT) == "":
            min_orig, max_orig = (a1, a2) if a1 < a2 else (a2, a1)
            a1c, a2c = reverse_compl(a1), reverse_compl(a2)
            min_compl, max_compl = (a1c, a2c) if a1c < a2c else (a2c, a1c)
            a1u, a2u = (min_orig, max_orig) if min_orig < min_compl else (min_compl, max_compl)
            uid = f"{chrom}:{bp}:{a1u}:{a2u}"
        else:
            uid = vid
        uid_col.append(uid)
    return uid_col


# Parse arguments --------------------------
parser = argparse.ArgumentParser(description="Constract universal variant IDs (see detailed description and example of use in the script file).")
parser.add_argument("--fname", help="Path to input file with chromosome, position, allele 1, allele 2 and variant ID columns.")
parser.add_argument("--chr", default="CHR", help="Chromosome column.")
parser.add_argument("--bp", default="BP", help="Position column.")
parser.add_argument("--a1", default="A1", help="Allele 1 column.")
parser.add_argument("--a2", default="A2", help="Allele 2 column.")
parser.add_argument("--id", default="ID", help="Variant ID column.")
parser.add_argument("--out", help="Output file name.")
args = parser.parse_args()


# Main -------------------------------------
col_names = [args.chr, args.bp, args.a1, args.a2, args.id] # the order should correspond to the order in STD_COL_NAMES
df = pd.read_csv(args.fname, delim_whitespace=True, usecols=col_names, dtype=str)
print(f"{df.shape[0]} variants loaded from {args.fname}")
df = std_format(df, col_names, STD_COL_NAMES)
uid_col = get_uid_col(df)
df["UID"] = uid_col

df[["ID", "UID"]].to_csv(args.out, sep='\t', index=False, header=False)
print(f"{args.out} saved.")

