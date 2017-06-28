# Examples:
# python sumstats_convert.py csv tests\case01.txt   tmp.txt  --force --auto
# python sumstats_convert.py mat tests\1234_ref.bim tmp.txt  --force --traits test
import pandas as pd
import numpy as np
from itertools import permutations
from scipy import stats
import scipy.io as sio
import os
import sys
import argparse
import six
from sumstats_convert_utils import *
import collections
from shutil import copyfile
import zipfile
import glob

def parse_args(args):
    parser = argparse.ArgumentParser(description="Convert summary statistics "
        "file into mat file through the intermediate csv file.")
    subparsers = parser.add_subparsers()
    parser_csv = subparsers.add_parser("csv", help="Create csv file with standard columns")
    parser_csv.add_argument("sumstats_file", type=str, help="Input file with summary statistics")
    parser_csv.add_argument("output_file", type=str, help="Output csv file.")

    # Generate parameters from describe_cname.
    # Keys in describe_cname will be used as standard columns names for the resulting csv file.
    for cname in sorted(cols._asdict()):
        parser_csv.add_argument("--{}".format(cname.lower()), default=None, type=str, help=describe_cname[cname])

    parser_csv.add_argument("--auto", action="store_true", default=False,
        help="Auto-detect column types based on a set of standard column names.")
    parser_csv.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunk to read the file.")
    parser_csv.add_argument("--force", action="store_true", default=False,
        help="Force overwrite target files if they exist.")
    parser_csv.add_argument("--head", default=0, type=int,
        help="How many header lines of the file to print out for visual inspection (0 to disable)")
    parser_csv.add_argument("--preview", default=0, type=int,
        help="How many chunks to output into the output (debug option to preview large files that take long time to parse)")
    parser_csv.add_argument("--sep", default='\s+', type=str, choices=[',', ';', '\t', ' '],
        help="Delimiter to use (',' ';' $' ' or $'\\t'). By default uses delim_whitespace option in pandas.read_table.")
    parser_csv.add_argument("--all-snp-info-23-and-me", default=None, type=str,
        help="all_snp_info file for summary stats in 23-and-me format")
    parser_csv.set_defaults(func=make_csv)

    parser_mat = subparsers.add_parser("mat", help="Create mat files that can "
        "be used as an input for cond/conj FDR and for CM3 model.")
    parser_mat.add_argument("ref_file", type=str,
        help="Tab-separated file with list of referense SNPs.")
    parser_mat.add_argument("csv_files", type=str, nargs='+',
        help="Tab-sepatated csv files.")
    parser_mat.add_argument("--mat-files", type=str, nargs='+', default=None,
        help="Output mat files.")
    parser_mat.add_argument("--traits", type=str, nargs='+', default=None,
        help="Trait names that will be used in mat files.")
    parser_mat.add_argument("--effect", default='BETA', type=str, choices=['BETA', 'OR', 'Z', 'LOGODDS'],
        help="Effect column.")
    parser_mat.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunk to read the file.")
    parser_mat.add_argument("--force", action="store_true", default=False,
        help="Force overwrite target files if they exist.") 
    parser_mat.set_defaults(func=make_mat)

    parser_lift = subparsers.add_parser("lift", help="Lift RS numbers to a newer version of SNPdb, and/or liftover chr:pos to another genomic build using UCSC chain files.")
    parser_lift.add_argument("sumstats_file", type=str, help="Input file with summary statistics")
    parser_lift.add_argument("output_file", type=str, help="Output csv file.")
    parser_lift.add_argument("--chain-file", default=None, type=str,
        help="Chain file to use for CHR:BP conversion between genomic builds")
    parser_lift.add_argument("--snp-chrpos", default=None, type=str,
        help="NCBI SNPChrPosOnRef file.")
    parser_lift.add_argument("--snp-history", default=None, type=str,
        help="NCBI SNPHistory file.")
    parser_lift.add_argument("--rs-merge-arch", default=None, type=str,
        help="NCBI RsMergeArch file.")
    parser_lift.add_argument("--force", action="store_true", default=False,
        help="Force overwrite target files if they exist.")
    parser_lift.add_argument("--keep-bad-snps", action="store_true", default=False,
        help="Keep SNPs with underined rs# number or CHR:POS location in the output file.")
    parser_lift.add_argument("--na-rep", default='NA', type=str, choices=['NA', ''],
        help="Missing data representation.")
    parser_lift.add_argument("--gzip", action="store_true", default=False,
        help="A flag indicating whether to compress the resulting file with gzip.")
    parser_lift.set_defaults(func=make_lift)

    parser_rs = subparsers.add_parser("rs", help="Augument summary statistic file with SNP RS number from reference file. "
        "Merging is done on chromosome and position. The output is written back into the original file. ")
    parser_rs.add_argument("ref_file", type=str,
        help="Tab-separated file with list of referense SNPs.")
    parser_rs.add_argument("csv_files", type=str, nargs='+',
        help="Tab-sepatated csv files.")
    parser_rs.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunk to read the file.")
    parser_rs.set_defaults(func=make_rs)

    parser_ls = subparsers.add_parser("ls", help="Display information about columns present in standard sumstat files.")
    parser_ls.add_argument("path", type=str,
        help="File or regular expresion of files to show the info.")
    parser_ls.add_argument("--snp", action="store_true", default=False,
        help="Include information about how many SNPs are there in the file. This is slow because we need to scan the entire file.")
    parser_ls.set_defaults(func=make_ls)

    parser_mattocsv = subparsers.add_parser("mat-to-csv", help="Convert matlab .mat file with logpvec and zvec into CSV files.")
    parser_mattocsv.add_argument("ref_file", type=str,
        help="Tab-separated file with list of referense SNPs.")
    parser_mattocsv.add_argument("mat_files", type=str, nargs='+',
        help="Input .mat files.")
    parser_mattocsv.add_argument("--csv-files", type=str, nargs='+', default=None,
        help="Output .csv files.")
    parser_mattocsv.add_argument("--force", action="store_true", default=False,
        help="Force overwrite target files if they exist.")
    parser_mattocsv.add_argument("--na-rep", default='NA', type=str, choices=['NA', ''],
        help="Missing data representation.")
    parser_mattocsv.add_argument("--gzip", action="store_true", default=False,
        help="A flag indicating whether to compress the resulting file with gzip.")

    parser_mattocsv.set_defaults(func=mat_to_csv)

    parser_ldsctomat = subparsers.add_parser("ldsc-to-mat", help="Convert .sumstats, .ldscore, .M, .M_5_50 and binary .annot files from LD score regression to .mat files.")
    parser_ldsctomat.add_argument("ref_file", type=str,
        help="Tab-separated file with list of referense SNPs.")
    parser_ldsctomat.add_argument("out_file", type=str,
        help="Output .mat file.")
    parser_ldsctomat.add_argument("--sumstats", type=str, default=None,
        help="Name of .sumstats.gz file")
    parser_ldsctomat.add_argument("--ldscore", type=str, default=None,
        help="Name of .ldscore.gz files, where symbol @ indicates chromosome index. Example: baseline.@.l2.ldscore.gz")
    parser_ldsctomat.add_argument("--annot", type=str, default=None,
        help="Name of .annot.gz files, where symbol @ indicates chromosome index. Example: baseline.@.annot.gz")
    parser_ldsctomat.add_argument("--M", type=str, default=None,
        help="Name of .M files, where symbol @ indicates chromosome index. Example: baseline.@.l2.M")
    parser_ldsctomat.add_argument("--M-5-50", type=str, default=None,
        help="Name of .M_5_50 files, where symbol @ indicates chromosome index. Example: baseline.@.l2.M_5_50")
    parser_ldsctomat.add_argument("--force", action="store_true", default=False,
        help="Force overwrite target files if they exist.")
    parser_ldsctomat.set_defaults(func=ldsc_to_mat)

    return parser.parse_args(args)

### =================================================================================
###                          Implementation for parser_csv
### ================================================================================= 
def set_clean_args_cnames(args):
    """
    Inspect column names in user args, and clean them according to sumstats_convert_utils.clean_header()
    Raises an exception if either before or after cleaning some column names are duplicated.
    """
    cnames = [x.lower() for x in cols._asdict()]
    args_cnames = [args[x] for x in cnames if args[x] is not None]
    if len(args_cnames) != len(set(args_cnames)):
        raise(ValueError('Duplicated argument: {}'.format(find_duplicates(args_cnames))))

    for cname in cnames:
        if args[cname] is not None:
            args[cname] = clean_header(args[cname])
    args_clean_cnames = [args[x] for x in cnames if args[x] is not None]
    if len(args_clean_cnames) != len(set(args_clean_cnames)):
        raise(ValueError('Cleaning rules yield duplicated argument: {}'.format(find_duplicates(args_clean_cnames))))

def set_clean_file_cnames(df):
    """
    Inspect column names in pd.DataFrame, and clean them according to sumstats_convert_utils.clean_header()
    Raises an exception if either before or after cleaning some column names are duplicated.
    """
    file_cnames = df.columns
    if len(file_cnames) != len(set(file_cnames)):
        raise(ValueError('Unable to process input file due to duplicated column names'))
    
    clean_file_cnames = [clean_header(x) for x in file_cnames]
    if len(clean_file_cnames) != len(set(clean_file_cnames)):
        raise(ValueError('Cleaning column names resulted in duplicated column names: {}'.format(clean_file_cnames)))

    df.columns = clean_file_cnames

def find_duplicates(values):
    return [item for item, count in collections.Counter(values).items() if count > 1]

def find_auto_cnames(args, clean_file_cnames):
    """
    Auto-detect column using a set of default columns names 
    """
    
    cnames = [x.lower() for x in cols._asdict()]
    user_args = [args[x] for x in cnames if args[x] is not None]

    for default_cname, cname in default_cnames.items():
        # Ignore default cname if it is explicitly provided by the user
        if (cname.lower() not in args) or args[cname.lower()]:
            continue

        # Ignore default cname if it is not present among file columns
        if clean_header(default_cname) not in clean_file_cnames:
            continue
        
        # Ignore default cname if user took the column for something else
        if clean_header(default_cname) in user_args:
            continue

        args[cname.lower()] = clean_header(default_cname)

def check_input_file(file):
    if not os.path.isfile(file):
        raise ValueError("Input file does not exist: {f}".format(f=file))

def check_output_file(file, force=False):
    # Delete target file if user specifies --force option
    if force:
        try:
            os.remove(file)
        except OSError:
            pass

    # Otherwise raise an error if target file already exists
    if os.path.isfile(file) and not force:
        raise ValueError("Output file already exists: {f}".format(f=file))

    # Create target folder if it doesn't exist
    output_dir = os.path.dirname(file)
    if output_dir and not os.path.isdir(output_dir): os.makedirs(output_dir)  # ensure that output folder exists

def make_csv(args):
    """
    Based on file with summary statistics creates a tab-delimited csv file with standard columns.
    """
    check_input_file(args.sumstats_file)
    if args.all_snp_info_23_and_me: check_input_file(args.all_snp_info_23_and_me)
    set_clean_args_cnames(vars(args))
    check_output_file(args.output_file, args.force)

    print('Reading summary statistics file {}...'.format(args.sumstats_file))
    if args.head > 0: print_header(args.sumstats_file, lines=args.head)

    reader = pd.read_table(args.sumstats_file, dtype=str, sep=args.sep, chunksize=args.chunksize)
    reader_23_and_me = pd.read_table(args.all_snp_info_23_and_me, dtype=str, sep=args.sep, chunksize=args.chunksize) if args.all_snp_info_23_and_me else None
    n_snps = 0
    with open(args.output_file, 'a') as out_f:
        for chunk_index, chunk in enumerate(reader):
            if reader_23_and_me:
                chunk = pd.concat([chunk, next(reader_23_and_me)], axis=1)
                chunk = chunk.loc[:, ~chunk.columns.duplicated()]

            original_file_cname = chunk.columns
            set_clean_file_cnames(chunk)

            if chunk_index == 0:  # First chunk => analyze column names
                if args.auto: find_auto_cnames(vars(args), chunk.columns)
            
                # Find the map from (cleaned) column name to a standard cname (as it will appear in the resulting file)
                cnames = [x.lower() for x in cols._asdict()]
                cname_map = {vars(args)[cname] : cname.upper() for cname in cnames if vars(args)[cname] is not None}

                # Report any missing columns that user requested to put in the resulting file
                cname_missing = [x for x in cname_map if x not in chunk.columns]
                if cname_missing: raise(ValueError('Columns {} are missing in the input file'.format(cname_missing)))

                # Describe the mapping between old and new column names that we are going to perform
                print('Interpret column names as follows:')
                for original in original_file_cname:
                    cname = cname_map.get(clean_header(original))
                    print("\t{o} : {d} ({e})".format(o=original, d=cname, e="Will be deleted" if not cname else describe_cname[cname]))
                if not cname_map: raise(ValueError('Arguments imply to delete all columns from the input file. Did you forget --auto flag?'))

            chunk.drop([x for x in chunk.columns if x not in cname_map], axis=1, inplace=True)
            chunk.rename(columns=cname_map, inplace=True)

            # Split CHR:POS column into two
            if cols.CHRPOS in chunk.columns:
                chunk[cols.CHR], chunk[cols.BP] = chunk[cols.CHRPOS].str.split(':', 1).str
                chunk.drop(cols.CHRPOS, axis=1, inplace=True)

            # Split A1/A2 column into two
            if cols.A1A2 in chunk.columns:
                chunk[cols.A1], chunk[cols.A2] = chunk[cols.A1A2].str.split('/', 1).str
                chunk.drop(cols.A1A2, axis=1, inplace=True)

            # Ensure standard labels in CHR column
            if cols.CHR in chunk.columns:
                chunk[cols.CHR].fillna(-9, inplace=True)
                chunk[cols.CHR] = format_chr(chunk[cols.CHR])

            # Ensure that alleles are coded as capital letters
            if cols.A1 in chunk.columns: chunk[cols.A1] = chunk[cols.A1].str.upper()
            if cols.A2 in chunk.columns: chunk[cols.A2] = chunk[cols.A2].str.upper()

            chunk = chunk.sort_index(axis=1)
            chunk.to_csv(out_f, index=False, header=(chunk_index==0), sep='\t', na_rep='NA')
            n_snps += len(chunk)
            print("{n} lines processed".format(n=(chunk_index+1)*args.chunksize))
            if args.preview and (chunk_index+1) >= args.preview:
                print('Abort reading input file due to --preview flag.')
                break

        print("{n} SNPs saved to {f}".format(n=n_snps, f=args.output_file))

### =================================================================================
###                          Implementation for parser_mat
### ================================================================================= 
def get_str_list_sign(str_list):
    return np.array([-1 if e[0]=='-' else 1 for e in str_list], dtype=np.int)

_base_complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
def _complement(seq):
    return "".join([_base_complement[b] for b in seq])

ALLELES = "AGCT"
COMPLEMENT = {''.join(s): _complement(s) for s in permutations(ALLELES, 2)}
AMBIGUOUS = [aa for aa, cc in COMPLEMENT.items() if aa == cc[::-1]]

def make_mat(args):
    """
    Takes csv files (created with the csv task of this script).
    Require columns: SNP, P, and one of the signed summary statistics columns (BETA, OR, Z, LOGODDS).
    Creates corresponding mat files which can be used as an input for the conditional fdr model.
    Only SNPs from the reference file are considered. Zscores of ambiguous SNPs are set to NA.
    """
    if args.mat_files is None:
        args.mat_files = []
        for csv_f in args.csv_files:
            mat_f = os.path.splitext(csv_f)[0] + ".mat"
            args.mat_files.append(mat_f)
    if args.traits is None:
        args.traits = []
        for csv_f in args.csv_files:
            args.traits.append('')
    if (len(args.csv_files) != len(args.mat_files) or
            len(args.csv_files) != len(args.traits)):
        raise ValueError(("Number of input csv files should be equal to the "
            "number of output mat files and number of traits!"))
    check_input_file(args.ref_file)
    for csv_f, mat_f in zip(args.csv_files, args.mat_files):
        check_input_file(csv_f)
        check_output_file(mat_f, args.force)

    # if signed_effect is true, take effect column as string to handle correctly
    # case of truncated numbers, e.g.: 0.00 and -0.00 should have different sign
    signed_effect = False if args.effect == cols.OR else True
    effect_col_dtype = str if signed_effect else np.float

    print('Reading reference file {}...'.format(args.ref_file))
    usecols = [cols.SNP, cols.A1, cols.A2]
    reader = pd.read_table(args.ref_file, sep='\t', usecols=usecols,
        chunksize=args.chunksize)
    ref_dict = {}
    for chunk in reader:
        gtypes = (chunk[cols.A1] + chunk[cols.A2]).apply(str.upper)
        #TODO?: add check whether some id is already in ref_dict
        ref_dict.update(dict(zip(chunk[cols.SNP], gtypes)))
    ref_dict = {i: (aa, COMPLEMENT[aa], aa[::-1], COMPLEMENT[aa[::-1]])
            for i, aa in ref_dict.items()}
    ref_snps = pd.read_table(args.ref_file, sep='\t', usecols=[cols.SNP], squeeze=True)
    #TODO?: add check whether ref_snps contains duplicates
    print("Reference dict contains {d} snps.".format(d=len(ref_dict)))

    for csv_f, mat_f, trait in zip(args.csv_files, args.mat_files, args.traits):
        print('Reading summary statistics file {}...'.format(csv_f))
        reader = pd.read_table(csv_f, sep='\t', usecols=[cols.A1, cols.A2, cols.SNP, cols.PVAL, args.effect],
                               chunksize=args.chunksize, dtype={args.effect:effect_col_dtype}, float_precision='high')
        df_out = None
        for i, chunk in enumerate(reader):
            chunk = chunk.loc[chunk[cols.SNP].isin(ref_dict),:]
            if chunk.empty: continue
            gtypes = (chunk[cols.A1] + chunk[cols.A2]).apply(str.upper)
            # index of SNPs that have the same alleles as indicated in reference
            ind = [gt in ref_dict[sid]
                for sid, gt in zip(chunk[cols.SNP], gtypes)]
            chunk = chunk.loc[ind,:]
            gtypes = gtypes[ind]
            log10pv = -np.log10(chunk[cols.PVAL].values)
            # not_ref_effect = [
            #   1 if effect allele in data == other allele in reference
            #   -1 if effect allele in data == effect allele in reference ]
            # So zscores with positive effects will be positive and zscores with 
            # negative effects will stay negative, since
            # stats.norm.ppf(chunk[cols.PVAL]*0.5) is always negetive (see zvect
            # calculation below).
            not_ref_effect = np.array([-1 if gt in ref_dict[sid][:2] else 1
                for sid, gt in zip(chunk[cols.SNP], gtypes)])
            #TODO: check proportion of positive and negative effects
            if signed_effect:
                # effect column has str type
                # -1 if effect starts with '-' else 1
                effect_sign = get_str_list_sign(chunk[args.effect])
            else:
                # effect column has np.float type
                # 1 if effect >=1 else -1
                if (chunk[args.effect] < 0).any():
                    raise ValueError("OR column contains negative values")
                effect_sign = np.sign(chunk[args.effect].values - 1)
                effect_sign[effect_sign == 0] = 1
            effect_sign *= not_ref_effect
            zvect = stats.norm.ppf(chunk[cols.PVAL].values*0.5)*effect_sign
            ind_ambiguous = [j for j,gt in enumerate(gtypes) if gt in AMBIGUOUS]
            # set zscore of ambiguous SNPs to nan
            zvect[ind_ambiguous] = np.nan
            #TODO: check whether output df contains duplicated rs-ids (warn)
            df = pd.DataFrame({"pvalue": log10pv, "zscore":zvect}, index=chunk[cols.SNP])
            if df_out is None: df_out = df
            else: df_out = df_out.append(df)
            print("{n} lines processed, {m} SNPs matched with reference file".format(n=(i+1)*args.chunksize, m=len(df_out)))
        if df_out.empty: raise(ValueError("No SNPs match after joining with reference data"))
        dup_index = df_out.index.duplicated(keep=False)
        if dup_index.any():
            print("Duplicated SNP ids detected:")
            print(df_out[dup_index])
            print("Keeping only the first occurance.")
        df_out = df_out[~df_out.index.duplicated(keep='first')]

        print('Writing .mat file...')
        df_ref_aligned = pd.DataFrame(columns=["pvalue", "zscore"], index=ref_snps)
        df_ref_aligned["pvalue"] = df_out["pvalue"]
        df_ref_aligned["zscore"] = df_out["zscore"]
        save_dict = {"logpvec"+trait: df_ref_aligned["pvalue"].values, "zvec"+trait: df_ref_aligned["zscore"].values}
        sio.savemat(mat_f, save_dict, format='5', do_compression=False,
            oned_as='column', appendmat=False)
        print("%s created" % mat_f)

### =================================================================================
###                          Implementation for parser_lift
### =================================================================================
def make_lift(args):
    """
    Lift RS numbers to a newer version of SNPdb, and/or
    liftover chr:pos to another genomic build using NCBI chain files.
    """
    from pyliftover import LiftOver
    from lift_rs_numbers import LiftRsNumbers

    check_input_file(args.sumstats_file)
    check_output_file(args.output_file, args.force)

    if args.chain_file is not None: check_input_file(args.chain_file)
    if args.snp_chrpos is not None: check_input_file(args.snp_chrpos)
    if args.snp_history is not None: check_input_file(args.snp_history)
    if args.rs_merge_arch is not None: check_input_file(args.rs_merge_arch)

    if (args.snp_history is not None) != (args.rs_merge_arch is not None):
        raise(ValueError('--snp-history and --rs-merge-arch must be used together'))

    print('Reading summary statistics file {}...'.format(args.sumstats_file))
    df = pd.read_table(args.sumstats_file, sep='\t')
    print('Done, {} markers found'.format(len(df)))

    lift_bp = None; lift_rs = None; snp_chrpos = None

    if (args.chain_file is not None) and (cols.CHR in df) and (cols.BP in df):
        print('Reading {}...'.format(args.chain_file))
        lift_bp = LiftOver(args.chain_file)

    if (args.snp_history is not None) and (cols.SNP in df):
        lift_rs = LiftRsNumbers(hist_file=args.snp_history, merge_file=args.rs_merge_arch)

    if (args.snp_chrpos is not None) and (cols.SNP in df) and (cols.CHR in df) and (cols.BP in df):
        print('Reading {}...'.format(args.snp_chrpos))
        snp_chrpos = pd.read_table(args.snp_chrpos, sep='\t', header=None, usecols=[0,1,2])
        snp_chrpos.columns=['snp_id','chr','pos'] #,'orien','neighbor_snp_list','isPAR']
        snp_chrpos.dropna(subset=['pos'],inplace=True)           # drop NA positions
        snp_chrpos['pos']=snp_chrpos['pos'].astype(np.int) + 1   # convert to integer unity-based positions
        snp_chrpos['chr'].replace({'X':23, 'Y':24, 'PAR': 25, 'MT': 26}, inplace=True)   # standarize chr labels
        snp_chrpos['chr']=snp_chrpos['chr'].astype(np.int)
        print('Done, found {} entries'.format(len(snp_chrpos)))
        # if snp_chrpos.duplicated(subset=['snp_id'], keep=False).sum() != 0:
        #   raise('Duplicated snp_id found in SNPChrPosOnRef file')

    indices_with_old_chrpos = range(len(df))  # indices with original chr:pos
    fixes = []
    if cols.SNP in df:
        if lift_rs is not None:
            # Fix1 brings forward SNP rs# numbers and set SNP rs# to None for SNPs found in SNPHistory table
            df[cols.SNP], stats = lift_rs.lift(df[cols.SNP])
            fixes.append('{} rs# numbers has changed based on RsMergeArch table'.format(stats['lifted']))
            print(stats)

        if (snp_chrpos is not None) and (cols.CHR in df) and (cols.BP in df):
            # Fix2 set chr:pos based on SNPChrPosOnRef table (only applies to SNPs with rs# number, not to other markers)
            df['SNP_LIFTED'] = [(int(x[2:]) if (x and x.startswith('rs') and x[2:].isdigit()) else -1) for x in df[cols.SNP]]
            df = pd.merge(df, snp_chrpos, how='left', left_on='SNP_LIFTED', right_on='snp_id')
            print('Warning: there are {} SNPs with a valid RS number, but not found in SNPChrPosOnRef'.format(((df['SNP_LIFTED'] != -1) & df['snp_id'].isnull()).sum()))

            idx = ((df['pos'] != df[cols.BP]) | (df['chr'] != df[cols.CHR])) & ~df['pos'].isnull() & ~df['chr'].isnull()
            df.ix[idx, cols.BP] = df.ix[idx, 'pos'].astype(int)
            df.ix[idx, cols.CHR] = df.ix[idx, 'chr'].astype(int)
            fixes.append('{} markers receive new CHR:POS based on SNPChrPosOnRef table'.format(idx.sum()))

            indices_with_old_chrpos = [i for (i, b) in enumerate(df['pos'].isnull() | df['chr'].isnull()) if b]
            df.drop(['SNP_LIFTED', 'snp_id', 'chr', 'pos'], axis=1, inplace=True)

    if (lift_bp is not None) and (cols.CHR in df) and (cols.BP in df):
        # Fix3 lifts chr:pos to another genomic build. This step is OFF by default, only applies if user provided chain file.
        # Note that lifting with pyliftover is rather slow, so we apply this step only to markers that are not in SNPChrPosOnRef table.
        print('Lift CHR:POS for {} SNPs to another genomic build...'.format(len(indices_with_old_chrpos)))
        unique = 0; multi = 0; failed = 0
        for i, index in enumerate(indices_with_old_chrpos):
            if (i+1) % 100 == 0: print('Finish {} SNPs'.format(i+1))
            chri = int(df.ix[index, cols.CHR]); bp = int(df.ix[index, cols.BP]); snp = df.ix[index, cols.SNP]
            lifted = lift_bp.convert_coordinate('chr{}'.format(chri), bp)
            if (lifted is None) or (len(lifted) == 0):
                #print('Unable to lift SNP {} at chr{}:{}, delete'.format(snp, chri, bp))
                df.ix[index, cols.CHR] = None
                df.ix[index, cols.BP] = None
                failed += 1
                continue
            if len(lifted) > 1:
                print('Warning: SNP {} at chr{}:{} lifts to multiple position, use first.'.format(snp, chri, bp))
                multi += 1
            if len(lifted) == 1:
                unique += 1

            df.ix[index, cols.CHR] = int(lifted[0][0][3:])
            df.ix[index, cols.BP] = lifted[0][1]
        print('Done, {} failed, {} unique, {} multi'.format(failed, unique, multi))
        fixes.append('{} markers receive new CHR:POS based on liftover chain files'.format(unique + multi))

    if not args.keep_bad_snps:
        if lift_rs and (cols.SNP in df):
            df_len = len(df)
            df.dropna(subset=[cols.SNP], inplace=True)                  # Fix4 due to SNP being deleted from dbSNP
            fixes.append("{n} markers were dropped based on SNPHistory table".format(n = df_len - len(df)))

        if (cols.CHR in df) and (cols.BP in df):
            if lift_bp:
                df_len = len(df)
                df.dropna(subset=[cols.CHR, cols.BP], inplace=True)     # Fix5, due to failed liftover across genomic builds
                fixes.append("{n} markers were dropped due to failed CHR:POS lift".format(n = df_len - len(df)))
            df[cols.CHR] = df[cols.CHR].astype(int)
            df[cols.BP] = df[cols.BP].astype(int)

    df.to_csv(args.output_file + ('.gz' if args.gzip else ''),
        index=False, header=True, sep='\t', na_rep=args.na_rep,
        compression='gzip' if args.gzip else None)

    print("{n} SNPs saved to {f}".format(n=len(df), f=args.output_file))
    print('Summary: \n\t{}'.format("\n\t".join(fixes)))

### =================================================================================
###                          Implementation for parser_rs
### =================================================================================
def make_rs(args):
    """
    Augument summary statistic file with SNP RS number from reference file.
    Merging is done on chromosome and position.
    The output is written back into the original file.
    """
    check_input_file(args.ref_file)
    for csv_f in args.csv_files:
        check_input_file(csv_f)
        check_output_file(csv_f + '.tmp', force=True)
        check_output_file(csv_f + '.tmp2', force=True)

    print('Reading reference file {}...'.format(args.ref_file))
    ref_file = pd.read_table(args.ref_file, sep='\t', usecols=[cols.SNP, cols.CHR, cols.BP])
    print("Reference dict contains {d} snps.".format(d=len(ref_file)))

    for csv_f in args.csv_files:
        print('Reading summary statistics file {}...'.format(csv_f))
        copyfile(csv_f, csv_f + ".tmp2")  # don't open the original file (otherwise it might be hard to remove it afterwards)
        reader = pd.read_table(csv_f + ".tmp2", sep='\t', chunksize=args.chunksize)
        n_snps = 0
        with open(csv_f + '.tmp', 'a') as out_f:
            for chunk_index, chunk in enumerate(reader):
                if cols.SNP in chunk: chunk.drop(cols.SNP, axis=1, inplace=True)
                chunk = pd.merge(chunk, ref_file, how='left', on=[cols.CHR, cols.BP])
                chunk = chunk.sort_index(axis=1)
                chunk.to_csv(out_f, index=False, header=(chunk_index==0), sep='\t', na_rep='NA')
                n_snps += len(chunk)
                print("{n} lines processed".format(n=(chunk_index+1)*args.chunksize))
        os.remove(csv_f)
        os.rename(csv_f + '.tmp', csv_f)
        check_output_file(csv_f + '.tmp2', force=True)  # try remove, ignore if fails
        print("{n} SNPs saved to {f}".format(n=n_snps, f=csv_f))

### =================================================================================
###                          Implementation for make_ls
### =================================================================================
def make_ls(args):
    ml = max([len(os.path.basename(file)) for file in glob.glob(args.path)])
    print('{f}\t{n}\t{c}'.format(f='file'.ljust(ml),n='#snp'.ljust(9),c='\t'.join(sorted(cols._asdict()))))
    num_snps = 'n/a'
    for file in glob.glob(args.path):
        if not os.path.isfile(file): continue
        if args.snp: num_snps = sum(1 for line in get_compression_and_open(file))-1
        for chunk in pd.read_table(file, sep='\t', chunksize=1):
            print('{f}\t{n}\t{c}'.format(f=os.path.basename(file).ljust(ml), c='\t'.join([str(x in chunk) for x in sorted(cols._asdict())]),n=str(num_snps).ljust(9)))
            break
    print('Columns description:')
    for cname in sorted(cols._asdict()):
        print('{c}\t{d}'.format(c=cname, d=describe_cname[cname]))

### =================================================================================
###                          Implementation for mat_to_csv
### =================================================================================
def mat_to_csv(args):
    check_input_file(args.ref_file)
    if args.csv_files is None:
        args.csv_files = []
        for mat_f in args.mat_files:
            csv_f = os.path.splitext(mat_f)[0] + ".csv"
            args.csv_files.append(csv_f)
    for mat_f, csv_f in zip(args.mat_files, args.csv_files):
        check_input_file(mat_f)
        check_output_file(csv_f, args.force)

    print('Reading reference file {}...'.format(args.ref_file))
    ref_file = pd.read_table(args.ref_file, sep='\t', usecols=[cols.SNP, cols.A1, cols.A2])
    print("Reference dict contains {d} snps.".format(d=len(ref_file)))

    for mat_f, csv_f in zip(args.mat_files, args.csv_files):
        print('Reading .mat file {}...'.format(mat_f))
        df = ref_file.copy()
        sumstats = sio.loadmat(mat_f)
        for key in sumstats.keys():
            if key.lower().startswith('zvec'): df['Z'] = sumstats[key]
            if key.lower().startswith('logpvec'): df['PVAL'] = np.power(10, -sumstats[key])
            if key.lower().startswith('nvec'): df['N'] = sumstats[key]
        df.to_csv(csv_f + '.gz' if args.gzip else '',
                  index=False, header=True, sep='\t', na_rep=args.na_rep,
                  compression='gzip' if args.gzip else None)

### =================================================================================
###                          Implementation for ldsc_to_mat
### =================================================================================
def ldsc_to_mat(args):
    check_input_file(args.ref_file)
    if args.sumstats: check_input_file(args.sumstats)
    for chri in range(1, 23):
        if args.ldscore: check_input_file(args.ldscore.replace('@', str(chri)))
        if args.annot: check_input_file(args.annot.replace('@', str(chri)))
        if args.M: check_input_file(args.M.replace('@', str(chri)))
        if args.M_5_50: check_input_file(args.M_5_50.replace('@', str(chri)))
    check_output_file(args.out_file, args.force)

    print('Reading reference file {}...'.format(args.ref_file))
    ref_file = pd.read_table(args.ref_file, sep='\t', usecols=[cols.SNP, cols.A1, cols.A2])
    print("Reference dict contains {d} snps.".format(d=len(ref_file)))

    save_dict = {}
    if args.sumstats:
        df_sumstats = pd.read_csv(args.sumstats, sep='\t')
        df_len = len(df_sumstats)
        print('Read {} SNPs from --sumstats file'.format(df_len))

        df_sumstats = pd.merge(ref_file, df_sumstats, how='left', on=cols.SNP)
        df_sumstats = df_sumstats.dropna(how='any')
        print('Removed {} SNPs not in --ref file'.format(df_len - len(df_sumstats)))
        df_len = len(df_sumstats)

        df_sumstats['ALLELES'] = df_sumstats.A1_x + df_sumstats.A2_x + df_sumstats.A1_y + df_sumstats.A2_y
        df_sumstats = df_sumstats[filter_alleles(df_sumstats['ALLELES'])].copy()
        print('Removed {} SNPs with alleles not matching --ref file'.format(df_len - len(df_sumstats)))
        print('{} SNPs remain'.format(len(df_sumstats)))

        df_z = df_sumstats.Z.copy()
        df_sumstats.Z = align_alleles(df_sumstats.Z, df_sumstats['ALLELES'])
        print('Flip z score for {} SNPs'.format((df_sumstats.Z != df_z).sum()))

        df_sumstats.drop(['A1_x', 'A2_x', 'A1_y', 'A2_y','ALLELES'], axis=1, inplace=True)
        df_sumstats = pd.merge(ref_file, df_sumstats, how='left', on=cols.SNP)
        save_dict['zvec'] = df_sumstats["Z"].values
        save_dict['nvec'] = df_sumstats["N"].values

    if args.ldscore:
        df_ldscore = pd.concat([pd.read_csv(args.ldscore.replace('@', str(chri)), delim_whitespace=True)  for chri in range(1, 23)])
        df_ldscore.drop([x for x in ['CHR', 'BP', 'CM', 'MAF'] if x in df_ldscore], inplace=True, axis=1)
        print('Shape of ldscore file: {shape}'.format(shape=df_ldscore.shape))
        df_ldscore = pd.merge(ref_file[['SNP']], df_ldscore, how='left', on='SNP')
        del df_ldscore['SNP']
        print('Shape of ldscore file after merge: {shape}'.format(shape=df_ldscore.shape))
        save_dict['annonames'] = list(df_ldscore.columns)
        save_dict['annomat'] = df_ldscore.values

    if args.annot:
        df_annot = pd.concat([pd.read_csv(args.annot.replace('@', str(chri)), delim_whitespace=True)  for chri in range(1, 23)])
        df_annot.drop([x for x in ['CHR', 'BP', 'CM', 'MAF'] if x in df_annot], inplace=True, axis=1)
        print('Shape of annots file: {shape}'.format(shape=df_annot.shape))
        df_annot = pd.merge(ref_file[['SNP']], df_annot, how='left', on='SNP')
        del df_annot['SNP']
        print('Shape of annots file after merge: {shape}'.format(shape=df_annot.shape))
        save_dict['annonames_bin'] = list(df_annot.columns)
        save_dict['annomat_bin'] = df_annot.values

    if args.M_5_50:
        m_5_50 = pd.concat([pd.read_csv(args.M_5_50.replace('@', str(chri)), delim_whitespace=True, header=None) for chri in range(1, 23)])
        m_5_50 = np.atleast_2d(m_5_50.sum().values)
        print('M_5_50={}'.format(m_5_50))
        save_dict['M_5_50']=m_5_50

    if args.M:
        m = pd.concat([pd.read_csv(args.M.replace('@', str(chri)), delim_whitespace=True, header=None) for chri in range(1, 23)])
        m = np.atleast_2d(m.sum().values)
        print('M={}'.format(m))
        save_dict['M']=m

    sio.savemat(args.out_file, save_dict, format='5', do_compression=False, oned_as='column', appendmat=False)
    print('Result written to {f}'.format(f=args.out_file))

### =================================================================================
###                                Main section
### ================================================================================= 
if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    args.func(args)
    print("Done")

    # Examples
    # (TBD: this section must be fixed)
    
    # csv task
    # python sumstats_convert.py csv ../data/pgc/adhd_y/daner_ADD9_0712_info08_snps_maf001.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_adhd_y_2m_info08_snp_maf001.csv --id SNP --effect OR --pval P --effectA A1 --otherA A2 --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/pgc/adhd_meta_20160919/archive/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_adhd_2016_2m.csv --id SNP --effect OR --pval P --effectA A1 --otherA A2 --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/23andME/archive/sumstats_23andMe_to_pgc.csv.gz ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/23andMe_adhd_2m.csv --signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/ssgac/EduYears_Main.txt.gz ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/ssgac_eduyears_2m.csv --id MarkerName --effect Beta --pval  Pval --effectA A1 --otherA A2 --signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/ssgac/SWB_Full.txt.gz ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/ssgac_swb_2m.csv --id MarkerName --effect Beta --pval  Pval --effectA A1 --otherA A2 --signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/ssgac/Neuroticism_Full.txt.gz ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/ssgac_neuro_2m.csv --id MarkerName --effect Beta --pval  Pval --effectA A1 --otherA A2 --signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/ssgac/DS_Full.txt.gz ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/ssgac_ds_2m.csv --id MarkerName --effect Beta --pval  Pval --effectA A1 --otherA A2 --signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/gpc/GPC-2.EXTRAVERSION.full.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/gpc_extraversion_2m.csv --id RSNUMBER --effect Effect_nw --pval PVALUE --effectA Allele1 --otherA Allele2 --signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/gpc/GPC-2.NEUROTICISM.full_header.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/gpc_neuroticism_2m.csv --id RSNUMBER --effect Effect_nw --pval PVALUE --effectA Allele1 --otherA Allele2 --signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/pgc/cross_disorder/pgc.cross.ADD4.2013-05_tabs.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_adhd_2013_2m.csv --effect or --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/pgc/cross_disorder/pgc.cross.BIP11.2013-05_tabs.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_bip_2013_2m.csv --effect or --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/pgc/cross_disorder/pgc.cross.AUT8.2013-05_tabs.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_aut_2013_2m.csv --effect or --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/pgc/cross_disorder/pgc.cross.MDD9.2013-05_tabs.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_mdd_2013_2m.csv --effect or --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/pgc/cross_disorder/pgc.cross.SCZ17.2013-05_tabs.txt ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_scz_2013_2m.csv --effect or --ref-id SNP --ref-a1 A1 --ref-a2 A2
    # python sumstats_convert.py csv ../data/ssgac/SSGAC_College_Rietveld2013_publicrelease.txt  ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/ssgac_college2013_2m.csv --id MarkerName --effect OR --pval Pvalue --effectA Effect_Allele --otherA Other_Allele --ref-id SNP --ref-a1 A1 --ref-a2 A2

    # mat task
    # python sumstats_convert.py mat ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_adhd_y_2m_info08_snp_maf001.csv --ref-id SNP --traits adhd_y
    # python sumstats_convert.py mat ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_adhd_2016_2m.csv --ref-id SNP --traits adhd_2016
    # python sumstats_convert.py mat ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/23andMe_adhd_2m.csv --ref-id SNP --traits adhd_23andMe
    # python sumstats_convert.py mat ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/gpc_extraversion_2m.csv ../tmp/gpc_neuroticism_2m.csv --ref-id SNP --traits extra_gpc neuro_gpc
    # python sumstats_convert.py mat ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/ssgac_eduyears_2m.csv ../tmp/ssgac_swb_2m.csv ../tmp/ssgac_neuro_2m.csv ../tmp/ssgac_ds_2m.csv --ref-id SNP --traits eduyears_ssgac swb_ssgac neuro_ssgac ds_ssgac
    # python sumstats_convert.py mat ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/pgc_adhd_2013_2m.csv ../tmp/pgc_bip_2013_2m.csv ../tmp/pgc_aut_2013_2m.csv ../tmp/pgc_mdd_2013_2m.csv ../tmp/pgc_scz_2013_2m.csv --ref-id SNP --traits adhd_2013 bip_2013 aut_2013 mdd_2013 scz_2013
    # python sumstats_convert.py mat ../projects/adhd/cond_fdr/2m_template/2558411_ref.bim ../tmp/ssgac_college2013_2m.csv --ref-id SNP --traits college2013_ssgac
