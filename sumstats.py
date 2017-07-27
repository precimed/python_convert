#!/usr/bin/env python
'''
(c) 2016 Oleksandr Frei and Alexey A. Shadrin
Various utilities for GWAS summary statistics.
'''

import pandas as pd
import numpy as np
from itertools import permutations
from scipy import stats
import scipy.io as sio
import os
import time, sys, traceback
import argparse
import six
from sumstats_utils import *
import collections
import re
from shutil import copyfile
import zipfile
import glob
import socket
import getpass

__version__ = '1.0.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* sumstats.py: utilities for GWAS summary statistics\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2016-2017 Oleksandr Frei and Alexey A. Shadrin\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

def parse_args(args):
    parser = argparse.ArgumentParser(description="A collection of various utilities for GWAS summary statistics.")
    subparsers = parser.add_subparsers()

    # 'csv' utility : load raw summary statistics file and convert it into a standardized format
    parser_csv = subparsers.add_parser("csv",
        help='Load raw summary statistics file and convert it into a standardized format: '
        'tab-separated file with standard column names, standard chromosome labels, NA label for missing data, etc. '
        'The conversion does not change the number of lines in the input files (e.g. no filtering is done on markers). '
        'Unrecognized columns are removed from the summary statistics file. '
        'The remaining utilities in sumstats.py work with summary statistics files in the standardized format.')

    parser_csv.add_argument("--sumstats", type=str, help="[required] Raw input file with summary statistics.")
    parser_csv.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_csv.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    # Generate parameters from describe_cname.
    # Keys in describe_cname will be used as standard columns names for the resulting csv file.
    for cname in sorted(cols._asdict()):
        parser_csv.add_argument("--{}".format(cname.lower()), default=None, type=str, help=describe_cname[cname])

    parser_csv.add_argument("--auto", action="store_true", default=False,
        help="Auto-detect column types based on a set of standard column names.")
    parser_csv.add_argument("--ignore", type=str, nargs='+',
        help="List of column names in the original file to ignore during auto-detection")
    parser_csv.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunk to read the file.")
    parser_csv.add_argument("--head", default=0, type=int,
        help="How many header lines of the file to print out for visual inspection (0 to disable)")
    parser_csv.add_argument("--preview", default=0, type=int,
        help="How many chunks to output into the output (debug option to preview large files that take long time to parse)")
    parser_csv.add_argument("--skip-validation", action="store_true", default=False,
        help="Skip validation of the resulting csv file")
    parser_csv.add_argument("--sep", default='\s+', type=str, choices=[',', ';', '\t', ' '],
        help="Delimiter to use (',' ';' $' ' or $'\\t'). By default uses delim_whitespace option in pandas.read_table.")
    parser_csv.add_argument("--all-snp-info-23-and-me", default=None, type=str,
        help="all_snp_info file for summary stats in 23-and-me format")
    parser_csv.set_defaults(func=make_csv)

    # 'qc' utility: miscellaneous quality control and filtering procedures
    parser_qc = subparsers.add_parser("qc",
        help="Miscellaneous quality control and filtering procedures")

    parser_qc.add_argument("--sumstats", type=str, help="[required] Input file with summary statistics in standardized format")
    parser_qc.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_qc.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_qc.add_argument("--exclude-ranges", type=str, nargs='+',
        help='Exclude SNPs in ranges of base pair position, for example MHC. '
        'The syntax is chr:from-to, for example 6:25000000-35000000. Multiple regions can be excluded.')
    parser_qc.add_argument("--dropna-cols", type=str, nargs='+',
        help='List of column names. SNPs with missing values in either of the columns will be excluded.')
    parser_qc.add_argument("--fix-dtype-cols", type=str, nargs='+',
        help='List of column names. Ensure appropriate data type for the columns (CHR, BP - int, PVAL - float, etc)')
    parser_qc.add_argument("--max-or", type=float, default=None,
        help='Filter SNPs with OR (odds ratio) exceeding specified threshold')
    parser_qc.set_defaults(func=make_qc)

    # 'mat' utility: load summary statistics into matlab format
    parser_mat = subparsers.add_parser("mat", help="Create mat files that can "
        "be used as an input for cond/conj FDR and for CM3 model. "
        "Takes csv files (created with the csv task of this script). "
        "Require columns: SNP, P, and one of the signed summary statistics columns (BETA, OR, Z, LOGODDS). "
        "Creates corresponding mat files which can be used as an input for the conditional fdr model. "
        "Only SNPs from the reference file are considered. Zscores of strand ambiguous SNPs are set to NA. "
        "To use CHR:POS for merging summary statistics with reference file consider 'rs' utility "
        "which auguments summary statistics with SNP column (first run 'sumstats.py rs ...', "
        "then feed the resulting file into sumstats.py mat ...)")

    parser_mat.add_argument("--sumstats", type=str, help="[required] Input file with summary statistics in standardized format.")
    parser_mat.add_argument("--ref", type=str, help="[required] Tab-separated file with list of referense SNPs.")
    parser_mat.add_argument("--out", type=str, help="[required] File to output the result. File should end with .mat extension.")
    parser_mat.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_mat.add_argument("--trait", type=str, default='',
        help="Trait name that will be used in mat file. Can be kept empty, in this case the variables will be named 'logpvec', 'zvec' and 'nvec'")
    parser_mat.add_argument("--effect", default=None, type=str, choices=['BETA', 'OR', 'Z', 'LOGODDS'],
        help="Effect column. Default is to auto-detect. In case if multiple effect columns are present in the input file"
        " a warning will be shown and the first column is taken according to the priority list: "
        "Z (highest priority), BETA, OR, LOGODDS (lowest priority).")
    parser_mat.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunk to read the file.")
    parser_mat.set_defaults(func=make_mat)

    # 'lift' utility: lift RS numbers to a newer version of SNPdb, and/or liftover chr:pos to another genomic build using UCSC chain files
    parser_lift = subparsers.add_parser("lift",
        help="Lift RS numbers to a newer version of SNPdb, "
        "and/or liftover chr:pos to another genomic build using UCSC chain files. "
        "WARNING: this utility may use excessive amount of memory (up and beyong 32 GB of RAM).")

    parser_lift.add_argument("--sumstats", type=str, help="[required] Input file with summary statistics in standardized format")
    parser_lift.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_lift.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_lift.add_argument("--chain-file", default=None, type=str,
        help="Chain file to use for CHR:BP conversion between genomic builds")
    parser_lift.add_argument("--snp-chrpos", default=None, type=str,
        help="NCBI SNPChrPosOnRef file.")
    parser_lift.add_argument("--snp-history", default=None, type=str,
        help="NCBI SNPHistory file.")
    parser_lift.add_argument("--rs-merge-arch", default=None, type=str,
        help="NCBI RsMergeArch file.")
    parser_lift.add_argument("--keep-bad-snps", action="store_true", default=False,
        help="Keep SNPs with undefined rs# number or CHR:POS location in the output file.")
    parser_lift.add_argument("--na-rep", default='NA', type=str, choices=['NA', ''],
        help="Missing data representation.")
    parser_lift.add_argument("--gzip", action="store_true", default=False,
        help="A flag indicating whether to compress the resulting file with gzip.")
    parser_lift.set_defaults(func=make_lift)

    # 'rs' utility: augument summary statistic file with SNP RS number from reference file
    parser_rs = subparsers.add_parser("rs",
        help="Augument summary statistic file with SNP RS number from reference file. "
        "Merging is done on chromosome and position.")

    parser_rs.add_argument("--sumstats", type=str, help="[required] Input file with summary statistics in standardized format")
    parser_rs.add_argument("--ref", type=str, help="[required] Tab-separated file with list of referense SNPs.")
    parser_rs.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_rs.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_rs.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunk to read the file.")
    parser_rs.set_defaults(func=make_rs)

    # 'ls' utility: display information about columns of a standardized summary statistics file
    parser_ls = subparsers.add_parser("ls",
        help="Report information about standard sumstat files, "
        "including the set of columns available, number of SNPs, etc.")

    parser_ls.add_argument("--path", type=str, help="[required] File or regular expresion of the files to include in the report.")
    parser_ls.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_ls.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_ls.add_argument("--snp", action="store_true", default=False,
        help="Include information about how many SNPs are there in the file. "
        "This option is slow because causes a complete scan of all input files.")
    parser_ls.set_defaults(func=make_ls)

    # 'mat-to-csv' utility: convert matlab .mat file with logpvec and zvec into CSV files
    parser_mattocsv = subparsers.add_parser("mat-to-csv",
        help="Convert matlab .mat file with logpvec, zvec and (optionally) nvec into CSV files.")

    parser_mattocsv.add_argument("--mat", type=str, help="[required] Input mat file.")
    parser_mattocsv.add_argument("--ref", type=str, help="[required] Tab-separated file with list of referense SNPs.")
    parser_mattocsv.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_mattocsv.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_mattocsv.add_argument("--na-rep", default='NA', type=str, choices=['NA', ''],
        help="Missing data representation.")
    parser_mattocsv.add_argument("--gzip", action="store_true", default=False,
        help="A flag indicating whether to compress the resulting file with gzip.")
    parser_mattocsv.set_defaults(func=mat_to_csv)

    # 'ldsc-to-mat' utility: convert data from LD score regression formats to .mat files
    parser_ldsctomat = subparsers.add_parser("ldsc-to-mat",
        help="Convert .sumstats, .ldscore, .M, .M_5_50 and "
        "binary .annot files from LD score regression to .mat files.")

    parser_ldsctomat.add_argument("--ref", type=str, help="[required] Tab-separated file with list of referense SNPs.")
    parser_ldsctomat.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_ldsctomat.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

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
    parser_ldsctomat.set_defaults(func=ldsc_to_mat)

    # 'diff-mat' utility: compare two .mat files with logpvec, zvec and nvec, and report the differences
    parser_diffmat = subparsers.add_parser("diff-mat",
        help="Compare two .mat files with logpvec, zvec and nvec, "
        "and report the differences.")

    parser_diffmat.add_argument("--mat1", type=str, default=None, help="[required] Name of the first .mat file")
    parser_diffmat.add_argument("--mat2", type=str, default=None, help="[required] Name of the second .mat file")
    parser_diffmat.add_argument("--ref", type=str, help="[required] Tab-separated file with list of referense SNPs.")
    parser_diffmat.add_argument("--out", type=str, help="[required] File to output the result.")
    parser_diffmat.add_argument("--force", action="store_true", default=False, help="Allow sumstats.py to overwrite output file if it exists.")

    parser_diffmat.add_argument("--sumstats", type=str, default=None,
        help="Optionally, the name of the source summary statistics file in standardized .csv format. "
        "Assuming that both .mat files originate from this file diff-mat will produce output file "
        "to help investigate where the differences came from.")
    parser_diffmat.set_defaults(func=diff_mat)

    return parser.parse_args(args)

### =================================================================================
###                          Implementation for parser_csv
### ================================================================================= 
def set_clean_args_cnames(args):
    """
    Inspect column names in user args, and clean them according to sumstats_utils.clean_header()
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
    Inspect column names in pd.DataFrame, and clean them according to sumstats_utils.clean_header()
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

        # Ignore default cname if user explicitly asked to ignore it
        if args['ignore'] and (clean_header(default_cname) in [clean_header(x) for x in args['ignore']]):
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

def make_csv(args, log):
    """
    Based on file with summary statistics creates a tab-delimited csv file with standard columns.
    """
    check_input_file(args.sumstats)
    if args.all_snp_info_23_and_me: check_input_file(args.all_snp_info_23_and_me)
    set_clean_args_cnames(vars(args))
    check_output_file(args.out, args.force)

    log.log('Reading summary statistics file {}...'.format(args.sumstats))
    if args.head > 0:
        log.log('File header:')
        header = get_header(args.sumstats, lines=args.head)
        for line in header: log.log(line)

    reader = pd.read_table(args.sumstats, dtype=str, sep=args.sep, chunksize=args.chunksize)
    reader_23_and_me = pd.read_table(args.all_snp_info_23_and_me, dtype=str, sep=args.sep, chunksize=args.chunksize) if args.all_snp_info_23_and_me else None
    n_snps = 0
    with open(args.out, 'a') as out_f:
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
                log.log('Interpret column names as follows:')
                for original in original_file_cname:
                    cname = cname_map.get(clean_header(original))
                    log.log("\t{o} : {d} ({e})".format(o=original, d=cname, e="Will be deleted" if not cname else describe_cname[cname]))
                if not cname_map: raise(ValueError('Arguments imply to delete all columns from the input file. Did you forget --auto flag?'))

                final_cols = set(cname_map.values())  # final list of columns in the resulting file
                if (cols.CHRPOS not in final_cols) and (cols.CHR not in final_cols): log.log('Warning: CHR column ({}) is not found'.format(describe_cname[cols.CHR]))
                if (cols.CHRPOS not in final_cols) and (cols.BP not in final_cols): log.log('Warning: BP column ({}) is not found'.format(describe_cname[cols.BP]))
                if cols.SNP not in final_cols: log.log('Warning: SNP column ({}) is not found'.format(describe_cname[cols.SNP]))
                if cols.PVAL not in final_cols: log.log('Warning: PVAL column ({}) is not found'.format(describe_cname[cols.PVAL]))
                if (cols.A1 not in final_cols) and (cols.A1A2 not in final_cols): log.log('Warning: A1 column ({}) is not found'.format(describe_cname[cols.A1]))
                if (cols.A2 not in final_cols) and (cols.A1A2 not in final_cols): log.log('Warning: A2 column ({}) is not found'.format(describe_cname[cols.A2]))
                effect_size_column_count = int(cols.Z in final_cols) + int(cols.OR in final_cols) + int(cols.BETA in final_cols) + int(cols.LOGODDS in final_cols)
                if effect_size_column_count == 0: log.log('Warning: None of the columns indicate effect direction: typically either BETA, OR, LOGODDS or Z column is expected')
                if effect_size_column_count > 1: log.log('Warning: Multiple columns indicate effect direction: typically only one of BETA, OR, LOGODDS and Z columns is expected')

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
            print("{f}: {n} lines processed".format(f=args.sumstats, n=(chunk_index+1)*args.chunksize))
            if args.preview and (chunk_index+1) >= args.preview:
                log.log('Abort reading input file due to --preview flag.')
                break

        if n_snps == 0: raise(ValueError('Input summary stats file appears to be empty.'))
        log.log("Done. {n} SNPs saved to {f}".format(n=n_snps, f=args.out))

    if not args.skip_validation:
        log.log('Validate the resulting file...')
        reader = pd.read_table(args.out, sep='\t', chunksize=args.chunksize)
        n_snps = 0
        for chunk_index, chunk in enumerate(reader):
            if chunk_index==0:
                log.log('Column types:\n\t' + '\n\t'.join([column + ':' + str(dtype) for (column, dtype) in zip(chunk.columns, chunk.dtypes)]))
            n_snps += len(chunk)
        log.log("Done. {n} SNPs read from {f}".format(n=n_snps, f=args.out))

### =================================================================================
###                          Implementation for parser_qc
### ================================================================================= 
def make_qc(args, log):
    check_input_file(args.sumstats)
    check_output_file(args.out, args.force)

    # Interpret --exclude-ranges input
    ChromosomeRange = collections.namedtuple('ChromosomeRange', ['chr', 'from_bp', 'to_bp'])
    exclude_ranges = []
    if args.exclude_ranges is not None:
        for exclude_range in args.exclude_ranges:
            try:
                range = ChromosomeRange._make([int(x) for x in exclude_range.replace(':', ' ').replace('-', ' ').split()[:3]])
            except Exception as e:
                raise(ValueError('Unable to interpret exclude range "{}", chr:from-to format is expected.'.format(exclude_range)))
            exclude_ranges.append(range)
            log.log('Exclude SNPs on chromosome {} from BP {} to {}'.format(range.chr, range.from_bp, range.to_bp))

    if args.dropna_cols is None: args.dropna_cols = []
    if args.fix_dtype_cols is None: args.fix_dtype_cols = []

    if args.exclude_ranges is not None:
        args.dropna_cols.extend(['CHR', 'BP'])
        args.fix_dtype_cols.extend(['CHR', 'BP'])
    if args.fix_dtype_cols is not None:
         for col in args.fix_dtype_cols:
            if cols_type_map[col] == int:
                args.dropna_cols.append(col)
    if args.max_or is not None:
        args.fix_dtype_cols.append('OR')

    log.log('Reading sumstats file {}...'.format(args.sumstats))
    sumstats = pd.read_table(args.sumstats, sep='\t', dtype=str)
    log.log("Sumstats file contains {d} markers.".format(d=len(sumstats)))

    if len(args.dropna_cols) > 0:
        sumstats_len = len(sumstats)
        sumstats.dropna(subset=args.dropna_cols, inplace=True)
        if len(sumstats) != sumstats_len: log.log('Drop {} markers because of missing values'.format(sumstats_len - len(sumstats)))
    for col in args.fix_dtype_cols:
        if col in sumstats:
            log.log('Set column {} dtype to {}'.format(col, cols_type_map[col]))
            if cols_type_map[col] == float or cols_type_map[col] == np.float64:
                sumstats[col] = pd.to_numeric(sumstats[col], errors='coerce')
                if col in args.dropna_cols:
                    sumstats_len = len(sumstats)
                    sumstats.dropna(subset=[col], inplace=True)
                    if len(sumstats) != sumstats_len: log.log('Drop {} markers after dtype conversion'.format(sumstats_len - len(sumstats)))
            else:
                sumstats[col] = sumstats[col].astype(cols_type_map[col])
    for range in exclude_ranges:
        idx = (sumstats[cols.CHR] == range.chr) & (sumstats[cols.BP] >= range.from_bp) & (sumstats[cols.BP] < range.to_bp)
        if idx.sum() > 0:
            log.log('Exclude {} SNPs in range {}:{}-{}'.format(idx.sum(), range.chr, range.from_bp, range.to_bp))
            sumstats = sumstats[~idx].copy()
    if args.max_or is not None and cols.OR in sumstats:
        sumstats_len = len(sumstats)
        sumstats.drop(sumstats.OR > args.max_or, inplace=True)
        if len(sumstats) != sumstats_len: log.log('Drop {} markers because OR exceeded threshold'.format(sumstats_len - len(sumstats)))
    sumstats.to_csv(args.out, index=False, header=True, sep='\t', na_rep='NA')
    log.log("{n} SNPs saved to {f}".format(n=len(sumstats), f=args.out))

### =================================================================================
###                          Implementation for parser_mat
### ================================================================================= 
def get_str_list_sign(str_list):
    return np.array([-1 if e[0]=='-' else 1 for e in list(str_list)], dtype=np.int)

_base_complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
def _complement(seq):
    return "".join([_base_complement[b] for b in seq])

ALLELES = "AGCT"
COMPLEMENT = {''.join(s): _complement(s) for s in permutations(ALLELES, 2)}
AMBIGUOUS = [aa for aa, cc in COMPLEMENT.items() if aa == cc[::-1]]

def make_mat(args, log):
    """
    Takes csv files (created with the csv task of this script).
    Require columns: SNP, P, and one of the signed summary statistics columns (BETA, OR, Z, LOGODDS).
    Creates corresponding mat files which can be used as an input for the conditional fdr model.
    Only SNPs from the reference file are considered. Zscores of strand ambiguous SNPs are set to NA.
    """
    check_input_file(args.ref)
    check_input_file(args.sumstats)
    check_output_file(args.out, args.force)

    columns = list(pd.read_table(args.sumstats, sep='\t', nrows=0).columns)
    log.log('Columns in {}: {}'.format(args.sumstats, columns))
    if args.effect is None:
        if cols.Z in columns: args.effect = cols.Z
        elif cols.BETA in columns: args.effect = cols.BETA
        elif cols.OR in columns: args.effect = cols.OR
        elif cols.LOGODDS in columns: args.effect = cols.LOGODDS
        else: raise(ValueError('Signed effect column is not detected in {}'.format(args.sumstats)))
        effect_size_column_count = np.sum([int(c in columns) for c in [cols.Z, cols.BETA, cols.OR, cols.LOGODDS]])
        if effect_size_column_count > 1: log.log('Warning: Multiple columns indicate effect direction')
        log.log('Use {} column as effect direction.'.format(args.effect))
    missing_columns = [c for c in [cols.A1, cols.A2, cols.SNP, cols.PVAL, args.effect] if c not in columns]
    if missing_columns: raise(ValueError('{} columns are missing'.format(missing_columns)))

    # if signed_effect is true, take effect column as string to handle correctly
    # case of truncated numbers, e.g.: 0.00 and -0.00 should have different sign
    signed_effect = False if args.effect == cols.OR else True
    effect_col_dtype = str if signed_effect else np.float

    log.log('Reading reference file {}...'.format(args.ref))
    usecols = [cols.SNP, cols.A1, cols.A2]
    reader = pd.read_table(args.ref, sep='\t', usecols=usecols,
        chunksize=args.chunksize)
    ref_dict = {}
    for chunk in reader:
        gtypes = (chunk[cols.A1] + chunk[cols.A2]).apply(str.upper)
        #TODO?: add check whether some id is already in ref_dict
        ref_dict.update(dict(zip(chunk[cols.SNP], gtypes)))
    ref_dict = {i: (aa, COMPLEMENT[aa], aa[::-1], COMPLEMENT[aa[::-1]])
            for i, aa in ref_dict.items()}
    ref_snps = pd.read_table(args.ref, sep='\t', usecols=[cols.SNP], squeeze=True)
    #TODO?: add check whether ref_snps contains duplicates
    log.log("Reference dict contains {d} snps.".format(d=len(ref_dict)))

    log.log('Reading summary statistics file {}...'.format(args.sumstats))
    reader = pd.read_table(args.sumstats, sep='\t', usecols=[cols.A1, cols.A2, cols.SNP, cols.PVAL, args.effect],
                            chunksize=args.chunksize, dtype={args.effect:effect_col_dtype}, float_precision='high')
    df_out = None
    for i, chunk in enumerate(reader):
        if i==0: log.log('Column types:\n\t' + '\n\t'.join([column + ':' + str(dtype) for (column, dtype) in zip(chunk.columns, chunk.dtypes)]))
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
            effect_sign = get_str_list_sign(chunk[args.effect].astype(str))
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
        print("{f}: {n} lines processed, {m} SNPs matched with reference file".format(f=args.sumstats, n=(i+1)*args.chunksize, m=len(df_out)))
    if df_out.empty: raise(ValueError("No SNPs match after joining with reference data"))
    dup_index = df_out.index.duplicated(keep=False)
    if dup_index.any():
        log.log("Duplicated SNP ids detected:")
        log.log(df_out[dup_index])
        log.log("Keeping only the first occurance.")
    df_out = df_out[~df_out.index.duplicated(keep='first')]

    log.log('Writing .mat file...')
    df_ref_aligned = pd.DataFrame(columns=["pvalue", "zscore"], index=ref_snps)
    df_ref_aligned["pvalue"] = df_out["pvalue"]
    df_ref_aligned["zscore"] = df_out["zscore"]
    save_dict = {"logpvec"+args.trait: df_ref_aligned["pvalue"].values, "zvec"+args.trait: df_ref_aligned["zscore"].values}
    sio.savemat(args.out, save_dict, format='5', do_compression=False,
        oned_as='column', appendmat=False)
    log.log("%s created" % args.out)

### =================================================================================
###                          Implementation for parser_lift
### =================================================================================
def make_lift(args, log):
    """
    Lift RS numbers to a newer version of SNPdb, and/or
    liftover chr:pos to another genomic build using NCBI chain files.
    """
    from pyliftover import LiftOver
    from lift_rs_numbers import LiftRsNumbers

    check_input_file(args.sumstats)
    check_output_file(args.out, args.force)

    if args.chain_file is not None: check_input_file(args.chain_file)
    if args.snp_chrpos is not None: check_input_file(args.snp_chrpos)
    if args.snp_history is not None: check_input_file(args.snp_history)
    if args.rs_merge_arch is not None: check_input_file(args.rs_merge_arch)

    if (args.snp_history is not None) != (args.rs_merge_arch is not None):
        raise(ValueError('--snp-history and --rs-merge-arch must be used together'))

    log.log('Reading summary statistics file {}...'.format(args.sumstats))
    df = pd.read_table(args.sumstats, sep='\t')
    log.log('Done, {} markers found'.format(len(df)))

    lift_bp = None; lift_rs = None; snp_chrpos = None

    if (args.chain_file is not None) and (cols.CHR in df) and (cols.BP in df):
        log.log('Reading {}...'.format(args.chain_file))
        lift_bp = LiftOver(args.chain_file)

    if (args.snp_history is not None) and (cols.SNP in df):
        lift_rs = LiftRsNumbers(hist_file=args.snp_history, merge_file=args.rs_merge_arch)

    if args.snp_chrpos is not None:
        log.log('Reading {}...'.format(args.snp_chrpos))
        snp_chrpos = pd.read_table(args.snp_chrpos, sep='\t', header=None, usecols=[0,1,2])
        snp_chrpos.columns=['snp_id','chr','pos'] #,'orien','neighbor_snp_list','isPAR']
        snp_chrpos.dropna(subset=['pos'],inplace=True)           # drop NA positions
        snp_chrpos['pos']=snp_chrpos['pos'].astype(np.int) + 1   # convert to integer unity-based positions
        snp_chrpos['chr'].replace({'X':23, 'Y':24, 'PAR': 25, 'MT': 26}, inplace=True)   # standarize chr labels
        snp_chrpos['chr']=snp_chrpos['chr'].astype(np.int)
        log.log('Done, found {} entries'.format(len(snp_chrpos)))
        # if snp_chrpos.duplicated(subset=['snp_id'], keep=False).sum() != 0:
        #   raise('Duplicated snp_id found in SNPChrPosOnRef file')

    indices_with_old_chrpos = range(len(df))  # indices with original chr:pos
    fixes = []
    if (cols.SNP in df) and (lift_rs is not None):
        # Fix1 brings forward SNP rs# numbers and set SNP rs# to None for SNPs found in SNPHistory table
        df[cols.SNP], stats = lift_rs.lift(df[cols.SNP])
        fixes.append('{} rs# numbers has changed based on RsMergeArch table'.format(stats['lifted']))
        log.log(stats)

    if (cols.SNP in df) and (snp_chrpos is not None):
        # Fix2 set chr:pos based on SNPChrPosOnRef table (only applies to SNPs with rs# number, not to other markers)
        df['SNP_LIFTED'] = [(int(x[2:]) if (x and x.startswith('rs') and x[2:].isdigit()) else -1) for x in df[cols.SNP]]
        df = pd.merge(df, snp_chrpos, how='left', left_on='SNP_LIFTED', right_on='snp_id')
        log.log('Warning: there are {} SNPs with a valid RS number, but not found in SNPChrPosOnRef'.format(((df['SNP_LIFTED'] != -1) & df['snp_id'].isnull()).sum()))

        if cols.CHR in df and cols.BP in df:
            idx = ((df['pos'] != df[cols.BP]) | (df['chr'] != df[cols.CHR])) & ~df['pos'].isnull() & ~df['chr'].isnull()
            df.ix[idx, cols.BP] = df.ix[idx, 'pos'].astype(int)
            df.ix[idx, cols.CHR] = df.ix[idx, 'chr'].astype(int)
            fixes.append('{} markers receive new CHR:POS based on SNPChrPosOnRef table'.format(idx.sum()))
        else:
            idx = ~df['pos'].isnull() & ~df['chr'].isnull()
            df[cols.BP] = np.nan; df[cols.CHR] = np.nan
            df.ix[idx, cols.BP] = df.ix[idx, 'pos'].astype(int)
            df.ix[idx, cols.CHR] = df.ix[idx, 'chr'].astype(int)
            fixes.append('{} markers receive CHR:POS based on SNPChrPosOnRef table'.format(idx.sum()))

        indices_with_old_chrpos = [i for (i, b) in enumerate(df['pos'].isnull() | df['chr'].isnull()) if b]
        df.drop(['SNP_LIFTED', 'snp_id', 'chr', 'pos'], axis=1, inplace=True)

    if (cols.CHR in df) and (cols.BP in df) and (cols.SNP not in df) and (snp_chrpos is not None):
        # Fix3 set SNP rs# based on SNPChrPosOnRef table based on CHR:POS
        df = pd.merge(df, snp_chrpos, how='left', left_on=[cols.CHR, cols.BP], right_on=['chr', 'pos'])
        df['snp_id'][df['snp_id'].isnull()] = -1
        df['snp_id'] = 'rs' + df['snp_id'].astype(int).astype(str)
        df['snp_id'][df['snp_id'] == 'rs-1'] = None
        df['SNP'] = df['snp_id']
        fixes.append('{} markers receive SNP rs# based on SNPChrPosOnRef table'.format((~df[cols.SNP].isnull()).sum()))
        df.drop(['snp_id', 'chr', 'pos'], axis=1, inplace=True)

    if lift_bp is not None:
        # Fix4 lifts chr:pos to another genomic build. This step is OFF by default, only applies if user provided chain file.
        # Note that lifting with pyliftover is rather slow, so we apply this step only to markers that are not in SNPChrPosOnRef table.
        log.log('Lift CHR:POS for {} SNPs to another genomic build...'.format(len(indices_with_old_chrpos)))
        unique = 0; multi = 0; failed = 0
        for i, index in enumerate(indices_with_old_chrpos):
            if (i+1) % 100 == 0: print('Finish {} SNPs'.format(i+1))
            chri = int(df.ix[index, cols.CHR]); bp = int(df.ix[index, cols.BP]); snp = df.ix[index, cols.SNP]
            lifted = lift_bp.convert_coordinate('chr{}'.format(chri), bp)
            if (lifted is None) or (len(lifted) == 0):
                #log.log('Unable to lift SNP {} at chr{}:{}, delete'.format(snp, chri, bp))
                df.ix[index, cols.CHR] = None
                df.ix[index, cols.BP] = None
                failed += 1
                continue
            if len(lifted) > 1:
                log.log('Warning: SNP {} at chr{}:{} lifts to multiple position, use first.'.format(snp, chri, bp))
                multi += 1
            if len(lifted) == 1:
                unique += 1

            df.ix[index, cols.CHR] = int(lifted[0][0][3:])
            df.ix[index, cols.BP] = lifted[0][1]
        log.log('Done, {} failed, {} unique, {} multi'.format(failed, unique, multi))
        fixes.append('{} markers receive new CHR:POS based on liftover chain files'.format(unique + multi))

    if not args.keep_bad_snps:
        if cols.SNP in df:
            df_len = len(df)
            df.dropna(subset=[cols.SNP], inplace=True)                  # Fix5 due to SNP being deleted from dbSNP
            if len(df) < df_len:
                fixes.append("{n} markers were dropped based on SNPHistory table and/or SNPChrPosOnRef table".format(n = df_len - len(df)))

        if (cols.CHR in df) and (cols.BP in df):
            df_len = len(df)
            df.dropna(subset=[cols.CHR, cols.BP], inplace=True)     # Fix6, due to failed liftover across genomic builds
            if len(df) < df_len:
                fixes.append("{n} markers were dropped due to missing CHR:POS information or due to failed CHR:POS lift".format(n = df_len - len(df)))
            df[cols.CHR] = df[cols.CHR].astype(int)
            df[cols.BP] = df[cols.BP].astype(int)

    df.to_csv(args.out + ('.gz' if args.gzip else ''),
        index=False, header=True, sep='\t', na_rep=args.na_rep,
        compression='gzip' if args.gzip else None)

    log.log("{n} SNPs saved to {f}".format(n=len(df), f=args.out))
    log.log('Summary: \n\t{}'.format("\n\t".join(fixes)))

### =================================================================================
###                          Implementation for parser_rs
### =================================================================================
def make_rs(args, log):
    """
    Augument summary statistic file with SNP RS number from reference file.
    Merging is done on chromosome and position.
    """
    check_input_file(args.ref)
    check_input_file(args.sumstats)
    check_output_file(args.out, force=args.force)

    log.log('Reading reference file {}...'.format(args.ref))
    ref_file = pd.read_table(args.ref, sep='\t', usecols=[cols.SNP, cols.CHR, cols.BP])
    log.log("Reference dict contains {d} snps.".format(d=len(ref_file)))

    log.log('Reading summary statistics file {}...'.format(args.sumstats))
    reader = pd.read_table(args.sumstats, dtype=str, sep='\t', chunksize=args.chunksize)
    n_snps = 0
    with open(args.out, 'a') as out_f:
        for chunk_index, chunk in enumerate(reader):
            if cols.SNP in chunk: chunk.drop(cols.SNP, axis=1, inplace=True)
            chunk.BP = chunk.BP.astype(int)
            chunk.CHR = chunk.CHR.astype(int)
            chunk = pd.merge(chunk, ref_file, how='left', on=[cols.CHR, cols.BP])
            chunk = chunk.sort_index(axis=1)
            chunk.to_csv(out_f, index=False, header=(chunk_index==0), sep='\t', na_rep='NA')
            n_snps += len(chunk)
            print("{f}: {n} lines processed".format(f=args.sumstats, n=(chunk_index+1)*args.chunksize))
    log.log("{n} SNPs saved to {f}".format(n=n_snps, f=args.out))

### =================================================================================
###                          Implementation for make_ls
### =================================================================================
def make_ls(args, log):
    ml = max([len(os.path.basename(file)) for file in glob.glob(args.path)])
    cols_list = [x for x in cols._asdict() if x not in ['A1A2', 'CHRPOS']]
    log.log('{f}\t{n}\t{c}'.format(f='file'.ljust(ml),n='#snp'.ljust(9),c='\t'.join([x.replace('NCONTROL', 'NCONT.') for x in cols_list])))
    num_snps = 'n/a'
    for file in glob.glob(args.path):
        if not os.path.isfile(file): continue
        if args.snp: num_snps = sum(1 for line in get_compression_and_open(file))-1
        for chunk in pd.read_table(file, sep='\t', chunksize=1):
            log.log('{f}\t{n}\t{c}'.format(f=os.path.basename(file).ljust(ml), c='\t'.join([('YES' if x in chunk else '-') for x in cols_list]),n=str(num_snps).ljust(9)))
            break
    log.log('Columns description:')
    for cname in sorted(cols._asdict()):
        log.log('{c}\t{d}'.format(c=cname, d=describe_cname[cname]))

### =================================================================================
###                          Implementation for mat_to_csv
### =================================================================================
def mat_to_csv(args, log):
    check_input_file(args.ref)
    check_input_file(args.mat)
    check_output_file(args.out, args.force)

    log.log('Reading reference file {}...'.format(args.ref))
    ref_file = pd.read_table(args.ref, sep='\t', usecols=[cols.SNP, cols.A1, cols.A2])
    log.log("Reference dict contains {d} snps.".format(d=len(ref_file)))

    log.log('Reading .mat file {}...'.format(args.mat))
    df = ref_file.copy()
    sumstats = sio.loadmat(args.mat)
    for key in sumstats.keys():
        if key.lower().startswith('zvec'):
            df['Z'] = sumstats[key]
            log.log('Found zvec, {} non-nan values.'.format(np.sum(~np.isnan(sumstats[key]))))
        if key.lower().startswith('logpvec'):
            df['PVAL'] = np.power(10, -sumstats[key])
            log.log('Found logpvec, {} non-nan values.'.format(np.sum(~np.isnan(sumstats[key]))))
        if key.lower().startswith('nvec'):
            df['N'] = sumstats[key]
            log.log('Found nvec, {} non-nan values.'.format(np.sum(~np.isnan(sumstats[key]))))

    df.to_csv(args.out,
        index=False, header=True, sep='\t', na_rep=args.na_rep,
        compression='gzip' if args.gzip else None)
    log.log('Result is written into {}'.format(args.out))

### =================================================================================
###                          Implementation for ldsc_to_mat
### =================================================================================
def ldsc_to_mat(args, log):
    check_input_file(args.ref)
    if args.sumstats: check_input_file(args.sumstats)
    for chri in range(1, 23):
        if args.ldscore: check_input_file(args.ldscore.replace('@', str(chri)))
        if args.annot: check_input_file(args.annot.replace('@', str(chri)))
        if args.M: check_input_file(args.M.replace('@', str(chri)))
        if args.M_5_50: check_input_file(args.M_5_50.replace('@', str(chri)))
    check_output_file(args.out, args.force)

    log.log('Reading reference file {}...'.format(args.ref))
    ref_file = pd.read_table(args.ref, sep='\t', usecols=[cols.SNP, cols.A1, cols.A2])
    log.log("Reference dict contains {d} snps.".format(d=len(ref_file)))

    save_dict = {}
    if args.sumstats:
        df_sumstats = pd.read_csv(args.sumstats, sep='\t')
        df_len = len(df_sumstats)
        log.log('Read {} SNPs from --sumstats file'.format(df_len))

        df_sumstats = pd.merge(ref_file, df_sumstats, how='left', on=cols.SNP)
        df_sumstats = df_sumstats.dropna(how='any')
        log.log('Removed {} SNPs not in --ref file'.format(df_len - len(df_sumstats)))
        df_len = len(df_sumstats)

        df_sumstats['ALLELES'] = df_sumstats.A1_x + df_sumstats.A2_x + df_sumstats.A1_y + df_sumstats.A2_y
        df_sumstats = df_sumstats[filter_alleles(df_sumstats['ALLELES'])].copy()
        log.log('Removed {} SNPs with alleles not matching --ref file'.format(df_len - len(df_sumstats)))
        log.log('{} SNPs remain'.format(len(df_sumstats)))

        df_z = df_sumstats.Z.copy()
        df_sumstats.Z = align_alleles(df_sumstats.Z, df_sumstats['ALLELES'])
        log.log('Flip z score for {} SNPs'.format((df_sumstats.Z != df_z).sum()))

        df_sumstats.drop(['A1_x', 'A2_x', 'A1_y', 'A2_y','ALLELES'], axis=1, inplace=True)
        df_sumstats = pd.merge(ref_file, df_sumstats, how='left', on=cols.SNP)
        save_dict['zvec'] = df_sumstats["Z"].values
        save_dict['nvec'] = df_sumstats["N"].values

    if args.ldscore:
        df_ldscore = pd.concat([pd.read_csv(args.ldscore.replace('@', str(chri)), delim_whitespace=True)  for chri in range(1, 23)])
        df_ldscore.drop([x for x in ['CHR', 'BP', 'CM', 'MAF'] if x in df_ldscore], inplace=True, axis=1)
        log.log('Shape of ldscore file: {shape}'.format(shape=df_ldscore.shape))
        df_ldscore = pd.merge(ref_file[['SNP']], df_ldscore, how='left', on='SNP')
        del df_ldscore['SNP']
        log.log('Shape of ldscore file after merge: {shape}'.format(shape=df_ldscore.shape))
        save_dict['annonames'] = list(df_ldscore.columns)
        save_dict['annomat'] = df_ldscore.values

    if args.annot:
        df_annot = pd.concat([pd.read_csv(args.annot.replace('@', str(chri)), delim_whitespace=True)  for chri in range(1, 23)])
        df_annot.drop([x for x in ['CHR', 'BP', 'CM', 'MAF'] if x in df_annot], inplace=True, axis=1)
        log.log('Shape of annots file: {shape}'.format(shape=df_annot.shape))
        df_annot = pd.merge(ref_file[['SNP']], df_annot, how='left', on='SNP')
        del df_annot['SNP']
        log.log('Shape of annots file after merge: {shape}'.format(shape=df_annot.shape))
        save_dict['annonames_bin'] = list(df_annot.columns)
        save_dict['annomat_bin'] = df_annot.values

    if args.M_5_50:
        m_5_50 = pd.concat([pd.read_csv(args.M_5_50.replace('@', str(chri)), delim_whitespace=True, header=None) for chri in range(1, 23)])
        m_5_50 = np.atleast_2d(m_5_50.sum().values)
        log.log('M_5_50={}'.format(m_5_50))
        save_dict['M_5_50']=m_5_50

    if args.M:
        m = pd.concat([pd.read_csv(args.M.replace('@', str(chri)), delim_whitespace=True, header=None) for chri in range(1, 23)])
        m = np.atleast_2d(m.sum().values)
        log.log('M={}'.format(m))
        save_dict['M']=m

    sio.savemat(args.out, save_dict, format='5', do_compression=False, oned_as='column', appendmat=False)
    log.log('Result written to {f}'.format(f=args.out))

### =================================================================================
###                          Implementation for diff_mat
### =================================================================================
def diff_mat(args, log):
    check_input_file(args.mat1)
    check_input_file(args.mat2)
    check_input_file(args.ref)
    if args.sumstats: check_input_file(args.sumstats)
    check_output_file(args.out, args.force)

    def read_mat_file(filename, log):
        log.log('Reading .mat file {}...'.format(filename))
        mat = sio.loadmat(filename)
        log.log('Found variables: {}'.format([x for x in mat.keys() if x not in ['__version__', '__header__', '__globals__']]))
        nvec = None; zvec = None; logpvec = None
        for key in mat.keys():
            if key.lower().startswith('zvec'): zvec = mat[key]
            if key.lower().startswith('logpvec'): logpvec = mat[key]
            if key.lower().startswith('nvec'): nvec = mat[key]
        return (zvec, logpvec, nvec)

    def compare_vectors(v1, v2,  log):
        s1 = v1 is not None
        s2 = v2 is not None
        log.log('\tIs present in the first  file? {}'.format('YES' if s1 else 'NO'))
        log.log('\tIs present in the second file? {}'.format('YES' if s2 else 'NO'))
        if not s1 or not s2: return (None, None)

        s1 = len(v1); s2 = len(v2)
        log.log('\tHave equal length in both files? {}'.format('YES, {}'.format(s1) if s1==s2 else 'NO, {} vs {}'.format(s1, s2)))
        if s1 != s2: return (None, None)

        s = all(np.isnan(v1) == np.isnan(v2))
        log.log('\tHave equal pattern of defined values? {}'.format('YES' if s else 'NO'))

        if not s:
            s1 = np.sum(np.isnan(v1)); s2 = np.sum(np.isnan(v2))
            log.log('\t\tHave equal number of undefined values? {}'.format('YES, {}'.format(s1) if s1==s2 else 'NO, {} vs {} - difference is {}'.format(s1, s2, abs(s1-s2))))

            s1 = np.sum(~np.isnan(v1) & np.isnan(v2))
            s2 = np.sum(np.isnan(v1) & ~np.isnan(v2))
            log.log('\t\tHow many values are defined in the first file but not in the second file? {}'.format(s1))
            log.log('\t\tHow many values are defined in the second file but not in the first file? {}'.format(s2))

        # Select values defined in both vectors
        idx_def = ~np.isnan(v1) & ~np.isnan(v2)
        def_v1 = v1[idx_def]
        def_v2 = v2[idx_def]
        idx_diff_sign = ((def_v1 <= 0) & (def_v2 > 0)) | ((def_v1 > 0) & (def_v2 <= 0))
        idx_diff_value = np.greater(abs(def_v1 - def_v2), 1e-5)
        log.log('\tAll defined values are equal in both files? {}'.format('YES' if all(def_v1 == def_v2) else 'NO, {} values are different'.format(np.sum(def_v1 != def_v2))))
        if not all(def_v1 == def_v2):
            s = np.sum(idx_diff_sign)
            log.log('\t\tMean(Std) for the first  file? {} ({})'.format(np.mean(def_v1), np.std(def_v1)))
            log.log('\t\tMean(Std) for the second file? {} ({})'.format(np.mean(def_v2), np.std(def_v2)))
            log.log('\t\tAll values have equal sign? {}'.format('YES' if s == 0 else 'NO, {} values have different sign'.format(s)))
            log.log('\t\tMaximum difference of absolute values? {}'.format(max(abs(abs(def_v1) - abs(def_v2)))))
            log.log('\t\tNumber of values with absolute difference above 1e-5? {}'.format(sum(idx_diff_value)))

        # Return index of SNPs that we consider to be different
        # First vector, idx_diff_nan_or_sign, indicate where nan pattern is different or sign is different
        # Second vector, idx_diff_nan_or_sign_or_value, indicates where nan pattern is different or sign is different or value is different
        idx_diff_nan_or_sign = (np.isnan(v1) != np.isnan(v2))
        idx_diff_nan_or_sign[idx_def] = idx_diff_sign
        idx_diff_nan_or_sign_or_value = (np.isnan(v1) != np.isnan(v2))
        idx_diff_nan_or_sign_or_value[idx_def] = idx_diff_value | idx_diff_sign
        return (idx_diff_nan_or_sign, idx_diff_nan_or_sign_or_value)

    [zvec1, logpvec1, nvec1] = read_mat_file(args.mat1, log)
    [zvec2, logpvec2, nvec2] = read_mat_file(args.mat2, log)

    log.log('{}:'.format('zvec'))
    (zvec_diff, _) = compare_vectors(zvec1, zvec2, log)

    log.log('{}:'.format('logpvec'))
    (_, logpvec_diff) = compare_vectors(logpvec1, logpvec2, log)

    log.log('{}:'.format('nvec'))
    (_, nvec_diff) = compare_vectors(nvec1, nvec2, log)

    diff = np.logical_or.reduce([diff for diff in [zvec_diff, logpvec_diff, nvec_diff] if diff is not None])
    log.log('Overall, {} markers appears to be different.'.format(np.sum(diff)))

    if np.sum(diff) == 0:
        open(args.out, 'a').close()
        return

    log.log('Reading reference file {}...'.format(args.ref))
    ref = pd.read_table(args.ref, sep='\t', usecols=[cols.SNP, cols.CHR, cols.BP])
    log.log("Reference dict contains {d} snps.".format(d=len(ref)))

    # Insert not-null vectors into the data frame
    vectors = {'zvec1':zvec1, 'zvec2':zvec2, 'logpvec1':logpvec1, 'logpvec2':logpvec2, 'nvec1':nvec1, 'nvec2':nvec2}
    for k in vectors:
        if (vectors[k] is not None) and (len(vectors[k]) == len(ref)):
            ref[k] = vectors[k]
    ref = ref.loc[[i for i, x in enumerate(diff) if x]]

    if args.sumstats:
        log.log('Reading sumstats file {}...'.format(args.sumstats))
        sumstats = pd.read_table(args.sumstats, sep='\t')
        log.log("Sumstats file contains {d} markers.".format(d=len(sumstats)))

        sumstats_cols = sumstats.columns
        if cols.SNP in sumstats_cols:
            sumstats.columns = [(x + '_onSNP' if x != 'SNP' else x) for x in sumstats_cols]
            ref = pd.merge(ref, sumstats, how='left', on='SNP')
        if cols.CHR in sumstats_cols and cols.BP in sumstats_cols:
            sumstats.columns = [(x + '_onCHRPOS' if x not in ['CHR', 'BP'] else x) for x in sumstats_cols]
            ref = pd.merge(ref, sumstats, how='left', on=['CHR', 'BP'])

    ref.to_csv(args.out, sep='\t', index=False, na_rep='NA')
    log.log('Result is written into {}'.format(args.out))

### =================================================================================
###                                Misc stuff and helpers
### =================================================================================
def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh):
        self.fh = fh
        self.log_fh = open(fh + '.log', 'w')

        # remove error file from previous run if it exists
        try:
            os.remove(fh + '.error')
        except OSError:
            pass

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        print(msg)
        self.log_fh.write(str(msg).rstrip() + '\n')

    def error(self, msg):
        '''
        Print to log file, error file and stdout with a single command.
        '''
        print(msg)
        self.log_fh.write(str(msg).rstrip() + '\n')
        with open(self.fh + '.error', 'w') as error_fh:
            error_fh.write(str(msg).rstrip() + '\n')

def wait_for_system_memory(log):
    # safety guard to ensure that at least 25% of system memory is available
    try:
        import psutil
        if psutil.virtual_memory().percent < 80:
            return
        log.log('Waiting for at least 20% of system memory to be available...')
        while psutil.virtual_memory().percent > 80:
            time.sleep(5)
    except ImportError:
        pass

### =================================================================================
###                                Main section
### ================================================================================= 
if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    if args.out is None:
        raise ValueError('--out is required.')

    log = Logger(args.out)
    start_time = time.time()

    try:
        defaults = vars(parse_args([sys.argv[1]]))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = MASTHEAD
        header += "Call: \n"
        header += './sumstats.py {} \\\n'.format(sys.argv[1])
        options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        wait_for_system_memory(log)
        log.log('Beginning analysis at {T} by {U}, host {H}'.format(T=time.ctime(), U=getpass.getuser(), H=socket.gethostname()))

        # run the analysis
        args.func(args, log)

    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.error( traceback.format_exc(ex) )
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
