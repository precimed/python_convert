import pandas as pd
import numpy as np
from itertools import permutations
from scipy import stats
import scipy.io as sio
import os
import sys
import argparse
from summary_stats_Utils import *

def get_str_list_sign(str_list):
    return np.array([-1 if e[0]=='-' else 1 for e in str_list], dtype=np.int)

_base_complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
def _complement(seq):
    return "".join([_base_complement[b] for b in seq])

ALLELES = "AGCT"
COMPLEMENT = {''.join(s): _complement(s) for s in permutations(ALLELES, 2)}
AMBIGUOUS = [aa for aa, cc in COMPLEMENT.items() if aa == cc[::-1]]

def parse_args(args):
    parser = argparse.ArgumentParser(description="Convert summary statistics "
        "file into mat file through the intermediate csv file.")
    subparsers = parser.add_subparsers()
    parser_csv = subparsers.add_parser("csv",
        help="Create csv file with three columns: snpid, pvalue, zscores.")
    parser_csv.add_argument("sumstats_file", type=str,
        help="Tab-separated input file.")
    parser_csv.add_argument("ref_file", type=str,
        help="Tab-separated file with list of referense SNPs.")
    parser_csv.add_argument("output_file", type=str, help="Output csv file.")
    parser_csv.add_argument("--id", default=None, type=str,
        help="SNP id column.")
    parser_csv.add_argument("--pval", default=None, type=str,
        help="P-value column.")
    parser_csv.add_argument("--effectA", default=None, type=str,
        help="Effect allele column.")
    parser_csv.add_argument("--otherA", default=None, type=str,
        help="Other allele column.")
    parser_csv.add_argument("--effect", default=None, type=str,
        help="Effect column.")
    parser_csv.add_argument("--signed-effect", action="store_true",
        default=False, help="Effect is signed.")
    parser_csv.add_argument("--ref-id", default="SNP", type=str,
        help="Id column in reference file.")
    parser_csv.add_argument("--ref-a1", default="A1", type=str,
        help="First allele column in reference file.")
    parser_csv.add_argument("--ref-a2", default="A2", type=str,
        help="Second allele column in reference file.")
    parser_csv.add_argument("--chunksize", default=100000, type=int,
        help="Size of chunck to read the file.")
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
    parser_mat.add_argument("--ref-id", default="SNP", type=str,
        help="Id column in reference file.")
    parser_mat.set_defaults(func=make_mat)

    parser_mat = subparsers.add_parser("test",
        help="Run test.")
    parser_mat.set_defaults(func=run_test)

    return parser.parse_args(args)


def make_csv(args):
    """ Based on file with summary statistics creates a tab-delimited csv file
    with 3 columns: snpid, pvalue, zscore.
    Only SNPs from the rederence file are considered. Zscores of ambiguous SNPs
    are set to NA.
    """
    if os.path.isfile(args.output_file):
        raise ValueError("Output file already exists!")

    cname_translation = find_column_name_translation(args.sumstats_file,
        snp=args.id, p=args.pval, a1=args.effectA, a2=args.otherA,
        odds_ratio=args.effect if not args.signed_effect else None,
        beta=args.effect if args.signed_effect else None)
    cname_description = {x: describe_cname[cname_translation[x]] for x in cname_translation if cname_translation[x] != 'UNKNOWN'}
    print('Interpreting column names from summary statistics file as follows:')
    print('\n'.join([x + ':\t' + cname_description[x] for x in cname_description]) + '\n')
    cname_skip = [x for x in cname_translation if cname_translation[x ] == 'UNKNOWN']
    if cname_skip: print('Skip the remaining columns ({}).'.format(', '.join(cname_skip)))
    cname = {v: k for k, v in iteritems(cname_translation)}  # inverse mapping (ignore case when key has multiple values)
    if not args.id: args.id = cname.get('SNP')
    if not args.pval: args.pval = cname.get('P')
    if not args.effectA: args.effectA = cname.get('A1')
    if not args.otherA: args.otherA = cname.get('A2')
    if not args.effect:
        args.effect = next(iter(filter(None, [cname.get('BETA'), cname.get('LOG_ODDS'), cname.get('Z'), cname.get('OR')])), None)
        args.signed_effect = False if 'OR' in cname else True

    print('Reading reference file {}...'.format(args.ref_file))
    usecols = [args.ref_id, args.ref_a1, args.ref_a2]
    reader = pd.read_table(args.ref_file, sep='\t', usecols=usecols,
        chunksize=args.chunksize)
    ref_dict = {}
    for chunk in reader:
        gtypes = (chunk[args.ref_a1] + chunk[args.ref_a2]).apply(str.upper)
        #TODO?: add check whether some id is already in ref_dict
        ref_dict.update(dict(zip(chunk[args.ref_id], gtypes)))
    ref_dict = {i: (aa, COMPLEMENT[aa], aa[::-1], COMPLEMENT[aa[::-1]])
            for i, aa in ref_dict.items()}

    print("Reference dict contains %d snps." % len(ref_dict))

    print('Reading summary statistics file {}...'.format(args.sumstats_file))
    out_columns = ["snpid", "pvalue", "zscore"]
    col_attr = ["id", "pval", "effectA", "otherA", "effect"]
    usecols = [getattr(args, c) for c in col_attr]
    # if signed_effect is true, take effect column as string to handle correctly
    # case of truncated numbers, e.g.: 0.00 and -0.00 should have different sign
    effect_col_dtype = str if args.signed_effect else np.float
    reader = pd.read_table(args.sumstats_file, sep='\t', usecols=usecols,
        chunksize=args.chunksize, dtype={args.effect:effect_col_dtype})
    n_snps = 0
    with open(args.output_file, 'a') as out_f:
        out_f.write("%s\n" % "\t".join(out_columns))
        for i, chunk in enumerate(reader):
            chunk = chunk.loc[chunk[args.id].isin(ref_dict),:]
            gtypes = (chunk[args.effectA] + chunk[args.otherA]).apply(str.upper)
            # index of SNPs that have the same alleles as indicated in reference
            ind = [gt in ref_dict[sid]
                for sid, gt in zip(chunk[args.id], gtypes)]
            chunk = chunk.loc[ind,:]
            gtypes = gtypes[ind]
            log10pv = -np.log10(chunk[args.pval].values)
            # not_ref_effect = [
            #   1 if effect allele in data == other allele in reference
            #   -1 if effect allele in data == effect allele in reference ]
            # So zscores with positive effects will be positive and zscores with 
            # negative effects will stay negative, since
            # stats.norm.ppf(chunk[args.pval]*0.5) is always negetive (see zvect
            # calculation below).
            not_ref_effect = np.array([-1 if gt in ref_dict[sid][:2] else 1
                for sid, gt in zip(chunk[args.id], gtypes)])
            #TODO: check proportion of positive and negative effects
            if args.signed_effect:
                # effect column has str type
                # -1 if effect starts with '-' else 1
                effect_sign = get_str_list_sign(chunk[args.effect])
            else:
                # effect column has np.float type
                # 1 if effect >=1 else -1
                if (chunk[args.effect] < 0).any():
                    raise ValueError("Effect column contains negative values, "
                        "but --signed-effect option was not used.")
                effect_sign = np.sign(chunk[args.effect].values - 1)
                effect_sign[effect_sign == 0] = 1
            effect_sign *= not_ref_effect
            zvect = stats.norm.ppf(chunk[args.pval].values*0.5)*effect_sign
            ind_ambiguous = [j for j,gt in enumerate(gtypes) if gt in AMBIGUOUS]
            # set zscore of ambiguous SNPs to nan
            zvect[ind_ambiguous] = np.nan
            #TODO: check whether output df contains duplicated rs-ids (warn)
            df = pd.DataFrame({"pvalue": log10pv, "zscore":zvect},
                index=chunk[args.id])
            df.to_csv(out_f, index=True, header=False, sep='\t', na_rep='NA')
            n_snps += len(df)
            print("{n} lines processed".format(n=(i+1)*args.chunksize))
        print("{n} SNPs saved to {f}".format(n=n_snps, f=args.output_file))


def make_mat(args):
    """ Takes csv files (created with the csv task of this script) with three
    columns: snpid, pvalue, zscore. Creates corresponding mat files which can be
    used as an input for the conditional fdr model.
    """
    if args.mat_files is None:
        args.mat_files = []
        for csv_f in args.csv_files:
            mat_f = os.path.splitext(csv_f)[0] + ".mat"
            args.mat_files.append(mat_f)
    if args.traits is None:
        args.traits = []
        for csv_f in args.csv_files:
            trait = os.path.splitext(os.path.basename(csv_f))[0]
            args.traits.append(trait)
    if (len(args.csv_files) != len(args.mat_files) or
            len(args.csv_files) != len(args.traits)):
        raise ValueError(("Number of input csv files should be equal to the "
            "number of output mat files and number of traits!"))
    for csv_f, mat_f in zip(args.csv_files, args.mat_files):
        if not os.path.isfile(csv_f):
            raise ValueError("%s input file not found!" % csv_f)
        if os.path.isfile(mat_f):
            raise ValueError("%s output file already exists!" % mat_f)

    ref_snps = pd.read_table(args.ref_file, sep='\t', usecols=[args.ref_id],
        squeeze=True)
    #TODO?: add check whether ref_snps contains duplicates
    print("Reference file contains %d snps" % len(ref_snps))

    for csv_f, mat_f, trait in zip(args.csv_files, args.mat_files, args.traits):
        df_out = pd.DataFrame(columns=["pvalue", "zscore"], index=ref_snps)
        #TODO: check whether input df contains duplicated rs-ids (error)
        df = pd.read_csv(csv_f, sep='\t', index_col=0)
        df_out["pvalue"] = df["pvalue"]
        df_out["zscore"] = df["zscore"]
        save_dict = {"logpvec_"+trait: df_out["pvalue"].values,
            "zvec_"+trait: df_out["zscore"].values}
        # df_out.to_csv(mat_f+'.tmp', index=True, header=False, sep='\t')
        sio.savemat(mat_f, save_dict, format='5', do_compression=False,
            oned_as='column', appendmat=False)
        print("%s created" % mat_f)


def run_test(args):
    raise NotImplementedError("Testing is not yet implemented.")

    out_csv = "test/output/test.csv"
    out_mat = "test/output/test.mat"
    for filename in [out_csv, out_mat]:
        try:
            os.remove(filename)
            print("Removing %s" % filename)
        except OSError:
            pass

    print("Running csv task.")
    csv_args = ["csv", "test/input/pgc.adhd.full.2012-10_tabs_reduced.txt",
        "test/input/2558411_ref.bim", out_csv, "--effect", "zscore",
        "--signed-effect", "--ref-id", "SNP", "--ref-a1", "A1",
        "--ref-a2", "A2"]
    args = parse_args(csv_args)
    args.func(args)

    print("Running mat task.")
    mat_args = ["mat", "test/input/2558411_ref.bim", out_csv, "--ref-id", "SNP",
        "--traits", "adhd", "--mat-files", out_mat]
    args = parse_args(mat_args)
    args.func(args)

    import scipy.io as sio
    ref_mat = "test/input/adhd.mat"
    print("Checking consistency of results with %s" % ref_mat)
    ref = sio.loadmat(ref_mat)
    test = sio.loadmat(out_mat)
    ref_p, ref_z = ref["logpvec_adhd"], ref["zvec_adhd"]
    test_p, test_z = test["logpvec_adhd"], test["zvec_adhd"]
    valid_ref_p, valid_ref_z = np.isfinite(ref_p), np.isfinite(ref_z)
    valid_test_p, valid_test_z = np.isfinite(test_p), np.isfinite(test_z)
    assert (valid_ref_p==valid_test_p).all(), "NAN pattern of pvals differes!"
    assert (valid_ref_z==valid_test_z).all(), "NAN pattern of zscores differes!"

    thresh_all = 1.e-2
    thresh_each = 1.e-3
    ref_p, ref_z = ref_p[valid_ref_p], ref_z[valid_ref_z]
    test_p, test_z = test_p[valid_test_p], test_z[valid_test_z]
    diff_p, diff_z = np.abs(ref_p - test_p), np.abs(ref_z - test_z)
    assert np.linalg.norm(diff_p,1)<thresh_all, "Pvals differ significantly!"
    assert max(diff_p)<thresh_each, "Some pvals differ significantly!"
    assert np.linalg.norm(diff_z,1)<thresh_all, "Zscores differ significantly!"
    assert max(diff_z)<thresh_each, "Some zscores differ significantly!"
    print("Testing successful!")


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    args.func(args)
    print("Done")

    # Examples
    
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
