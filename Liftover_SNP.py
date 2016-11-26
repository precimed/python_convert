import os, re
import numpy as np
import pandas as pd
import argparse
import logging
import random
from distutils.version import StrictVersion

from pyliftover import LiftOver

Intro = r'''
Lifting SNPs 'rs' number and genomic position across different builds.

Lift 'rs' number:
----------------
    1. Based on NCBI SNP merge and history files: RsMergeArch and SNPHistory
    2. SNPs with ID coded as 'CHR:POS' will be left untouched
    Details: http://genome.sph.umich.edu/wiki/LiftOver

Lift genomic position:
----------------------
    1. Only performed when the CHR and POS columns are specified
    2. Based on UCSC build converting files:
        hg17ToHg19.over.chain.gz
        hg18ToHg19.over.chain.gz
    3. It is possible that one position in old build maps to multiple
    positions in new build. In such case,
        The SNP is removed from cleaned data
        The matching with highest score were used in dup_*.txt
    4. It is also possible that no new position matches the old
        The SNP is removed from the cleaned data
        The old position was kept in miss_*.txt

The results after running the script will be stored in the following files, where * is the name of the original file:
    - lifted_* is the main result. It contains the original input plus new columns with the lifted data (SNP, CHR, POS)

    - dup_* report duplicated rs and/or CHR:BP entries after lifting
    - miss_* report CHR:BP entries that couldn't be lifted by LiftOver tool with given chain file
    - multi_* report CHR:BP entries that lift into multiple locations with given chain file

    Reports:
    - lift_pos_result reports the results of lifting BP
    - lift_rs_result reports the results of lifting RS# numbers
    - summary_lift_pos reports summary of lifting BP
    - summary_lift_rs reports summary of lifting RS# numbers
'''

def myopen(fn):
    import gzip
    try:
        h = gzip.open(fn)
        ln = h.read(2) # read arbitrary bytes so check if @param fn is a gzipped file
    except:
        # cannot read in gzip format
        return open(fn)
    h.close()
    return gzip.open(fn)

def read_rs_history(histFile):
    RS_HISTORY = set() # store rs

    logging.info("Reading '{}' file...".format(histFile))
    for ln in myopen(histFile):
        fd = ln.strip().split('\t')
        if ln.find('re-activ') < 0:
            RS_HISTORY.add(fd[0])

    return RS_HISTORY


def read_rs_merge(mergFile):
    RS_MERGE = dict() # high_rs -> (lower_rs, current_rs)

    logging.info("Reading '{}' file...".format(mergFile))
    for ln in myopen(mergFile):
        fd = ln.strip().split('\t')
        h, l = fd[0], fd[1]
        c = fd[6]
        RS_MERGE[h] = (l, c)

    return RS_MERGE


# Returns list of tuples, where each tuple is (rs#, build, chr, pos)
def fetch_snps(snp_ids, verbose=False):
    from Bio import Entrez
    Entrez.email = "oleksandr.frei@gmail.com"

    def pull_var(v, line):
        return [x for x in line if x.startswith(v)][0].replace(v, '')

    def parse_snp(snp_info):

        snp = snp_info.split('\n')

        rsId = snp[0].split(" | ")[0]
        lineset = [x.split(' | ') for x in snp if x.startswith('CTG')]
        if len(lineset) == 0:
            return None

        try:
            build = pull_var("assembly=", lineset[0])
            chr = pull_var("chr=", lineset[0])
            pos = pull_var("chr-pos=", lineset[0])
        except:
            return None

        return rsId, build, chr, pos

    logging.info('Querying dbSNP for {} SNPs...'.format(len(snp_ids)))
    response = Entrez.efetch(db='SNP', id=','.join(snp_ids), rettype='flt', retmode='flt').read()
    logging.info('Done')
    if verbose:
        print response

    snp_infos = []
    for snp_info in filter(None, response.split('\n\n')):
        snp_infos.append(parse_snp(snp_info))
    return snp_infos

def lift_rs(rsvec, RS_HISTORY, RS_MERGE):
    RS_LIFTED = rsvec.copy(); nsnps = len(rsvec)
    RS_idx = np.empty((nsnps,), dtype='|S10')
    logging.info("Lifting rs# numbers for n={} SNPs...".format(nsnps))
    for i in xrange(nsnps):
        rs = rsvec[i]
        if (i+1) % 200000 == 0:
            logging.info("{} SNPs done".format(i+1))
        if rs not in RS_MERGE:
            RS_LIFTED[i] = rs; RS_idx[i] = 'unchanged'
            continue
        while True:
            if rs in RS_MERGE:
                rsLow, rsCurrent = RS_MERGE[rs]
                if rsCurrent not in RS_HISTORY and rsCurrent != '':
                    RS_LIFTED[i] = rsCurrent; RS_idx[i] = 'lifted'
                    break
                else:
                    rs = rsLow
            else:
                RS_LIFTED[i] = rs; RS_idx[i] = 'unlifted'
                break
    logging.info("Lifting rs# numbers is finished.")
    return RS_LIFTED, RS_idx
        
def lift_pos(posvec, chrvec, chainFile):
    logging.info("Lifting genomic positions...")
    nsnps = len(posvec)
    posvec = posvec -1;
    pos_lifted = np.empty((nsnps,), dtype='int32')
    chr_lifted = np.empty((nsnps,), dtype='int32')
    pos_indi = np.empty((nsnps,), dtype='|S10')
    dup_indi = np.empty((nsnps,), dtype='bool'); dup_indi.fill(False)
    lift = LiftOver(chainFile)
    for i in xrange(nsnps):
        if (i+1) % 200000 == 0:
            logging.info("{} SNPs done".format(i+1))
        pos = posvec[i]; chr = 'chr%d' % (chrvec[i],)
        tmp = lift.convert_coordinate(chr, pos)
        if not tmp:
            pos_lifted[i] = pos; pos_indi[i] = 'miss'; chr_lifted[i]=chrvec[i]
        elif len(tmp) > 1:
            pos_lifted[i] = tmp[0][1]; 
            chr_lifted[i] = re.sub('chr', '', tmp[0][0])
            pos_indi[i] = 'multi'
        else:
            pos_lifted[i] = tmp[0][1] 
            chr_lifted[i] = re.sub('chr', '', tmp[0][0])
            if pos == tmp[0][1]:
                pos_indi[i] = 'unchanged'
            else:
                pos_indi[i] = 'lifted'
    return pos_lifted+1, pos_indi, chr_lifted
                
def trim_ch_rs (sum_dat, snpCol, chrCol, with_pos):
    nsnps = sum_dat.shape[0]
    logging.info("Parsing {} rs# numbers from the input file...".format(nsnps))
    chrnum_vec = np.empty((nsnps, ), dtype='int')
    rsvec_num = []; rsPattern = re.compile(r'rs[0-9]*')
    rsidx = np.empty((nsnps, ), dtype='bool'); rsidx.fill(False)
    for i in xrange(nsnps):
        if (i+1) % 200000 == 0:
            logging.info("{} SNPs done".format(i+1))
        rs = sum_dat.loc[:,snpCol][i]
        if with_pos:
            chr = sum_dat.loc[:,chrCol][i]
        if rsPattern.match(rs):
            rsidx[i] = True
            rsvec_num.append(re.sub('rs', '', rs))
        if with_pos:
            chrnum_vec[i] = int(re.sub('[chrCHR]', '', str(chr)))
    rsvec_num = np.array(rsvec_num)
    return rsvec_num, rsidx, chrnum_vec

def try_find_build(rs, pos):

    snps_info = fetch_snps(rs)
    #snps_info = [('rs3737728', 'GRCh38.p2', '1', '1086035'), ('rs3934834', 'GRCh38.p2', '1', '1070426'), ('rs9651273', 'GRCh38.p2', '1', '1096160')]

    logging.info("Loading liftover chain files...")
    lift38_19 = LiftOver('dataset/hg38ToHg19.over.chain.gz')
    lift19_18 = LiftOver('dataset/hg19ToHg18.over.chain.gz')
    lift19_17 = LiftOver('dataset/hg19ToHg17.over.chain.gz')
    logging.info("Done")

    for (rsId, build, true_chr, pos_hg38), source_pos in zip(snps_info, pos):
        try:
            if build != 'GRCh38.p2':  # assume a specific build we get from Entrez.efetch(db='SNP')
                continue
            source_pos -= 1
            pos_hg19 = lift38_19.convert_coordinate('chr{}'.format(true_chr), int(pos_hg38) - 1)[0][1]
            pos_hg18 = lift19_18.convert_coordinate('chr{}'.format(true_chr), pos_hg19)[0][1]
            pos_hg17 = lift19_17.convert_coordinate('chr{}'.format(true_chr), pos_hg19)[0][1]
            print "build={} {} chr{} source={} hg38={}{} hg19={}{} hg18={}{} hg17={}{}".format(
                build, rsId, true_chr, source_pos,
                pos_hg38, '*' if pos_hg38==source_pos else '',
                pos_hg19, '*' if pos_hg19==source_pos else '',
                pos_hg18, '*' if pos_hg18==source_pos else '',
                pos_hg17, '*' if pos_hg17==source_pos else '')
        except:
            pass

def lift_over(sumFile, outDir, histFile, mergFile, chainFile,
              snpCol, chrCol, posCol, bim=False, reffile="", find_build=False):
    logging.info("Reading input file '{}'...".format(sumFile))
    sum_dat = pd.read_csv(sumFile, sep=' +|\t', engine='python')
    logging.info("Done. Columns are: {}".format(", ".join(sum_dat.columns)))
    if bim:
        logging.info("Setting new column names based on BIM format")
        sum_dat.columns = ['CHR', 'SNP', 'GP', 'POS', 'A1', 'A2']
        snpCol='SNP'; chrCol='CHR'; posCol='POS'

    if snpCol not in sum_dat.columns:
        raise ValueError("Input file does not have {} column".format(snpCol))

    with_pos = chrCol is not None and chrCol != '-' and posCol is not None and posCol != '-'
    with_ref = reffile != None and reffile != ""

    if with_pos:
        if chrCol not in sum_dat.columns:
            raise ValueError("Input file does not have {} column".format(chrCol))
        if posCol not in sum_dat.columns:
            raise ValueError("Input file does not have {} column".format(posCol))

    if find_build:
        sample_size = 60
        sample = random.sample(xrange(sum_dat.shape[0]), sample_size)
        sum_dat_sample = sum_dat.ix[sample, :]
        sum_dat_sample = sum_dat_sample.sort_values(chrCol)
        sum_dat_sample.reset_index(inplace=True)
        try_find_build(sum_dat_sample[snpCol].as_matrix(), sum_dat_sample[posCol].as_matrix())
        return

    logging.info("Checking if there are duplicates by rs# number in the input file... ")
    duplicated = sum_dat.duplicated(snpCol, keep=False)
    if any(duplicated):
        logging.warning("{} duplicated rs# numbers were found in the input file".format(sum(duplicated)))
    else:
        logging.info("No duplicated rs# numbers were found in the input file")

    rsvec_num, rsidx, chrnum_vec = trim_ch_rs(sum_dat, snpCol, chrCol, with_pos)

    RS_HISTORY = read_rs_history(histFile) if isinstance(histFile, basestring) else histFile
    RS_MERGE = read_rs_merge(mergFile) if isinstance(mergFile, basestring) else mergFile

    lifted_rs, lift_rs_indi = lift_rs(rsvec_num, RS_HISTORY, RS_MERGE)
    summary_lift_rs(sum_dat.loc[:, snpCol][rsidx], lifted_rs, lift_rs_indi, 
            outDir)
    sum_dat.loc[:,'new_ID'] = sum_dat.loc[:, snpCol].copy()
    sum_dat.loc[rsidx, 'new_ID'] = np.array(['rs%s' % (s,) for s in lifted_rs])
    # TO-DO: It is better to flag all duplicate but the function with
    #   keep = False doesnt work ! So the first of multiple occurance will
    #   sneak into the clean dataset
    # ofrei: keep = False seems to work well for me with pandas 0.18.0
    sum_dat_dup_idx = sum_dat.duplicated(subset = ('new_ID',), keep=False)
    if with_pos:
        lifted_pos, lift_pos_indi ,chr_lifted = lift_pos(sum_dat[posCol],
                chrnum_vec, chainFile)
        summary_lift_pos(sum_dat.loc[:,snpCol], sum_dat.loc[:,chrCol], 
                sum_dat.loc[:, posCol], lifted_pos, lift_pos_indi, outDir)
        sum_dat.loc[:,'new_pos'] = lifted_pos
        sum_dat.loc[:,'new_chr'] = chr_lifted.astype('int')
        sum_dat.loc[:, 'postag'] = np.array(['%s:%s' % (str(c), str(p)) for
            c, p in zip(chr_lifted, lifted_pos)])
        sum_dat_dup_idx2 = sum_dat.duplicated(subset = ('postag'), keep=False)
        sum_dat_dup_idx = np.logical_or(sum_dat_dup_idx, sum_dat_dup_idx2)
        sum_dat_miss_idx = lift_pos_indi=='miss' 
        sum_dat_multi_idx = lift_pos_indi=='multi'
        if np.sum(sum_dat_multi_idx) > 0:
            sum_dat_multi = sum_dat.ix[sum_dat_multi_idx, :]
            multi_file = os.path.join(outDir, 'multi_%s' % (os.path.basename(sumFile),))
            sum_dat_multi.to_csv(multi_file, index=False, sep='\t')
            logging.info("Created {} file with {} entries".format(multi_file, sum_dat_multi.shape[0]))

    elif with_ref:
        logging.info("Lifting with BIM reference file...")
        refbim = pd.read_csv(reffile,delimiter='\t', header=None)
        refbim.columns = ['CHR', 'SNP', 'GP', 'POS', 'A1', 'A2']
        if any(duplicated):
            logging.warning("(!!!) pandas.merge has not been tested on how it merges duplicated entries (!!!!)")
        tmp = pd.merge(sum_dat, refbim, left_on='new_ID', right_on='SNP',
                how='left')
        sum_dat_miss_idx = np.isnan(tmp.POS)
        sum_dat.loc[~sum_dat_miss_idx, 'new_pos'] = tmp.POS[
                ~sum_dat_miss_idx].astype('int')
        sum_dat.loc[~sum_dat_miss_idx, chrCol] = tmp.CHR[
                ~sum_dat_miss_idx].astype('int')
        sum_dat_multi_idx = np.empty((sum_dat.shape[0],), dtype='bool')
        sum_dat_multi_idx.fill(False)
    else:
        sum_dat_miss_idx = np.empty((sum_dat.shape[0],), dtype='bool')
        sum_dat_miss_idx.fill(False)
        sum_dat_multi_idx = np.empty((sum_dat.shape[0],), dtype='bool')
        sum_dat_multi_idx.fill(False)
    if np.sum(sum_dat_dup_idx) > 0:
        sum_dat_dup = sum_dat.ix[sum_dat_dup_idx, :]
        dup_file = os.path.join(outDir, 'dup_%s' % (os.path.basename(sumFile),))
        sum_dat_dup.to_csv(dup_file, index=False, sep='\t')
        logging.info("Created {} file with {} entries".format(dup_file, sum_dat_dup.shape[0]))
    if np.sum(sum_dat_miss_idx) > 0:
        sum_dat_miss = sum_dat.ix[sum_dat_miss_idx, :]
        miss_file = os.path.join(outDir, 'miss_%s' % (os.path.basename(sumFile),))
        sum_dat_miss.to_csv(miss_file, index=False, sep='\t')
        logging.info("Created {} file with {} entries".format(miss_file, sum_dat_miss.shape[0]))
    rm_idx = (sum_dat_dup_idx.astype('int') + 
                sum_dat_miss_idx.astype('int') + 
                sum_dat_multi_idx.astype('int')) >=1
    if np.sum(rm_idx) > 0:
        sum_dat = sum_dat.ix[~rm_idx,:]
        logging.warning(
            "{0} entries removed from the input because of duplication, "
            "misses or multiple mappings between builds".format(np.sum(rm_idx)))
    if bim:
        logging.info("Updating plink files...")
        update_plinkfiles(outDir, sumFile, sum_dat, snpCol)
    else:
        result_file = os.path.join(outDir, 'lifted_%s' % (os.path.basename(sumFile),))
        logging.info("Saving the result to {}...".format(result_file))
        sum_dat.to_csv(result_file, index=False, sep='\t')
        logging.info("Done.")

def update_plinkfiles(outDir, sumFile, sum_dat, snpCol):
    tmp_ex = os.path.join(outDir, 'tmp_extract.txt')
    sum_dat.to_csv(tmp_ex, index=False, sep='\t', columns=(snpCol,))
    bf = re.sub('.bim','', sumFile)
    plink_cmd = r'''plink --bfile %s --extract %s --make-bed \
            --out %s''' % (bf, tmp_ex, sumFile+'.tmp')
    os.system(plink_cmd)
    tmp_bim = pd.read_csv(sumFile+'.tmp.bim',delimiter='\t', header=None)
    tmp_bim.columns = ['oCHR', 'oSNP', 'oGP', 'oPOS', 'oA1', 'oA2']
    tmp = pd.merge(tmp_bim, sum_dat, left_on='oSNP', right_on=snpCol,
            how='left')
    miss_idx = np.isnan(tmp.loc[:,'new_pos'])
    assert np.sum(miss_idx) == 0
    tmp_pos_file = os.path.join(outDir, 'tmp_update_pos.txt')
    tmp.to_csv(tmp_pos_file, index=False, sep='\t', 
            columns=('oSNP', 'new_pos'), header=None)
    plink_cmd = r'''plink --bfile %s --update-map %s 2 --make-bed \
                        --out %s''' % (sumFile+'.tmp',tmp_pos_file,
                                sumFile+'.tmp2')
    os.system(plink_cmd)
    tmp_bim2 = pd.read_csv(sumFile+'.tmp2.bim',sep='\t', header=None)
    tmp_bim2.columns = ['oCHR', 'oSNP', 'oGP', 'oPOS', 'oA1', 'oA2']
    tmp = pd.merge(tmp_bim2, sum_dat, left_on='oSNP', right_on=snpCol,
            how='left')
    tmp.to_csv(sumFile+'.tmp2.bim', sep='\t', header=None, index=False,
            columns=('oCHR', 'new_ID', 'oGP', 'oPOS', 'oA1', 'oA2'))
    os.system("mv %s %s " % (sumFile+'.tmp2.bim', bf+'_lifted.bim'))
    os.system("mv %s %s " % (sumFile+'.tmp2.bed', bf+'_lifted.bed'))
    os.system("mv %s %s " % (sumFile+'.tmp2.fam', bf+'_lifted.fam'))
    os.system("rm %s " % (sumFile+'.tmp.bim',))
    os.system("rm %s " % (sumFile+'.tmp.fam',))
    os.system("rm %s " % (sumFile+'.tmp.bed',))
    os.system("rm %s " % (tmp_pos_file,))
    os.system("rm %s " % (tmp_ex,))

def summary_lift_rs(orig_rs, new_rs, indivec, outDir):
    li_idx = indivec == 'lifted'
    ul_idx = indivec == 'unlifted'
    uc_idx = indivec == 'unchanged'
    orig_rs = np.array(orig_rs)
    summary_file = os.path.join(outDir, 'summary_lift_rs.txt')
    logging.info("Saving lift summary to '{}'...".format(summary_file))
    with open(summary_file, 'w') as f:
        f.write('Total number of SNPs with "rs" number: %d\n' % (
            len(orig_rs,)))
        f.write('\t Total number of SNPs with "rs" lifted: %d\n' % (
            np.sum(li_idx,)))
        f.write('\t Total number of SNPs with "rs" unchanged: %d\n' % (
            np.sum(uc_idx,)))
        f.write('\t Total number of SNPs with "rs" cant lift: %d\n' % (
            np.sum(ul_idx,)))

    results_file = os.path.join(outDir, 'lift_rs_result.txt')
    logging.info("Saving lifted SNPs to '{}'...".format(summary_file))
    with open(results_file, 'w') as f:
        f.write('ORI_RS\tNEW_RS\tSTATUS\n')
        for i in xrange(len(orig_rs)):
            f.write('%s\trs%s\t%s\n' % (orig_rs[i], new_rs[i], indivec[i]))

def summary_lift_pos(orig_snp, chrvec, posvec, new_posvec, indivec, outDir):
    li_idx = indivec == 'lifted'
    uc_idx = indivec == 'unchanged'
    miss_idx = indivec == 'miss'
    multi_idx = indivec == 'multi'
    orig_snp = np.array(orig_snp); chrvec = np.array(chrvec)
    posvec = np.array(posvec); new_posvec = np.array(new_posvec)
    with open(os.path.join(outDir, 'summary_lift_pos.txt'), 'w') as f:
        f.write('Total number of SNPs: %d\n' % ( len(orig_snp,)))
        f.write('\t Total number of SNPs lifted: %d\n' % (np.sum(li_idx,)))
        f.write('\t Total number of SNPs unchanged: %d\n' % (np.sum(uc_idx,)))
        f.write('\t Total number of SNPs missed: %d\n' % (np.sum(miss_idx,)))
        f.write('\t Total number of SNPs with multiple locations in new build: %d\n' % ( np.sum(multi_idx,)))
    with open (os.path.join(outDir, 'lift_pos_result.txt'), 'w') as f:
        f.write('SNP\tCHR\t\ORI_POS\tNEW_POS\tSTATUS\n')
        for i in xrange(len(orig_snp)):
            f.write('%s\t%s\t%d\t%d\%s\n' % (orig_snp[i], str(chrvec[i]), 
                posvec[i], new_posvec[i], indivec[i]))

if __name__ == "__main__":
    import time
    import warnings
    warnings.simplefilter('ignore')
    tsts = time.time()
    parser = argparse.ArgumentParser(prog="Liftover_SNPs",
            formatter_class=argparse.RawTextHelpFormatter,
            description=Intro)

    parser.add_argument('input_file', type=str, help='Path of the input SNPs file')

    parser.add_argument('-s', '--snp', type=str, required=True, help='The name of the SNP field in the input file', dest='snp_column')
    parser.add_argument('-c', '--chr', type=str, help='The name of the Chromosome field name in the input file', default='-', dest='chr_column')
    parser.add_argument('-p', '--pos', type=str, help='The name of the BP field in input file', default='-', dest='pos_column')

    parser.add_argument(      '--output-folder', type=str, default='.', help='Output directory')
    parser.add_argument(      '--history-file', type=str, default='dataset/SNPHistory.bcp.gz',help='NCBI SNP build history file')
    parser.add_argument(      '--merge-file', type=str, default='dataset/RsMergeArch.bcp.gz', help='NCBI SNP merge file')
    parser.add_argument(      '--chain-file', type=str, default='dataset/hg18ToHg19.over.chain.gz', help='UCSC chain file')
    parser.add_argument(      '--find-build', action='store_true', help='Attempt to detect the build of the input file', default=False)

    parser.add_argument(       '--bim', action='store_true', help='(experimental option) update PLINT fileset bim file', default=False)
    parser.add_argument(       '--ref', type=str, help='(experimental option) Reference bim file', default='')

    parser.add_argument('-v',  '--verbose', action="store_true", help="increase output verbosity")

    args = parser.parse_args()

    logging_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging_level)

    if StrictVersion(pd.__version__) < StrictVersion("0.17.0"):
        logging.warning("Old pandas version detected. "
                        "Liftover_SNP script hasn not been tested with pandas={}", pd.__version__)

    # Check correctness of user-provided arguments
    if (args.chr_column == '-') != (args.pos_column == '-'):
        raise ValueError("Arguments --chr-column and --pos-column must be provided together")

    if args.chr_column == '-' and args.find_build:
        raise ValueError("Unable to find build without CHR:POS information")

    if not os.access(args.output_folder, os.F_OK):
        logging.warning("Output directory {} not exists, making one for you".format(args.output_folder))
        os.makedirs(args.output_folder)

    lift_over(args.input_file, args.output_folder, args.history_file, args.merge_file, args.chain_file,
              args.snp_column, args.chr_column, args.pos_column, args.bim, args.ref, args.find_build)

    logging.info('Finish at {}'.format(time.ctime()))
    ted = time.time()
    logging.info('Time taken {} mins {} sec'.format((ted-tsts)//60, np.round(ted-tsts) % 60))
