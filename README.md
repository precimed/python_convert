# A collection of various utilities for GWAS summary statistics.

## sumstats.py

sumstats.py is a collection of utilities that work with GWAS summary stats.
``csv`` utility reads raw summary statistics files
and convert them into a standardized format:
tab-separated file with standard
column names, standard chromosome labels,
NA label for missing data, etc.
``qc`` utility perform a set of highly customizable quality control procedures.
``mat`` utility re-saves summary stats in MATLAB format for cond/conj pleiotropy analysis.
``lift`` utility can lift genomic corredinats across genomic builds, and SNP rs numbers to a newer versions of SNPdb.

Some of the steps require additional data. Examples can be found [here](http://norment.myftp.org:8080/python_convert/).

```
usage: sumstats.py [-h]
                   {csv,qc,mat,lift,clump,rs,ls,mat-to-csv,ldsc-to-mat,frq-to-mat,ref-to-mat,ldsum,diff-mat} ...

A collection of various utilities for GWAS summary statistics.

positional arguments:
  {csv,qc,zscore,mat,lift,clump,rs,ls,mat-to-csv,ldsc-to-mat,frq-to-mat,ref-to-mat,ldsum,diff-mat,neff}
    csv                 Load raw summary statistics file and convert it into a
                        standardized format: tab-separated file with standard
                        column names, standard chromosome labels, NA label for
                        missing data, etc. The conversion does not change the
                        number of lines in the input files (e.g. no filtering
                        is done on markers). Unrecognized columns are removed
                        from the summary statistics file. The remaining
                        utilities in sumstats.py work with summary statistics
                        files in the standardized format.
    qc                  Miscellaneous quality control and filtering procedures
    zscore              Calculate z-score from p-value column and effect size
                        column
    mat                 Create mat files that can be used as an input for
                        cond/conj FDR and for CM3 model. Takes csv files
                        (created with the csv task of this script). Require
                        columns: SNP, P, and one of the signed summary
                        statistics columns (BETA, OR, Z, LOGODDS). Creates
                        corresponding mat files which can be used as an input
                        for the conditional fdr model. Only SNPs from the
                        reference file are considered. Zscores of strand
                        ambiguous SNPs are set to NA. To use CHR:POS for
                        merging summary statistics with reference file
                        consider 'rs' utility which auguments summary
                        statistics with SNP column (first run 'sumstats.py rs
                        ...', then feed the resulting file into sumstats.py
                        mat ...)
    lift                Lift RS numbers to a newer version of SNPdb, and/or
                        liftover chr:pos to another genomic build using UCSC
                        chain files. WARNING: this utility may use excessive
                        amount of memory (up and beyong 32 GB of RAM).
    clump               Perform LD-based clumping of summary stats. This works
                        similar to FUMA snp2gene functionality
                        (http://fuma.ctglab.nl/tutorial#snp2gene). Step 1. Re-
                        save summary stats into one file for each chromosome.
                        Step 2a Use 'plink --clump' to find independent
                        significant SNPs (default r2=0.6) Step 2b Use 'plink
                        --clump' to find lead SNPs, by clumping independent
                        significant SNPs (default r2=0.1) Step 3. Use 'plink
                        --ld' to find genomic loci around each independent
                        significant SNP (default r2=0.6) Step 4. Merge
                        together genomic loci which are closer than certain
                        threshold (250 KB) Step 5. Merge together genomic loci
                        that fall into exclusion regions, such as MHC Step 6.
                        Output genomic loci report, indicating lead SNPs for
                        each loci Step 7. Output candidate SNP report
    rs                  Augument summary statistic file with SNP RS number
                        from reference file. Merging is done on chromosome and
                        position. If SNP column already exists in --sumstats
                        file, it will be overwritten.
    ls                  Report information about standard sumstat files,
                        including the set of columns available, number of
                        SNPs, etc.
    mat-to-csv          Convert matlab .mat file with logpvec, zvec and
                        (optionally) nvec into CSV files.
    ldsc-to-mat         Convert .sumstats, .ldscore, .M, .M_5_50 and binary
                        .annot files from LD score regression to .mat files.
    frq-to-mat          Convert .frq files plink from .mat files.
    ref-to-mat          Convert reference files to .mat files.
    ldsum               convert plink .ld.gz files (pairwise ld r2) to ld
                        scores
    diff-mat            Compare two .mat files with logpvec, zvec and nvec,
                        and report the differences.
    neff                generate N column from NCASE and NCONTROL, as 4 / (1 /
                        NCASE + 1 / NCONTROL)  

optional arguments:
  -h, --help            show this help message and exit
```

For more information about each command call ``sumstats.py <command> --help``.

Examples:
```
python $(python_convert)/sumstats.py csv --sumstats scz2.snp.results.txt.gz --out PGC_SCZ_2014.csv --force --auto --head 5 --chr hg19chrc
python $(python_convert)/sumstats.py mat --sumstats PGC_SCZ_2014.csv --out PGC_SCZ_2014.mat --ref 2558411_ref.bim --force
```

Further examples can be found in [GWAS_SUMSTAT/Makefile](https://github.com/precimed/GWAS_SUMSTAT/blob/master/Makefile).

## sumstats.py clump

``clump`` utility determine 
  - independent significant SNPs
  - lead SNPs
  - genomic loci
  - candidate SNPs
using the same logic as FUMA's snp2gene. An example:

```
python sumstats.py clump \
	--clump-field FDR \
	--force  \
	--plink /home/oleksandr/plink/plink \
	--sumstats cond0p01_BIP_vs_COG/result.mat.csv \
	--bfile-chr /full/path/to/ref_1kG_phase3_EUR/chr@ \
	--exclude-ranges ['6:25119106-33854733', '8:7200000-12500000'] \
	--clump-p1 0.01 \
	--out cond0p01_BIP_vs_COG/result.clump
```

Here the input file ``results.mat.csv`` was converted from cond/conj FDR results using this script:

```
# Usage:
# python fdrmat2csv.py result.mat /space/syn03/1/data/GWAS/SUMSTAT/misc/9279485_ref.bim

import pandas as pd
import scipy.io as sio
import sys
import numpy as np
if __name__ == '__main__':
    sumstats = sio.loadmat(sys.argv[1])
    ref=pd.read_table(ys.argv[2], delim_whitespace=True)
    ref['FDR']=sumstats['fdrmat']
    ref[['CHR', 'SNP', 'BP', 'A1', 'A2', 'FDR']].to_csv(sys.argv[1] + '.csv', index=False, sep='\t')
```

## make_ld_matrix

Make LD matrix from reference data. Output either in matlab format or as dense lower triangular text file.
To run this tool you need to download reference data from http://ctg.cncr.nl/software/magma (for example g1000_eur).
Example:
```
python make_ld_matrix.py --ref 2558411_ref.bim --bfile g1000_eur --ld_window_r2 0.1 --savemat ldmat_p1.mat
```
For more info see [make_ld_matrix](./make_ld_matrix/README.md).
