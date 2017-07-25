# A collection of various utilities for GWAS summary statistics.

## sumstats.py

sumstats.py utility allows you to read raw summary statistics files
and convert them into a standardized format:
tab-separated file with standard
column names, standard chromosome labels,
NA label for missing data, etc.
This file can be further 
loaded into matlab format,
processed by various quality control procedures,
aligned to a set of reference markers,
lifted across genomic builds or versions of SNPdb.

```
usage: sumstats.py [-h]
                   {csv,qc,mat,lift,rs,ls,mat-to-csv,ldsc-to-mat,diff-mat} ...

A collection of various utilities for GWAS summary statistics.

positional arguments:
  {csv,qc,mat,lift,rs,ls,mat-to-csv,ldsc-to-mat,diff-mat}
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
    rs                  Augument summary statistic file with SNP RS number
                        from reference file. Merging is done on chromosome and
                        position.
    ls                  Report information about standard sumstat files,
                        including the set of columns available, number of
                        SNPs, etc.
    mat-to-csv          Convert matlab .mat file with logpvec, zvec and
                        (optionally) nvec into CSV files.
    ldsc-to-mat         Convert .sumstats, .ldscore, .M, .M_5_50 and binary
                        .annot files from LD score regression to .mat files.
    diff-mat            Compare two .mat files with logpvec, zvec and nvec,
                        and report the differences.

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

## make_ld_matrix

Make LD matrix from reference data. Output either in matlab format or as dense lower triangular text file.
To run this tool you need to download reference data from http://ctg.cncr.nl/software/magma (for example g1000_eur).
Example:
```
python make_ld_matrix.py --ref 2558411_ref.bim --bfile g1000_eur --ld_window_r2 0.1 --savemat ldmat_p1.mat
```
For more info see [make_ld_matrix](./make_ld_matrix/README.md).
