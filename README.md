# Python tools to work with summary statistics


## sumStats2ref.py 
Python tool to load in raw summary statistics and align them to 2.5M or 9M SNPs template.

```
python sumStats2ref.py --help
usage: Preprocess Summary stats [-h] [-F F] [-Ref REF] [-T T] [-O O]
                                [--forceID] [--snpCol SNPCOL] [--pCol PCOL]
                                [--effACol EFFACOL] [--othACol OTHACOL]
                                [--effCol EFFCOL] [--orCol ORCOL]
                                [--posCol POSCOL] [--chrCol CHRCOL]

Preprocess summary stats for matlab

optional arguments:
  -h, --help         show this help message and exit
  -F F               Summary stats file (default: None)
  -Ref REF           Reference file (default: None)
  -T T               Trait Name (default: None)
  -O O               Output DIR (default: .)
  --forceID          Force using SNP ID other than position (default: False)
  --snpCol SNPCOL    SNP ID field (default: SNP)
  --pCol PCOL        P value field (default: P)
  --effACol EFFACOL  Effective allele field (default: None)
  --othACol OTHACOL  The other allele field (default: None)
  --effCol EFFCOL    Effect size field (default: None)
  --orCol ORCOL      Odds ratio field (default: None)
  --posCol POSCOL    Genomic position field (default: None)
  --chrCol CHRCOL    Chromosome field (default: None)
==============================================================================
```

Conversion on 2.5M SNPs template

```
python sumStats2ref.py -F "CHARGE_general_cognitive_function_summary_results" -Ref 2558411_ref.bim -T COG_CHARGE -O COG_CHARGE --snpCol MarkerName --pCol P --effACol Effect_Allele --othACol Non_eff_Allele --effCol Beta 
```

Conversion on 9M SNPs template
```
python sumStats2ref.py -F "CHARGE_general_cognitive_function_summary_results" -Ref 9279485_ref.bim -T COG_CHARGE -O COG_CHARGE_9m --snpCol MarkerName --pCol P --effACol Effect_Allele --othACol Non_eff_Allele --effCol Beta 
```


## sumstats_convert (by @interCM )

A python tool to convert files with summary statistics into mat files that can
be used as an input for cond/conj FDR as well as for CM3 model.
The process has 2 stages. At the first stage an intermediate csv file is created
which contains 3 columns: snpid, pvalue, zscores. At the second stage this csv
file is used to construct a mat file with p-values and z scores restricted and
ordered according to the provided reference file. Produced mat file can be used
as an input for cond/conj FDR and CM3 models.

Compared to `sumStats2ref.py`, the `sumstats_convert` script has lower memory usage.

### Example of usage

1st stage. Create intermediate csv file:
```
python sumstats_convert.py csv \
../data/ssgac/EduYears_Main.txt.gz \
../projects/adhd/cond_fdr/2m_template/2558411_ref.bim \
../tmp/ssgac_eduyears_2m.csv \
--id MarkerName --effect Beta --pval  Pval --effectA A1 --otherA A2 \
--signed-effect --ref-id SNP --ref-a1 A1 --ref-a2 A2
```

2nd stage. Create mat file using csv file created at the first stage:
```
python sumstats_convert.py mat \
../projects/adhd/cond_fdr/2m_template/2558411_ref.bim \
../tmp/ssgac_eduyears_2m.csv \
--ref-id SNP --traits eduyears_ssgac
```

For more details about available arguments use:
```
python sumstats_convert.py --help
```


## make_ld_matrix.py

Make LD matrix from reference data. Output either in matlab format or as dense lower triangular text file.
To run this tool you need to download reference data from http://ctg.cncr.nl/software/magma (for example g1000_eur).
Example:
```
python make_ld_matrix.py --ref 2558411_ref.bim --bfile g1000_eur --ld_window_r2 0.1 --savemat ldmat_p1.mat
```
For more info run `python make_ld_matrix.py --help`.
```
