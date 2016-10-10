# python_convert
A python tool to load in raw summary statistics

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

Conversion on 2.5M SNPs template `/home/oleksandr/2558411_ref.bim`

```
python sumStats2ref.py -F "CHARGE_general_cognitive_function_summary_results" -Ref 2558411_ref.bim -T COG_CHARGE -O COG_CHARGE --snpCol MarkerName --pCol P --effACol Effect_Allele --othACol Non_eff_Allele --effCol Beta 
```

Conversion on 9M SNPs template `/home/oleksandr/9279485_ref.bim`
```
python sumStats2ref.py -F "CHARGE_general_cognitive_function_summary_results" -Ref 9279485_ref.bim -T COG_CHARGE -O COG_CHARGE_9m --snpCol MarkerName --pCol P --effACol Effect_Allele --othACol Non_eff_Allele --effCol Beta 
```
