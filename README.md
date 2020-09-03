# ABCD_IQ_Polygenic_Score
Notebook containing the analysis used for '[Polygenic Score of Intelligence is More Predictive of Crystallized than Fluid Performance Among Children](https://www.biorxiv.org/content/10.1101/637512v2)' R. Loughnan et. al.

The notebook assumes that ABCD data has been downloaded from [NIMH](https://nda.nih.gov/abcd), genetics have been pre-processed as described in the manscript. This should result in four files:
1. An `.Rds` file downloaded from [NIMH](https://nda.nih.gov/abcd) (`nda2.0.1.Rds`)
2. A file contining the IQ and Educational Attainment polygenic scores created using [PRSice](https://www.prsice.info) (`merged_prs_scores.tsv`)
3. A genetic ancestry file created using [fastStructure](https://rajanil.github.io/fastStructure/) (`Genetic_Ancestry_Factors.4.txt`)
4. A genetic PC file created using [PLINK](https://www.cog-genomics.org/plink/) (`plink2.eigenvec`)

