# MissR_geneflow
 
# How to cite:
Yue Shi, Kristen L. Bouska, Garrett J. McKinney, William Dokai, Andrew Bartels, Megan V. McPhee, Wesley A. Larson. Gene flow influences the genomic architecture of local adaptation in six riverine fish species. *Molecular Ecology* (2021). https://www.biorxiv.org/content/10.1101/2021.05.18.444736v1

# Data

Input or intermediate data for various analyses for the six study species. We use the following code for species names.

  - **bhmw** -> Bullhead Minnow
  - **blgl** -> Bluegill
  - **fwdm** -> Freshwater Drum
  - **cncf** -> Channel Catfish
  - **gzsd** -> Gizzard Shad
  - **ersn** -> Emerald Shiner

`./data/popmap/`: popmap filter (after filtering) for each species.

`./data/outlier_neutral/`: #chrom (tag) and #pos (tag.pos)
  - `\*fstol.list`: *F<sub>ST</sub>* outliers
  - `\*geaol.list`: GEA outliers
  - `\*fst_gea_ol.list`: outlier SNPs (union of *F<sub>ST</sub>* outliers and GEA outliers)
  - `\*neutral.list`: neutral SNPs
  - `\*outlier_summary.txt`: identified outliers using each of the 7 methods. Note: bs: Bayescan; arl: Arlequin; ofk: OutFLANK; bf: Bayenv2.  

`./data/alignment/`: list of aligned SNPs (after filtering) for each species.

**Note**:

vcf files (after filtering) for each species along with corresponding genepop files were deposited in Dyrad.

# Scripts

`./scripts/0_read_processing/`:
 - 00_process_radtags.slurm
 - 01_clone_filter.slurm
 - 02_ustacks.slurm
 - 03_cstacks.slurm
 - 04_sstacks.slurm
 - 05_stacks_TGP.slurm
 
`./scripts/1_snp_filtering/`:
 - 11_hdplot.r
 - 12_vcf_keep_highest_MAF.py
 - 13_countHetsMissing_genepop_sample-ncode3.pl
 


