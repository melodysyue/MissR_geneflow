# MissR_geneflow
 
# How to cite:
Yue Shi, Kristen L. Bouska, Garrett J. McKinney, William Dokai, Andrew Bartels, Megan V. McPhee, Wesley A. Larson. Gene flow influences the genomic architecture of local adaptation in six riverine fish species. *Molecular Ecology* (2021). https://www.biorxiv.org/content/10.1101/2021.05.18.444736v1

# Data

These data are in the '/data/' directory, which includes necessary input data for various analyses for the six study species. We use the following code for species names.
  - **bhmw** -> Bullhead Minnow
  - **blgl** -> Bluegill
  - **fwdm** -> Freshwater Drum
  - **cncf** -> Channel Catfish
  - **gzsd** -> Gizzard Shad
  - **ersn** -> Emerald Shiner

'./data/popmap/': popmap filter (after filtering) for each species.

'./data/outlier_neutral/': list of outlier SNPs and neutral SNPs for each species, along with identified outliers with each method. 

'./data/alignment/': list of aligned SNPs for each species.

./data/inversions/': list of SNPs within each identified putative inversions.


The following data were deposited in Dyrad.

vcf files (after filtering) for each species along with corresponding genepop files.


# Scripts

**./0_read_processing/**:
 - 00_process_radtags.slurm
 - 01_clone_filter.slurm
 - 02_ustacks.slurm
 - 03_cstacks.slurm
 - 04_sstacks.slurm
 - 05_stacks_TGP.slurm
 
**./1_snp_filtering/**:
 - 11_hdplot.r
 - 12_vcf_keep_highest_MAF.py
 - 13_countHetsMissing_genepop_sample-ncode3.pl
 


