library(bigsnpr)
#########################################
###Conduct LD thinning on neutral SNPs###
#########################################

rm(list=ls())

#plink <- download_plink2("./plink/")

plink <- "./plink/plink2.exe"
bedfile <- './plink/neutral/ersn.neutral.bed'
popmap <- read.table("/data/popmap/popmap_filtered_ersn.txt",header=T)

###check the data
(obj.bed <- bed(bedfile))

###Get the input ready
# Read from bed/bim/fam, it will create new files.
snp_readBed(bedfile) #only need to do this once
#this creates two new output: X.rds and X.bk
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("./plink/neutral/ersn.neutral.rds") 
#load the "bigSNP" object


# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")

###Get aliases for useful slots
G <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
CHR <- parse_number(CHR) #bigsnpr only takes chrs as integers
POS <- obj.bigSNP$map$physical.pos
# Check and make sure there is no missing data
big_counts(G, ind.col = 1:10)

G2 <- snp_fastImputeSimple(G) #impute missing data
big_counts(G2, ind.col = 1:10)


svd <- snp_autoSVD(G2, CHR, POS, 
                   min.mac=3,
                   max.iter=10,
                   ncores=nb_cores(),
                   roll.size=0) #set roll.size=0 otherwise it will say it is too big. so annoying!


###Output a list of SNPs after removing outliers
table(obj.bigSNP$fam$family.ID == popmap$sample) #make sure order is the same

subset=attr(svd, "subset")

thinned <- obj.bigSNP$map %>% 
  rownames_to_column("index") %>% 
  filter(index %in% subset) %>% 
  select(chr=chromosome, pos=physical.pos)

write.table(thinned, "./data/outlier_neutral/ersn.neutral_thin.list", quote=F, row.names = F)
