library(vegan)
library(psych) #pairs.panels
library(qvalue)
library(usdm) #VIF
library(tidyverse)
library(RColorBrewer)
library(adegenet)

rm(list=ls())

###################
###genetic input###
###################
#gen <- read.table("./rda_input/bhmw.rda.geno",header=FALSE)
loc <- read.table("./rowID_tag_tagpos/bhmw.chrom.pos.rowID.list",header=FALSE)
popmap <- read.table("/data/popmap/popmap_filtered_bhmw.txt",header=TRUE)

colnames(loc) <- c("id","chrom","pos")
colnames(gen)=loc$id
rownames(gen)=popmap$sample
gen[1:6,1:6]
sum(is.na(gen)) #check missingness
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp)) # No NAs
gen.imp[1:6,1:6]


###############
###env input###
###############
###no need to check multicollinearity 
env.pcs <- read.table("missR_env_PC12.txt", header=TRUE) 
colnames(env.pcs)=c("env_PC1", "env_PC2")
env.pcs$pop <- rownames(env.pcs)

ind.env <- left_join(popmap, env.pcs, by="pop") %>% 
  column_to_rownames("sample") %>% 
  dplyr::select(-pop)


#########
###RDA###
#########

dim(gen.imp)
dim(ind.env)
table(rownames(gen.imp)==rownames(ind.env))

rda.output <- rda(gen.imp ~., data=ind.env, scale=T)
rda.output
summary(rda.output)

###check model significance
anova.cca(rda.output, parallel=4, permutation=100) #significance test for full model
anova.cca(rda.output, parallel=4, by="axis", permutations=100) #significance test for each RDA axis

####################################################
###Identify Candidate SNPs under local adaptation###
####################################################
load.rda <- scores(rda.output, choices=c(1:2), display="species")  # Species scores for the first 2 constrained axes

### loci with extreme loading on the significant RDA axes only.
###this one assume the SNP distribution is normal.

##make sure the SNP loading distribution is normal ish.
hist(load.rda[,1], main="Loadings on RDA1") #this is normal distribution
hist(load.rda[,2], main="Loadings on RDA2")

##create a function to extract loci with extreme loadings
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 
ncand <- length(cand1) + length(cand2)
#ncand <- length(cand1) #if RDA axis2 is not significant

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")
#colnames(cand1) <- c("axis","snp","loading") #if RDA axis 2 is not significant


cand <- rbind(cand1, cand2)
#cand <- cand1 #if RDA axis 2 is not significant

cand$snp <- as.character(cand$snp)

length(cand$snp[duplicated(cand$snp)]) 
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections


###identify the env_PC that have the maximum correlation for each candidate SNP ###
##Note: in some cases, correlations may be strong for multiple variables 
##but here, for the sake of simplicity, we will just use maximum.

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,6] <- names(which.max(abs(bar[c(4,5)]))) # gives the variable
  cand[i,7] <- max(abs(bar[c(4,5)]))              # gives the correlation
}

colnames(cand)[6] <- "predictor"
colnames(cand)[7] <- "correlation"

head(cand)
table(cand$predictor) 
summary(cand)

#write.table(cand, "./rda_prda_output/bhmw.rda.ind.snp.candidates.txt",quote=F, row.names = F)




