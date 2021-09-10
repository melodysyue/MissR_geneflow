rm(list=ls())
library(vcfR)
library(tidyverse)
library(stringr)
library(here)

#reference
#McKinney et al. Molecular Ecology Resources (2016). https://doi.org/10.1111/1755-0998.12613

exampleLoci<-read.vcfR("ersn.noclone.G50.l90.G60.l70.G70.l50.recode.vcf")

#HDplot function
HDPlot<-function(vcfData){
  #set up results table
  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=10))
  colnames(HDplotTable)<-c("Chrom","Pos","depth_a","depth_b","ratio","num_hets","num_samples","het_perc","std","z")
  
  #get genotypes from vcf file
  genos<-extract.gt(vcfData, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  
  #get allele reads from vcf file
  reads<-extract.gt(vcfData, element = "AD", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  
  #replace . with 0
  reads<-gsub("\\.","0,0",reads)
  
  alleleReads<-apply(reads,2,function(x) str_split_fixed(x,",",2))
  alleleReads_1<-alleleReads[1:dim(reads)[1],]
  alleleReads_2<-alleleReads[dim(reads)[1]+1:dim(reads)[1],]
  #convert to numeric format
  alleleReads_1<-apply(alleleReads_1,2, function(x) as.numeric(x))
  alleleReads_2<-apply(alleleReads_2,2, function(x) as.numeric(x))
  #subset to heterozygous genotypes
  #make genotype matrix where heterozygotes are 1 and other genotypes are 0
  hetMatrix<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=0,'./.'=0,'0/1'=1,'1/0'=1))
  #multiply read count matrices by heterozygote matrix to get read counts for heterozygous genotypes
  alleleReads_1_het<-alleleReads_1*hetMatrix
  alleleReads_2_het<-alleleReads_2*hetMatrix
  #rows are loci and columns are samples
  #sum reads per allele per locus for heterozygous samples
  A_reads<-apply(alleleReads_1_het,1,sum)
  B_reads<-apply(alleleReads_2_het,1,sum)
  totalReads<-A_reads+B_reads
  ratio<-A_reads/totalReads
  std<-sqrt(totalReads*0.5*0.5)
  z<- -(totalReads/2-A_reads)/std
  #get percent heterozygosity for each locus
  numHets<-apply(hetMatrix,1,sum)
  hetPerc<-numHets/dim(hetMatrix)[2]
  
  #assign results to HDplotTable
  HDplotTable$Chrom<-vcfData@fix[,1]
  HDplotTable$Pos<-vcfData@fix[,2]
  HDplotTable$Locus_ID
  HDplotTable$depth_a<-A_reads
  HDplotTable$depth_b<-B_reads
  HDplotTable$ratio<-ratio
  HDplotTable$num_hets<-numHets
  HDplotTable$num_samples<-dim(hetMatrix)[2]
  HDplotTable$het_perc<-hetPerc
  HDplotTable$std<-std
  HDplotTable$z<-z
  return(HDplotTable)
}


example_output<-HDPlot(exampleLoci)

#plot results
#deviation plot
hd_deviation <- 
  ggplot()+
  geom_point(data=example_output,aes(x=het_perc,y=z),alpha=0.5)+
  scale_y_continuous(minor_breaks = seq(-100, 100, 5))+
  scale_x_continuous(breaks=seq(0,1,0.1),
                     minor_breaks=seq(0,1,0.05))
#ratio plot
hd_ratio <- 
  ggplot()+geom_point(data=example_output,aes(x=het_perc,y=ratio),alpha=0.5)+xlim(0,1)

write.csv(example_output, file="ersn_hdplot.csv")
ggsave(plot=hd_deviation,"ersn_hd_deviation.pdf", width=12, height=8)
ggsave(plot=hd_ratio,"ersn_hd_ratio.pdf", width=12, height=8)

#Create a blacklist of loci flagged as possible paralogs to remove using vcftools - MAKE SURE TO CHANGE THE het_per AND z PARAMETERS BELOW BASED ON YOUR OUTPUT

#Based on the output plots, we use the following threshold, het 0.55, z (-8 to 8)
zh=8
zl=-8
h=0.55

#plot it
example_output_class <- example_output %>% 
  mutate(class=ifelse(z <= zh & z >= zl & het_perc <= h, "singleton", "non_singleton"))

#deviation plot
hd_deviation <- 
  example_output_class %>% 
  drop_na() %>% 
  ggplot()+
  geom_point(aes(x=het_perc,y=z, color=class), alpha=0.5)+
  xlim(0,1)+
  scale_y_continuous(minor_breaks = seq(-50, 50, 5))+
  scale_x_continuous(breaks=seq(0,1,0.1),
                     minor_breaks=seq(0,1,0.05))
#ratio plot
hd_ratio <- 
  example_output_class %>% 
  drop_na() %>% 
  ggplot()+geom_point(aes(x=het_perc,y=ratio, color=class),alpha=0.5)+xlim(0,1)

write.csv(example_output_class, file="ersn_hdplot.csv")
ggsave(plot=hd_deviation,"ersn_hd_deviation.pdf", width=12, height=8)
ggsave(plot=hd_ratio,"ersn_hd_ratio.pdf", width=12, height=8)


###create a blacklist
high_het=example_output[which(example_output$het_perc>h),]
no_high_het=example_output[which(example_output$het_perc<h),]
no_hi_z=no_high_het[which(no_high_het$z>zh),]
no_low_z=no_high_het[which(no_high_het$z<zl),]
no_na=no_high_het[which(no_high_het$z=="NaN"),]
blacklist=rbind(high_het,no_hi_z,no_low_z, no_na)
write.csv(blacklist, file="ersn_hdplot._filtered.csv")

head(blacklist) #rowname was counted as a column. 
write.table(blacklist[,c("Chrom","Pos")],file="ersn_hdplot_positions_filtered.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

#Double check with the distribution
table(example_output_class$class)
table(is.na(example_output_class$class))
