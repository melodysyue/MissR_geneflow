library(tidyverse)
library(ggsignif) #geom_signif
library(scales)
library(patchwork)
library(ggpubr)
library(grid)

rm(list=ls())

aln <- read.table("./data/alignment/ersn.bwa.aln.stats.txt",header=T)
dxy <- read.table("./data/dxy/ersn.dxy.txt",header=T)
ol <- read.table("./data/outlier_neutral/ersn.fst_gea_ol.list",header=T)
ol <- ol %>% 
  mutate(snp = paste0("tig_", chrom, ":", pos))

aln.stat <- left_join(aln, dxy, by="locus")
c=9
ld <- read.table(paste0("./data/ld/ersn.chr9.ld"),header=T)
name="Emerald Shiner Chromosome 9"

###filter ild.stat to only the specified chromosome, and order by SNP position
aln.stat <- aln.stat %>% 
  mutate(tag_id=paste0("tig_", tag, ":", tag.pos)) %>% 
  mutate(ol = ifelse(tag %in% ol$chrom, "Outliers", "Background")) %>% 
  filter(chr==c) %>% 
  arrange(chr.pos) %>% 
  rownames_to_column("rowID") %>% 
  select(rowID, locus, tag, tag.pos, tag_id, chr, chr.pos, chr.cumpos, ol, Fstp, Ho, dxy_avg)

aln.stat$ol <- as.factor(aln.stat$ol)
levels(aln.stat$ol)
aln.stat$ol <- factor(aln.stat$ol, levels=c('Outliers', "Background"))

###Get LD information
###format ld data. 
#Get position information
aln <- aln %>% 
  select(.,locus:chr.cumpos, maf) %>% 
  mutate(snp=paste0("tig_",tag,":",tag.pos))

#Fill in positions in the LD table
ld.pos <- left_join(ld, aln, by=c("SNP_A"="snp")) %>% 
  dplyr::select(chr_A=chr, chr.pos_A=chr.pos, SNP_A, CHR_B, BP_B, SNP_B, R2)
ld.pos <- left_join(ld.pos, aln, by=c("SNP_B"="snp")) %>% 
  dplyr::select(chr_A, chr.pos_A, SNP_A, chr_B=chr, chr.pos_B=chr.pos, SNP_B,R2)

#filter loci with maf >=0.01 for LD calculation
maf <- aln %>% 
  filter(maf>=0.01) %>% 
  pull(snp)

ld.pos <- ld.pos %>% 
  mutate(dist_kb=abs(chr.pos_A-chr.pos_B)/1000) %>% 
  filter(SNP_A %in% maf & SNP_B %in%maf)

ld.pos <- ld.pos %>% 
  mutate(pairs = ifelse(SNP_A %in% ol$snp & SNP_B %in% ol$snp, "Outliers", "Background"))
ld.pos$pairs <- as.factor(ld.pos$pairs)
levels(ld.pos$pairs)
ld.pos$pairs <- factor(ld.pos$pairs, levels=c("Outliers", "Background"))

###outliers vs background

p_fstp <- aln.stat %>% 
  ggplot(aes(x=ol, y=Fstp, fill=ol))+
  geom_boxplot(alpha=0.5)+
  ylab("Fstp")+
  scale_fill_manual(values=c("red", "dimgray"))+
  theme_classic(base_size=20)+
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=1),
        legend.position = "bottom",
        legend.title = element_blank())

p_fstp


p_ho <- aln.stat %>% 
  ggplot(aes(x=ol, y=Ho, fill=ol))+
  geom_boxplot(alpha=0.5)+
  ylab("Ho")+
  scale_fill_manual(values=c("red", "dimgray"))+
  theme_classic(base_size=20)+
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=1))

p_ho

p_dxy <-aln.stat %>% 
  ggplot(aes(x=ol, y=dxy_avg, fill=ol))+
  geom_boxplot(alpha=0.5)+
  ylab("Dxy")+
  scale_fill_manual(values=c("red", "dimgray"))+
  theme_classic(base_size=20)+
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=1))

p_dxy

ld_box <- ld.pos %>% 
  ggplot(aes(x=pairs, y=R2, fill=pairs))+
  geom_boxplot(alpha=0.5)+
  ylab(expression(r^2))+
  scale_fill_manual(values=c("red", "dimgray"))+
  theme_classic(base_size=20)+
  theme(axis.title.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour="black",fill=NA, size=1))

ld_box

#legend
leg <- get_legend(p_fstp)
pdf("./Fig6_legend.pdf", width = 12, height = 1)
as_ggplot(leg)
dev.off()

pdf("./Fig6_ersn_chr9.pdf", width = 12, height = 4)
p <- p_fstp + p_ho + p_dxy + ld_box +
  plot_layout(nrow=1, guides="collect") & theme(legend.position = "none")

p + plot_annotation(
  title=name,
  theme=theme(plot.title = element_text(size=20))
)
dev.off()


