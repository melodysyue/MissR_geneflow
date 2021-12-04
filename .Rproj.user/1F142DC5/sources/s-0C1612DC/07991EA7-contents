library(tidyverse)
library(ggrepel)

rm(list=ls())

aln.stat <- read.table("./data/alignment/ersn.aln.hmm.islands.stats.txt",header=T)
aln.stat <- aln.stat %>% 
  mutate(chr.mbp = chr.cumpos/1e6) %>% 
  group_by(win) %>% 
  mutate(win.start = min(chr.mbp),
         win.end = max(chr.mbp))

ol.snps <- read.table("./data/outlier_neutral/ersn.fst_gea_ol.list", header=T) %>% 
  as.data.frame()
ild <- read.table("./data/hmm_inversions/ersn.aln.hmm.islands.target.list", header=T) #comment this out for gizzard shad

#only for emerald shiner
inv6 <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds3.chr6.ol.txt", header=T)
inv9 <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds1.chr9.ol.txt", header=T)
inv19 <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds2.chr19.ol.txt", header=T)
inv <- rbind(inv6, inv9, inv19)
inv <- left_join(inv, aln.stat, by=c("tag", "tag.pos")) %>% 
  select(chr, chr.pos, chr.mbp)
inv.plot <- inv %>% 
  group_by(chr) %>% 
  summarize(inv.min=min(chr.mbp), inv.max=max(chr.mbp)) %>% 
  as.data.frame()


species="Emerald Shiner"
#fst=sprintf("%.4f", round(0.0720,4)) #for bullhead minnow
#fst=sprintf("%.4f", round(0.0303,4)) #for bluegill
#fst=sprintf("%.4f", round(0.0050,4)) #for freshwater drum
#fst=sprintf("%.4f", round(0.0025,4)) #for channel catfish
#fst=sprintf("%.4f", round(0.0024,4)) #for gizzard shad
fst=sprintf("%.4f", round(0.0003,4)) #for emerald shiner



#######################
##labeling parameters##
#######################
###Order chromosomes###
aln.stat$chr <- as.factor(aln.stat$chr)
nCHR=length(levels(aln.stat$chr))
levels(aln.stat$chr)
#comment the following out for gizzard shad
aln.stat$chr <- factor(aln.stat$chr, levels=c(seq(1,nCHR-1,1),"Unplaced")) 
levels(aln.stat$chr)[levels(aln.stat$chr)=="Unplaced"] <- "Un"
levels(aln.stat$chr)


#calculate the difference between the maximum and minimum cumulative base pair 
#position for each chromosome and divide it by two to get the middle of each chromosome.
axis.set <- aln.stat %>% 
  group_by(chr) %>% 
  summarize(center = (max(chr.mbp) + min(chr.mbp)) / 2)

ylim=max(aln.stat$Fstp)+0.1 # have some space on the y axis
##replace nagative Fst value to a very small value close to 0 for plotting reasons
aln.stat <- aln.stat %>% 
  mutate(Fst_plot=if_else(Fstp<0, 0.001, Fstp))

###plot it

p <- ggplot(aln.stat) +
  geom_point(aes(x = chr.mbp, y = Fst_plot, color = chr), size=2) +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits=c(0, ylim)) +
  scale_color_manual(values = rep(c("dimgray","darkgrey"), nCHR)) +
  scale_size_continuous(range = c(1,4)) +
  labs(x = NULL, 
       y = "Fstp",
       title=paste0(species, " (Global Fstp = ", fst, ")")) + 
  theme_minimal(base_size = 15) +
  theme( 
    legend.position = "none",
    panel.border = element_rect(colour="black",fill=NA, size=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 12, vjust = 0.5))+
  
  #add points for genomic islands of high differentiation, comment this chunk for gizzard shad.
  geom_point(data=subset(aln.stat, win %in% ild$x), 
             aes(x=(win.start+win.end)/2, y=ylim*0.9), col="navyblue", size=3.5, alpha=0.5)+
  
  #add points for putative inversions. Only run this for emerald shiner. 
  geom_rect(data=inv.plot,
            aes(xmin=inv.min, xmax=inv.max, ymin=ylim*0.95, ymax=ylim*0.99), fill="purple", alpha=0.8)+
  
  # Add highlighted points for outlier SNPs
  geom_point(data=subset(aln.stat, tag %in% ol.snps$chrom), 
             aes(x = chr.mbp, y = Fst_plot), color="tomato", alpha=0.8, size=3.5)

p
pdf("./Fig2_F.pdf", width = 12, height=3)
p
dev.off()




