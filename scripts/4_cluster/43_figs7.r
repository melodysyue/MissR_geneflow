library(tidyverse)
library(IRanges)
library(ggsignif)
library(patchwork)

rm(list=ls())

species <- "Emerald Shiner"
aln <- read.table("./data/alignment/ersn.aln.hmm.islands.stats.txt", header=T)
ild <- read.table("./data/hmm_inversions/ersn.aln.hmm.islands.target.list", header=T)

win <- 100000
c <- 9
aln.c <- aln %>% 
  filter(chr==c)
min <- min(aln.c$chr.pos)
max <- max(aln.c$chr.pos)

starts <- seq(from=1, to=max, by=win)
ends <- starts+win-1
ends[ends>max] <- max


result <- NULL
#loop through each window and count the number of SNPs in each window

for (i in 1:length(starts)){
  start <- starts[i]
  end <- ends[i]
  n.snp <- aln.c %>% 
    filter(chr.pos>=start & chr.pos<=end) %>% 
    nrow()
  
  temp <- c(index=i, start=start, end=end,n.snp=n.snp)
  result <- rbind(result,temp) %>% 
    as.data.frame()
}

result <- result %>% 
  mutate(mid=(start+end)/2)

#find which windows overlap with HMM islands
aln.ild <- aln.c %>% 
  filter(win %in% ild$x) %>% 
  group_by(win) %>% 
  summarize(start=min(chr.pos), end=max(chr.pos))


query <- IRanges(result$start, result$end)
subject <- IRanges(aln.ild$start, aln.ild$end)

overlaps <- findOverlaps(query,subject)

ylim <- max(result$n.snp)


p1 <- result %>% 
  ggplot(aes(x=mid,y=n.snp))+
  geom_line()+
  geom_point(data=subset(result, index %in% overlaps@from), 
             aes(x=mid, y=ylim*0.95), col="navyblue", size=3.5, alpha=0.5)+
  theme_bw(base_size = 20)+
  xlab("Position")+
  ylab("Number of SNPs")+
  ggtitle(paste0(species, " Chromosome ", c))+
  scale_x_continuous(breaks=c(0,2e7,4e7),
                     labels=c("0", "20Mb", "40Mb"))

result.type <- result %>% 
  mutate(type = case_when(
    index %in% overlaps@from ~ "Island", 
    TRUE ~ "Background")) %>% 
  filter(n.snp>0)
result.type$type <- factor(result.type$type, levels = c("Island", "Background"))

means <- aggregate(n.snp ~ type, result.type, mean)
means$n.snp <- round(means$n.snp,2)
  
p2 <- result.type %>% 
  ggplot(aes(x=type, y=n.snp))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("Island", "Background")),
              map_signif_level = TRUE)+
  theme_bw(base_size = 20)+
  xlab("")+
  ylab("Number of SNPs")+
  geom_text(data=means, aes(label=n.snp, y=ylim*0.95), size=5)


pdf("FigS7_ersn_chr9.pdf", width = 12, height = 3)
p1+p2+
  plot_layout(widths = c(3,1))
dev.off()