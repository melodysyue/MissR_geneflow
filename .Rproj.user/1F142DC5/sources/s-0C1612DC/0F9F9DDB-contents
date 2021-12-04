library(tidyverse)
library(factoextra)
library(ggpubr)


rm(list=ls())
#############
##env data###
#############
pop.env <- read.csv("./data/env/missR_env_bypop_reduced.csv",header=TRUE)
rownames(pop.env) <- pop.env$X
pop.env <- pop.env[,-1]
rownames(pop.env) <- c("Pool 4", "Pool 8", "Pool 13", "Pool 26", "Open River", "La Grange")
colnames(pop.env) <- c("Chlf MED","DO MED", "pH MED", "Temp MED", "TN MED", "TP MED", "Turb MED", "Ice MED", "Snow MED",
                       "Chlf IQR", "Temp IQR", "TN IQR", "TP IQR", "Turb IQR", "Ice IQR", "Snow IQR",
                       "DRG", "Flow", "Flow 7MIN", "Flow Max")

pca.r <- prcomp(pop.env, scale=TRUE, center=TRUE)


#Figure1C: Env PCA biplot
pdf("Fig1C_env_pca_biplot.pdf", width = 12, height = 9)
fviz_pca_biplot(pca.r,  #This reference line corresponds to the expected value if the contribution where uniform.
                axes=c(1,2),
                col.ind = "dodgerblue",# Individuals color
                col.var = "dimgray", # Variables color
                pointsize=5,
                repel=TRUE,
                labelsize=10)+ 
  labs(title="", x="PC1 (60.5%)", y="PC2 (33.1%)")+
  theme_bw(base_size = 25)

dev.off()

#SuppFigure1: Env PCs
#see how many PCs to retain
ev <- pca.r$sdev^2

#Kaiser-Guttman criterion and broken stick;
evplot = function(ev) {
  # Broken stick model (MacArthur 1957)
  n = length(ev)
  bsm = data.frame(j=seq(1:n), p=0)
  bsm$p[1] = 1/n
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p = 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))
  barplot(ev, main="Kaiser-Guttmann criterion", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="Broken stick model", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

pdf("SuppFig1A.pdf", width = 12, height=6)
evplot(ev)
dev.off()

p2 <- fviz_contrib(pca.r, choice = "var",axes=1)+coord_flip()+
  labs(title = "Contribution of variables to PC1", x="")+
  theme_bw(base_size=15)


p3 <- fviz_contrib(pca.r, choice = "var",axes=2)+coord_flip()+
  labs(title = "Contribution of variables to PC2", x="")+
  theme_bw(base_size=15)

pdf("SuppFig1B.pdf", width = 12, height=9)
ggarrange(p2,p3)
dev.off()


#extract environmental PCs
pca.r$x[,1:2]
plot(pca.r$x[,1])
plot(pca.r$x[,2])
write.table(pca.r$x[,1:2], "missR_env_PC12.txt", quote=F, row.names = T)



