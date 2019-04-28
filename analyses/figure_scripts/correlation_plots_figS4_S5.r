library(psych)
library(parcor)
library(fpc)
library(gplots)
library(rprojroot)

# You will need to replace with the path to the data directory on your machine
setwd(find_root("metabochem.Rproj"))
setwd("analyses/data/")
train.dif.s <- read.csv("significance_results/healthstate_anova_wsig_control_training_serum_nonpara.txt")
train.dif.p <- read.csv("significance_results/healthstate_anova_wsig_control_training_plasma_nonpara.txt")

load("training_set_tq_nonormalize.rda")
d$Set <- paste0(rep("Training", nrow(d)), d$Health_State)
d <- d[, c(1:5, 136, 6:135)]
d_serum <- d[d$Organ == "Serum", ]
d_plasma <- d[d$Organ == "Plasma", ]
health_df_serum <- d_serum
health_df_plasma <- d_plasma

test.dif.s <- read.csv("significance_results/healthstate_anova_wsig_control_test_serum_nonpara.txt")
test.dif.p <- read.csv("significance_results/healthstate_anova_wsig_control_test_plasma_nonpara.txt")

load("test_set_tq_nonormalize.rda")
d$Set <- paste0(rep("Test", nrow(d)), d$Health_State)
d <- d[, c(1:5, 136, 6:135)]
d_serum <- d[d$Organ == "Serum", ]
d_plasma <- d[d$Organ == "Plasma", ]
health_df_serum.test <- d_serum
health_df_plasma.test <- d_plasma

# cluster based on spearmen and partial correlations-------------

single.met.list <- list(health_df_serum, health_df_serum.test,
                        health_df_plasma, health_df_plasma.test)

head(single.met.list[[4]])

# regress intensity values on the covariates
for(i in seq_along(single.met.list)) {
  for(j in 7:ncol(single.met.list[[i]])) {
    lm.fit <- lm(single.met.list[[i]][, j] ~ Smoking_Status + Gender, data = single.met.list[[i]])
    single.met.list[[i]][, j] <- lm.fit$residuals
  }
}

names(single.met.list) <- c("serum_training", "serum_test",
                            "plasma_training", "plasma_test")

dif.list <- list(train.dif.s, train.dif.s, train.dif.p, train.dif.p)

hclust.average <- function(x) hclust(x, method="average")
palette <- colorRampPalette(c("blue","yellow"))(n=300)

for (i in 1:3) {
  x <- single.met.list[[i]][, -seq(1, 6)]
  partial.cor <- pls.net(as.matrix(x), ncomp = 30)
  partial.cor <- partial.cor$pcor

  partial.cor <- partial.cor[c(which(dif.list[[i]]$adj_pvalues < .075)),
                             c(which(dif.list[[i]]$adj_pvalues < .075))]

  colnames(partial.cor) <- colnames(single.met.list[[i]])[which(dif.list[[i]]$adj_pvalues < .075) + 6]
  rownames(partial.cor) <- colnames(single.met.list[[i]])[which(dif.list[[i]]$adj_pvalues < .075) + 6]

  d <- as.dist(1 - abs(partial.cor))
  ks <- seq(2, nrow(as.matrix(d)) - 1, by = 1)
  hc <- hclust(d, method = "average")
  ASW <- sapply(ks, FUN=function(k) {
    fpc::cluster.stats(d, cutree(hc,k))$avg.silwidth
  })
  nClus <- ks[which.max(ASW)]
  tiff(paste0("",
             names(single.met.list)[i], "par_asw.tiff"), res = 300,
       width = 7, height = 7, units = "in")
  plot(y=ASW, x=ks, type="l", ylab="Average Silhouette Width", xlab="Number of Clusters",
       xaxt = "n")
  axis(side = 1, at = ks)
  dev.off()

  plot(hclust(d, method ="average"), main="Average Linkage with Correlation-Based Distance", xlab="", sub ="")

  tiff(paste0("",
             names(single.met.list)[i], "par_cor.tiff"),
       res = 300, width = 7, height = 7, units = "in")
  heatmap.2(data.matrix(d), trace="none", density.info="none", col=palette, mar=c(10,10),cexCol=1, cexRow=1,
            Rowv = T, Colv=T, keysize = 1,hclustfun = hclust.average)
  dev.off()

  d <- as.dist(1 - abs(cor(single.met.list[[i]][, c(which(dif.list[[i]]$adj_pvalues < .075) + 6)],
                           method = "spearman")))
  ks <- seq(2, nrow(as.matrix(d)) - 1, by = 1)
  hc <- hclust(d, method = "average")
  ASW <- sapply(ks, FUN=function(k) {
    fpc::cluster.stats(d, cutree(hc,k))$avg.silwidth
  })
  nClus <- ks[which.max(ASW)]
  tiff(paste0("",
             names(single.met.list)[i], "spe_asw.tiff"), res = 300,
       width = 7, height = 7, units = "in")
  plot(y=ASW, x=ks, type="l", ylab="Average Silhouette Width", xlab="Number of Clusters",
       xaxt = "n")
  axis(side = 1, at = ks)
  dev.off()

  tiff(paste0("",
             names(single.met.list)[i], "spe_cor.tiff"),
       res = 300, width = 7, height = 7, units = "in")
  heatmap.2(data.matrix(d), trace="none", density.info="none", col=palette, mar=c(10,10),cexCol=1, cexRow=1,
            Rowv = T, Colv=T, keysize = 1,hclustfun = hclust.average)
  dev.off()
}

