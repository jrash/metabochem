## ----setup, include=FALSE------------------------------------------------
library(car)
library(chemmodlab)
library(magrittr)
library(dplyr)
library(glmnet)
library(caret)
library(doParallel)
library(pROC)
library(gplots)
library(ChemmineR)
library(rprojroot)

rm(list = ls())

##################
# Get new signficant metabolites
##################
setwd(find_root("metabochem.Rproj"))
setwd("analyses/data")

# Normalization and log transformation

sdfset <- read.SDFset("common_test_training_molecule-v2.sdf")

# read output from website and label activity

# test for significance with only compounds that have chemical structures
f1 <- read.delim(file = "sample_factors_training.txt", header = T, row.names = 1)
f2 <- read.csv(file = "sample_metabolites_training_excol_fix.csv", header = T, row.names = 1)

# removing the observations with no factor levels
f1 <- na.omit(f1)
rownames(f1) <- f1$Sample.name

# create patient identifier column
f1$Subject_name <- gsub(" Plasma| Serum", "", f1$Subject_name, perl = T)

f1 <- f1[, -2]
colnames(f1) <- c("Patient","Organ", "Health_State", "Smoking_Status", "Gender")
head(f1)

f2 <- t(f2)

#---------Get original metabolite names

# fix compound names so that they are the same
# format as the file provided by Melaine
ids <- sdfid(sdfset)
ids.new <- gsub(" |,","_",ids, perl = T)
colnames(f2) <- gsub(" |,","_", colnames(f2), perl = T)

# 130 compounds provided
sum(ids.new %in% colnames(f2))
dim(f2)

f2 <- f2[, colnames(f2) %in% ids.new]
orig.mets <- colnames(f2)[colnames(f2) %in% ids.new]

# Processing the Training Set
# ------------------------------------------------------------------------
train.dif.s <- 
  read.csv("significance_results/healthstate_anova_wsig_control_training_serum_nonpara.txt")
train.dif.p <- 
  read.csv("significance_results/healthstate_anova_wsig_control_training_plasma_nonpara.txt")

load("training_set_tq_normalize.rda")
health_df <- d
health_df$Patient <- NULL


# ori.serum.mets.05 <- ori.serum.mets.05[-4]
serum.met.idx <- which(train.dif.s$adj_pvalues < .075) + 4
plasma.met.idx <- which(train.dif.p$adj_pvalues < .075) + 4

# make gender binary
health_df$Gender <- ifelse(health_df$Gender == "F", 0, 1)
health_df$Health_State <- ifelse(health_df$Health_State == "Healthy", 0, 1)
health_df$Smoking_Status <- ifelse(health_df$Smoking_Status == "Former", 0, 1)
health_df$Patient <- NULL

#----only analyze the samples collected from Serum
health_df_serum <- health_df[health_df$Organ == "Serum", ]
health_df_serum_full <- health_df_serum 
health_df_serum_full$Organ <- NULL

# only keep significant metabolites
colnames(health_df_serum)[serum.met.idx]
health_df_serum <- health_df_serum[c(1:4, serum.met.idx)]

health_df_serum$Organ <- NULL

#----only analyze the samples collected from Plasma
health_df_plasma <- health_df[health_df$Organ == "Plasma", ]
health_df_plasma_full <- health_df_plasma 
health_df_plasma_full$Organ <- NULL

# only keep significant metabolites
colnames(health_df_plasma)[plasma.met.idx]
health_df_plasma <- health_df_plasma[c(1:4, plasma.met.idx)]

health_df_plasma$Organ <- NULL

health_df$Organ <- ifelse(health_df$Organ == "Serum", 0, 1)

health_df <- health_df[, c(2,1,3:ncol(health_df))]

fing.file <- 
  "metabolics_fingerprint.csv"

# Processing the Test Set
# ------------------------------------------------------------------------

test.dif.s <- read.csv("significance_results/healthstate_anova_wsig_control_test_serum_nonpara.txt")
test.dif.p <- read.csv("significance_results/plasma/healthstate_anova_wsig_control_test_plasma_nonpara.txt")

load("test_set_tq_normalize.rda")

health_df.test <- d
health_df.test$Patient <- NULL

# make gender binary
health_df.test$Gender <- ifelse(health_df.test$Gender == "F", 0, 1)
health_df.test$Health_State <- ifelse(health_df.test$Health_State == "Healthy", 0, 1)
health_df.test$Smoking_Status <- ifelse(health_df.test$Smoking_Status == "Former", 0, 1)
health_df.test$Patient <- NULL

#----only analyze the samples collected from Serum
health_df_serum.test <- health_df.test[health_df.test$Organ == "Serum", ]
health_df_serum_full.test <- health_df_serum.test
health_df_serum_full.test$Organ <- NULL

# only keep significant metabolites
colnames(health_df_serum.test)[serum.met.idx]
health_df_serum.test <- health_df_serum.test[c(1:4, serum.met.idx)]

health_df_serum.test$Organ <- NULL

#----only analyze the samples collected from Plasma
health_df_plasma.test <- health_df.test[health_df.test$Organ == "Plasma", ]
health_df_plasma_full.test <- health_df_plasma.test 
health_df_plasma_full.test$Organ <- NULL

# only keep significant metabolites
colnames(health_df_plasma.test)[plasma.met.idx]
health_df_plasma.test <- health_df_plasma.test[c(1:4, plasma.met.idx)]

health_df_plasma.test$Organ <- NULL

health_df.test$Organ <- ifelse(health_df.test$Organ == "Serum", 0, 1)

health_df.test <- health_df.test[, c(2,1,3:ncol(health_df.test))]

# ----clusters only significant metabolites in plasma or serum
sumClusterIntensitiesPS <- function(f, linkage, clus.df, sample){
  fing <- read.csv(f, row.names = 2)
  mets <- colnames(clus.df)[4:length(colnames(clus.df))]
  fing <- fing[rownames(fing) %in% mets, ] 
  fing$row.ID = NULL
  
  d_fing <- dist(fing, method = "binary")
  ks <- seq(2, nrow(fing) - 1, by = 1)
  hc <- hclust(d_fing, method = linkage)
  ASW <- sapply(ks, FUN=function(k) {
    fpc::cluster.stats(d_fing, cutree(hc,k))$avg.silwidth
  })
  nClus <- ks[which.max(ASW)]
  
  tiff(paste0("",
              sample, "_structure_asw.tiff"), res = 300,
       width = 7, height = 7, units = "in")
  plot(y=ASW, x=ks, type="l", ylab="Average Silhouette Width", xlab="Number of Clusters",
       xaxt = "n")
  axis(side = 1, at = ks)
  dev.off()
  
  kclus <- cutree(hc,nClus)
  
  group_met <- clus.df[-c(1:3)]
  group_met <- t(group_met)
  kclus <- as.data.frame(kclus)
  
  group_met <- merge(kclus, group_met, by.x = "row.names", by.y = "row.names")
  row.names(group_met) <- group_met$Row.names
  group_met$Row.names <- NULL
  
  # you want to sum the raw intensity values for the clustered metabolites
  # not the log transformed values, that would be equivalent to multiplying
  # the raw intensity values, much more sensitive to noise if take the product
  # dont know which is better really, but summing the raw values seems to make
  # the most sense.  This is equvalent to summing the raw values, because
  # each patient metabolite uses the same normalization constant
  group_met[, -1] <- 2^group_met[, -1]
  group_met_des <- group_met %>% 
    group_by(kclus) %>% 
    summarise_each(funs(sum))
  group_met_des[, -1] <- log(group_met_des[, -1], base = 2)
  
  group_met_des <- as.data.frame(group_met_des)
  
  # Summed the *raw* metabolite intensities
  group_met_des <- t(group_met_des[-1])
  
  health.state <- clus.df$Health_State
  
  mets.clus.intens <- clus.df[, c(2, 3)]
  mets.clus.intens <- merge(mets.clus.intens, group_met_des, by.x= "row.names", by.y = "row.names")
  row.names(mets.clus.intens) <- mets.clus.intens$Row.names
  mets.clus.intens$Row.names <- NULL
  
  list(mets.clus.intens, kclus)
}

# serum fingerprint distance heatmaps

fing <- 
  read.csv("metabolics_fingerprint.csv",
           row.names = 2)
mets <- colnames(health_df_serum)[4:length(colnames(health_df_serum))]
fing.serum <- fing[rownames(fing) %in% mets, ] 
fing.serum$row.ID = NULL
d_fing <- dist(fing.serum, method = "binary")

d_fing_mat <- as.matrix(d_fing)
matmat <- data.matrix(d_fing_mat < .3)
matmat <- apply(matmat, 2, as.numeric)
heatmap.2(matmat, trace = "none", cexCol =.75)

heatmap.2(data.matrix(d_fing))
hclust.average <- function(x) hclust(x, method="average")
palette <- colorRampPalette(c("blue","yellow"))(n=300)

tiff("serum_structure_heatmap.tiff", res = 300,
     width = 7, height = 7, units = "in")
heatmap.2(data.matrix(d_fing), trace="none", density.info="none", col=palette, mar=c(10,10),cexCol=1, cexRow=1, 
          Rowv = T, Colv=T, keysize = 1,hclustfun = hclust.average)
dev.off()

mets.clus.serum.train <- sumClusterIntensitiesPS(f = fing.file,
                                                 linkage = "average",
                                                 clus.df = health_df_serum,
                                                 sample = "serum")[[1]]

mets.clus.serum.train <- cbind(health_df_serum$Health_State, mets.clus.serum.train)

mets.clus.serum.test <- sumClusterIntensitiesPS(f = fing.file,
                                                linkage = "average",
                                                clus.df = health_df_serum.test,
                                                sample = "serum")[[1]]
mets.clus.serum.test <- cbind(health_df_serum.test$Health_State, mets.clus.serum.test)
kclus.serum <- sumClusterIntensitiesPS(f = fing.file,
                                       linkage = "average",
                                       clus.df = health_df_serum.test,
                                       sample = "serum")[[2]]

# plasma fingerprint distance heatmaps

fing <- read.csv("metabolics_fingerprint.csv",
                 row.names = 2)
mets <- colnames(health_df_plasma)[4:length(colnames(health_df_plasma))]
fing.serum <- fing[rownames(fing) %in% mets, ] 
fing.serum$row.ID = NULL
d_fing <- dist(fing.serum, method = "binary")

tiff("plasma_structure_heatmap.tiff", res = 300,
     width = 7, height = 7, units = "in")
heatmap.2(data.matrix(d_fing), trace="none", density.info="none", col=palette, mar=c(10,10),cexCol=1, cexRow=1, 
          Rowv = T, Colv=T, keysize = 1,hclustfun = hclust.average)
dev.off()

mets.clus.plasma.train <- sumClusterIntensitiesPS(f = fing.file,
                                                  linkage = "average",
                                                  clus.df = health_df_plasma,
                                                  sample = "plasma")[[1]]

mets.clus.plasma.train <- cbind(health_df_plasma$Health_State, mets.clus.plasma.train)

mets.clus.plasma.test<- sumClusterIntensitiesPS(f = fing.file,
                                                linkage = "average",
                                                clus.df = health_df_plasma.test,
                                                sample = "plasma")[[1]]
mets.clus.plasma.test <- cbind(health_df_plasma.test$Health_State, mets.clus.plasma.test)
kclus.plasma <- sumClusterIntensitiesPS(f = fing.file,
                                        linkage = "average",
                                        clus.df = health_df_plasma.test,
                                        sample = "plasma")[[2]]

# Save the metabolite data for other ML models

save(health_df_serum_full, health_df_plasma_full,
     health_df_serum_full.test, health_df_plasma_full.test, 
     serum.met.idx, plasma.met.idx, 
     kclus.plasma, kclus.serum,
     file = "ML_data.RDATA")