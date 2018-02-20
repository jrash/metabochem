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

rm(list = ls())

## ----process, echo = T---------------------------------------------------
train.dif.s <- 
  read.csv("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/significant_metabolites/serum/healthstate_anova_wsig_control_training_serum_nonpara.txt")
train.dif.p <- 
  read.csv("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/significant_metabolites/plasma/healthstate_anova_wsig_control_training_plasma_nonpara.txt")
train.dif.s$significance_adj <- train.dif.s$adj_pvalues < .075
train.dif.p$significance_adj <- train.dif.p$adj_pvalues < .075

load("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/training_v_test/training_set_tq_normalize.rda")
health_df <- d
health_df$Patient <- NULL

# make gender binary
health_df$Gender <- ifelse(health_df$Gender == "F", 0, 1)
health_df$Health_State <- ifelse(health_df$Health_State == "Healthy", 0, 1)
health_df$Smoking_Status <- ifelse(health_df$Smoking_Status == "Former", 0, 1)
health_df$Patient <- NULL

#----only analyze the samples collected from Serum
health_df_serum <- health_df[health_df$Organ == "Serum", ]

# only keep significant metabolites
colnames(health_df_serum)[which(train.dif.s$adj_pvalues < .075) + 4]
health_df_serum <- health_df_serum[c(1:4, which(train.dif.s$adj_pvalues < .075) + 4)]

health_df_serum$Organ <- NULL

#----only analyze the samples collected from Plasma
health_df_plasma <- health_df[health_df$Organ == "Plasma", ]

# only keep significant metabolites
colnames(health_df_plasma)[which(train.dif.p$adj_pvalues < .075) + 4]
health_df_plasma <- health_df_plasma[c(1:4, which(train.dif.p$adj_pvalues < .075) + 4)]

health_df_plasma$Organ <- NULL

health_df$Organ <- ifelse(health_df$Organ == "Serum", 0, 1)

# ------------------------------------------------------------------------
health_df <- health_df[, c(2,1,3:ncol(health_df))]
head(health_df[, 1:10])

# ----cluster-------------------------------------------------------------
sumClusterIntensities <- function(f, linkage, clus.df){
  fing <- read.csv(f, row.names = 2)
  mets <- colnames(clus.df)[5:length(colnames(clus.df))]
  fing <- fing[rownames(fing) %in% mets, ] 
  fing$row.ID = NULL
  
  d_fing <- dist(fing, method = "binary")
  ks <- seq(2, nrow(fing) - 1, by = 10)
  hc <- hclust(d_fing, method = linkage)
  ASW <- sapply(ks, FUN=function(k) {
    fpc::cluster.stats(d_fing, cutree(hc,k))$avg.silwidth
  })
  nClus <- ks[which.max(ASW)]
  plot(y=ASW,x=ks,type="l", ylab="Average Silhouette Width", xlab="Clusters")
  kclus <- cutree(hc,nClus)
  
  group_met <- clus.df[-c(1:4)]
  group_met <- t(group_met)
  kclus <- as.data.frame(kclus)
  
  group_met <- merge(kclus, group_met, by.x = "row.names", by.y = "row.names")
  row.names(group_met) <- group_met$Row.names
  group_met$Row.names <- NULL
  
  # you want to sum the raw intensity values for the clustered metabolites
  # not the log transformed values, that would be equivalent to multiplying
  # the raw intensity values, much more sensitive to noise if take the product
  # dont know which is better really, but summing the raw values seems to make
  # the most sense.  I think this is equvalent to summing the raw values, because
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
  
  mets.clus.intens <- clus.df[, c(2:4)]
  mets.clus.intens <- merge(mets.clus.intens, group_met_des, by.x= "row.names", by.y = "row.names")
  row.names(mets.clus.intens) <- mets.clus.intens$Row.names
  mets.clus.intens$Row.names <- NULL
  
  mets.clus.intens
}

fing.file <- 
  "C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/metabolics_fingerprint.csv"
mets.clus.intens <- sumClusterIntensities(f = fing.file, linkage = "average", clus.df = health_df)
head(mets.clus.intens[1:10])

mets.intens <- health_df[, -1]
# mets.intens <- data.matrix(mets.intens)
# mets.clus.intens <- data.matrix(mets.clus.intens)
des.lengths <- c(ncol(mets.intens), ncol(mets.clus.intens))
big.df <- data.frame(cbind(health_df$Health_State, mets.intens, mets.clus.intens))

# Processing the Test Set

test.dif.s <- read.csv("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/significant_metabolites/serum/healthstate_anova_wsig_control_test_serum_nonpara.txt")
test.dif.p <- read.csv("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/significant_metabolites/plasma/healthstate_anova_wsig_control_test_plasma_nonpara.txt")

load("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/training_v_test/test_set_tq_normalize.rda")

health_df.test <- d
health_df.test$Patient <- NULL

# make gender binary
health_df.test$Gender <- ifelse(health_df.test$Gender == "F", 0, 1)
health_df.test$Health_State <- ifelse(health_df.test$Health_State == "Healthy", 0, 1)
health_df.test$Smoking_Status <- ifelse(health_df.test$Smoking_Status == "Former", 0, 1)
health_df.test$Patient <- NULL

#----only analyze the samples collected from Serum
health_df_serum.test <- health_df.test[health_df.test$Organ == "Serum", ]

# only keep significant metabolites
colnames(health_df_serum.test)[serum.met.idx]
health_df_serum.test <- health_df_serum.test[c(1:4, serum.met.idx)]

health_df_serum.test$Organ <- NULL

#----only analyze the samples collected from Plasma
health_df_plasma.test <- health_df.test[health_df.test$Organ == "Plasma", ]

# only keep significant metabolites
colnames(health_df_plasma.test)[plasma.met.idx]
health_df_plasma.test <- health_df_plasma.test[c(1:4, plasma.met.idx)]

health_df_plasma.test$Organ <- NULL

health_df.test$Organ <- ifelse(health_df.test$Organ == "Serum", 0, 1)

# ------------------------------------------------------------------------
health_df.test <- health_df.test[, c(2,1,3:ncol(health_df.test))]
head(health_df.test[, 1:10])

mets.clus.intens.test <- sumClusterIntensities(f = fing.file, linkage = "average",
                                               clus.df = health_df.test)
head(mets.clus.intens[1:10])
mets.intens.test <- health_df.test[, -1]
des.lengths <- c(ncol(mets.intens), ncol(mets.clus.intens))
big.df.test <- data.frame(cbind(health_df.test$Health_State, mets.intens.test, mets.clus.intens.test))

###########################
#-------------group metabolites by chemical structure threshold instead
###########################

# ----cluster only significant metabolites in plasma or serum
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
  
  tiff(paste0("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/",
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
  # the most sense.  I think this is equvalent to summing the raw values, because
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
  read.csv("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/metabolics_fingerprint.csv",
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

tiff("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/serum_structure_heatmap.tiff", res = 300,
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

fing <- read.csv("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/metabolics_fingerprint.csv",
                 row.names = 2)
mets <- colnames(health_df_plasma)[4:length(colnames(health_df_plasma))]
fing.serum <- fing[rownames(fing) %in% mets, ] 
fing.serum$row.ID = NULL
d_fing <- dist(fing.serum, method = "binary")

d_fing_mat <- as.matrix(d_fing)
matmat <- data.matrix(d_fing_mat < .3)
matmat <- apply(matmat, 2, as.numeric)
gplots::heatmap.2(matmat, trace = "none", cexCol =.75)


tiff("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/plasma_structure_heatmap.tiff", res = 300,
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


# fascinating if find natural cluster structure, then the 
# clusters of metabolites all have the same fold change direction
# nvm wasnt searching well enough

# cystine glutamic and aspartic acid are clustered together, but have
# different direction of fold change. If you sum these guys, you are going to blow
# it. But a majority vote approach might work well.

#############
#---------Build multimetabolite clusters
#############

# Serum -------------------------------

# remove correlated columns
head(health_df_serum)

single.met.list <- list(health_df_plasma, health_df_serum,
                        health_df_plasma.test, health_df_serum.test)

# regress intensity values on the covariates
for(i in seq_along(single.met.list)) {
  for(j in 4:ncol(single.met.list[[i]])) {
    lm.fit <- lm(single.met.list[[i]][, j] ~ Smoking_Status + Gender, data = single.met.list[[i]])
    single.met.list[[i]][, j] <- lm.fit$residuals
  }
}

health_df_serum <- single.met.list[[2]][, -c(2, 3)]
health_df_serum.test <- single.met.list[[4]][, -c(2, 3)]

cor.cols <- findCorrelation(cor(data.matrix(health_df_serum)), cutoff = .9)
cor.cols

cv.preds <- matrix(ncol = length(2:ncol(health_df_serum)), nrow = nrow(health_df_serum))
cv.probs <- matrix(ncol = length(2:ncol(health_df_serum)), nrow = nrow(health_df_serum))
test.preds <- matrix(ncol = length(2:ncol(health_df_serum)), nrow = nrow(health_df_serum.test))
test.probs <- matrix(ncol = length(2:ncol(health_df_serum)), nrow = nrow(health_df_serum.test))

serum.perf.df <- data.frame(matrix(nrow = ncol(health_df_serum) + 8, ncol = 9))
colnames(serum.perf.df) <- c("model", "error.cv", "sens.cv", "spec.cv", "AUC.cv",
                             "error.test", "sens.test", "spec.test", "AUC.test")
y.cv <- health_df_serum$Health_State
y.test <- health_df_serum.test$Health_State

for (i in 2:ncol(health_df_serum)) {
  
  # LOOCV
  serum.perf.df[i-1, 1] <- colnames(health_df_serum)[i]
  cv.pred <- vector(length = nrow(health_df_serum))
  cv.prob <- vector(length = nrow(health_df_serum))
  for (j in 1:nrow(health_df_serum)) {
    glm.fit <- glm(Health_State ~ ., data = health_df_serum[-j, c(1, i)], family="binomial")
    cv.probs[j, i-1] <- predict(glm.fit, health_df_serum[j, c(1, i)], type ="response")
  }
  cv.thresh <- coords(roc(y.cv, cv.probs[, i-1]), "best", ret = "threshold")
  cv.pred <- ifelse(cv.probs[, i-1] > cv.thresh, 1, 0)
  table(cv.pred, y.cv)
  
  # error rate
  serum.perf.df[i-1, 2] <- mean(cv.pred == y.cv)
  
  # sensitivity.test
  idx <- y.cv == 1
  serum.perf.df[i-1, 3] <- mean(y.cv[idx] == cv.pred[idx])
  
  # specificity.test
  idx <- y.cv == 0
  serum.perf.df[i-1, 4] <- mean(y.cv[idx] == cv.pred[idx])
  
  # auc
  if (i == 2) {
    tiff("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/serum_auc.tiff", res = 300,
         width = 7, height = 7, units = "in")
    plot(roc(y.cv, cv.probs[, i-1]), col = "blue", lty = "dashed", print.auc = T,
         print.auc.x = 0.25, print.auc.y = 0.05, xlim = c(1,0))
  }
  serum.perf.df[i-1, 5] <- as.numeric(auc(y.cv, cv.probs[, i-1]))
  
  cv.preds[, i-1] <- cv.pred
  
  # Test Set
  
  set.seed(835)
  glm.fit <- glm(Health_State ~ ., data = health_df_serum[, c(1, i)], family="binomial")
  test.probs[, i-1] <- predict(glm.fit, health_df_serum.test[, c(1, i)], type ="response")
  glm.pred <- ifelse(test.probs[, i-1] > cv.thresh, 1, 0)
  
  #error rate
  serum.perf.df[i-1, 6] <- mean(glm.pred == y.test)
  
  # sensitivity.test
  idx <- y.test == 1
  serum.perf.df[i-1, 7] <- mean(y.test[idx] == glm.pred[idx])
  
  # specificity.test
  idx <- y.test == 0
  serum.perf.df[i-1, 8] <- mean(y.test[idx] == glm.pred[idx])
  
  # auc
  if (i == 2) {
    plot(roc(y.test, test.probs[, i-1]), col = "black", print.auc = T,
         print.auc.x = 0.25, print.auc.y = 0.1, xlim = c(1,0), add = T)
  }
  serum.perf.df[i-1, 9] <- as.numeric(auc(y.test, test.probs[, i-1]))
  
  test.preds[, i-1] <- glm.pred
}

serum.perf.df
colnames(cv.preds) <- colnames(health_df_serum)[-1]
colnames(cv.probs) <- colnames(health_df_serum)[-1]
cv.pred.multi.maj.vote <- vector(length = nrow(health_df_serum))
cv.probs.multi <- vector(length = nrow(health_df_serum))

serum.perf.df[ncol(health_df_serum), 1] <- "G+A_maj_vote"
serum.perf.df[ncol(health_df_serum) + 1, 1] <- "G+A_avg_prob"

for (j in 1:nrow(health_df_serum)) {
  cv.pred.multi.maj.vote[j] <- mean(cv.preds[j, c(1, 2)])
  cv.probs.multi[j] <- mean(cv.probs[j, c(1, 2)])
}
cv.pred.multi.maj.vote <- ifelse(cv.pred.multi.maj.vote >= .5, 1, 0)
cv.thresh <- coords(roc(y.cv, cv.probs.multi), "best", ret = "threshold")
cv.pred.multi.avg.prob <- ifelse(cv.probs.multi > cv.thresh, 1, 0)

#error rate
serum.perf.df[ncol(health_df_serum), 2] <- mean(cv.pred.multi.maj.vote == y.cv)
serum.perf.df[ncol(health_df_serum) + 1, 2] <- mean(cv.pred.multi.avg.prob == y.cv)

# sensitivity.test
idx <- y.cv == 1
serum.perf.df[ncol(health_df_serum), 3] <- mean(y.cv[idx] == cv.pred.multi.maj.vote[idx])
serum.perf.df[ncol(health_df_serum) + 1, 3] <- mean(y.cv[idx] == cv.pred.multi.avg.prob[idx])

# specificity.test
idx <- y.cv == 0
serum.perf.df[ncol(health_df_serum), 4] <- mean(y.cv[idx] == cv.pred.multi.maj.vote[idx])
serum.perf.df[ncol(health_df_serum) + 1, 4] <- mean(y.cv[idx] == cv.pred.multi.avg.prob[idx])

# auc
plot(roc(y.cv, cv.probs.multi), add = TRUE, col = "orange", lty = "dashed", print.auc = T,
     print.auc.x = .25, print.auc.y = .15)
serum.perf.df[ncol(health_df_serum), 5] <- as.numeric(auc(y.cv, cv.pred.multi.maj.vote))
serum.perf.df[ncol(health_df_serum) + 1, 5] <- as.numeric(auc(y.cv, cv.probs.multi))


colnames(test.preds) <- colnames(health_df_serum)[-1]
colnames(test.probs) <- colnames(health_df_serum)[-1]
test.pred.multi <- vector(length = nrow(health_df_serum.test))
test.probs.multi <- vector(length = nrow(health_df_serum.test))
for (j in 1:nrow(health_df_serum)) {
  test.pred.multi[j] <- mean(test.preds[j, c(1, 2, 3)])
  test.probs.multi[j] <- mean(test.probs[j, c(1, 2, 3)])
}
test.pred.multi.maj.vote <- ifelse(test.pred.multi >= .5, 1, 0)
test.pred.multi.avg.prob <- ifelse(test.probs.multi > cv.thresh, 1, 0)

#error rate
serum.perf.df[ncol(health_df_serum), 6] <- mean(test.pred.multi.maj.vote == y.test)
serum.perf.df[ncol(health_df_serum) + 1, 6] <- mean(test.pred.multi.avg.prob == y.test)

# sensitivity.test
idx <- y.test == 1
serum.perf.df[ncol(health_df_serum), 7] <- mean(y.test[idx] == test.pred.multi.maj.vote[idx])
serum.perf.df[ncol(health_df_serum) + 1, 7] <- mean(y.test[idx] == test.pred.multi.avg.prob[idx])

# specificity.test
idx <- y.test == 0
serum.perf.df[ncol(health_df_serum), 8] <- mean(y.test[idx] == test.pred.multi.maj.vote[idx])
serum.perf.df[ncol(health_df_serum) + 1, 8] <- mean(y.test[idx] == test.pred.multi.avg.prob[idx])

# auc
plot(roc(y.test, test.probs.multi), add = TRUE, col = "red", print.auc = T,
     print.auc.x = .25, print.auc.y = .2)
dev.off()
serum.perf.df[ncol(health_df_serum), 9] <- as.numeric(auc(y.test, test.pred.multi.maj.vote))
serum.perf.df[ncol(health_df_serum) + 1, 9] <- as.numeric(auc(y.test, test.probs.multi))

# make your crazy roc curves you crazy bastard---------------

# Single metabolite------------------------

# LOOCV------------------------------
i <- 2
cv.prob <- vector(length = nrow(health_df_serum))
# LOOCV
for (j in 1:nrow(health_df_serum)) {
  glm.fit <- glm(Health_State ~ ., data = health_df_serum[-j, c(1, i)], family="binomial")
  cv.prob[j]  <- predict(glm.fit, health_df_serum[j, c(1, i)], type ="response")
}

# make ROC
y.cv <- health_df_serum$Health_State
plot(roc(y.cv, cv.prob), col = "grey", print.auc = T)

# External Validation--------------------------------------
set.seed(835)
glm.fit <- glm(Health_State ~ ., data = health_df_serum[, c(1, i)], family="binomial")
test.prob <- predict(glm.fit, health_df_serum.test[, c(1, i)], type ="response")
glm.pred <- ifelse(test.prob > cv.thresh, 1, 0)

# test roc
y.test <- health_df_serum.test$Health_State
plot(roc(y.test, test.prob), add = TRUE, col = "blue", lty = "dashed")

# Multi metabolite------------------------

colnames(cv.preds) <- colnames(health_df_serum)[-1]
colnames(cv.probs) <- colnames(health_df_serum)[-1]
cv.probs.multi <- vector(length = nrow(health_df_serum))

for (j in 1:nrow(health_df_serum)) {
  cv.probs.multi[j] <- mean(cv.probs[j, c(1, 2, 3)])
}
cv.thresh <- coords(roc(y.cv, cv.probs.multi), "best", ret = "threshold")
cv.pred.multi.avg.prob <- ifelse(cv.probs.multi > cv.thresh, 1, 0)

plot(auc(y.cv, cv.probs.multi))

serum.perf.df

kclus.serum
thresh.serum

# Now try clusters -------------------------

# remove correlated columns
clus.met.list <- list(mets.clus.plasma.train, mets.clus.serum.train, 
                      mets.clus.plasma.test, mets.clus.serum.test)

# regress intensity values on the covariates
for(i in seq_along(clus.met.list)) {
  colnames(clus.met.list[[i]])[1] <- "Health_State"
  for(j in 4:ncol(clus.met.list[[i]])) {
    lm.fit <- lm(clus.met.list[[i]][, j] ~ Smoking_Status + Gender, data = clus.met.list[[i]])
    clus.met.list[[i]][, j] <- lm.fit$residuals
  }
}

mets.clus.serum.train <- clus.met.list[[2]][, -c(2, 3)]
mets.clus.serum.test <- clus.met.list[[4]][, -c(2, 3)]

head(mets.clus.serum.train)

cor.cols <- findCorrelation(cor(data.matrix(mets.clus.serum.train)), cutoff = .9)
cor.cols

cv.preds <- matrix(ncol = length(2:ncol(mets.clus.serum.train)), nrow = nrow(mets.clus.serum.train))
cv.probs <- matrix(ncol = length(2:ncol(mets.clus.serum.train)), nrow = nrow(mets.clus.serum.train))
test.preds <- matrix(ncol = length(2:ncol(mets.clus.serum.train)), nrow = nrow(mets.clus.serum.test))
test.probs <- matrix(ncol = length(2:ncol(mets.clus.serum.train)), nrow = nrow(mets.clus.serum.test))

y.cv <- mets.clus.serum.train$Health_State
y.test <- mets.clus.serum.test$Health_State

for (i in 2:ncol(mets.clus.serum.train)) {
  
  # LOOCV
  serum.perf.df[i+ncol(health_df_serum), 1] <- colnames(mets.clus.serum.train)[i]
  cv.pred <- vector(length = nrow(mets.clus.serum.train))
  cv.prob <- vector(length = nrow(mets.clus.serum.train))
  for (j in 1:nrow(mets.clus.serum.train)) {
    glm.fit <- glm(Health_State ~ ., data = mets.clus.serum.train[-j, c(1, i)], family="binomial")
    cv.probs[j, i-1] <- predict(glm.fit, mets.clus.serum.train[j, c(1, i)], type ="response")
  }
  
  cv.thresh <- coords(roc(y.cv, cv.probs[, i-1]), "best", ret = "threshold")
  cv.pred <- ifelse(cv.probs[, i-1] > cv.thresh, 1, 0)
  table(cv.pred, y.cv)
  
  # error rate
  serum.perf.df[i+ncol(health_df_serum), 2] <- mean(cv.pred == y.cv)
  
  # sensitivity.test
  idx <- y.cv == 1
  serum.perf.df[i+ncol(health_df_serum), 3] <- mean(y.cv[idx] == cv.pred[idx])
  
  # specificity.test
  idx <- y.cv == 0
  serum.perf.df[i+ncol(health_df_serum), 4] <- mean(y.cv[idx] == cv.pred[idx])
  
  # auc
  # plot(roc(y.cv, cv.probs[, i-1]))
  serum.perf.df[i+ncol(health_df_serum), 5] <- as.numeric(auc(y.cv, cv.probs[, i-1]))
  
  cv.preds[, i-1] <- cv.pred
  
  # Test Set
  
  set.seed(835)
  glm.fit <- glm(Health_State ~ ., data = mets.clus.serum.train[, c(1, i)], family="binomial")
  test.probs[, i-1] <- predict(glm.fit, mets.clus.serum.test[, c(1, i)], type ="response")
  glm.pred <- ifelse(test.probs[, i-1] > cv.thresh, 1, 0)
  
  #error rate
  serum.perf.df[i+ncol(health_df_serum), 6] <- mean(glm.pred == y.test)
  
  # sensitivity.test
  idx <- y.test == 1
  serum.perf.df[i+ncol(health_df_serum), 7] <- mean(y.test[idx] == glm.pred[idx])
  
  # specificity.test
  idx <- y.test == 0
  serum.perf.df[i+ncol(health_df_serum), 8] <- mean(y.test[idx] == glm.pred[idx])
  
  # auc
  # plot(roc(y.test, test.probs[, i-1]))
  serum.perf.df[i+ncol(health_df_serum), 9] <- as.numeric(auc(y.test, test.probs[, i-1]))
  
  test.preds[, i-1] <- glm.pred
}

serum.perf.df

write.csv(serum.perf.df, "serum_performances.csv")

# Plasma ------------------------

# remove correlated columns

health_df_plasma <- single.met.list[[1]][, -c(2, 3)]
health_df_plasma.test <- single.met.list[[3]][, -c(2, 3)]

cv.preds <- matrix(ncol = length(2:ncol(health_df_plasma)), nrow = nrow(health_df_plasma))
cv.probs <- matrix(ncol = length(2:ncol(health_df_plasma)), nrow = nrow(health_df_plasma))
test.preds <- matrix(ncol = length(2:ncol(health_df_plasma)), nrow = nrow(health_df_plasma.test))
test.probs <- matrix(ncol = length(2:ncol(health_df_plasma)), nrow = nrow(health_df_plasma.test))

plasma.perf.df <- data.frame(matrix(nrow = ncol(health_df_plasma) + 12, ncol = 9))
colnames(plasma.perf.df) <- c("model", "error.cv", "sens.cv", "spec.cv", "AUC.cv",
                              "error.test", "sens.test", "spec.test", "AUC.test")
y.cv <- health_df_plasma$Health_State
y.test <- health_df_plasma.test$Health_State

i <- 14
colnames(health_df_plasma)

for (i in 2:ncol(health_df_plasma)) {
  
  # LOOCV
  plasma.perf.df[i-1, 1] <- colnames(health_df_plasma)[i]
  cv.pred <- vector(length = nrow(health_df_plasma))
  cv.prob <- vector(length = nrow(health_df_plasma))
  for (j in 1:nrow(health_df_plasma)) {
    glm.fit <- glm(Health_State ~ ., data = health_df_plasma[-j, c(1, i)], family="binomial")
    cv.probs[j, i-1] <- predict(glm.fit, health_df_plasma[j, c(1, i)], type ="response")
  }
  cv.thresh <- coords(roc(y.cv, cv.probs[, i-1]), "best", ret = "threshold")
  cv.pred <- ifelse(cv.probs[, i-1] > cv.thresh, 1, 0)
  table(cv.pred, y.cv)
  
  # error rate
  plasma.perf.df[i-1, 2] <- mean(cv.pred == y.cv)
  
  # sensitivity.test
  idx <- y.cv == 1
  plasma.perf.df[i-1, 3] <- mean(y.cv[idx] == cv.pred[idx])
  
  # specificity.test
  idx <- y.cv == 0
  plasma.perf.df[i-1, 4] <- mean(y.cv[idx] == cv.pred[idx])
  
  # auc
  if (i == 14) {
    tiff("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/qsar_patient_profiles/plasma_auc.tiff", res = 300,
         width = 7, height = 7, units = "in")
    plot(roc(y.cv, cv.probs[, i-1]), col = "grey", lty = "dashed", print.auc = T,
         print.auc.x = 0.25, print.auc.y = 0.025, xlim = c(1,0))
  }
  plasma.perf.df[i-1, 5] <- as.numeric(auc(y.cv, cv.probs[, i-1]))
  
  cv.preds[, i-1] <- cv.pred
  
  # Test Set
  
  set.seed(835)
  glm.fit <- glm(Health_State ~ ., data = health_df_plasma[, c(1, i)], family="binomial")
  test.probs[, i-1] <- predict(glm.fit, health_df_plasma.test[, c(1, i)], type ="response")
  glm.pred <- ifelse(test.probs[, i-1] > cv.thresh, 1, 0)
  
  #error rate
  plasma.perf.df[i-1, 6] <- mean(glm.pred == y.test)
  
  # sensitivity.test
  idx <- y.test == 1
  plasma.perf.df[i-1, 7] <- mean(y.test[idx] == glm.pred[idx])
  
  # specificity.test
  idx <- y.test == 0
  plasma.perf.df[i-1, 8] <- mean(y.test[idx] == glm.pred[idx])
  
  # auc
  if (i == 14) {
    plot(roc(y.test, test.probs[, i-1]), col = "black", print.auc = T,
         print.auc.x = 0.25, print.auc.y = 0.075, xlim = c(1,0), add = T)
  }
  plasma.perf.df[i-1, 9] <- as.numeric(auc(y.test, test.probs[, i-1]))
  
  test.preds[, i-1] <- glm.pred
  
  
}
colnames(cv.preds) <- colnames(health_df_plasma)[-1]
colnames(cv.probs) <- colnames(health_df_plasma)[-1]
colnames(test.preds) <- colnames(health_df_plasma)[-1]
colnames(test.probs) <- colnames(health_df_plasma)[-1]


plasma.perf.df

kclus.plasma

i <- 4

# cv.preds <- matrix(ncol = length(1:max(kclus.plasma)), nrow = nrow(health_df_plasma))
# cv.probs <- matrix(ncol = length(1:max(kclus.plasma)), nrow = nrow(health_df_plasma))
# test.preds <- matrix(ncol = length(1:max(kclus.plasma)), nrow = nrow(health_df_plasma.test))
# test.probs <- matrix(ncol = length(1:max(kclus.plasma)), nrow = nrow(health_df_plasma.test))

for (i in 1:max(kclus.plasma)) {

  cv.pred.multi.maj.vote <- vector(length = nrow(health_df_plasma))
  cv.probs.multi <- vector(length = nrow(health_df_plasma))
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 1] <- paste0("V", i, "_maj_vote")
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 1] <- paste0("V", i, "_avg_prob")
  for (j in 1:nrow(health_df_plasma)) {
    cv.pred.multi.maj.vote[j] <- mean(cv.preds[j, which(kclus.plasma == i)])
    cv.probs.multi[j] <- mean(cv.probs[j, which(kclus.plasma == i)])
  }
  cv.pred.multi.maj.vote <- ifelse(cv.pred.multi.maj.vote >= .5, 1, 0)
  cv.thresh <- coords(roc(y.cv, cv.probs.multi), "best", ret = "threshold")
  cv.pred.multi.avg.prob <- ifelse(cv.probs.multi > cv.thresh, 1, 0)
  
  #error rate
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 2] <- mean(cv.pred.multi.maj.vote == y.cv)
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 2] <- mean(cv.pred.multi.avg.prob == y.cv)
  
  # sensitivity.test
  idx <- y.cv == 1
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 3] <- mean(y.cv[idx] == cv.pred.multi.maj.vote[idx])
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 3] <- mean(y.cv[idx] == cv.pred.multi.avg.prob[idx])
  
  # specificity.test
  idx <- y.cv == 0
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 4] <- mean(y.cv[idx] == cv.pred.multi.maj.vote[idx])
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 4] <- mean(y.cv[idx] == cv.pred.multi.avg.prob[idx])
  
  # auc
  # plot(roc(y.cv, cv.probs.multi))
  if (i == 4) {
    plot(roc(y.cv, cv.probs.multi), col = "orange", lty = "solid", print.auc = T,
         print.auc.x = 0.25, print.auc.y = 0.125, xlim = c(1,0))
  }
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 5] <- as.numeric(auc(y.cv, cv.pred.multi.maj.vote))
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 5] <- as.numeric(auc(y.cv, cv.probs.multi))
  
  # Test set
  

  test.pred.multi <- vector(length = nrow(health_df_plasma.test))
  test.probs.multi <- vector(length = nrow(health_df_plasma.test))
  for (j in 1:nrow(health_df_plasma.test)) {
    test.pred.multi[j] <- mean(test.preds[j, which(kclus.plasma == i)])
    test.probs.multi[j] <- mean(test.probs[j, which(kclus.plasma == i)])
  }
  test.pred.multi.maj.vote <- ifelse(test.pred.multi >= .5, 1, 0)
  test.pred.multi.avg.prob <- ifelse(test.probs.multi > cv.thresh, 1, 0)
  
  #error rate
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 6] <- mean(test.pred.multi.maj.vote == y.test)
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 6] <- mean(test.pred.multi.avg.prob == y.test)
  
  # sensitivity.test
  idx <- y.test == 1
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 7] <- mean(y.test[idx] == test.pred.multi.maj.vote[idx])
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 7] <- mean(y.test[idx] == test.pred.multi.avg.prob[idx])
  
  # specificity.test
  idx <- y.test == 0
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 8] <- mean(y.test[idx] == test.pred.multi.maj.vote[idx])
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 8] <- mean(y.test[idx] == test.pred.multi.avg.prob[idx])
  
  # auc
  if (i == 4) {
    plot(roc(y.test, test.pred.multi.maj.vote), col = "blue", lty = "dashed", print.auc = T,
         print.auc.x = 0.25, print.auc.y = 0.175, xlim = c(1,0), add = T)
  }
  dev.off()
  plasma.perf.df[ncol(health_df_plasma) + 2*i-2, 9] <- as.numeric(auc(y.test, test.pred.multi.maj.vote))
  plasma.perf.df[ncol(health_df_plasma) + 2*i-1, 9] <- as.numeric(auc(y.test, test.probs.multi))
}


plasma.perf.df

kclus.plasma
thresh.plasma

# Now try clusters -------------------------

# remove correlated columns

mets.clus.plasma.train <- clus.met.list[[1]][, -c(2, 3)]
mets.clus.plasma.test <- clus.met.list[[3]][, -c(2, 3)]

head(mets.clus.plasma.train)

cor.cols <- findCorrelation(cor(data.matrix(mets.clus.plasma.train)), cutoff = .9)
cor.cols

cv.preds <- matrix(ncol = length(2:ncol(mets.clus.plasma.train)), nrow = nrow(mets.clus.plasma.train))
cv.probs <- matrix(ncol = length(2:ncol(mets.clus.plasma.train)), nrow = nrow(mets.clus.plasma.train))
test.preds <- matrix(ncol = length(2:ncol(mets.clus.plasma.train)), nrow = nrow(mets.clus.plasma.test))
test.probs <- matrix(ncol = length(2:ncol(mets.clus.plasma.train)), nrow = nrow(mets.clus.plasma.test))

y.cv <- mets.clus.plasma.train$Health_State
y.test <- mets.clus.plasma.test$Health_State


cur.row <- ncol(health_df_plasma) + 2*max(kclus.plasma) - 2

for (i in 2:ncol(mets.clus.plasma.train)) {
  
  # LOOCV
  plasma.perf.df[i+cur.row, 1] <- colnames(mets.clus.plasma.train)[i]
  cv.pred <- vector(length = nrow(mets.clus.plasma.train))
  cv.prob <- vector(length = nrow(mets.clus.plasma.train))
  for (j in 1:nrow(mets.clus.plasma.train)) {
    glm.fit <- glm(Health_State ~ ., data = mets.clus.plasma.train[-j, c(1, i)], family="binomial")
    cv.probs[j, i-1] <- predict(glm.fit, mets.clus.plasma.train[j, c(1, i)], type ="response")
  }
  cv.thresh <- coords(roc(y.cv, cv.probs[, i-1]), "best", ret = "threshold")
  cv.pred <- ifelse(cv.probs[, i-1] > cv.thresh, 1, 0)
  table(cv.pred, y.cv)
  
  # error rate
  plasma.perf.df[i+cur.row, 2] <- mean(cv.pred == y.cv)
  
  # sensitivity.test
  idx <- y.cv == 1
  plasma.perf.df[i+cur.row, 3] <- mean(y.cv[idx] == cv.pred[idx])
  
  # specificity.test
  idx <- y.cv == 0
  plasma.perf.df[i+cur.row, 4] <- mean(y.cv[idx] == cv.pred[idx])
  
  # auc
  # if (i == 2) {
  #   plot(roc(y.cv, cv.probs[, i-1]), col = "orange", lty = "dashed", print.auc = T,
  #        print.auc.x = 0.25, print.auc.y = 0.125, xlim = c(1,0))
  # }
  plasma.perf.df[i+cur.row, 5] <- as.numeric(auc(y.cv, cv.probs[, i-1]))
  
  cv.preds[, i-1] <- cv.pred
  
  # Test Set
  
  set.seed(835)
  glm.fit <- glm(Health_State ~ ., data = mets.clus.plasma.train[, c(1, i)], family="binomial")
  test.probs[, i-1] <- predict(glm.fit, mets.clus.plasma.test[, c(1, i)], type ="response")
  glm.pred <- ifelse(test.probs[, i-1] > cv.thresh, 1, 0)
  
  #error rate
  plasma.perf.df[i+cur.row, 6] <- mean(glm.pred == y.test)
  
  # sensitivity.test
  idx <- y.test == 1
  plasma.perf.df[i+cur.row, 7] <- mean(y.test[idx] == glm.pred[idx])
  
  # specificity.test
  idx <- y.test == 0
  plasma.perf.df[i+cur.row, 8] <- mean(y.test[idx] == glm.pred[idx])
  
  auc
  plasma.perf.df[i+cur.row, 9] <- as.numeric(auc(y.test, test.probs[, i-1]))
  test.preds[, i-1] <- glm.pred
}

plasma.perf.df

write.csv(plasma.perf.df, "plasma_performances.csv")

kclus.plasma
