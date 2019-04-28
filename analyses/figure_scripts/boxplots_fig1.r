library(ggplot2)
library(rprojroot)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
setwd(find_root("metabochem.Rproj"))
setwd("analyses/data/")

train.dif.s <- read.csv("healthstate_anova_wsig_control_training_serum_nonpara.txt")
train.dif.p <- read.csv("healthstate_anova_wsig_control_training_plasma_nonpara.txt")

load("training_set_tq_nonormalize.rda")
d$Set <- paste0(rep("Training", nrow(d)), d$Health_State)
d <- d[, c(1:5, 136, 6:135)]
d_serum <- d[d$Organ == "Serum", ]
d_plasma <- d[d$Organ == "Plasma", ]

train.sig.serum <- d_serum[, c(3, which(train.dif.s$adj_pvalues < .075) + 6)]
head(train.sig.serum)

train.sig.plasma <- d_plasma[, c(3, which(train.dif.p$adj_pvalues < .075) + 6)]
dim(train.sig.plasma)

test.dif.s <- read.csv("healthstate_anova_wsig_control_test_serum_nonpara.txt")
test.dif.p <- read.csv("healthstate_anova_wsig_control_test_plasma_nonpara.txt")

load("test_set_tq_normalize.rda")
d$Set <- paste0(rep("Training", nrow(d)), d$Health_State)
d <- d[, c(1:5, 136, 6:135)]
d_serum <- d[d$Organ == "Serum", ]
d_plasma <- d[d$Organ == "Plasma", ]


sum(test.dif.s$significance_adj)

test.sig.serum <- d_serum[, c(3, which(train.dif.s$adj_pvalues < .075) + 6)]
dim(test.sig.serum)

test.sig.plasma <- d_plasma[, c(3, which(train.dif.p$adj_pvalues < .075) + 6)]
head(test.sig.plasma)

# insert covariate adjustment here

# single.met.list <- list(health_df_serum, test.sig.serum,
#                         health_df_plasma, test.sig.plasma)
# 
# head(single.met.list[[4]])
# 
# # regress intensity values on the covariates
# for(i in seq_along(single.met.list)) {
#   for(j in 7:ncol(single.met.list[[i]])) {
#     lm.fit <- lm(single.met.list[[i]][, j] ~ Smoking_Status + Gender, data = single.met.list[[i]])
#     single.met.list[[i]][, j] <- lm.fit$residuals
#   }
# }



train.sig.serum$Set <- rep("TRAIN", nrow(train.sig.serum))
test.sig.serum$Set <- rep("TEST", nrow(test.sig.serum))

sig.serum <- rbind(train.sig.serum, test.sig.serum)

head(sig.serum)

dfmelt<-melt(sig.serum, measure.vars = 2:6)
head(dfmelt)
levels(dfmelt$variable) <- sapply(sub("_", " ", levels(dfmelt$variable)), simpleCap)
dfmelt$Set <- as.factor(dfmelt$Set)
dfmelt$Set <- factor(dfmelt$Set, levels = c("TRAIN", "TEST"))

# show the p-values

train.dif.s[train.dif.s$adj_pvalues < .075, ]
test.dif.s[train.dif.s$adj_pvalues < .075, ]

train.dif.p[train.dif.p$adj_pvalues < .075, ]
test.dif.p[train.dif.p$adj_pvalues < .075, ]


p <- ggplot(dfmelt, aes(x=Health_State, y=value,fill=Health_State))+
  geom_boxplot()+
  facet_wrap(variable~Set, nrow=2, scales = "fixed")+
  labs(x="Health State")+
  theme(axis.text.x=element_text(angle=45, hjust = 1),
        strip.text.x = element_text(color = "Black", size = 12, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"))
print(p)


# Make list of variable names to loop over.
var_list = combn(names(iris)[1:3], 2, simplify=FALSE)

sig.serum$Set <- factor(sig.serum$Set, levels = c("TRAIN", "TEST"))

ggplot(sig.serum[, c(1, 7, 2)], aes_string(x = "Health_State", y = colnames(sig.serum)[2],
       fill = "Health_State"))+
  geom_boxplot()+
  facet_grid(.~Set) +
  labs(x="")+
  labs(y="Log Base 2 Intensity")+
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
        axis.text.y=element_text(size = 18),
        strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"),
        axis.title.y = element_text(size = 18))

i <- 2

tiff(paste0("../analyses/serum_",
            colnames(sig.serum)[i],".tiff"), res = 300,
     width = 7, height = 7, units = "in")
p <- ggplot(sig.serum[, c(1, 7, i)], aes_string(x = "Health_State", y = colnames(sig.serum)[i]))+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
  facet_grid(.~Set) +
  labs(x="")+
  labs(y="Log Base 2 Intensity")+
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
        axis.text.y=element_text(size = 18),
        strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"),
        axis.title.y = element_text(size = 18),
        legend.position="none")

print(p)
dev.off()

for (i in 3:5) {
  tiff(paste0("../analyses/serum_",
              colnames(sig.serum)[i],".tiff"), res = 300,
              width = 7, height = 7, units = "in")
  p <- ggplot(sig.serum[, c(1, 7, i)], aes_string(x = "Health_State", y = colnames(sig.serum)[i]))+
    geom_boxplot(outlier.alpha = 0)+
    geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
    facet_grid(.~Set) +
    labs(x="")+
    labs(y="")+
    theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
          axis.text.y=element_text(size = 18),
          strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
          strip.text.y = element_text(color = "Black"),
          strip.background = element_rect(colour="black", fill="#CCCCFF"),
          axis.title.y = element_text(size = 18),
          legend.position="none")

  print(p)
  dev.off()
}

i <- 6

tiff(paste0("../analyses/serum_",
            colnames(sig.serum)[i],".tiff"), res = 300,
     width = 7, height = 7, units = "in")
p <- ggplot(sig.serum[, c(1, 7, i)], aes_string(x = "Health_State", y = colnames(sig.serum)[i]))+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
  facet_grid(.~Set) +
  labs(x="")+
  labs(y="")+
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
        axis.text.y=element_text(size = 18),
        strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"),
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.position="top")

print(p)
dev.off()


train.sig.plasma$Set <- rep("TRAIN", nrow(train.sig.plasma))
test.sig.plasma$Set <- rep("TEST", nrow(test.sig.plasma))

sig.plasma <- rbind(train.sig.plasma, test.sig.plasma)

head(sig.plasma)

# Make list of variable names to loop over.
sig.plasma$Set <- factor(sig.plasma$Set, levels = c("TRAIN", "TEST"))

colnames(sig.plasma) <- sub("-|-", "", colnames(sig.plasma))
colnames(sig.plasma) <- sub("-|-", "", colnames(sig.plasma))
colnames(sig.plasma) <- sub("[0-9]", "", colnames(sig.plasma))

i <- 2

tiff(paste0("../analyses/plasma_",
            colnames(sig.plasma)[i],".tiff"), res = 300,
     width = 7, height = 7, units = "in")
p <- ggplot(sig.plasma[, c(1, 12, i)], aes_string(x = "Health_State", y = colnames(sig.plasma)[i]))+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
  facet_grid(.~Set) +
  labs(x="")+
  labs(y="Log Base 2 Intensity")+
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
        axis.text.y=element_text(size = 18),
        strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"),
        axis.title.y = element_text(size = 18),
        legend.position="none")

print(p)
dev.off()


for (i in 3:10) {
  tiff(paste0("../analyses/plasma_",
              colnames(sig.plasma)[i],".tiff"), res = 300,
       width = 7, height = 7, units = "in")
  p <- ggplot(sig.plasma[, c(1, 12, i)], aes_string(x = "Health_State", y = colnames(sig.plasma)[i]))+
    geom_boxplot(outlier.alpha = 0)+
    geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
    facet_grid(.~Set) +
    labs(x="")+
    labs(y="")+
    theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
          axis.text.y=element_text(size = 18),
          strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
          strip.text.y = element_text(color = "Black"),
          strip.background = element_rect(colour="black", fill="#CCCCFF"),
          axis.title.y = element_text(size = 18),
          legend.position="none")

  print(p)
  dev.off()
}

i <- 8

tiff(paste0("../analyses/plasma_",
            colnames(sig.plasma)[i],".tiff"), res = 300,
     width = 7, height = 7, units = "in")
p <- ggplot(sig.plasma[, c(1, 12, i)], aes_string(x = "Health_State", y = colnames(sig.plasma)[i]))+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
  facet_grid(.~Set) +
  labs(x="")+
  labs(y="Log Base 2 Intensity")+
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
        axis.text.y=element_text(size = 18),
        strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"),
        axis.title.y = element_text(size = 18),
        legend.position="none")

print(p)
dev.off()


i <- 7

tiff(paste0("../analyses/plasma_",
            colnames(sig.plasma)[i],".tiff"), res = 300,
     width = 8, height = 7, units = "in")
p <- ggplot(sig.plasma[, c(1, 12, i)], aes_string(x = "Health_State", y = colnames(sig.plasma)[i]))+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
  facet_grid(.~Set) +
  labs(x="")+
  labs(y="")+
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
        axis.text.y=element_text(size = 18),
        strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"),
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size=18))

print(p)
dev.off()


i <- 11

tiff(paste0("../analyses/plasma_",
            colnames(sig.plasma)[i],".tiff"), res = 300,
     width = 8, height = 7, units = "in")
p <- ggplot(sig.plasma[, c(1, 12, i)], aes_string(x = "Health_State", y = colnames(sig.plasma)[i]))+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = .2), aes(color = Health_State)) +
  facet_grid(.~Set) +
  labs(x="")+
  labs(y="")+
  theme(axis.text.x=element_text(angle=45, hjust = 1, size = 18),
        axis.text.y=element_text(size = 18),
        strip.text.x = element_text(color = "Black", size = 18, face = "bold"),
        strip.text.y = element_text(color = "Black"),
        strip.background = element_rect(colour="black", fill="#CCCCFF"),
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size=18))

print(p)
dev.off()


boxplot(test.sig.serum$xylose~test.sig.serum$Health_State, main = "xylose")
stripchart(d_serum$xylose ~ d_serum$Health_State, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

boxplot(d_serum$glutamic_acid~d_serum$Health_State, main = "glutamate")
stripchart(d_serum$glutamic_acid~d_serum$Health_State, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

boxplot(d_serum$aspartic_acid~d_serum$Health_State, main = "aspartate")
stripchart(d_serum$aspartic_acid~d_serum$Health_State, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')






boxplot(d_serum$xylose~d_serum$Health_State, main = "xylose")
stripchart(d_serum$xylose ~ d_serum$Health_State, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

boxplot(d_serum$glutamic_acid~d_serum$Health_State, main = "glutamate")
stripchart(d_serum$glutamic_acid~d_serum$Health_State, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

boxplot(d_serum$aspartic_acid~d_serum$Health_State, main = "aspartate")
stripchart(d_serum$aspartic_acid~d_serum$Health_State, vertical = TRUE,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
