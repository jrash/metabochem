setwd("C:/Users/jrash/ubuntu_share/Google Drive/FourchesLab/metabolimics/significant_metabolites")

train.dif.s <- read.csv("serum/healthstate_anova_wsig_control_training_serum_nonpara.txt")
train.dif.p <- read.csv("plasma/healthstate_anova_wsig_control_training_plasma_nonpara.txt")

test.dif.s <- read.csv("serum/healthstate_anova_wsig_control_test_serum_nonpara.txt")
test.dif.p <- read.csv("plasma/healthstate_anova_wsig_control_test_plasma_nonpara.txt")

# show the p-values

train.sig.s <- train.dif.s$adj_pvalues < .075
test.sig.s <- test.dif.s$adj_pvalues < .075

train.sig.p <- train.dif.p$adj_pvalues < .075
test.sig.p <- test.dif.p$adj_pvalues < .075


train.clus.s <- train.dif.s$variables %in% c("glutamic_acid",
                                           "aspartic_acid", "N-acetylglutamate",
                                           "glutamine", "asparagine", "cystine")
sum(train.clus.s)
table(train.sig.s, train.clus.s)

p.values.fish <- vector(length = 4)
names(p.values.fish) <- c("train.s", "test.s", "train.p", "test.p")

p.values.fish[1] <- fisher.test(train.sig.s, train.clus.s, alternative = "greater")$p.value

test.clus.s <- test.dif.s$variables %in% c("glutamic_acid",
                                             "aspartic_acid", "N-acetylglutamate",
                                             "glutamine", "asparagine", "cystine")
sum(test.clus.s)
table(test.sig.s, test.clus.s)
p.values.fish[2] <- fisher.test(test.sig.s, test.clus.s, alternative = "greater")$p.value

train.clus.p <- train.dif.p$variables %in% c("glutamic_acid",
                                             "aspartic_acid", "N-acetylglutamate",
                                             "glutamine", "asparagine", "cystine")
sum(train.clus.p)
table(train.sig.p, train.clus.p)
p.values.fish[3] <- fisher.test(train.sig.p, train.clus.p, alternative = "greater")$p.value

test.clus.p <- test.dif.p$variables %in% c("glutamic_acid",
                                           "aspartic_acid", "N-acetylglutamate",
                                           "glutamine", "asparagine", "cystine")
sum(test.clus.p)
table(test.sig.p, test.clus.p)
p.values.fish[4] <- fisher.test(test.sig.p, test.clus.p, alternative = "greater")$p.value


p.adjust(p.values.fish, method = "BH") < .05
p.adjust(p.values.fish, method = "BH") < .01
p.adjust(p.values.fish, method = "BH") < .001

clus.idx <- train.dif.s$variables %in% c("glutamic_acid",
                             "aspartic_acid", "N-acetylglutamate",
                             "glutamine", "asparagine", "cystine")
train.dif.s[clus.idx, ][5, ]
test.dif.s[clus.idx, ][5, ]
train.dif.p[clus.idx, ][5, ]
test.dif.p[clus.idx, ][5, ]

#------------Pathway analysis enrichment

p.values.fish <- vector(length = 8)

#----Serum

# Alanine, aspartate and glutamate metabolism

num.in.path <- c(9, 4, 6, 5, 4, 7, 7, 12)
num.in.clus <- c(4, 2, 0, 0, 2, 1, 2, 5)

for (i in 1:8) {
  path.mem <- rep(0, 130)
  path.mem[1:num.in.path[i]] <- rep(1, num.in.path[i])
  
  clus.mem <- rep(0, 130)
  if (num.in.clus[i] != 0) {
    clus.mem[1:num.in.clus[i]] <- rep(1, num.in.clus[i])
    p.values.fish[i] <- fisher.test(path.mem, clus.mem, alternative = "greater")$p.value
  } else {
    p.values.fish[i] <- 1
  }
  
}

p.adjust(p.values.fish, method = "BH") < .05
p.adjust(p.values.fish, method = "BH") < .01
p.adjust(p.values.fish, method = "BH") < .001



