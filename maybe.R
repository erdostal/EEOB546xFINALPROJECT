setwd("~/Box Sync/Classes/BioDataSkills_BCB546X_Spr17/EEOB546xFINALPROJECT")
library(geomorph)
library(reshape2)
library(ggplot2)
library(contrast)
#input and structure data
cluster_table <- t(head(as.matrix(read.table("table_clusters.txt", header = T)),1000))

#Table 2 (two-way ANOVA test, as done in the paper)
new_table <- cluster_table
row.names(new_table) <- c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5")
melted <- melt(new_table)
Y <- melted$value
Species <- as.factor(melted$Var1)
Cluster <- as.factor(melted$Var2)
linear_model <- lm(Y~Cluster*Species)
anova(linear_model) #same as in paper




###Figure 3 (with colors, but no BH correction) 
#differences between this figure and Simon's figure is highlighted on Simon_Fig3_discrepency_noBH.png

cluster_table <- as.data.frame(t(cluster_table))
cluster_table$A1ave <- rowMeans(subset(cluster_table, select = c(A1_155_, A1_073_, A1_097_), na.rm = TRUE))
cluster_table$A2ave <- rowMeans(subset(cluster_table, select = c(A2_255_, A2_034_, A2_044_, A2_099_, A2_101_), na.rm = TRUE))
cluster_table$D5ave <- rowMeans(subset(cluster_table, select = c(D5_002_, D5_031_, D5_004_, D5_053_, D5_ggg_), na.rm = TRUE))

#Constrasts
Contrast.A1.A2 <- contrast(linear_model, list(Cluster = levels(Cluster), Species = "A1"), list(Cluster = levels(Cluster), Species = "A2"))
Contrast.A1.D5 <- contrast(linear_model, list(Cluster = levels(Cluster), Species = "A1"), list(Cluster = levels(Cluster), Species = "D5"))
Contrast.A2.D5 <- contrast(linear_model, list(Cluster = levels(Cluster), Species = "A2"), list(Cluster = levels(Cluster), Species = "D5"))


#A1A2 Differences
cluster_table$A1A2.direction <- ifelse((Contrast.A1.A2$Pvalue < 0.05), 1,0)
cluster_table$A1A2.direction.BH <- ifelse((p.adjust(Contrast.A1.A2$Pvalue, method = "BH") < 0.05), 1,0)
cluster_table$A1A2.sign <- sign(Contrast.A1.A2$testStat)
cluster_table$A1A2.significant <- as.factor(cluster_table$A1A2.sign * cluster_table$A1A2.direction)
cluster_table$A1A2.significant.BH <- as.factor(cluster_table$A1A2.sign * cluster_table$A1A2.direction.BH)

#A1D5 Differences
cluster_table$A1D5.direction <- ifelse((Contrast.A1.D5$Pvalue < 0.05), 1,0)
cluster_table$A1D5.direction.BH <- ifelse((p.adjust(Contrast.A1.D5$Pvalue, method = "BH") < 0.05), 1,0)
cluster_table$A1D5.sign <- sign(Contrast.A1.D5$testStat)
cluster_table$A1D5.significant <- as.factor(cluster_table$A1D5.sign * cluster_table$A1D5.direction)
cluster_table$A1D5.significant.BH <- as.factor(cluster_table$A1D5.sign * cluster_table$A1D5.direction.BH)

#A2D5 Differences
cluster_table$A2D5.direction <- ifelse((Contrast.A2.D5$Pvalue < 0.05), 1,0)
cluster_table$A2D5.direction.BH <- ifelse((p.adjust(Contrast.A2.D5$Pvalue, method = "BH") < 0.05), 1,0)
cluster_table$A2D5.sign <- sign(Contrast.A2.D5$testStat)
cluster_table$A2D5.significant <- as.factor(cluster_table$A2D5.sign * cluster_table$A2D5.direction)
cluster_table$A2D5.significant.BH <- as.factor(cluster_table$A2D5.sign * cluster_table$A2D5.direction.BH)


####Table 3 (not same numbers as paper, but same methods (excpet I didn't do BH correction))
A1A2.clusters <- sum(cluster_table$A1A2.significant != 0)
A1D5.clusters <- sum(cluster_table$A1D5.significant != 0)
A2D5.clusters <- sum(cluster_table$A2D5.significant != 0)
A1A2.clusters
A1D5.clusters
A2D5.clusters


A1A2.clusters.BH <- sum(cluster_table$A1A2.significant.BH != 0)
A1D5.clusters.BH <- sum(cluster_table$A1D5.significant.BH != 0)
A2D5.clusters.BH <- sum(cluster_table$A2D5.significant.BH != 0)
A1A2.clusters.BH
A1D5.clusters.BH
A2D5.clusters.BH




#Plotting Figure 3 (With color, with and without BH correction)
ggplot(A1ave~A2ave, data=cluster_table, mapping = aes(x = A1ave, y = A2ave, color = A1A2.significant)) + geom_point() + geom_abline(intercept = 0, slope = 1) + xlim(0, 3500) + ylim(0,3500)
ggplot(A1ave~A2ave, data=cluster_table, mapping = aes(x = A1ave, y = A2ave, color = A1A2.significant.BH)) + geom_point() + geom_abline(intercept = 0, slope = 1) + xlim(0, 3500) + ylim(0,3500)

ggplot(A1ave~D5ave, data=cluster_table, mapping = aes(x = A1ave, y = D5ave, color = A1D5.significant)) + geom_point() + geom_abline(intercept = 0, slope = 1) + xlim(0, 3500) + ylim(0,3500)
ggplot(A1ave~D5ave, data=cluster_table, mapping = aes(x = A1ave, y = D5ave, color = A1D5.significant.BH)) + geom_point() + geom_abline(intercept = 0, slope = 1) + xlim(0, 3500) + ylim(0,3500)

ggplot(A2ave~D5ave, data=cluster_table, mapping = aes(x = A2ave, y = D5ave, color = A2D5.significant)) + geom_point() + geom_abline(intercept = 0, slope = 1) + xlim(0, 3500) + ylim(0,3500)
ggplot(A2ave~D5ave, data=cluster_table, mapping = aes(x = A2ave, y = D5ave, color = A2D5.significant.BH)) + geom_point() + geom_abline(intercept = 0, slope = 1) + xlim(0, 3500) + ylim(0,3500)





########################################################
########### What Simon Should Have Done  ###############
########################################################



#factors for manova tests
AD <- as.factor(c("A", "A", "A", "A", "A", "A", "A", "D", "D", "D", "D", "A", "D"))
A1A2D5 <- as.factor(c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5"))
A1.ancestral <- as.factor(c("A", "Ancestral", "A", "A", "A", "A", "A", "D", "D", "D", "D", "A", "D"))


#manova for 2 species (A and D)
two.species <- procD.lm(cluster_table~AD)
two.z <- two.species$aov.table[[1,6]]
two.rand.mean <- mean(two.species$random.SS)
two.rand.sd <- sd(two.species$random.SS)
two.species #significant 


#manova for 3 species (A1, A2, and D5)
A1A2D5.species <- procD.lm(cluster_table~A1A2D5)
A1A2D5.z <- A1A2D5.species$aov.table[[1,6]]
A1A2D5.rand.mean <- mean(A1A2D5.species$random.SS)
A1A2D5.rand.sd <- sd(A1A2D5.species$random.SS)
A1A2D5.species #significant


#manova for 2 species, and wild progenitor (D5, A, and A1-73 as wild progenitor)
A1.ancestral.species <- procD.lm(cluster_table~A1.ancestral)
A1.ancestral.z <- A1.ancestral.species$aov.table[[1,6]]
A1.ancestral.rand.mean <- mean(A1.ancestral.species$random.SS)
A1.ancestral.rand.sd <- sd(A1.ancestral.species$random.SS)
A1.ancestral.species #significant

#####Z-tests
#can be compared to Z distribution. See Adams & Collyer, 2016 Evolution. "On the comparison of the strengthof morphological integration acrossmorphometric datasets "
A1A2D5.z.test <- (((two.z-two.rand.mean)-(A1A2D5.z - A1A2D5.rand.mean)) / sqrt((two.rand.sd)^2 + (A1A2D5.rand.sd)^2))
A1.ancestral.z.test <- (((two.z-two.rand.mean)-(A1.ancestral.z - A1.ancestral.rand.mean)) / sqrt((two.rand.sd)^2 + (A1.ancestral.rand.sd)^2))


######results of z-test
A1A2D5.z.test
#p-value of x-score for 0.48 is .6844; for z-score of 0.49 is .6879
A1.ancestral.z.test
#p-value for Z-score of 0.57 is 0.7175; p-value for a-score of 0.58 is 0.7190 

