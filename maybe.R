setwd("~/Box Sync/Classes/BioDataSkills_BCB546X_Spr17/EEOB546xFINALPROJECT")
library(geomorph)
library(reshape2)
library(ggplot2)
library(contrast)
#input and structure data
cluster_matrix <- t(head(as.matrix(read.table("table_clusters.txt", header = T)),1000))

#Table 2 (two-way ANOVA test, as done in the paper)
new_table <- cluster_matrix
row.names(new_table) <- c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5")
melted <- melt(new_table)
Y <- melted$value
Species <- as.factor(melted$Var1)
Cluster <- as.factor(melted$Var2)
linear_model <- lm(Y~Cluster*Species)
anova(linear_model) #same as in paper




###Figure 3 (with colors, but no BH correction) 
#differences between this figure and Simon's figure is highlighted on Simon_Fig3_discrepency_noBH.png

cluster_table <- as.data.frame(t(cluster_matrix))
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


####Table 3 (not same numbers as paper, but same methods 

#without BH correction
A1A2.clusters <- sum(cluster_table$A1A2.significant != 0)
A1D5.clusters <- sum(cluster_table$A1D5.significant != 0)
A2D5.clusters <- sum(cluster_table$A2D5.significant != 0)
A1A2.clusters
A1D5.clusters
A2D5.clusters

#with BH correction
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
A1A2D5 <- as.factor(c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5"))
AD <- as.factor(c("A", "A", "A", "A", "A", "A", "A", "D", "D", "D", "D", "A", "D"))
Awild <- as.factor(c("A", "Awild", "A", "A", "A", "A", "A", "D", "D", "D", "D", "A", "D"))
PermManova.3species <- advanced.procD.lm(cluster_matrix~A1A2D5, ~ 1, groups = ~A1A2D5)
PermManova.2species <- advanced.procD.lm(cluster_matrix~AD, ~ 1, groups = ~AD)
PermManova.commonwild <- advanced.procD.lm(cluster_matrix~Awild, ~ 1, groups = ~Awild)
PermManova.3species$P.means.dist
PermManova.2species$P.means.dist
PermManova.commonwild$P.means.dist


