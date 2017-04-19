setwd("~/Box Sync/Classes/BioDataSkills_BCB546X_Spr17/EEOB546xFINALPROJECT")
library(geomorph)
library(reshape2)
library(ggplot2)
library(contrast)
#input and structure data
cluster_table <- t(head(as.matrix(read.table("table_clusters.txt", header = T)),1000))

#two-way anova test (as done in the paper)
new_table = cluster_table
row.names(new_table) <- c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5")
melted <- melt(new_table)
Y <- melted$value
Species <- as.factor(melted$Var1)
Cluster <- as.factor(melted$Var2)
linear_model <- lm(Y~Cluster*Species)
anova(linear_model) #same as in paper



#Figure 4 (and maybe 3)
new_table <- as.data.frame(t(cluster_table))
new_table$A1ave <- rowMeans(subset(new_table, select = c(A1_155_, A1_073_, A1_097_), na.rm = TRUE))
new_table$A2ave <- rowMeans(subset(new_table, select = c(A2_255_, A2_034_, A2_044_, A2_099_, A2_101_), na.rm = TRUE))

A.averages <- subset(new_table, select = c(A1ave, A2ave))
A.melted <- melt(t(A.averages))
Y.data <- A.melted$value
A.species <- as.factor(A.melted$Var1)
A.clusters <- as.factor(A.melted$Var2)
GLM_clusters <- lm(Y.data~A.species*A.clusters)

#GLM_clusters <- glm(Y~Cluster*Species)
contrasts <- contrast(GLM_clusters, list(A.clusters = levels(A.clusters), A.species = "A2ave"), list(A.clusters = levels(A.clusters), A.species = "A1ave"))
Pvalue <- contrasts$Pvalue
test <- ifelse((Pvalue < 0.05),1,0)
sum(test)

ggplot(A1ave~A2ave, data=new_table, mapping = aes(x = A1ave, y = A2ave, color=significant)) + geom_point()



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



##results of z-test
A1A2D5.z.test
#p-value of x-score for 0.48 is .6844; for z-score of 0.49 is .6879
A1.ancestral.z.test
#p-value for Z-score of 0.57 is 0.7175; p-value for a-score of 0.58 is 0.7190 
