library(lazyeval)
library(ggplot2)
#Package we used to replicate PCA (CAN NOT PRESENT THIS. CONFLICTS WITH PUBLICATION) 
library(factoextra)

#############################################
sessionInfo()

##############################################

# import file from comparative analysis table
cluster_table <- read.table("table_clusters.txt", header = T)
cluster_table <- as.data.frame(cluster_table)

#sort table based on total number of clusters across species
#cluster_table$sum <- rowSums(cluster_table)
#cluster_table <- cluster_table[order(-cluster_table$sum),]
#cluster_table <- subset(cluster_table, select = -c(sum))
cluster_table <- head(cluster_table,1000)

#creating table of A clusters (check final table, there was an A2 not in line with the rest of "A"s)
Acluster_table <- cluster_table[,1:8]
#creating table of D clusters
Dcluster_table <- cluster_table[,9:13]

################### no normalization PCA ###################
#Use
#transpose data, accessions now row names
cluster_table <- t(cluster_table)

#Variable that runs principle component analysis on cluster_table
cluster.pca <- prcomp(cluster_table, scale = TRUE, na.rm = TRUE)
cluster.eig <- get_eigenvalue(cluster.pca)

# barplot(cluster.eig[, 2], names.arg=1:nrow(cluster.eig), main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="steelblue")

# fviz_screeplot(cluster.pca, ncp=10, choice="eigenvalue")

#X is an item within cluster.pc from previous command
ind.coord <- cluster.pca$x


#subsetting data, taking first 1000 columns
subsections <- as.data.frame(cluster_table[,1:2])
subsections$sub <-c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5")


subsections <- subsections[,-c(1:2)]
subfac <- as.factor(subsections)

png("cotton_outgroup.PCA.direct.annot.new.png", 1000, 1000, pointsize=20)
fviz_pca_ind(cluster.pca, habillage=subfac) + theme_minimal()
dev.off()


#Recreate Table 2
install.packages("reshape2")
library(reshape2)
new_table = cluster_table
row.names(new_table) <- c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5")
melted <- melt(new_table)
Y <- melted$value
Species <- as.factor(melted$Var1)
Cluster <- as.factor(melted$Var2)
anova(lm(Y~Cluster*Species))

#Recreate Figure 4 #still needs work 
new_table <- as.data.frame(t(new_table))
linear_model <- glm(A1~A2, data=new_table)
ggplot(new_table, aes(x=A1, y=A2)) + geom_point(size=2) + geom_abline(intercept=0, slope=1)
