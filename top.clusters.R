library(ggplot2)
#Package Simon used with original data
library(contrast)
#Package Corrinne has used to replicate similar data 
library(factoextra)

#############################################
sessionInfo()

##############################################

# import file from comparative analysis table
cluster_table <- read.csv("test_data.csv", header = T)

#creating table of A clusters (check final table, there was an A2 not in line with the rest of "A"s)
Acluster_table <- cluster_table[,1:8]
#creating table of D clusters
Dcluster_table <- cluster_table[,9:13]

################### no normalization PCA ###################
#transpose data, accessions now row names
cluster_table <- t(cluster_table)

cluster.pca <- prcomp(cluster_table, scale = TRUE)
cluster.eig <- get_eigenvalue(cluster.pca)

# barplot(cluster.eig[, 2], names.arg=1:nrow(cluster.eig), main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="steelblue")

# fviz_screeplot(cluster.pca, ncp=10, choice="eigenvalue")

ind.coord <- cluster.pca$x

subsections <- as.data.frame(cluster_table[,1:2])
subsections$sub <-c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "D5", "kirkii", "kokia")

subsections$CL0001 <- NULL
subsections$CL0002 <- NULL

subfac <- as.factor(subsections[,1])

png("cotton_outgroup.PCA.direct.annot.png", 1000, 1000, pointsize=20)
fviz_pca_ind(cluster.pca, habillage=subfac) + theme_minimal()
dev.off()


########### DESeq PCA ############

library('DESeq2')

cluster_table <- read.table("annotated.counts", header = T, row.names=1, sep="\t")

samples <- data.frame(row.names = colnames(cluster_table), condition = c("A1_155", "A1_073", "A1_097", "A2_255", "A2_034", "A2_044", "A2_099",  "A2_101", "D5_002", "D5_031", "D5_004", "D5_053","D5ggg", "kirkii", "kokia"))

data_count <- DESeqDataSetFromMatrix(countData = cluster_table,colData = samples,design = ~ condition)

results <- DESeq(data_count)

# Plot dispersions, turn off first and third lines to visualize in line
png("OG_dispersion.annot.png", 1000, 1000, pointsize=20)
plotDispEsts(results, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld_results <- rlogTransformation(results)

# Generate PCA

library(devtools)
# devtools::install_github("hoesler/rwantshue")
library('rwantshue')

scheme <- iwanthue()
scheme$hex(5)

png("cotton_outgroup.PCA.DESeq.annot.png", 1000, 1000, pointsize=20)

OGPCA <- plotPCA(rld_results, intgroup="condition", returnData=TRUE)
OGPCA$group <-c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "D5", "kirkii", "kokia")
percentVar <- round(100*attr(OGPCA, "percentVar"))

ggplot(OGPCA, aes(PC1, PC2, color=condition, shape=group)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() + scale_colour_manual(values=scheme$hex(21))

dev.off()

########### characterize composition ###########

annot_clust <- read.table("annotated.counts", header = T, row.names=1, sep="\t")
annot_clust$cluster <- NULL

# 9.5 multiplier represents # reads (x) * 95nt/read * 1 kb/1000nt * 100% = # reads * 9.5 = # Kb in entire genome for that class 
Kbamount <- data.frame(annot_clust[1], apply(annot_clust[2:16], 2, function (x) x*9.5))
KBsum <- aggregate(. ~Lineage, data=Kbamount, FUN=sum)

KBsum$A1 <- rowMeans(KBsum[,2:4])
KBsum$A2 <- rowMeans(KBsum[,5:9])
KBsum$D5 <- rowMeans(KBsum[,10:14])

KBsum$A1min <- apply(KBsum[,2:4], 1, min)
KBsum$A2min <- apply(KBsum[,5:9], 1, min)
KBsum$D5min <- apply(KBsum[,10:14], 1, min)
KBsum$kirkiimin <- KBsum$kirkii
KBsum$kokia_min <- KBsum$kokia_

KBsum$A1max <- apply(KBsum[,2:4], 1, max)
KBsum$A2max <- apply(KBsum[,5:9], 1, max)
KBsum$D5max <- apply(KBsum[,10:14], 1, max)
KBsum$kirkiimax <- KBsum$kirkii
KBsum$kokia_max <- KBsum$kokia_

KBsum <- KBsum[,-(2:14)]
KBm <- melt(KBsum[,-(7:16)])
min <- c(KBsum$kirkiimin, KBsum$kokia_min, KBsum$A1min, KBsum$A2min, KBsum$D5min)
max <- c(KBsum$kirkiimax, KBsum$kokia_max, KBsum$A1max, KBsum$A2max, KBsum$D5max)
KBm$min <- min
KBm$max <- max

limits <- aes(ymax=KBm$max, ymin=KBm$min)

dodge <- position_dodge(width=0.9)

png("Figure_TE.amounts.png", 1000, 1000, pointsize=20)

ggplot(KBm, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + scale_y_log10() + scale_fill_manual(breaks=c("kirkii", "kokia_", "A1", "A2", "D5"), values=c("blue3", "green3", "orchid", "orchid4", "slategrey")) + geom_errorbar(limits, position = dodge)
dev.off()

sum(KBsum$kirkii)/1000
[1] 110.3615
sum(KBsum$kokia_)/1000
[1] 109.4685


#### kirkii vs kokia clusters ####

png("Figure_TE.comparisons.png", 1000, 1000, pointsize=20)

p1 <- ggplot(annot_clust, aes(x=kokia_, y=kirkii, shape=signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1) + scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3")) +  scale_x_continuous(expand = c(0, 0), limits=c(0,750)) + scale_y_continuous(expand = c(0, 0), limits=c(0,750))
p2 <- ggplot(annot_clust[annot_clust$Lineage == "LTR", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1)+ scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,400)) + scale_y_continuous(expand = c(0, 0), limits=c(0,400))
p3 <- ggplot(annot_clust[annot_clust$Lineage == "LTR/Gypsy", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1)+ scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,750)) + scale_y_continuous(expand = c(0, 0), limits=c(0,750))
p4 <- ggplot(annot_clust[annot_clust$Lineage == "LTR/Copia", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1)+ scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,150)) + scale_y_continuous(expand = c(0, 0), limits=c(0,150))

grid.arrange(p1,p2,p3,p4, ncol=2)

dev.off()

