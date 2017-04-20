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

<<<<<<< HEAD
#Recreate Figure 4 #still needs work 
new_table <- as.data.frame(t(new_table))
linear_model <- glm(A1~A2, data=new_table)
ggplot(new_table, aes(x=A1, y=A2)) + geom_point(size=2) + geom_abline(intercept=0, slope=1)
=======
#Recreate Figure 3 #still needs work 
new_table <- as.data.frame(t(new_table))
linear_model <- glm(A1~A2, data=new_table)
ggplot(new_table, aes(x=A1, y=A2)) + geom_point(size=2) + geom_abline(intercept=0, slope=1)

#Recreate Figure 4, not pretty

# import file from comparative analysis table and take top 1000 clusters, as per SRB paper
clust <- read.table("table_clusters.txt", header = T, sep="\t")
clust <- head(clust, 1000)

#ordination analysis

#transform data and make the euclidean distance matrix
mydata <- t(clust)
d <- dist(mydata, method = "euclidean")

#SRB performs classical MDS (principal coordinate analysis) using cmdscale, minimally
#he mentions this package specifically, which I why I don't know where vegan comes in
cmdfit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]

png("SRB.Fig4.png", 5000, 5000, pointsize=12, res=600)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(x, y, labels = row.names(mydata), cex=.7)
points(x, y, pch=19)
dev.off()

# Simon also claimed to use vegan/ordihull/ordispider to do this analysis
# ordihull/spider appear to rely on vegan output
# there are not enough details to determine if he used vegan for the ordination analysis
# and how he would have conducted it: community dissimilarities, eigenvector methods, detrended correspondence analysis, etc
# not to mention there are tons of other options
# AFAIK ordihull/spider appears to make the output of vegan more pretty

########### characterize composition ###########
##Corrinne's code for Corrinne's file
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

##### BEGIN EMMA'S CODE FOR FIGURE 2
#Read averaged cluster counts per category
KBsum <- read.csv("f.csv", header = TRUE)

#Load reshape in order to melt data frame and create min/max columns
library(reshape2)
min <- c(KBsum$A1min, KBsum$A2min,KBsum$D5min)
max <- c(KBsum$A1max, KBsum$A2max, KBsum$D5max)
KBm <- melt(KBsum)
KBm$min <- min
KBm$max <- max

#Load appropriate packages for visualzing data
library(ggplot2)
library(scales)
library(gridExtra)
#Turn min/max into limits for range bars
limits <- aes(ymax=KBm$max, ymin=KBm$min)

#Keep width of bars the same even though data size will change
dodge <- position_dodge(width=0.9)

#Create space for image 
png("Figure_TE.amounts.png", 1000, 1000, pointsize=20)
#List of colors to use in graph
cols <- c("red","blue","green","yellow","orange","orchid","orange","darkgreen","orchid4")
#List of parameters for image
ggplot(KBm, aes(x=X, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + scale_y_log10() + scale_fill_manual(breaks=c("A1", "A2", "D5"), values = cols) + geom_errorbar(limits, position = dodge) 




#### kirkii vs kokia clusters ####

png("Figure_TE.comparisons.png", 1000, 1000, pointsize=20)

p1 <- ggplot(annot_clust, aes(x=, y=kirkii, shape=signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1) + scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3")) +  scale_x_continuous(expand = c(0, 0), limits=c(0,750)) + scale_y_continuous(expand = c(0, 0), limits=c(0,750))
p2 <- ggplot(annot_clust[annot_clust$Lineage == "LTR", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1)+ scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,400)) + scale_y_continuous(expand = c(0, 0), limits=c(0,400))
p3 <- ggplot(annot_clust[annot_clust$Lineage == "LTR/Gypsy", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1)+ scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,750)) + scale_y_continuous(expand = c(0, 0), limits=c(0,750))
p4 <- ggplot(annot_clust[annot_clust$Lineage == "LTR/Copia", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, size=1)+ scale_color_manual(breaks=c("positive", "negative"), values=c("blue3", "green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,150)) + scale_y_continuous(expand = c(0, 0), limits=c(0,150))

grid.arrange(p1,p2,p3,p4, ncol=2)

dev.off()












>>>>>>> e96254db25f62eb526e9e086b783720381db962b
