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
#Acluster_table <- cluster_table[,1:8]
#creating table of D clusters
#Dcluster_table <- cluster_table[,9:13]

################### no normalization PCA ###################
#Use
#transpose data, accessions now row names
cluster_table <- t(cluster_table)

#Variable that runs principle component analysis on cluster_table
cluster.pca <- prcomp(cluster_table, scale = TRUE, na.rm = TRUE)
#cluster.eig <- get_eigenvalue(cluster.pca)

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

library(reshape2)
new_table = cluster_table
row.names(new_table) <- c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "A2", "D5")
melted <- melt(new_table)
Y <- melted$value
Species <- as.factor(melted$Var1)
Cluster <- as.factor(melted$Var2)
anova(lm(Y~Cluster*Species))

#Recreate Figure 2
library(ggplot2)
library(reshape2)
#Read in data frame
stable <- read.table("f.csv", sep = ",", header =TRUE)
#Create new rows that are averages of each species
stable$A1 <- rowMeans(stable[,2:4])
stable$A2 <- rowMeans(stable[,5:9])
stable$D5 <- rowMeans(stable[,10:14])
#Create standard deviation function
std <- function(x) sd(x)/sqrt(length(x))
#Find standard deviation of each species
stable$A1err <- apply(stable[,2:4], 1, std)
stable$A2err <- apply(stable[,5:9], 1, std)
stable$D5err <- apply(stable[,10:14], 1, std)
#Create Standard Deviation minimum
stable$A1min <- stable$A1 - stable$A1err
stable$A2min <- stable$A2 - stable$A2err
stable$D5min <- stable$D5 - stable$D5err
#And Standard Deviation maximum
stable$A1max <- stable$A1 + stable$A1err
stable$A2max <- stable$A2 + stable$A2err
stable$D5max <- stable$D5 + stable$D5err
#Remove categorizations of TE that we don't care about
stable <- stable[,-(2:14)]
#Melt dataframe in order to plot properly
melted <- melt(stable[,(1:4)])
min <- c(stable$A1min, stable$A2min, stable$D5min)
max <- c(stable$A1max, stable$A2max, stable$D5max)
melted$min <- min
melted$max <- max
#Create standard error bars
limits <- aes(ymax=melted$max, ymin=melted$min)
#Maintain width of bars
dodge <- position_dodge(width=0.9)
#Create image space
png("Figure_TE.amounts.png", 7500, 5000, pointsize=12, res=600)
#Plot within image space
ggplot(melted, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + geom_errorbar(limits, position = dodge, width=0.3) + labs(y="Mbp/1C", x="") + theme_set(theme_grey(base_size=12)) + scale_fill_hue(l=40)
dev.off()

#Figure 3
library(ggplot2)
new_table <- as.data.frame(t(cluster_table))
new_table$A1ave <- rowMeans(subset(new_table, select = c(A1_155_, A1_073_, A1_097_), na.rm = TRUE))
new_table$A2ave <- rowMeans(subset(new_table, select = c(A2_255_, A2_034_, A2_044_, A2_099_, A2_101_), na.rm = TRUE))
ggplot(A1ave~A2ave, data=new_table, mapping = aes(x = A1ave, y = A2ave)) + geom_point() + geom_abline(intercept=0, slope=1)


png("cotton.ordination.png", 5000, 5000, pointsize=12, res=600)
ggplot(cmdpoints, aes(x=V1, y=V2, color=genome, shape=genome)) + geom_point(size=2) + xlab("PCoA component 1") + ylab("PCoA component 2") + scale_color_manual(values=c("A1"="orchid", "A2"="orchid4", "D5"="slategrey","kirkii"="blue3", "kokia"="green3"))+ geom_text_repel(aes(label=species))### Recreate Figure 4 Simon's way
mydata <- t(annot_clust[,-c(1)])
d <- dist(mydata, method = "euclidean")

cmdfit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim

x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric        MDS",        type="n")
text(x, y, labels = row.names(mydata), cex=.7)


#Recreate Figure 4, not pretty ##Need to fix to recreate actual figure, not the way it should have been done 

#ordination analysis

#make the euclidean distance matrix
d <- dist(cluster_table, method = "euclidean")

#SRB performs classical MDS (principal coordinate analysis) using cmdscale, minimally
#he mentions this package specifically, which I why I don't know where vegan comes in
cmdfit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]

png("SRB.Fig4.png", 5000, 5000, pointsize=12, res=600)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(x, y, labels = row.names(cluster_table), cex=.7)
points(x, y, pch=19)
dev.off()

# Simon also claimed to use vegan/ordihull/ordispider to do this analysis
# ordihull/spider appear to rely on vegan output
# there are not enough details to determine if he used vegan for the ordination analysis
# and how he would have conducted it: community dissimilarities, eigenvector methods, detrended correspondence analysis, etc
# not to mention there are tons of other options
# AFAIK ordihull/spider appears to make the output of vegan more pretty

