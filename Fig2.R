##################################################################
put the below in the clipboard, or load in the file :)

Lineage	A1_155	A1_073	A1_097	A2_255	A2_034	A2_044	A2_099	A2_101	D5_002	D5_031	D5_004	D5_053	D5_ggg
RLC	36.6985	36.3945	34.732	38.304	38.342	36.2235	37.5535	36.2995	34.827	34.0575	38.342	35.853	35.055
RLG	458.204	451.6395	477.033	553.983	528.8365	497.762	482.8185	544.958	109.003	107.7965	121.9325	103.5405	123.424
RLX	500.8685	503.766	505.4855	524.4095	526.8415	514.2825	518.206	573.705	184.9935	169.0905	187.093	158.194	165.0625
RXX	5.0635	5.434	4.237	4.541	4.883	4.8355	4.8355	4.921	3.116	3.1635	3.0115	3.325	3.097
TXX	57.057	41.99	61.807	47.158	45.714	60.23	51.623	45.486	49.951	38.76	35.34	28.7755	29.8015


##################################################################

library(ggplot2)
library(reshape2)
stable <- read.table("clipboard", sep = "\t", header =TRUE)

stable$A1 <- rowMeans(stable[,2:4])
stable$A2 <- rowMeans(stable[,5:9])
stable$D5 <- rowMeans(stable[,10:14])

std <- function(x) sd(x)/sqrt(length(x))

stable$A1err <- apply(stable[,2:4], 1, std)
stable$A2err <- apply(stable[,5:9], 1, std)
stable$D5err <- apply(stable[,10:14], 1, std)

stable$A1min <- stable$A1 - stable$A1err
stable$A2min <- stable$A2 - stable$A2err
stable$D5min <- stable$D5 - stable$D5err

stable$A1max <- stable$A1 + stable$A1err
stable$A2max <- stable$A2 + stable$A2err
stable$D5max <- stable$D5 + stable$D5err

stable <- stable[,-(2:14)]

melted <- melt(stable[,(1:4)])
min <- c(stable$A1min, stable$A2min, stable$D5min)
max <- c(stable$A1max, stable$A2max, stable$D5max)
melted$min <- min
melted$max <- max

limits <- aes(ymax=melted$max, ymin=melted$min)

dodge <- position_dodge(width=0.9)

png("Figure_TE.amounts.png", 7500, 5000, pointsize=12, res=600)

ggplot(melted, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + geom_errorbar(limits, position = dodge, width=0.3) + labs(y="Mbp/1C", x="") + theme_set(theme_grey(base_size=12)) + scale_fill_hue(l=40)
dev.off()
