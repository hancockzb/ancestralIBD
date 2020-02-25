library(ggplot2)
library(ape)
library(bio3d)
#Randomly sampled SNPs
snps<-read.dna("ts_45_mig001_1.fas",format="fasta",as.character=T)
snps_subset <- snps[, sample(ncol(snps), 3000)]
snps_sub <- as.fasta(snps_subset)
write.fasta(snps_sub, file = "mig001_45_1.fas")

#Merging the different TMRCA simulations from fastsimcoal2
###Example, does not include all files
PA_1m_Center <- read.table("PA_1m_Center.txt", header=TRUE, sep=",", check.names = FALSE)
PA_1m_Center <- data.frame(PA_1m_Center[1:12], stack(PA_1m_Center[1:ncol(PA_1m_Center)]))
PA_1m_Center <- subset(PA_1m_Center, select = c(values, ind))
PA_1m_Center$DemeSampled <- "Center"
PA_1m_Ends <- read.table("PA_1m_Ends.txt", header=TRUE, sep=",",check.names = FALSE)
PA_1m_Ends <- data.frame(PA_1m_Ends[1:12], stack(PA_1m_Ends[1:ncol(PA_1m_Ends)]))
PA_1m_Ends <- subset(PA_1m_Ends, select = c(values, ind))
PA_1m_Ends$DemeSampled <- "Ends"
Island_1m <- read.table("Island_1m.txt", header=TRUE, sep=",", check.names = FALSE)
Island_1m <- data.frame(Island_1m[1:12], stack(Island_1m[1:ncol(Island_1m)]))
Island_1m <- subset(Island_1m, select = c(values, ind))
Island_1m$DemeSampled <- "Island"
PA_1m_all <- rbind(PA_1m_Center, PA_1m_Ends, Island_1m)
PA_01m_all <- rbind(PA_1m_Center, PA_1m_Ends, Island_1m)
PA_001m_all <- rbind(PA_1m_Center, PA_1m_Ends, Island_1m)
PA_1m_all$MigrationRate <- "0.1"
PA_01m_all$MigrationRate <- "0.01"
PA_001m_all$MigrationRate <- "0.001"
All <- rbind(PA_1m_all, PA_01m_all, PA_001m_all)

#Pairwise divergence times
div_p1p10 <- read.table("div_09_mig001_1.txt", header=TRUE, sep=",")
div_p5p6 <- read.table("div_45_mig001_1.txt", header = TRUE, sep = ",")
div_p1p10$Deme_sampled <- "End"
div_p5p6$Deme_sampled <- "Center"
div_mig1 <- rbind(div_p1p10, div_p5p6)
div_mig1$Migration_rate <- "0.1"
div_mig01 <- rbind(div_p1p10, div_p5p6)
div_mig01$Migration_rate <- "0.01"
div_mig001 <- rbind(div_p1p10, div_p5p6)
div_mig001$Migration_rate <- "0.001"
div_50 <- rbind(div_mig1, div_mig01, div_mig001)
div_50$Ratio <- "50"
div_50$Divergence_time <- div_50$Pairwise_divergence/(2*1e-7)
div_50$TMRCA_TD <- (div_50$Divergence_time - 50000)/50000
div_50$Ne <- (div_50$Pairwise_divergence - (2*50000*1e-7))/(4*1e-7)
div_25 <- rbind(div_mig1, div_mig01, div_mig001)
div_25$Ratio <- "25"
div_25$Divergence_time <- div_25$Pairwise_divergence/(2*1e-7)
div_25$TMRCA_TD <- (div_25$Divergence_time - 25000)/25000
div_25$Ne <- (div_25$Pairwise_divergence - (2*25000*1e-7))/(4*1e-7)
div_10 <- rbind(div_mig1, div_mig01, div_mig001)
div_10$Ratio <- "10"
div_10$Divergence_time <- div_10$Pairwise_divergence/(2*1e-7)
div_10$TMRCA_TD <- (div_10$Divergence_time - 10000)/10000
div_10$Ne <- (div_10$Pairwise_divergence - (2*10000*1e-7))/(4*1e-7)
div_5 <- rbind(div_mig1, div_mig01, div_mig001)
div_5$Ratio <- "5"
div_5$Divergence_time <- div_5$Pairwise_divergence/(2*1e-7)
div_5$TMRCA_TD <- (div_5$Divergence_time - 5000)/5000
div_5$Ne <- (div_5$Pairwise_divergence - (2*5000*1e-7))/(4*1e-7)
div_1 <- rbind(div_mig1, div_mig01, div_mig001)
div_1$Ratio <- "1"
div_1$Divergence_time <- div_1$Pairwise_divergence/(2*1e-7)
div_1$TMRCA_TD <- (div_1$Divergence_time - 1000)/1000
div_1$Ne <- (div_1$Pairwise_divergence - (2*1000*1e-7))/(4*1e-7)
div_all <- rbind(div_50, div_25, div_10, div_5, div_1)
div_all$Ratio <- factor(div_all$Ratio, levels=c(50, 25, 10, 5, 1), labels = c("50", "25", "10", "5", "1"))

#Load trace files into R, subset, analyze
log_09_50 <- read.table("09_1.txt", header=TRUE, sep=",")
data1 <- subset(log_09_50, select = -c(state1, state2, state3))
data2 <- data.frame(data1[1:3], stack(data1[1:ncol(data1)]))
data_09_50 <- subset(data2, select = -c(X0.1, X0.01, X0.001))
data_09_50$DemeSampled <- "Ends"
data_09_50$TD_ND <- "1"
log_45_50 <- read.table("45_1.txt", header=TRUE, sep=",")
data_1 <- subset(log_45_50, select = -c(state1, state2, state3))
data_2 <- data.frame(data_1[1:3], stack(data_1[1:ncol(data_1)]))
data_45_50 <- subset(data_2, select = -c(X0.1, X0.01, X0.001))
data_45_50$DemeSampled <- "Center"
data_45_50$TD_ND <- "1"
data_50 <- rbind(data_09_50, data_45_50)
data_25<- rbind(data_09_50, data_45_50)
data_10 <- rbind(data_09_50, data_45_50)
data_5 <- rbind(data_09_50, data_45_50)
data_1 <- rbind(data_09_50, data_45_50)
data_50$values <- data_50$values/1e-7
data_50$TMRCA_TD <- (data_50$values - 50000)/50000
data_25$values <- data_25$values/1e-7
data_25$TMRCA_TD <- (data_25$values - 25000)/25000
data_10$values <- data_10$values/1e-7
data_10$TMRCA_TD <- (data_10$values - 10000)/10000
data_5$values <- data_5$values/1e-7
data_5$TMRCA_TD <- (data_5$values - 5000)/5000
data_1$values <- data_1$values/1e-7
data_1$TMRCA_TD <- (data_1$values - 1000)/1000
data_all <- rbind(data_50, data_25, data_10, data_5, data_1)
data_all$TD_ND <- factor(data_all$TD_ND, levels=c(50, 25, 10, 5, 1), labels = c("50", "25", "10", "5", "1"))

#Visualize the data
ggplot(data_all,aes(x = ind, y = TMRCA_TD, fill = DemeSampled)) + geom_violin() + stat_summary(fun.y = median, geom = "point", size = 2, shape = 23, position = position_dodge(width = 0.9)) + facet_grid("TD_ND", scales = "free") + geom_hline(yintercept = 0, linetype = "dashed") + xlab("Migration rate") + ylab("(TMRCA-TD)/TD")
ggplot(data_50) + geom_density(aes( x = TMRCA_TD, fill = DemeSampled), alpha = 0.5) + facet_grid("ind", scales = "free") + xlim(-1, 1) + geom_vline(xintercept = 0, linetype = "dashed")
ggplot(data_25) + geom_density(aes( x = TMRCA_TD, fill = DemeSampled), alpha = 0.5) + facet_grid("ind", scales = "free") + xlim(-0.5, 8) + geom_vline(xintercept = 0, linetype = "dashed")
ggplot(data_10) + geom_density(aes( x = TMRCA_TD, fill = DemeSampled), alpha = 0.5) + facet_grid("ind", scales = "free") + xlim(-0.5, 10) + geom_vline(xintercept = 0, linetype = "dashed")
ggplot(data_5) + geom_density(aes( x = TMRCA_TD, fill = DemeSampled), alpha = 0.5) + facet_grid("ind", scales = "free") + xlim(-0.5, 15) + geom_vline(xintercept = 0, linetype = "dashed")
ggplot(data_1) + geom_density(aes( x = TMRCA_TD, fill = DemeSampled), alpha = 0.5) + facet_grid("ind", scales = "free") + xlim(0, 50) + geom_vline(xintercept = 0, linetype = "dashed")

#ANOVAs
model_50 <- aov(TMRCA_TD~DemeSampled*ind, data = data_50)
model_25 <- aov(TMRCA_TD~DemeSampled*ind, data = data_25)
model_10 <- aov(TMRCA_TD~DemeSampled*ind, data = data_10)
model_5 <- aov(TMRCA_TD~DemeSampled*ind, data = data_5)
model_1 <- aov(TMRCA_TD~DemeSampled*ind, data = data_1)
