#PCA
#Import raw datafile and change to data frame
raw <- read.delim2("~/stage/wetransfer_platinum-kidney-data_2026-03-05_1053/sampleData_R-ODAF_shiny_cisplatin.txt")
raw_df <- as.data.frame(raw)


##Cisplatin
#Make seperate data frame for cisplatin-treated cells
raw_cis <- raw_df[, c("Cisplatin_5um_48h_R1", "Cisplatin_5um_72h_R3", "Cisplatin_5um_48h_R3", "Cisplatin_5um_24h_R3", "Cisplatin_5um_48h_R2", "Cisplatin_5um_24h_R2",
                         "Cisplatin_5um_72h_R1", "Cisplatin_5um_72h_R2")]

#Make data frame numeric
raw_cis[raw_cish] <- lapply(raw_cis[raw_cish], as.numeric)

#Make new/empty data fram
raw_cish <- c("Cisplatin_5um_48h_R1", "Cisplatin_5um_72h_R3", "Cisplatin_5um_48h_R3", "Cisplatin_5um_24h_R3", "Cisplatin_5um_48h_R2", "Cisplatin_5um_24h_R2",
              "Cisplatin_5um_72h_R1", "Cisplatin_5um_72h_R2")


#Make new column average over all samples
raw_cis$mean_cis <- rowMeans(raw_cis[raw_cish], na.rm = TRUE)




##Carboplatin
#Make seperate data frames for carboplatin-treated cells
raw_car <- raw_df[, c("Carboplatin_60um_48h_R1", "Carboplatin_60um_72h_R3", "Carboplatin_60um_48h_R3", "Carboplatin_60um_24h_R3", "Carboplatin_60um_48h_R2", "Carboplatin_60um_24h_R2", "Carboplatin_60um_72h_R1", "Carboplatin_60um_72h_R2", "Carboplatin_60um_24h_R1")]
raw_carh <- c("Carboplatin_60um_48h_R1", "Carboplatin_60um_72h_R3", "Carboplatin_60um_48h_R3", "Carboplatin_60um_24h_R3", "Carboplatin_60um_48h_R2", "Carboplatin_60um_24h_R2",
              "Carboplatin_60um_72h_R1", "Carboplatin_60um_72h_R2", "Carboplatin_60um_24h_R1")
#Make numeric
raw_car[raw_carh] <-lapply(raw_car[raw_carh], as.numeric)

#Make new column average over all samples
raw_car$mean_car <- rowMeans(raw_car[raw_carh], na.rm = TRUE)


#Oxaliplatin
#Make seperate data frames for oxaliplatin-treated cells
raw_ox <- raw_df[, c("Oxaliplatin_15um_48h_R1", "Oxaliplatin_15um_72h_R3", "Oxaliplatin_15um_48h_R3", "Oxaliplatin_15um_24h_R3", "Oxaliplatin_15um_48h_R2", "Oxaliplatin_15um_24h_R2",
                     "Oxaliplatin_15um_72h_R1", "Oxaliplatin_15um_72h_R2", "Oxaliplatin_15um_24h_R1")]
raw_oxh <- c("Oxaliplatin_15um_48h_R1", "Oxaliplatin_15um_72h_R3", "Oxaliplatin_15um_48h_R3", "Oxaliplatin_15um_24h_R3", "Oxaliplatin_15um_48h_R2", "Oxaliplatin_15um_24h_R2",
             "Oxaliplatin_15um_72h_R1", "Oxaliplatin_15um_72h_R2", "Oxaliplatin_15um_24h_R1")

#Make data frame numeric
raw_ox[raw_oxh] <- lapply(raw_ox[raw_oxh], as.numeric)

#Make new column average over all samples
raw_ox$mean_ox <- rowMeans(raw_ox[raw_oxh], na.rm = TRUE)


#Make new data frame with mean columns of three Pt-based compounds
df_mean <- data.frame(
  mean_cis = raw_cis$mean_cis,
  mean_car =raw_car$mean_car,
  mean_ox = raw_ox$mean_ox
)


#Switch rows and columns
df_mean1 <- as.data.frame(t(df_mean))


#Check class
class(df_mean1)


#Transform data frame to matrix
mat_mean <- as.matrix(df_mean1)


#Remove the 0 values
mat_mean_clean <- mat_mean[, apply(mat_mean, 2, function(x) sd(x, na.rm = TRUE) > 0)]


#Change column names
rownames(mat_mean_clean) <- c("cisplatin", "carboplatin", "oxaliplatin")


#pca  
pca <- prcomp(mat_mean_clean, scale. = TRUE)
pc_eigenvalues <- (pca$sdev^2) / (sum(pca$sdev^2)) 
per_var <- pc_eigenvalues * 100


#pca plot
plot(pca$x[,1], pca$x[,2],
     xlab = "PC1 (58.1%)",
     ylab = "PC2 (41.9%)",
     pch = 19,
     cex = 2,
     main = "PCA of mean sample values per Pt-based compound")
text(pca$x[,1], pca$x[,2],
     labels = rownames(mat_mean_clean),
     cex = 0.8,
     pos = 3)
abline(h=0, lty = "dotted", col = "gray", lwd = 2)
abline(v=0, lty = "dotted", col = "gray", lwd = 2)
