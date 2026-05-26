#PCA
#1. Import normalized transcriptome count datafile and change to data frame
raw <- read.delim2("data/count_matrix.txt")
raw_df <- as.data.frame(raw)


#Cisplatin
#2. Make seperate data frame for cisplatin-treated cells
raw_cis <- raw_df[, c("Cisplatin_5um_48h_R1", "Cisplatin_5um_72h_R3", "Cisplatin_5um_48h_R3", "Cisplatin_5um_24h_R3", "Cisplatin_5um_48h_R2", 
                      "Cisplatin_5um_24h_R2", "Cisplatin_5um_72h_R1", "Cisplatin_5um_72h_R2")]

#3. Make data frame numeric
raw_cis[] <- lapply(raw_cis, as.numeric)

#4. Make new column taking the average over all replicates and time points
raw_cis$mean_cis <- rowMeans(raw_cis, na.rm = TRUE)


#Carboplatin
#5. Make seperate data frame for carboplatin-treated cells
raw_car <- raw_df[, c("Carboplatin_60um_48h_R1", "Carboplatin_60um_72h_R3", "Carboplatin_60um_48h_R3", "Carboplatin_60um_24h_R3", "Carboplatin_60um_48h_R2", 
                      "Carboplatin_60um_24h_R2", "Carboplatin_60um_72h_R1", "Carboplatin_60um_72h_R2", "Carboplatin_60um_24h_R1")]

#6. Make data frame numeric
raw_car[] <-lapply(raw_car, as.numeric)

#7. Make new column taking the average over all replicates and time points
raw_car$mean_car <- rowMeans(raw_car, na.rm = TRUE)


#Oxaliplatin
#8. Make seperate data frame for oxaliplatin-treated cells
raw_ox <- raw_df[, c("Oxaliplatin_15um_48h_R1", "Oxaliplatin_15um_72h_R3", "Oxaliplatin_15um_48h_R3", "Oxaliplatin_15um_24h_R3", "Oxaliplatin_15um_48h_R2", 
                     "Oxaliplatin_15um_24h_R2", "Oxaliplatin_15um_72h_R1", "Oxaliplatin_15um_72h_R2", "Oxaliplatin_15um_24h_R1")]

#9. Make data frame numeric
raw_ox[] <- lapply(raw_ox, as.numeric)

#10. Make new column taking the average over all replicates and time points
raw_ox$mean_ox <- rowMeans(raw_ox, na.rm = TRUE)


#11. Make new data frame with average columns of three Pt-based compounds
df_mean <- data.frame(
  mean_cis = raw_cis$mean_cis,
  mean_car =raw_car$mean_car,
  mean_ox = raw_ox$mean_ox
)


#12. Switch rows for columns as samples should be the columns and genes rows
df_mean1 <- as.data.frame(t(df_mean))


#13. Transform data frame to matrix because pca can only be calculated for a matrix 
mat_mean <- as.matrix(df_mean1)


#14. Remove the zero values
mat_mean_clean <- mat_mean[, apply(mat_mean, 2, function(x) sd(x, na.rm = TRUE) > 0)]


#15. Change sample names
rownames(mat_mean_clean) <- c("cisplatin", "carboplatin", "oxaliplatin")


#16. Perform pca and calculate eigenvalues
pca <- prcomp(mat_mean_clean, scale. = TRUE)
pc_eigenvalues <- (pca$sdev^2) / (sum(pca$sdev^2)) 
per_var <- pc_eigenvalues * 100


#17. Create biplot of pca
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
