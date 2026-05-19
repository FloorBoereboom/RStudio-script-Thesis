#Merge three Pt-based DEG tables into a single shared DEG dataset
#Import complete DEGs tables obtained from R-ODAF Shiny app
cis <- read.delim("~/stage/wetransfer_platinum-kidney-data_2026-03-05_1053/DS1_Cisplatin_FULL_DESeq2_results_FDR0.01_20260324_105312.txt")
car <- read.delim ("~/stage/wetransfer_platinum-kidney-data_2026-03-05_1053/DS1_Carboplatin_FULL_DESeq2_results_FDR0.01_20260324_111113.txt")
ox <- read.delim("~/stage/wetransfer_platinum-kidney-data_2026-03-05_1053/DS1_Oxaliplatin_FULL_DESeq2_results_FDR0.01_20260324_112033.txt")


#Merge three DEG tables by shared gene identifiers
merge1 <- merge(cis, car, by = "GeneID")
merge2 <- merge(merge1, ox, by = "GeneID")


#Retrieve number of significant DEGs
sig1 <- merge2[merge2$padj.x < 0.01, ]
sig2 <- sig1[sig1$padj.y < 0.01, ]
sig3 <- sig2[sig2$padj < 0.01, ]


#Change names of columns
colnames(merge2) <- c("GeneID", "GeneSymbol.cis", "baseMean.cis", "Log2FC.cis", "IfcSE.cis", "stat.cis", "pvalue.cis", "padj.cis", "GeneSymbol.car", "baseMean.car", "Log2FC.car", "IfcSE.car", "stat.car", "pvalue.car", "padj.car", "GeneSymbol.ox", "baseMean.ox", "Log2FC.ox", "IfcSE.ox", "stat.ox", "pvalue.ox", "padj.ox")


#Export table 
write.table(merge2, file = "~/stage/R-Studio/merge.txt", sep = ",", row.names = FALSE, quote = FALSE)
