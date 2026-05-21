#Merge three Pt-based DEG tables into a single shared DEG dataset
#1. Import complete DEGs tables obtained from R-ODAF Shiny app
cis <- read.delim("data/cisplatin.txt")
car <- read.delim ("data/carboplatin.txt")
ox <- read.delim("data/oxaliplatin.txt")


#2. Merge three DEG tables by shared gene identifiers
merge1 <- merge(cis, car, by = "GeneID")
merge2 <- merge(merge1, ox, by = "GeneID")

#3. Change names of columns
colnames(merge2) <- c("GeneID", "GeneSymbol.cis", "baseMean.cis", "Log2FC.cis", "IfcSE.cis", "stat.cis", "pvalue.cis", "padj.cis", "GeneSymbol.car", 
                      "baseMean.car", "Log2FC.car", "IfcSE.car", "stat.car", "pvalue.car", "padj.car", "GeneSymbol.ox", "baseMean.ox", "Log2FC.ox", 
                      "IfcSE.ox", "stat.ox", "pvalue.ox", "padj.ox")

#4. Retrieve number of significant DEGs
sig1 <- merge2[merge2$padj.cis < 0.01, ]
sig2 <- sig1[sig1$padj.car < 0.01, ]
sig3 <- sig2[sig2$padj.ox < 0.01, ]

#5. Export table to file
write.table(merge2, file = "~/stage/R-Studio/merge.txt", sep = ",", row.names = FALSE, quote = FALSE)
