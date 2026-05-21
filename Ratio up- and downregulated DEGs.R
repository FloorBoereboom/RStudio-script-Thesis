#1. Import the three Pt-based compound significant DEGs tables (padj<0.01)
cis <- read.delim("data/significant_DEGs_cisplatin.txt")
car <- read.delim ("data/significant_DEGs_carboplatin.txt")
ox <- read.delim("data/significant_DEGs_oxaliplatin.txt")


# Cisplatin
# 2. Calculate the percentage of significant downregulated DEGs 
cisup <- cis[cis$log2FoldChange > 0,]
cisdown <- cis[cis$log2FoldChange < 0, ]
percentage_cis <- nrow(cisdown) / nrow(cis) * 100


# Carboplatin
#3. Calculate the percentage of significant downregulated DEGs
carup <- car[car$log2FoldChange > 0, ]
cardown <- car[car$log2FoldChange < 0, ]
percentage_car <- nrow(cardown) / nrow(car) * 100

# Oxaliplatin
#4. Calculate the percentage of significant downregulated DEGs
oxup <- ox[ox$log2FoldChange > 0,]
oxdown <- ox[ox$log2FoldChange < 0,]
percentage_ox <- nrow(oxdown) / nrow(ox) * 100

