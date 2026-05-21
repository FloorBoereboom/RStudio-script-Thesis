# Retreive number of significant DEGs overlapping between two Pt-based compounds
#1. Import complete DEGs tables (including non-significant values)
cis <- read.delim("data/cisplatin_full_DEGs.txt") 
car <- read.delim ("data/carboplatin_full_DEGs.txt")
ox <- read.delim("data/oxaliplatin_full_DEGs.txt")


#Retreive overlap of DEGs between compounds cisplatin and carboplatin
merge1 <- merge(cis, car, by = "GeneID")
sig1 <- merge1[merge1$padj.cis < 0.01, ]
sig2 <- sig1[sig1$padj.car < 0.01, ]


#Retreive overlap of DEGs between compounds cisplatin and oxaliplatin
merge2 <- merge(cis, ox, by = "GeneID")
sig1 <- merge2[merge2$padj.cis < 0.01, ]
sig2 <- sig1[sig1$padj.car < 0.01, ]

#Retreive overlap of DEGs between compounds carboplatin and oxaliplatin
merge3 <- merge(car, ox, by = "GeneID")
sig1 <- merge3[merge3$padj.car < 0.01, ]
sig2 <- sig1[sig1$padj.ox < 0.01, ]
