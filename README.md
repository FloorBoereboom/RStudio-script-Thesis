This repository supports the bachelor thesis developed by Floor Boereboom under the supervision of Marvin Martens. 

Transcriptomic RNA-seq data were obtained from proximal tubule epithelial cells (ciPTEC) exposed to cisplatin, carboplatin and oxaliplatin in a previous study by Barnes et al. (2026). 
These data served as input for downstream analysis and consisted of a normalized transcriptome count matrix and corresponding metadata table per Pt-based compound. 
The datasets were preprocessed using the R-ODAF Shiny app which generated DEGs for each Pt-based compound.

Analysis performed:
1. Merged DEGs: The complete DEG tables (including non-significant padj values) from the three Pt-based compounds were merged together to obtain a shared DEG table.
2. Overlapping DEGs: The significant DEG tables (padj<0.01) were merged to obtain the number of overlapping DEGs between two Pt-based compounds.
3. PCA: PCA was performed on the normalized transcriptome count matrix to explore variance between the three compounds.
4. Figure PCA: The resulting PCA plot
5. Ratio up- and downregulated genes: The percentage of downregulated genes in the significant DEGs datasets (padj<0.01) was calculated for each Pt-based compound.
