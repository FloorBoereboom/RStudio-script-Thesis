This repository supports the bachelor thesis developed by Floor Boereboom under the supervision of Marvin Martens. 

Transcriptomic RNA-seq data were obtained from proximal tubule epithelial cells (ciPTEC) exposed to Pt-based compounds in a previous study by Barnes et al. (2026). 
These data served as input for downstream analysis and consisted of a normalized transcriptome count matrix and corresponding metadata table. 
The datasets were preprocessed using the R-ODAF Shiny app which generated DEGs for each Pt-based compound.

Analysis performed:
1. Merged DEGs: DEG tables from the three Pt-based compounds were merged together to obtain a shared DEG table.
2. Overlapping DEGs: The number of DEGs overlapping between two Pt-based compounds was retreived
3. PCA: PCA was performed to explore variance in the dataset
4. Figure PCA: The resulting PCA plot
5. Ratio up- and downregulated genes: The percentage of downregulated genes in the significant DEGs dataset was calculated for each Pt-based compound
