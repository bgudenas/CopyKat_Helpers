library(pheatmap)
colos = colorRampPalette(colors = c("blue4","blue","white","red","red4") )(100)

copy_fils = list.files("./TenX/Copykat", full.names = TRUE, recursive = TRUE, pattern = "CNA_results.txt")
preds_fils = list.files("./TenX/Copykat", full.names = TRUE, recursive = TRUE, pattern = "prediction.txt")

for ( i in 8:length(copy_fils)) {
samp = unlist(lapply(stringr::str_split(copy_fils[i], "\\/" ), "[[", 4 ))
outplot = paste0("./TenX/Copykat/", samp, "/", samp, "_CNA_heatmap.png")
if (!file.exists(outplot)){
  
CNA = read.table(copy_fils[i], header = TRUE)
print(i)
if (ncol(CNA) > 20000 ){next}
preds = read.table(preds_fils[i], header = TRUE)

Cmat = CNA[ ,4:ncol(CNA)]
rownames(Cmat) = paste0("row", seq(nrow(Cmat)))
#chroms = data.frame(row.names = paste0("row", seq(nrow(Cmat))), "Chromosome" = as.factor(CNA$chrom ))
chroms = CNA$chrom
rm(CNA)

Cmat = scale(Cmat)
Cmat[Cmat > 4 ] = 4
Cmat[Cmat < -4 ] = -4

gaps_chroms = c()
for ( j in 1:(length(chroms)-1)){
  if (chroms[j] != chroms[j+1]){
    gaps_chroms = c(gaps_chroms, j) 
  }
}

#Cmat[rowVars(Cmat) > 0.1, ]

png(outplot, width = 7, height = 10, units = "in", res = 120)
pheatmap(Cmat,
         color = colos,
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = "none",
         gaps_row = gaps_chroms,
        annotation_row = data.frame(row.names = rownames(Cmat), "Chromosome" = as.factor(chroms)),
        annotation_col =  data.frame(row.names = preds$cell.names, "Ploidy" = preds$copykat.pred ),
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         cluster_rows = FALSE,
        main = paste0(samp, " [Cells = ", ncol(Cmat), "]")
        )

dev.off()

  }
}
