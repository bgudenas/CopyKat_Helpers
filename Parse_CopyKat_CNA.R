
## Loops through all copykat outs and creates list containing a CNA x Cell matrix and ploidy predictions
Plot_sc_CNAs = function( copykat_results_dir  ){
  #copykat_results_dir = "./Smart_Seq/CopyKat/"
  fil_paths = list.files(copykat_results_dir, recursive = TRUE, full.names = TRUE, pattern = "CNA_results.txt")
  
  CNA_mat = matrix()
  for ( i in fil_paths ){
    df = read.table(i,  sep = "\t", header = TRUE)
    rownames(df) = paste0(df$chrom, "_", df$chrompos )
    df = df[ ,-c(1:3)]
   colnames(df) = stringr::str_remove(colnames(df), "^X" )
    CNA_mat = cbind(CNA_mat, df )
  }
  
  fil_paths = list.files(copykat_results_dir, recursive = TRUE, full.names = TRUE, pattern = "prediction.txt")
  
  pred_mat = c()
  for ( i in fil_paths ){
    df = read.table(i,  sep = "\t", header = TRUE)
    pred_mat = rbind(pred_mat, df)
  }
  output = list("CNA_mat" = CNA_mat[ ,-1], "Predictions" = pred_mat )
  return(output)
}


## creates a gene by cell matrix and populates gene,cell with copy number at that loci
gene_CNAs = function(CNA, map, ID = "external_gene_name"){
# map = readRDS("../Annots/Annotables/hg38.rds")
 # CNA = read.table("./Data/test_copykat_CNA_results.txt", header = TRUE)#
  
map = map[!duplicated(map$ensembl_gene_id), ]
map = map[!duplicated(map$external_gene_name), ]
map = map[map$chromosome_name %in% CNA$chrom, ]

genes = map[ ,ID]
geneMat = matrix(nrow= length(genes), ncol = (ncol(CNA)-3), data = 0)
colnames(geneMat) = colnames(CNA)[4:ncol(CNA)]
rownames(geneMat) = genes

for ( gene in 1:nrow(geneMat)){
  for (cell in 1:ncol(geneMat)){
    geneMat[gene, cell] = get_gene_CN(gene_index = gene, cell_index = cell )
  }
}
return(geneMat)
}


## Function: Get Copy number of bin closest to gene[gene_index] of cell[cell_index]
get_gene_CN = function( gene_index, cell_index ){
chrom = map$chromosome_name[gene_index]
start_pos = map$start_position[gene_index]
end_pos = map$end_position[gene_index]


Cmap = CNA$chrompos[CNA$chrom == chrom] 
evec = Cmap - end_pos
svec = Cmap - start_pos

ebin = which.min(abs(evec))
sbin = which.min(abs(svec))


if (ebin == sbin){
  gCN = CNA[CNA$chrom == chrom ,cell_index][ebin]
}
else if (ebin != sbin ){
  sCN = CNA[CNA$chrom == chrom ,cell_index][sbin]
  eCN = CNA[CNA$chrom == chrom ,cell_index][ebin]
  gCN = mean(c(sCN, eCN))
}
return(gCN)
}

#get_gene_CN(gene_index = sample(1:nrow(map), 1), 
#            cell_index = sample(4:ncol(CNA), 1))



#gmat = gene_CNAs(CNA = CNA, map = map )
