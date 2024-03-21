#This code is to produce gene activity and gene expression plots starting from the processed PVH or ARC seurat objects.
#Original code by: Ji-Dung Luo

#Load libraries
library(Seurat)
library(pheatmap)
library(ggplot2)

#Load seurat object - In this case using the PVH object
seu <- readRDS(file = "seuOBJ_snRNA_snATAC_allSample_PVH_merge_subClust.rds")

#Extract metadata
meta <- seu@meta.data

#For viewing cell numbers
table(meta$cluster[meta$tech == "ATAC"])
table(meta$cluster[meta$tech == "RNA"])


count <- GetAssayData(seu,assay = "RNA",slot = "data") #Extract gene counts
count <- as.matrix(count) 
clust_atac <- unique(meta$cluster[meta$tech == "ATAC"]) #Find cells from scATACseq in each cluster
clust_rna <- unique(meta$cluster[meta$tech == "RNA"]) #Find cells from scRNAseq in each cluster
clust_share <- clust_rna[clust_rna %in% clust_atac] #Find overlap

mark_ge <- c("Agtr1a","Sv2c","Gm49732","Bdnf","Pde3a","Gm39185","Esr2","Ligp","Plekhg1","Ebf3","Cpa6","Htr1d","Gm29771") #Genes you wish to plot

mark_ge <- mark_ge[mark_ge %in% rownames(count)] #Extract counts for desired gene markers

#
calRes_list <- list()
#
#For atac
cell_atac <- rownames(meta)[meta$tech == "ATAC"]
meta_atac <- meta[meta$tech=="ATAC",]
#

calPlot_list <- lapply(c(0:24),function(clust,ge,mat,meta){ #Apply function to each cluster from 0 to 24
  cell <- rownames(meta)[meta$cluster==clust]  #Extract cell names for the current cluster from the metadata
  exp <- mat[ge,cell]  #Extract expression data for the desired genes and cells from the matrix
  exp <- rowSums(exp) #Calculate the total expression for each gene across cells in the cluster
  return(exp)},mark_ge,count,meta_atac)

names(calPlot_list) <- c(0:24) #Set names as clusters
calPlot_dat <- do.call(cbind,calPlot_list) #Combine the elements of the list into a single data frame
calPlot_dat <- calPlot_dat[!matrixStats::rowSds(calPlot_dat) == 0,] #Remove rows with zero stdev from the data frame
calRes_list[["ATAC"]] <- calPlot_dat #Store results

pheatmap::pheatmap(calPlot_dat,scale = "row",cluster_rows = FALSE,cluster_cols = FALSE) #Plot gene activity


# Process gene expression
cell_rna <- rownames(meta)[meta$tech == "RNA"]
meta_rna <- meta[meta$tech=="RNA",]
#
calPlot_list <- lapply(mark_ge,function(ge,mat,clust,meta){
  exp <- mat[ge,] #Extract expression data for the current gene of interest
  exp_list <- lapply(clust,function(x,exp,meta2){ 
    cell_sel <- rownames(meta2)[meta2$cluster==x]  #Select cells corresponding to the current cluster
    ex <- exp[cell_sel]
 #Calculate median expression of gene and proportion of cells expressing the gene of interest
    to <- length(exp) 
    pos <- length(ex[ex > 0])
    med <- mean(ex)

    res <- setNames(c(med,(pos/to)*100),c("exp","prop")) #Store results
    return(res)},exp,meta)
  names(exp_list) <- clust
  exp_mat <- do.call(rbind,exp_list)
  exp_dat <- data.frame(gene=ge,cluster=rownames(exp_mat),stringsAsFactors = FALSE)
  exp_dat <- cbind(exp_dat,exp_mat)
  return(exp_dat)},count,clust_share,meta_rna)
names(calPlot_list) <- mark_ge
calPlot_dat <- do.call(rbind,calPlot_list)
calPlot_dat$cluster <- factor(calPlot_dat$cluster,levels = c(0:24))
calRes_list[["RNA"]] <- calPlot_dat

#Plot gene expression
dat <- calRes_list$RNA
dat <- dat[!dat$cluster %in% c(9,15,20,21,23,24),]
dat <- dat[dat$gene %in% rownames(calRes_list$ATAC),]
dat$gene <- factor(dat$gene,levels = rev(rownames(calRes_list$ATAC)))
ggplot(dat,aes(x=gene,y=cluster,size=prop,color=exp))+geom_point()+
    labs(x="",y="Clusters",color="avg_exp",size="Percent")+
    scale_color_gradient(low="white",high = "darkblue")+
    theme_classic()+coord_flip()



