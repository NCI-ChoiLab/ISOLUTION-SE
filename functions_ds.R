library(ggplot2)
library(dplyr)
library(stringr)
library(stringi)
library(reshape2)
library(Matrix)
library(Hmisc)
library(wesanderson)
library(magrittr)
library(geomtextpath)
load('eqtl/eQTL.Rdata')
expr_list_eqtl = list()
expr_list_eqtl_files = list.files('eqtl/eQTL_expression/')
for (f in expr_list_eqtl_files){
  celltype = str_split_fixed(f,'\\.',2)[,1]
  expr_list_eqtl[[celltype]] = readRDS(paste0('eqtl/eQTL_expression/',f))
  rm(celltype)
  rm(f)
}
rm(expr_list_eqtl_files)

ct_name <- readRDS("eqtl/ct_name.rds")
SCT_list = list()
for (ct in ct_name) {
  rds <- readRDS(file = paste0("eqtl/SCT_data/",ct,".rds"))
  rds -> SCT_list[[ct]]
}

# #####################
# # downsampling
# names(SCT_list)
# for (i in 1:41) {
#   mtx <- SCT_list[[i]]
#   if (ncol(mtx) > 3000) {
#     idx <- sample(ncol(mtx), 3000, replace = FALSE)
#     mtx_ds <- mtx[,idx]
#     SCT_list[[i]] <- mtx_ds
#   }
# }
# for (ct in ct_name) {
#   rds <- SCT_list[[ct]]
#   saveRDS(rds, file = paste0("eqtl/SCT_data/",ct,".rds"))
# }
# ###########################
eQTL_plot_pub = function(celltype,rs,gene){
  exprs <- expr_list_eqtl[[celltype]][,-c(1:4)]
  rownames(exprs) <- expr_list_eqtl[[celltype]]$gene_id
  samples <- colnames(exprs)
  rs2 = snp_info_eqtl[which(snp_info_eqtl$snp == rs),'variant_id']
  gene_sub = gene_info_eqtl[which(gene_info_eqtl$phenotype_name == gene), 'phenotype_id']
  datainput2 <- data.frame(snp = genotypes_eqtl[samples,rs2],
                           expression = as.numeric(exprs[gene_sub,]))
  
  datainput2$snp_ra <- "0|0"
  datainput2$snp_ra[which(datainput2$snp == "01")] <- "1|1"
  datainput2$snp_ra[which(datainput2$snp == "02")] <- "0|1"
  df <- as.vector(table(datainput2$snp_ra))
  datainput2$snp_ra <- paste0("0|0\n", "n=", df[1])
  datainput2$snp_ra[which(datainput2$snp == "01")] <- paste0("1|1\n", "n=", df[3])
  datainput2$snp_ra[which(datainput2$snp == "02")] <- paste0("0|1\n", "n=", df[2])
  p2 <- ggplot(datainput2, aes(x=snp_ra, y=expression, color=snp_ra)) +
    geom_violin(fill = NA, trim=FALSE,linewidth = 0.25)+
    geom_boxplot(aes(middle=median(expression)), width=0.05, linewidth = 0.25, outlier.shape = NA) +
    geom_smooth(mapping = aes(x = snp_ra, y = expression, group = 1),formula = y~x, color = "gray",
                method='lm', size = .5, se =TRUE,fill = alpha("gray", .5) ) +
    labs(title="",x=rs, y = paste0("Normalized expression of ", gene))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(celltype, "\n",gene , ": ",rs, ": ", paste(rs2, "b38", sep = "_")))+
    theme(text = element_text(family = 'Arial'),
          plot.title = element_text(size = 16, color = 'black', face = 'plain'),
          axis.text.x = element_text(color = "black", size = 12,face = "plain"),
          axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 14, face = "plain"),
          axis.title.y = element_text(color = "black", size = 14, face = "plain"),
          panel.background = element_rect(fill= 'transparent'),
          plot.background = element_rect(fill = 'transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none")
  return(p2)
}


graph_expression = function(g){
  loi <- lapply(SCT_list, function(x) x[g, ])
  df_final = data.frame()
  for (i in 1:length(ct_name)){
    df = data.frame(loi[ct_name[i]])
    colnames(df) = 'Expression'
    df$Celltype = ct_name[i]
    rownames(df) = NULL
    df_final = rbind(df_final,df)
    rm(df)
  }
  
  listofcelltypes <- c('Alveolar Transitional Cells','AT1', 'AT2', 'Club', 'Goblet', 'Basal', 'Multiciliated', 'Neuroendocrine', 'Secretory Transitional Cells',
                       'Alveolar Mph', 'Proliferating Mph', 'NK Cells', 'Proliferating NK', 'CD8+ T Cells', 'CD4+ T Cells', 'Proliferating T', 'Interstitial Mph Perivascular', 'Monocyte-derived Mph', 'Classical Monocytes', 'Non-classical Monocytes', 'DC1', 'DC2', 'Migratory DCs', 'Mast Cells', 'B Cells', 'Plasma Cells', 'Plasmacytoid DCs', 
                       'Lymphatic EC Mature', 'Lymphatic EC Proliferating', 'EC Arterial', 'EC Venous Pulmonary', 'EC Venous Systemic', 'EC Aerocyte Capillary', 'EC General Capillary', 
                       'Adventitial Fibroblasts', 'Smooth Muscle', 'Subpleural Fibroblasts', 'Peribronchial Fibroblasts', 'Alveolar Fibroblasts', 'Mesothelium', 'Myofibroblasts')
  
  cols <- c('darkseagreen1' ,'lawngreen', 'green2','seagreen4', 'olivedrab', 'green1', 'palegreen2', 'limegreen', 'darkgreen', 
            'lightblue', 'royalblue1', 'slategray3', 'cyan1','steelblue', 'aquamarine4', 'cyan3', 'lightblue', 'darkturquoise', 'cornflowerblue', 'dodgerblue4', 'aquamarine', 'cadetblue', 'slateblue1','cornflowerblue','cadetblue3','lightsteelblue3','cyan2',
            'purple1','orchid1','mediumpurple3','mediumorchid2','darkviolet','hotpink','maroon',
            'orange', 'goldenrod2', 'gold', 'tan', 'darkorange3', 'sandybrown', 'lightsalmon') 
  
  df_final$Celltype = factor(df_final$Celltype, levels = listofcelltypes)
  
  ggplot(df_final, aes(x=Celltype, y=Expression, color=Celltype)) +
    geom_violin(aes(fill = Celltype), trim=T, scale = "width") +
    stat_summary(fun = median, geom='point', size = 2, colour = "black", shape = 95)+ 
    labs(y = paste0("Normalized expression of ", g),
         x = 'Cell type') +
    ggtitle(g)+
    theme(text = element_text(family = 'Arial'),
          plot.title = element_text(size = 16, color = 'black', face = 'plain'),
          axis.text.x = element_text(color = "black", size = 12,face = "plain", angle = 75, hjust = 1),
          axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 14, face = "plain"),
          axis.title.y = element_text(color = "black", size = 14, face = "plain"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none") + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
}




