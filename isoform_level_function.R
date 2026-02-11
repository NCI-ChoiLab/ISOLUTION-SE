library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
library(reshape2)
library(Matrix)
library(Hmisc)
library(wesanderson)
library(magrittr)
library(ggtranscript)
library(rtracklayer)
library(geomtextpath)
load("data/NCI_gtf.RData")
load("data/count_list.RData")
load("data/celltypes.RData")
snp_info <- readRDS("data/snp_info.rds")
Nominal_combined <- readRDS("data/Nominal_combined.rds")
Isoform_info <- readRDS("data/Isoform_info.rds")
genotypes <- readRDS("data/genotypes.rds")
df_transcript_catalog <- readRDS(file = "data/df_transcript_catalog.rds")
transcript_list <- readRDS(file = "data/transcript_list.rds")
# save data
# for (celltype in names(expr_list)) {
#   rds <- expr_list[[celltype]]
#   saveRDS(rds, file = paste0(celltype,".expr.rds"))
# }
# for (celltype in celltypes) {
#   celltype <- gsub(" ","_", celltype)
#   mtx <- norm_mtx_sum[,grepl(paste0(celltype,"-"), colnames(norm_mtx_sum))]
#   saveRDS(mtx, file = paste0(celltype,".mtx.rds"))
# }

# load expression data list
expr_files <- list.files("data/",pattern = ".expr.rds")
celltypes_expr <- str_split_fixed(expr_files, "\\.", 3)[,1]
expr_list <- list()
for (i in 1:length(expr_files)) {
  expr <- readRDS(file = paste0("data/", expr_files[i]))
  expr_list[[celltypes_expr[i]]] <- expr
}

# load normalized expression matrix
mtx_files <- list.files("data/",pattern = ".mtx.rds")
norm_mtx_sum <- NULL
for (i in 1:length(mtx_files)) {
  mtx <- readRDS(file = paste0("data/", mtx_files[i]))
  norm_mtx_sum <- cbind(norm_mtx_sum, mtx)
}

rm(mtx)
rm(expr)

cat.palette = c( "Known"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "fusion"="#ffc125", "Genomic"="#969696", "Antisense"="#66C2A4")
isoform_dist_plot <- function(gene_id, celltype){
  idx <- grepl(gene_id, Isoform_info$annot_gene_id)
  
  isoform_list <- Isoform_info$transcript_name_unique[idx]
  
  expr <- norm_mtx_sum[isoform_list,]
  group <- as.character(str_split_fixed(colnames(expr), "-", 2)[,1])
  
  df_long <- melt(as.data.frame(as.matrix(expr)))
  df_long$Celltype <- str_split_fixed(df_long$variable, "-", 2)[,1]
  df_long$Isoform <- rep(rownames(expr), ncol(expr))
  df_long_sub <- subset(df_long, Celltype == celltype)
  expr_sub <- expr[,which(as.character(group) == celltype)]
  ranks <- names(sort(rowSums(expr_sub),decreasing = TRUE))
  df_long_sub$Isoform <- factor(df_long_sub$Isoform, levels = ranks)
  if (length(ranks) > 37 ) {
    df_long_sub <- subset(df_long_sub, Isoform %in% ranks[1:37])
  }
  
  p <- ggplot(df_long_sub, aes(fill = Celltype, y = value, x = Isoform)) +
    geom_bar(
      position = "dodge", width = 0.8, stat = "summary", fun = "mean",
      color = "black", linewidth = .5
    ) + ylab("Normalized expression")+ggtitle(paste("Isoform composition within", celltype))+
    stat_summary(
      fun.data = mean_sdl, geom = "errorbar", color = "black",
      position = position_dodge(0.8), width = 0.2, linewidth = 0.5
    ) +
    geom_point(
      position = position_jitterdodge(0.5, dodge.width = .8),
      alpha = 0.5
    )+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(p)
}

plot_top_isoforms_by_ct <- function(gene_id, ntop = 4){
  idx <- grepl(gene_id, Isoform_info$annot_gene_id)
  
  isoform_list <- Isoform_info$transcript_name_unique[idx]
  
  expr <- norm_mtx_sum[isoform_list,]
  group <- as.character(str_split_fixed(colnames(expr), "-", 2)[,1])
  
  df_long <- melt(as.data.frame(as.matrix(expr)))
  df_long$Celltype <- str_split_fixed(df_long$variable, "-", 2)[,1]
  df_long$Isoform <- rep(rownames(expr), ncol(expr))
  ranks <- names(sort(rowSums(expr),decreasing = TRUE))
  Isoforms <- head(ranks, ntop)
  df_long_sub <- subset(df_long, Isoform %in% Isoforms)
  df_long_sub$Celltype <- factor(df_long_sub$Celltype, levels = gsub(" ", "_",celltypes))
  p <- ggplot(df_long_sub, aes(fill = Celltype, y = value, x = Celltype)) +
    geom_bar(
      position = "dodge", width = 0.8, stat = "summary", fun = "mean",
      color = "black", linewidth = .8
    ) + scale_fill_manual(values = cellcolors) + 
    facet_wrap(~Isoform, ncol = 1, scales = "free") +
    stat_summary(
      fun.data = mean_sdl, geom = "errorbar", color = "black",
      position = position_dodge(0.8), width = 0.2, linewidth = 0.8
    ) + theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                         legend.position = "none")+
    geom_point(
      position = position_jitterdodge(0.5, dodge.width = .8),
      alpha = 0.8
    )+ylab("Normalized expression")+xlab("")+
    ggtitle(paste0("Top ", ntop, " isoforms within ",gene_id , " across cell types"))
  return(p)
}

Search_transcript_name2 <- function(x){ # transcript id
  name <- Isoform_info$transcript_name_unique[which(Isoform_info$annot_transcript_id %in% x)]
  return(name)
}

isoQTL_plot_pub <- function(celltype,
                            rs,
                            transcript,
                            return_count = FALSE
){
  transcript_name <- Search_transcript_name2(transcript)
  counts <- count_list[[celltype]][,-c(1:4)]
  rownames(counts) <- count_list[[celltype]]$trascript_id
  exprs <- expr_list[[celltype]][,-c(1:4)]
  rownames(exprs) <- expr_list[[celltype]]$trascript_id
  samples <- colnames(counts)
  
  # count matrix
  datainput1 <- data.frame(snp = genotypes[samples,rs],
                           isoform = as.numeric(counts[transcript,]))
  datainput2 <- data.frame(snp = genotypes[samples,rs],
                           isoform = as.numeric(exprs[transcript,]))
  datainput1$snp_ra <- "0|0"
  datainput1$snp_ra[which(datainput1$snp == "01")] <- "1|1"
  datainput1$snp_ra[which(datainput1$snp == "02")] <- "0|1"
  df <- as.vector(table(datainput1$snp_ra))
  datainput1$snp_ra <- paste0("0|0\n", "n=", df[1])
  datainput1$snp_ra[which(datainput1$snp == "01")] <- paste0("1|1\n", "n=", df[3])
  datainput1$snp_ra[which(datainput1$snp == "02")] <- paste0("0|1\n", "n=", df[2])
  p1 <- ggplot(datainput1, aes(x=snp_ra, y=isoform, color=snp_ra)) +
    geom_violin(trim=FALSE,linewidth = 0.75)+
    geom_boxplot(aes(middle=mean(isoform)),
                 width=0.1, linewidth = 0.75)+
    geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
                method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
    labs(title="",x=rs, y = paste0("Count of ", transcript))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(celltype, "\n",transcript_name , ":",  rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
    theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
          axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 16, face = "plain"),
          axis.title.y = element_text(color = "black", size = 16, face = "plain"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none")
  datainput2$snp_ra <- "0|0"
  datainput2$snp_ra[which(datainput2$snp == "01")] <- "1|1"
  datainput2$snp_ra[which(datainput2$snp == "02")] <- "0|1"
  df <- as.vector(table(datainput2$snp_ra))
  datainput2$snp_ra <- paste0("0|0\n", "n=", df[1])
  datainput2$snp_ra[which(datainput2$snp == "01")] <- paste0("1|1\n", "n=", df[3])
  datainput2$snp_ra[which(datainput2$snp == "02")] <- paste0("0|1\n", "n=", df[2])
  p2 <- ggplot(datainput2, aes(x=snp_ra, y=isoform, color=snp_ra)) +
    geom_violin(trim=FALSE,linewidth = 0.75)+
    geom_boxplot(aes(middle=mean(isoform)),
                 width=0.1, linewidth = 0.75)+
    geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
                method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
    labs(title="",x=rs, y = paste0("Normalized expression of ", transcript))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(celltype, "\n",transcript_name , ": ",rs, ": ", paste(snp_info[rs,1], snp_info[rs,3], snp_info[rs,5], snp_info[rs,4], "b38", sep = "_")))+
    theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
          axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 16, face = "plain"),
          axis.title.y = element_text(color = "black", size = 16, face = "plain"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none")
  
  if(return_count){
    return(list(p1,p2))
  }else{
    return(p2)
  }
  
}

plot_isoform_structure <- function(gene_id_of_interest){
  exons <- gtf %>% 
    filter(
      !is.na(gene_name), 
      gene_id == gene_id_of_interest
    ) %>% 
    select(
      seqnames,
      start,
      end,
      strand,
      type,
      gene_name,
      transcript_name
      #transcript_biotype
    ) %>% 
    filter(type == "exon")
  exons$transcript_name_unique <- exons$transcript_name
  exons$transcript_name_unique[grepl("TALON", exons$transcript_name)] <- paste(unique(exons$gene_name), exons$transcript_name[grepl("TALON", exons$transcript_name)], sep = "-")
  exons <- left_join(exons, df_transcript_catalog, by = "transcript_name")
  exons$transcript_catalog <- factor(exons$transcript_catalog, levels = 
                                       c( "Known", "ISM", "NIC", 
                                          "NNC", "fusion", "Genomic", "Antisense"))
  base_plot <-  exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_name
    )) +
    geom_range(
      aes(fill = transcript_catalog)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_name"),
      aes(strand = strand),
      arrow.min.intron.length = 3500
    )
  p <- base_plot + theme_bw() + 
    ggtitle(paste(unique(exons$gene_name), " gene")) +  
    xlab(unique(exons$seqnames)) + ylab("Isoforms") + 
    scale_fill_manual(values = cat.palette) +
    theme(legend.position = "top", legend.spacing.x  =  unit(0.01, 'cm'), legend.key.size = unit(0.25, 'cm')) +
    theme(axis.text.x = element_text(color = "black", size = 12,face = "plain"),
          axis.text.y = element_text(color = "black", size = 12,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 16, face = "plain"),
          axis.title.y = element_text(color = "black", size = 16, face = "plain"),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1))
  gene_gtf <- gtf %>% 
    filter(
      !is.na(gene_name), 
      gene_id == gene_id_of_interest
    )
  rst <- list(p = p, gtf = gene_gtf, iso_n = length(unique(exons$transcript_name)))
  return(rst)
}


