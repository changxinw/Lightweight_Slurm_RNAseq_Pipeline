library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(openxlsx)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)

make_pseudobulk <- function(expr_mtx, group_mtx){
    genes = rownames(expr_mtx)
    group_mtx$cells = rownames(group_mtx)
    group_mtx$GroupN = paste0(group_mtx[, 1], '@', group_mtx[, 2], recycle0 = TRUE)
    groups = unique(group_mtx$GroupN)
    res_list = list()
    for (i in groups){
        cells = group_mtx[group_mtx$GroupN==i, "cells"]
        res_list[[i]] = rowSums(expr_mtx[, cells])
    }
    res_df = do.call("cbind", res_list)
    rownames(res_df) = genes
    colnames(res_df) = groups
    return(res_df)
}

plot_giggle = function(motif_df, top_n = 15, title = ""){
    motif_df = na.omit(motif_df)
    if (nrow(motif_df) == 0) {return(NULL)}
    motif_df = motif_df[order(-motif_df$giggle_score), ]
    motif_df$factor = factor(motif_df$factor, levels = rev(motif_df$factor))
    motif_df = motif_df[1:top_n, ]
    p <- ggplot(data=motif_df, aes(x=factor, y=giggle_score)) +
         labs(x = "TFs", y = "Overlap Score", title = title) +
         geom_bar(stat="identity") +
         coord_flip() 
#     print(motif_df$fold.enrichment)
    return(p)
}
# plot_motif(enriched.motifs)
plot_motifs = function(motif_df, top_n = 15, title = ""){
#     print(motif_df)
    motif_df = na.omit(motif_df)
    if (nrow(motif_df) == 0) {return(NULL)}
    motif_df = motif_df[order(motif_df$pvalue), ]
    motif_df$motif.name = factor(motif_df$motif.name, levels = rev(motif_df$motif.name))
    motif_df[motif_df$pvalue==0, "pvalue"] = min(motif_df$pvalue[motif_df$pvalue!=0]) / 2
    motif_df = motif_df[1:top_n, ]
    motif_df$logP = -log10(motif_df$pvalue)
    p <- ggplot(data=motif_df, aes(x=motif.name, y=logP, fill = fold.enrichment)) +
         labs(x = "Motifs", y = "-log10(p value)", title = title) +
         geom_bar(stat="identity") +
         coord_flip()
#     print(motif_df$fold.enrichment)
    return(p)
}

pseudoBulkDEPeaks <- function(expr_df, meta_df, contrast, output, obj, up_q=0.01, down_q=0.01, design= ~ Batch + Group){
    prefix <- paste0(output, "/", paste0(contrast[1], "_vs_", contrast[2]))
    print(prefix)
    tt = paste(contrast[1], "v.s.", contrast[2])

    meta_df <- meta_df[meta_df$Group %in% contrast, ]
    meta_df$Group <- factor(meta_df$Group, levels = contrast)
    expr_df <- expr_df[, rownames(meta_df)]

    dds <- DESeqDataSetFromMatrix(countData=expr_df, colData=meta_df, design= design)
    dds <- DESeq(dds)

    contrast <- c("Group", contrast)
    res_de <- as.data.frame(results(dds, contrast=contrast))
    res_de$Peaks <- rownames(res_de)
    res_de = na.omit(res_de)
    res_de = res_de[order(res_de$padj), ]
    write.csv(res_de, paste0(prefix, "_DEseq2_de_table.csv"))
    
    res_list <- list()
    res_list[["de_table"]] = res_de
    res_list[["up_peaks"]] <- rownames(na.omit(res_de) %>% dplyr::filter(padj<=up_q, log2FoldChange>=1))
    if (length(res_list[["up_peaks"]])<10) {res_list[["up_peaks"]] <- rownames(na.omit(res_de) %>% dplyr::filter(pvalue<=up_q, log2FoldChange>=1))}
    res_list[["down_peaks"]] <- rownames(na.omit(res_de) %>% dplyr::filter(padj<=down_q, log2FoldChange<=-1))
    if (length(res_list[["down_peaks"]])<10) {res_list[["down_peaks"]] <- rownames(na.omit(res_de) %>% dplyr::filter(pvalue<=down_q, log2FoldChange<=-1))}
    if (length(res_list[["up_peaks"]]) > 0){
        res_list[["up_motifs"]] <- FindMotifs(object = obj, features = res_list[["up_peaks"]])
        res_list[["up_tfs"]] <- RunGiggle(res_list[["up_peaks"]], output = output)
        p = plot_motifs(res_list[["up_motifs"]], title = paste(tt, "Up Motifs"))
        ggsave(paste0(prefix, "_de_upmotifs.pdf"), p, width = 5, height = 6)
        p = plot_giggle(res_list[["up_tfs"]], title = paste(tt, "Up TFs"))
        ggsave(paste0(prefix, "_de_uptfs.pdf"), p, width = 5, height = 6)
    }
    if (length(res_list[["down_peaks"]]) > 0){
        res_list[["down_motifs"]] <- FindMotifs(object = obj, features = res_list[["down_peaks"]])
        res_list[["down_tfs"]] <- RunGiggle(res_list[["down_peaks"]], output = output)
        p = plot_motifs(res_list[["down_motifs"]], title = paste(tt, "Down Motifs"))
        ggsave(paste0(prefix, "_de_downmotifs.pdf"), p, width = 5, height = 6)
        p = plot_giggle(res_list[["down_tfs"]], title = paste(tt, "Down TFs"))
        ggsave(paste0(prefix, "_de_downtfs.pdf"), p, width = 5, height = 6)
    }
    openxlsx::write.xlsx(res_list, paste0(prefix, "_de_result.xlsx"))
    return(res_list)
}

RunGiggle <- function(inputPeaks, organism = "GRCm38", output = "./"){
  gIndex = paste0("/datacommons/ydiaolab/genome_ref/giggle.all/giggle.", organism)
  gPath = "/hpc/group/jilab/changxin/miniconda3/envs/MAESTRO/bin/giggle"
  antFile = "/datacommons/ydiaolab/genome_ref/giggle.all/CistromeDB.sample.annotation.txt"
  outputBed = paste0(output, "giggle.bed")
  write.table(t(as.data.frame(strsplit(inputPeaks, "-"))), outputBed, quote = F, row.names = F, col.names = F, sep = "\t")
  cmd <- paste0("sort --buffer-size 2G -k1,1 -k2,2n -k3,3n ", outputBed, " | bgzip -c > ", outputBed, ".gz")
  system(cmd)
  cmd <- paste0(gPath, " search -i ", gIndex, " -q ", outputBed, ".gz -s > ", outputBed, ".result.xls")
#   print(cmd)
  system(cmd)
  resultDf <- read.table(paste0(outputBed, ".result.xls"), sep="\t", row.names=NULL, comment.char="", stringsAsFactors =  FALSE)
  resultDf <- resultDf[,-9]
  colnames(resultDf) <- c("file", "file_size", "overlaps", "odds_ratio", "fishers_two_tail", "fishers_left_tail", "fishers_right_tail", "combo_score")
  resultDf <- resultDf[resultDf$overlaps>0,]
  if(organism == "GRCh38"){
    rownames(resultDf) <- sapply(strsplit(resultDf$file, "human/"), function(x) return(gsub("_5foldPeak.bed.gz", "", x[2])))
  }
  if(organism == "GRCm38"){
    rownames(resultDf) <- sapply(strsplit(resultDf$file, "mouse/"), function(x) return(gsub("_5foldPeak.bed.gz", "", x[2])))
  }
  resultDf <- resultDf[,c("file_size", "overlaps", "combo_score")]
  antFile <- read.csv(antFile, sep="\t", row.names=1, stringsAsFactors = FALSE)  
  targetDf <- merge(resultDf, antFile, by.x=0, by.y=0)
  colnames(targetDf) <- c("sample_id", "sample_peak_number", "overlap_peak_number", "giggle_score", "GSM_id", "species", "factor", "factor_type", "cell_line", "cell_type", "tissue", "disease")                          
  targetDf$biological_resource <- apply(targetDf, 1, function(x) return(paste0(x[9:11], collapse=";")))
  targetDf <- targetDf[, c("sample_id", "GSM_id", "species", "factor", "factor_type", "biological_resource", "giggle_score", "sample_peak_number", "overlap_peak_number")]
  targetDf_tf <- targetDf[targetDf$factor_type=='tf',]
  targetDf_tf <- targetDf_tf[order(-targetDf_tf$giggle_score), ]
  targetDf_tf <- targetDf_tf[!duplicated(targetDf_tf$factor), ]

  cmd <- paste0("rm ", outputBed, "*")
  system(cmd)
  return(targetDf_tf)
}
                                        
                                        
### VocanoPlot with the output of DEseq2
library(ggplot2)
library(dplyr)
library(ggrepel)
library(Seurat)

VolcanoPlot <- function(data, top_n=10, fc=1, adj.p=0.05, marked_genes=c(), tt=""){
    xl <- -fc
    xr <- fc
    yp <- -log10(adj.p)
    ###label differentially expressed genes
    dif.data <- na.omit(data) %>%
      mutate(logP = -log10(padj)) %>%
      mutate(color = ifelse(log2FoldChange > xr & logP > yp,
                            yes = "Case", no = ifelse(log2FoldChange < xl & logP > yp, yes = "Control",  no = "none")))
    ###extract top differentially expressed genes
    up <- arrange(subset(dif.data, color == "Case"), desc(log2FoldChange))
    down <- arrange(subset(dif.data, color == "Control"), log2FoldChange)
    if(nrow(up) >= top_n && nrow(down) >= top_n){
      top_labelled <- rbind.data.frame(up[1:top_n,], down[1:top_n,])
    }else if (nrow(up) < top_n && nrow(down) >= top_n){
      top_labelled <- rbind.data.frame(up, down[1:top_n,])
    }else if (nrow(up) >= top_n && nrow(down) < top_n){
      top_labelled <- rbind.data.frame(up[1:top_n,], down)
    }else{
      top_labelled <- rbind.data.frame(up, down)
    }
    marked_df <- dif.data[dif.data[, "SYMBOL"] %in% marked_genes, ]
    if (nrow(marked_df) > 0){
        top_labelled <- rbind.data.frame(top_labelled, marked_df)
        top_labelled <- top_labelled[!duplicated(top_labelled), ]
    }
    
    ###volcano plot
    p <- ggplot(dif.data, aes(x = log2FoldChange, y = logP)) +
      geom_point(aes(color = factor(color)), size = 1.55, alpha = 0.8, na.rm = TRUE) + # add gene points
      theme_bw(base_size = 16) + # clean up theme
      theme(legend.position = "none") + # remove legend
      ggtitle(tt) +
      xlab(expression(log[2]("FC"))) + # x-axis label
      ylab(expression(-log[10]("adj.P.Val"))) + # y-axis label
      geom_vline(xintercept = xl, colour = "black", linetype = "dashed") + # add line at 0
      geom_vline(xintercept = xr, colour = "black", linetype = "dashed") + # add line at 0
      geom_hline(yintercept = yp, colour = "black", linetype = "dashed") + # p(0.05) = 1.3
      scale_color_manual(values = c("Case" = "#E64B35",
                                    "Control" = "#3182bd",
                                    "none" = "#636363")) + # change colors
      scale_y_continuous(trans = "log1p")+  # Scaled Y-axis with log1p function
      geom_label_repel(data = top_labelled,
                       aes(label = SYMBOL), fill = "white",
                       fontface = 'bold',box.padding = 0.25, color = 'black',
                       label.size = 0.15,point.padding = 0.5,segment.color = 'gold', max.overlaps = Inf)
#     ggsave(output, p, width=5, height=5)
    return(p)
}
                                        
                                        
pseudoBulkDE <- function(expr_df, meta_df, contrast, output, logfc=1, pcut=0.05, pvalue = FALSE, group = "Group", species = "human", up_q=0.05, down_q=0.05, design= ~ Batch + Group, rGene = FALSE, marked_genes = c(), ...){
    if (!file.exists(output)) {dir.create(output)}
    prefix <- paste0(output, "/", paste0(contrast[1], "_vs_", contrast[2]))
    print(prefix)
    tt = paste(contrast[1], "v.s.", contrast[2])
    
    meta_df <- meta_df[meta_df[, group] %in% contrast, , drop = FALSE]
    meta_df[,group] <- factor(meta_df[, group], levels = contrast)
    meta_df[,group] <- relevel(meta_df[, group], ref = contrast[2])
    expr_df <- expr_df[, rownames(meta_df)]
    
    dds <- DESeqDataSetFromMatrix(countData=expr_df, colData=meta_df, design= design)
    dds <- DESeq(dds) 
    
    contrast <- c(group, contrast)
    res_de <- as.data.frame(results(dds, contrast=contrast))
    res_de = res_de[order(res_de$padj, -abs(res_de$log2FoldChange)), ]
    write.csv(res_de, paste0(prefix, "_DEseq2_de_table.csv"))
    res_de$SYMBOL <- rownames(res_de)
    p <- VolcanoPlot(as.data.frame(res_de), fc = logfc, marked_genes = marked_genes, tt = tt)
    ggsave(paste0(prefix, "_Volcano_deseq_result.pdf"), p, width=6, height=6)
    
    ### remove the gene in OE or KO experiment
    if (rGene!=FALSE){
        res_de = res_de[rownames(res_de) != rGene, ]
    }

    gene_list <- list()
    if (pvalue) {
        gene_list[["up"]] <- rownames(na.omit(res_de) %>% dplyr::filter(pvalue<=pcut, log2FoldChange>=logfc))
        gene_list[["down"]] <- rownames(na.omit(res_de) %>% dplyr::filter(pvalue<=pcut, log2FoldChange<=-logfc))
    } else {
        gene_list[["up"]] <- rownames(na.omit(res_de) %>% dplyr::filter(padj<=pcut, log2FoldChange>=logfc))
        gene_list[["down"]] <- rownames(na.omit(res_de) %>% dplyr::filter(padj<=pcut, log2FoldChange<=-logfc))
    }
    
    if (species == "human"){
        orgdb = "org.Hs.eg.db"
    } else {
        orgdb = "org.Mm.eg.db"
    }
    
    enrich_list = list()
    enrich_list[["de_table"]] = res_de[order(res_de$padj), ]
    enrich_list[["up_genes"]] = gene_list[["up"]]
    enrich_list[["down_genes"]] = gene_list[["down"]]
    de_result_list = list()
    de_result_list[["contrast"]] <- contrast
    de_result_list[["dds"]] <- dds
    de_result_list[["de_table"]] = res_de[order(res_de$padj), ]
    if (length(gene_list[["up"]])>10){
        up_enrich <- enrichGO(gene = gene_list[["up"]], OrgDb = orgdb, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff=up_q, universe=rownames(res_de), ...)
    #     if (nrow(up_enrich>0)){
    #         up_enrich@result[, "EnrichRatio"] <- (as.numeric(sapply(strsplit(up_enrich@result$GeneRatio, "/"), function(i) i[1])) * as.numeric(sapply(strsplit(up_enrich@result$BgRatio, "/"), function(i) i[2]))) / (as.numeric(sapply(strsplit(up_enrich@result$GeneRatio, "/"), function(i) i[2])) * as.numeric(sapply(strsplit(up_enrich@result$BgRatio, "/"), function(i) i[1])))
    #         up_enrich@result <- up_enrich@result[up_enrich@result$qvalue <= 0.05, ]
    #         up_enrich@result <- up_enrich@result[order(-up_enrich@result[, "EnrichRatio"]), ]
    #     }
        if (nrow(up_enrich)!=0){
            p <- dotplot(up_enrich, color="qvalue", x="GeneRatio/BgRatio", showCategory = 30)
            ggsave(paste0(prefix, "_Dotplot_enrichment_upgenes_deseq_result.pdf"), p, width=10, height=10)
        }
        enrich_list[["up_enrich"]] = up_enrich@result
        de_result_list[["up_enrich"]] <- up_enrich
    }
    
    if (length(gene_list[["down"]])>10){
        down_enrich <- enrichGO(gene = gene_list[["down"]], OrgDb = orgdb, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff=down_q, universe=rownames(res_de), ...)
    #     if (nrow(down_enrich)>0){
    #         down_enrich@result[, "EnrichRatio"] <- (as.numeric(sapply(strsplit(down_enrich@result$GeneRatio, "/"), function(i) i[1])) * as.numeric(sapply(strsplit(down_enrich@result$BgRatio, "/"), function(i) i[2]))) / (as.numeric(sapply(strsplit(down_enrich@result$GeneRatio, "/"), function(i) i[2])) * as.numeric(sapply(strsplit(down_enrich@result$BgRatio, "/"), function(i) i[1])))
    #         down_enrich@result <- down_enrich@result[down_enrich@result$qvalue <= 0.05, ]
    #         down_enrich@result <- down_enrich@result[order(-down_enrich@result[, "EnrichRatio"]), ]
    #     }
        if (nrow(down_enrich)!=0){
            p <- dotplot(down_enrich, color="qvalue", x="GeneRatio/BgRatio", showCategory = 30)
            ggsave(paste0(prefix, "_Dotplot_enrichment_downgenes_deseq_result.pdf"), p, width=10, height=10)
        }
        enrich_list[["down_enrich"]] = down_enrich@result
        de_result_list[["down_enrich"]] <- down_enrich
    }
    ### GSEA enrichment analysis
    gene_vec = na.omit(res_de)$log2FoldChange
    names(gene_vec) = rownames(na.omit(res_de))
    gsea_enrich = gseGO(sort(gene_vec, decreasing = TRUE), OrgDb = orgdb, keyType = "SYMBOL")
    if (nrow(gsea_enrich)!=0){
            p <- dotplot(gsea_enrich, x="NES", showCategory = 30)
            ggsave(paste0(prefix, "_Dotplot_enrichment_gsea_deseq_result.pdf"), p, width=10, height=10)
        }
    enrich_list[["gsea_enrich"]] = gsea_enrich@result
    de_result_list[["gsea_enrich"]] <- gsea_enrich
    ### save result
    saveRDS(de_result_list, paste0(prefix, "_de_result_list.rds"))
    openxlsx::write.xlsx(enrich_list, paste0(prefix, "_de_result.xlsx"))
    return(de_result_list)
}


### enrichment analysis of cell markers
compute_enrich = function(genes, species = "mm10", ...){
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(clusterProfiler)

    if (species == "human"){
        orgdb = "org.Hs.eg.db"
    } else {
        orgdb = "org.Mm.eg.db"
    }

    res_enrich <- enrichGO(gene = genes,
                           OrgDb = orgdb,
                           keyType = "SYMBOL",
                           ont = "BP",
                           pAdjustMethod = "fdr",
                           pvalueCutoff = 1,
                           qvalueCutoff=0.05,
                           ...)
    return(res_enrich)
}


### enrichment analysis of cell markers
compute_gsea = function(genes, species = "mm10", ...){
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(clusterProfiler)

    if (species == "hg38"){
        orgdb = "org.Hs.eg.db"
    } else {
        orgdb = "org.Mm.eg.db"
    }
    gsea_enrich = gseGO(sort(genes, decreasing = TRUE), OrgDb = orgdb, keyType = "SYMBOL")
    return(gsea_enrich)
}

plot_lisa = function(lisa_path, top_n = 15, title = ""){
    lisa_df = read.csv(lisa_path)
    lisa_df = lisa_df[, c(1, 2)]
    colnames(lisa_df) = c("TFs", "Samples")
    lisa_df$Num = sapply(lisa_df$Samples, function(x){
        -log10(as.numeric(strsplit(x, ";")[[1]][2]))
    })
    lisa_df = lisa_df[1:top_n, ]
    lisa_df$TFs = factor(lisa_df$TFs, levels = rev(lisa_df$TFs))
    p <- ggplot(data=lisa_df, aes(x=TFs, y=Num)) +
         labs(x = "TFs", y = "-log10P", title = title) +
         geom_bar(stat="identity") +
         coord_flip()
    return(p)
}
