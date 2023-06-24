# Code for manuscript: Unraveling the Heterogeneity of Tumor-Infiltrating Myeloid Cells in Immune Checkpoint Blockade: A Single-Cell Pan-cancer analysis

library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggtext)
library(patchwork)

cdir <- "Source_Data/Myeloid_Figures/"
dir.create(cdir)
setwd("Source_Data/")

########################################################################################################
########################################################################################################
# Summary
data <- NULL
data <- readRDS("Myeloid_Pan_Cancer_Post_Filtered.RDS")

x <- NULL
x <- data.frame(unique(data@meta.data[,c("Dataset_ID","Cancer","Sample_ID","Patient_ID","Treatment_Group","Treatment","Drug_Name","Response")]))
row.names(x) <- NULL
x <- x[order(x$Cancer, decreasing = F),]
xlsx::write.xlsx(x, "PanCancer_Myeloid_Table_1.xlsx", row.names = F, sheetName = "Table_1")

View(table(data@meta.data[,c("Dataset_ID","Cancer","Treatment_Group","Treatment","Drug_Name","Response")]))

plotx$Overall_Statistics <- paste(plotx$Count, "(",plotx$Overall_Proportion,"%)", sep = "")
plotx$Overall_Statistics <- paste(plotx$Label, plotx$Overall_Statistics, sep = "\n") # plotx$Cancer,
plotx$Intra_Statistics <- paste(plotx$Count, "(",plotx$Intra_Proportion,"%)", sep = "")
plotx$Intra_Statistics <- paste(plotx$Label, plotx$Intra_Statistics, sep = "\n") # plotx$Cancer, 
plotx <- plotx[order(plotx$Cancer, decreasing = F),]
plotx$Index <- 1:nrow(plotx)
saveRDS(plotx, "Myeloid_Figure1A_Summary_Study_Plotx.RDS")

plotx <- readRDS("Myeloid_Figure1A_Summary_Study_Plotx.RDS")
plotx <- plotx %>% add_row(Index = nrow(plotx)+1, Count = 0)
plotx$Overall_Statistics[which(plotx$Index == nrow(plotx))] <- ""

p <- NULL
p <- ggplot(plotx, aes(Index, Count, color = Cancer, label = Overall_Statistics)) +
  geom_segment(aes(x = Index, xend = Index, y = 0, yend = Count), size = 1.2) +
  geom_rect(aes(xmin = 1, xmax = nrow(plotx), ymin = -50, ymax = 0), fill = "grey97", color = "grey97") + 
  geom_point(aes(size = Count)) +
  rcartocolor::scale_color_carto_d(palette = "Prism") +
  scale_size(range = c(1, 8), limits = c(0, max(plotx$Count)), guide = "none") +
  # scale_color_manual(values = cancercols) +
  coord_polar() +
  theme_void()+
  geom_text_repel(aes(label = Overall_Statistics, y = Count+20), size = 3, max.overlaps  = Inf)

somePDFPath = paste(cdir, "Myeloid_Figure1A_Summary_Study.pdf", sep = "")
pdf(file=somePDFPath, width=10, height=8,pointsize=12)
print(p)
dev.off()

plotx <- readRDS("Myeloid_Figure1A_Summary_Study_Plotx.RDS")
plotx$Index <- as.numeric(factor(plotx$Label, levels = sort(unique(plotx$Label), decreasing = T)))

plotx <- plotx %>% 
  add_row(Index = max(plotx$Index)+1, Count = 0, Cancer = unique(plotx$Cancer)) %>% 
  filter(!is.na(Cancer))

p2 <- NULL
p2 <- ggplot(plotx, aes(Index, Intra_Proportion, color = Cancer, label = Intra_Statistics)) +
  geom_segment(aes(x = Index, xend = Index, y = 0, yend = Intra_Proportion), size = 1.2) +
  geom_rect(aes(xmin = 1, xmax = max(plotx$Index), ymin = -20, ymax = 0), fill = "grey97", color = "grey97") + 
  geom_point(aes(size = Count)) +
  rcartocolor::scale_color_carto_d(palette = "Prism", guide = "none") +
  scale_size(range = c(1, 8), limits = c(0, max(plotx$Intra_Proportion)), guide = "none") +
  # scale_color_manual(values = cancercols) +
  theme_void()+
  facet_wrap(~ Cancer, nrow = 2) +
  coord_polar()+
  geom_text_repel(aes(label = Intra_Statistics, y = Intra_Proportion+0.3), size = 3, max.overlaps  = Inf)

somePDFPath = paste(cdir, "Myeloid_Figure1A_Summary_Study_By_Cancer.pdf", sep = "")
pdf(file=somePDFPath, width=10, height=8,pointsize=12)
print(p2)
dev.off()

################################################################################
################################################################################
# UMAP
data <- NULL
data <- readRDS("Myeloid_Pan_Cancer_Post_Filtered.RDS")
data$Batch <- paste("Batch0",as.numeric(as.factor(data$Dataset_ID)), sep = "")

somePDFPath = paste(cdir, "Figure1B_UMAP_Batch_Effects.pdf", sep = "")
pdf(file=somePDFPath, width=10, height=8,pointsize=12)
DimPlot(data, group.by = "Batch", cols = color_conditions$colorful)+NoAxes()+NoGrid()
data <- RunUMAP(data, reduction = "pca", dims = 1:30)
DimPlot(data, group.by = "Batch", cols = color_conditions$colorful)+NoAxes()+NoGrid()
dev.off()

ctcols <- readRDS("Settings/Myeloid_Celltype_Cols_Pan_Cancer.RDS")

somePDFPath = paste(cdir, "Figure1B_Myeloid_UMAP.pdf", sep = "")
pdf(file=somePDFPath, width=10, height=6,pointsize=12)
DimPlot(data, group.by = "Manual_Annot", cols = ctcols)+NoAxes()+NoGrid()
DimPlot(data, group.by = "Manual_Annot", split.by = "Treatment_Group", ncol = 2, cols = ctcols)+NoAxes()+NoGrid()
DimPlot(subset(data, subset = Treatment_Group == "Post" & Response %in% c("Responder","Non-responder")), group.by = "Manual_Annot", split.by = "Response", ncol = 2, cols = ctcols)+NoAxes()+NoGrid()
dev.off()

########################################################################################################
########################################################################################################
# Proportion

ccols <- readRDS("Settings/Myeloid_Celltype_Cols_Pan_Cancer.RDS")
temp <- data.frame(Manual_Annot = sort(names(ccols)))
temp$Manual_Annot <- factor(temp$Manual_Annot)
temp$Manual_Annot_Index <- as.numeric(temp$Manual_Annot)
temp <- temp[order(temp$Manual_Annot_Index, decreasing = F),]

data <- NULL
data <- readRDS("Myeloid_Pan_Cancer_Post_Filtered.RDS")
data$IID <- paste(data$Cancer, data$Manual_Annot, data$Treatment_Group, sep = "|")
prop <- NULL
prop <- data.frame(table(data$IID))
colnames(prop) <- c("IID","Count")
prop$Cancer <- gsub("(.*?)\\|(.*?)\\|(.*)","\\1",prop$IID, ignore.case = T)
prop$Manual_Annot <- gsub("(.*?)\\|(.*?)\\|(.*)","\\2",prop$IID, ignore.case = T)
prop$Treatment_Group <- gsub("(.*?)\\|(.*?)\\|(.*)","\\3",prop$IID, ignore.case = T)
prop$Group <- "Myeloid"
prop$ID <- paste(prop$Cancer,prop$Treatment_Group,sep = "|")
prop <- split(prop, prop$ID)
prop <- lapply(prop, function(x){
  x$Proportion <- x$Count/sum(x$Count)
  return(x)
})
prop <- do.call(rbind.data.frame, prop)

ccols
prop$IIID <- paste(prop$Group, "(", prop$Treatment_Group, ")", sep = "")
prop$IIID <- factor(prop$IIID, levels = c("Myeloid(Pre)","Myeloid(Post)"))
prop$Treatment_Group <- factor(prop$Treatment_Group, levels = c("Pre","Post"))
prop$Manual_Annot_Index <- temp[match(prop$Manual_Annot, temp$Manual_Annot),"Manual_Annot_Index"]
saveRDS(prop,"Myeloid_Figure1B_Pre_Post_Prop.RDS")

prop <- readRDS("Myeloid_Figure1B_Pre_Post_Prop.RDS")

somePDFPath = paste(cdir, "Myeloid_Figure1B_Pre_Post_Barplot.pdf", sep = "")
pdf(file=somePDFPath, width=14, height=6,pointsize=12)
ggplot(prop, aes(Cancer, Proportion, fill = Manual_Annot, label = Manual_Annot_Index))+
  geom_bar(position = "stack", stat = "identity")+
  facet_wrap(~IIID, ncol = 4)+
  # geom_text_repel(size = 4, position = position_stack(vjust = 0), lineheight = 0.7, segment.alpha = 0.3)+
  scale_fill_manual(values = ccols)+theme_linedraw(base_size = 20)+
  RotatedAxis()+theme(legend.position = "bottom", legend.title= element_blank()) +
  guides(color=guide_legend(ncol=2))+xlab("")
dev.off()

####### Post R/NR Bar Plot #######
data <- subset(data, subset = Treatment_Group == "Post")
data <- subset(data, subset = Response %in% c("Responder","Non-responder"))
data$IID <- paste(data$Cancer, data$Manual_Annot, data$Response, sep = "|")

prop <- data.frame(table(data$IID))
colnames(prop) <- c("IID","Count")
prop$Cancer <- gsub("(.*?)\\|(.*?)\\|(.*)","\\1",prop$IID, ignore.case = T)
prop$Manual_Annot <- gsub("(.*?)\\|(.*?)\\|(.*)","\\2",prop$IID, ignore.case = T)
prop$Response <- gsub("(.*?)\\|(.*?)\\|(.*)","\\3",prop$IID, ignore.case = T)
prop$Group <- "Myeloid"
prop$ID <- paste(prop$Cancer,prop$Response,sep = "|")
prop <- split(prop, prop$ID)
prop <- lapply(prop, function(x){
  x$Proportion <- x$Count/sum(x$Count)
  return(x)
})
prop <- do.call(rbind.data.frame, prop)

ccols
prop$IIID <- paste(prop$Group, "(", prop$Response, ")", sep = "")
prop$IIID <- factor(prop$IIID, levels = c("Myeloid(Responder)","Myeloid(Non-responder)"))
prop$Manual_Annot_Index <- temp[match(prop$Manual_Annot, temp$Manual_Annot),"Manual_Annot_Index"]
saveRDS(prop,"Myeloid_Figure1C_Post_RNR_Prop.RDS")

prop <- readRDS("Myeloid_Figure1C_Post_RNR_Prop.RDS")

somePDFPath = paste(cdir, "Myeloid_Figure1C_Post_RNR_Barplot.pdf", sep = "")
pdf(file=somePDFPath, width=14, height=6,pointsize=12)
ggplot(prop, aes(Cancer, Proportion, fill = Manual_Annot, label = Manual_Annot_Index))+
  geom_bar(position = "stack", stat = "identity")+
  facet_wrap(~IIID, ncol = 4)+
  # geom_text_repel(size = 4, position = position_stack(vjust = 0), lineheight = 0.7, segment.alpha = 0.3)+
  scale_fill_manual(values = ccols)+theme_linedraw(base_size = 20)+
  RotatedAxis()+theme(legend.position = "bottom", legend.title= element_blank()) +
  guides(color=guide_legend(ncol=2))+xlab("")
dev.off()

################################################################################
################################################################################
# Find pseudotime pathways (GSEA)

data <- NULL
data <- readRDS("Myeloid_Pan_Cancer_Post_Filtered.RDS")

files <- NULL
files <- list.files("Pseudotime_Myeloid/", pattern = "RDS", full.names = T)

for(i in 1:length(files)){
  print(i)
  cname <- NULL
  cname <- gsub(".*Myeloid_State_DEGs_(.*).RDS","\\1",files[i], ignore.case = T)
  m <- NULL
  m <- readRDS(files[i])
  x <- NULL
  x <- subset(data, cells = row.names(m))
  x$State <- m[match(colnames(x), row.names(m)),"State"]
  Idents(x) <- "State"
  cmarkers <- NULL
  cmarkers <- FindAllMarkers(x, min.pct = 0.25)
  cmarkers <- cmarkers[which(cmarkers$p_val_adj < 0.05),]
  saveRDS(cmarkers, paste("Myeloid_State_DEGs_",cname, ".RDS", sep = ""))
}

files <- NULL
files <- list.files(pattern = "Myeloid_State_DEGs_")
files <- files[grep("Top200", files, ignore.case = T, invert = T)]

library(msigdbr)
library(fgsea)
msigref <- msigdbr(species = "Homo sapiens", category = "H")
msigref <- msigref %>% split(x = .$gene_symbol, f = .$gs_name)

goall <- NULL
keggall <- NULL
gseaall <- NULL

for(i in 1:length(files)){
  print(i)
  x <- NULL
  x <- readRDS(files[i])
  cname <- NULL
  cname <- gsub(".*Myeloid_State_DEGs_(.*)_CytoTraceRoot_Pseudotime_(.*).RDS","\\1_\\2",files[i], ignore.case = T)
  
  cclust <- NULL
  cclust <- sort(unique(x$cluster))
  table(x$cluster)
  
  for(j in 1:length(cclust)){
    print(j)
    r <- NULL
    r <- x[which(x$cluster == cclust[j]),]
    r <- r[grep("RPL[0-9]+|^RPS|^MT-", r$gene, ignore.case = T, invert = T),] # 
    r <- r[order(r$avg_log2FC, decreasing = T),]
    
    crank <- NULL
    crank <- r$avg_log2FC
    names(crank) <- r$gene
    Rhallmark <- NULL
    Rhallmark <- fgseaMultilevel(msigref, stats = crank, nPermSimple = 1000, eps = 0)
    if(nrow(Rhallmark) > 0){
      gseaall <- rbind(gseaall, data.frame(Group = cname, State = cclust[j], Rhallmark))
    }
    
    r <- r[which(r$avg_log2FC > 0),]
    ctop <- NULL
    ctop <- data.frame(gene=unique(r$gene)[1:ifelse(nrow(r)>200,200,nrow(r))])
    entrez_out <- NULL
    entrez_out <- AnnotationDbi::select(org.Hs.eg.db, keys = ctop$gene, keytype = 'SYMBOL', columns = 'ENTREZID')
    ctop$ENTREZID <- entrez_out[match(ctop$gene, entrez_out$SYMBOL),"ENTREZID"]
    ctop <- ctop[which(!is.na(ctop$ENTREZID)),]
    
    Rkegga <- NULL
    Rkegga <- kegga(unique(ctop$ENTREZID), species="Hs")
    Rkegga <- Rkegga[order(Rkegga$P.DE, decreasing = F),]
    keggall <- rbind(keggall, data.frame(Group = cname, State = cclust[j], Rkegga))
    
    Rgoana <- NULL
    Rgoana <- goana(unique(ctop$ENTREZID), species="Hs")
    Rgoana <- Rgoana[order(Rgoana$P.DE, decreasing = F),]
    goall <- rbind(goall, data.frame(Group = cname, State = cclust[j], Rgoana))
    
  }
  
}

saveRDS(goall, paste("Myeloid_GO_Terms_Top200_State_DEGs_CytoTraceRoot_Pseudotime.RDS", sep = ""))
saveRDS(keggall, paste("Myeloid_KEGG_Top200_State_DEGs_CytoTraceRoot_Pseudotime.RDS", sep = ""))
saveRDS(gseaall, paste("Myeloid_GSEA_Top200_State_DEGs_CytoTraceRoot_Pseudotime.RDS", sep = ""))

################################################################################
################################################################################
# Heatmap regulons (R vs NR different ones)
files <- NULL
files <- list.files(path = "Pseudotime_Myeloid/", pattern = "Sig_Genes", full.names = T)
files <- files[grep("Responder", files, ignore.case = T)]

library(ComplexHeatmap)
library(rcartocolor)

cname <- NULL
cname <- gsub(".*\\/.*Mono2_(.*)_Sig_Genes_Pseudotime_(.*).RDS","\\1",files, ignore.case = T)
files <- split(files, cname)

current <- NULL

for(i in 1:length(files)){
  print(i)
  subx <- NULL
  for(j in 1:length(files[[i]])){
    x <- NULL
    x <- readRDS(files[[i]][[j]])
    x$BH <- p.adjust(x$pval, method = "BH")
    x$NegLog10FDR <- -log10(x$padjust)
    x$Label <- paste(x$Celltype, gsub("","",x$ID), sep = ":")
    x$Label <- gsub("(.*)(CD4|CD8)_Post_(.*)","\\1Post(\\3)",x$Label)
    x <- x[which(x$padjust < 0.01),]
    subx <- rbind(subx, x)
  }
  current <- rbind(current, subx)
}

saveRDS(current, "Figure3D_Regulon_Summary_Statistics_Selected_Celltypes.RDS")

current <- readRDS("Figure3D_Regulon_Summary_Statistics_Selected_Celltypes.RDS")

cgroups <- NULL
cgroups <- split(unique(current[,c("Label","Celltype")])$Label, unique(current[,c("Label","Celltype")])$Celltype)
cgenes <- NULL
for(j in 1:length(cgroups)){
  subx <- NULL
  subx <- current[which(current$Label %in% cgroups[[j]]),]
  subx <- dcast(subx, Label ~ Marker, value.var = "Estimate")
  row.names(subx) <- subx$Label
  subx <- subx[,-1]
  subx[is.na(subx)] <- 0
  
  cthreshold <- NULL
  cthreshold <- sign(subx)
  cthreshold <- which(!colSums(cthreshold) %in% c(-2,2))
  cgenes <- c(cgenes, colnames(subx)[cthreshold])
}

plotx <- NULL
plotx <- current # [which(current$Marker %in% unique(cgenes)),]
cgenesestimate <- NULL
cgenesestimate <- plotx[which(abs(plotx$Estimate) > 0.5),]
cgenes <- unique(c(cgenes, cgenesestimate$Marker))
# plotx <- plotx[which(plotx$Marker %in% unique(cgenes$Marker)),]

plotx <- dcast(plotx, Label ~ Marker, value.var = "Estimate")
row.names(plotx) <- plotx$Label
plotx <- plotx[,-1]
plotx[is.na(plotx)] <- 0
Heatmap(plotx[,cgenes], col = rev(color_conditions$Spectral), split = gsub("(.*):.*","\\1",row.names(plotx)))
scaledplotx <- NULL
scaledplotx <- scale(plotx)
scaledplotx <- t(scale(t(scaledplotx)))
Heatmap(scaledplotx[,cgenes], col = rev(color_conditions$Spectral), split = gsub("(.*):.*","\\1",row.names(scaledplotx)))

ctcols <- readRDS("Settings/Final_Celltype_Cols_Pan_Cancer.RDS")
ctcols <- ctcols[which(names(ctcols) %in% unique(current$Celltype))]

ctcols <- NULL
ctcols <- randomcoloR::distinctColorPalette(k = length(unique(current$Celltype)))
names(ctcols) <- unique(current$Celltype)

ca = rowAnnotation(show_legend = c(T),
                   Cell_Type = gsub("(.*):.*","\\1",row.names((scaledplotx)), ignore.case = T),
                   col = list(
                     Cell_Type = ctcols[gsub("(.*):.*","\\1",row.names((scaledplotx)), ignore.case = T)]))

somePDFPath = paste(cdir, "Figure3D_Regulon_Summary_Statistics_Selected_Celltypes.pdf", sep = "")
pdf(file=somePDFPath, width=18, height=4,pointsize=12)
ht <- Heatmap((scaledplotx[,cgenes]),
              column_title = "Top Pseudotime Differentiated TFs",
              column_title_gp = gpar(face = 2),
              name = "Estimate",
              # top_annotation = ca,
              left_annotation = ca,
              row_split = gsub("(.*):.*","\\1",row.names(scaledplotx[,cgenes])),
              column_names_rot = 45,
              col = rev(color_conditions$RedYellowBlue),
              border = "white",
              row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 6),
              row_title_gp = gpar(fontsize = 0),
              row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
              show_column_dend = TRUE,
              show_row_dend = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE)
draw(ht, padding = unit(c(2, 30, 2, 2), "mm"))
dev.off()

################################################################################
################################################################################
# Common regulon

allsummary <- readRDS("Myeloid_Figure3_All_Groups_Regulon_Activity_Difference_Celltypes.RDS")
allsummary$Denominator_Activity <- NULL
allsummary$Numerator_Activity <- NULL
allsummary$ID <- paste(allsummary$Cancer, allsummary$Regulon, allsummary$Manual_Annot, allsummary$Comparison, sep = "_")

files <- NULL
files <- list.files(pattern = "Myeloid_Figure3_.*AUC_")
for(i in 1:length(files)){
  x <- NULL
  x <- readRDS(files[i])
  x$Regulon <- gsub("(.*) \\(.*\\)$","\\1",x$Regulon, ignore.case = T)
  x$Regulon <- gsub("_extended","",x$Regulon, ignore.case = T)
  if(length(which(x$Response %in% c("Post(NR)","Post(R)"))) > 0){
    x$Comparison <- "Post(R)_vs_Post(NR)"
    x$Category <- x$Response
  }else{
    x$Comparison <- "Post_vs_Pre"
    x$Category <- x$Treatment_Group
  }
  x$ID <- paste(x$Cancer, x$Regulon, x$Manual_Annot, x$Comparison, sep = "_")
  subx <- NULL
  subx <- x[which(x$Category %in% c("Pre","Post(NR)")),]
  allsummary[which(allsummary$ID %in% subx$ID),"Denominator_Activity"] <- subx[match(allsummary[which(allsummary$ID %in% subx$ID),"ID"], subx$ID),"AUC"]
  subx <- NULL
  subx <- x[which(x$Category %in% c("Post","Post(R)")),]
  allsummary[which(allsummary$ID %in% subx$ID),"Numerator_Activity"] <- subx[match(allsummary[which(allsummary$ID %in% subx$ID),"ID"], subx$ID),"AUC"]
}
allsummary[is.na(allsummary)] <- 0
saveRDS(allsummary, "Myeloid_Figure3_All_Groups_Regulon_Activity_Difference_Celltypes.RDS")

allsummary <- readRDS("Myeloid_Figure3_All_Groups_Regulon_Activity_Difference_Celltypes.RDS")

allsummary <- allsummary[order(allsummary$AD, decreasing = T),]
allsummary$Regulon <- as.character(allsummary$Regulon)
groups <- unique(allsummary$Comparison)
# allstats <- NULL

for(i in 1:length(groups)){
  x <- NULL
  x <- allsummary[which(allsummary$Comparison == groups[i]),]
  
  cgenes <- NULL
  cgenes <- unique(x[which(abs(x$AD) > 0.1 & abs(x$log2FC_Expression) > 0.25),"Regulon"])
  x <- x[which(x$Regulon %in% cgenes),]
  
  ################################################
  # up/down in all celltypes and all cancers
  ccombis <- NULL
  ccombis <- unique(paste(x$Cancer, x$Manual_Annot, sep = "|"))
  #  | (x$AD < 0 & x$log2FC_Expression < 1)
  csum <- NULL
  csum <- data.frame(sort(table(x[which(x$AD > 0 & x$AD > 0.2 & x$log2FC_Expression > 1),"Regulon"])))
  csum <- csum[which(csum$Freq > length(ccombis)/5),]
  # csum <- csum[order(csum$Freq, decreasing = T),]
  cregs <- NULL
  if(length(csum) > 0){
    cregs <- sort(unique(csum$Var1))
  }
  csum <- NULL
  csum <- data.frame(sort(table(x[which(x$AD < 0 & abs(x$AD) > 0.2 & x$log2FC_Expression < 1),"Regulon"])))
  csum <- csum[which(csum$Freq > length(ccombis)/5),]
  # csum <- csum[order(csum$Freq, decreasing = T),]
  if(length(csum) > 0){
    if(!is.null(cregs)){
      cregs <- unique(c(cregs, sort(unique(csum$Var1))))
    }else{
      cregs <- unique(sort(unique(csum$Var1)))
    }
  }
  
  if(length(cregs) > 0){
  pdf(paste(cdir,"Figure3B_Top_Common_Regulons/PanCancer_Top_Common/Myeloid_Figure3B_",groups[i],"_Top_Common_Regulons_Across_Cancers_Across_Celltypes.pdf", sep = ""), width = 18, height = 16, pointsize = 12)
  for(j in 1:length(cregs)){
    plotx <- NULL
    plotx <- x[which(x$Regulon == cregs[j]),]
    plotx <- melt(plotx, id.vars = c("ID","Comparison","Cancer","Manual_Annot"), measure.vars = c("Denominator_Activity","Numerator_Activity","Denominator_Expression","Numerator_Expression"))
    plotx$Group <- gsub(".*(Expression|Activity)","\\1",plotx$variable, ignore.case = T)
    plotx$Category <- as.character(plotx$variable)
    plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Pre"
    plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post"
    plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Post(NR)"
    plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post(R)"
    plotx$Category <- factor(plotx$Category, levels = c("Pre","Post","Post(NR)","Post(R)"))
    
    p <- NULL
    p <- ggplot(plotx, aes(Category, value)) +
      # geom_boxplot(aes(x = rep(c(-3, 3), each = 10), group = group), fill = 'steelblue') +
      geom_point(aes(fill = Manual_Annot),color = "black", shape = 21, size = 5, alpha = 0.8) +
      geom_line(aes(group = ID, color = Manual_Annot))+
      stat_compare_means(method = "t.test", paired = T)+
      facet_wrap(~Cancer+Group, scales = "free")+
      # scale_shape_manual(values = 1:length(unique(plotx$Cancer)))+
      theme_classic(base_size = 15)+scale_color_manual(values = ctcols)+
      scale_fill_manual(values = ctcols)+
      theme(plot.margin = unit(c(2,2,2,2), "cm"),
            axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.3, "cm"),
            strip.text.x = element_text(size = 15),
            plot.title = element_text(size = 20, face=2, hjust = 0.5),
            strip.background = element_rect(colour="white", fill="white"),
            legend.position = "right")+
      xlab("")+ylab("Expression/Activity")+
      guides(fill=guide_legend(title="Cell Type"), color=guide_legend(title="Cell Type"))+
      ggtitle(paste(unique(plotx$Comparison),":",cregs[j], sep = "")) # +
      # geom_text(data = ctest,  mapping = aes(x = -Inf, y = -Inf, label = Label),
      #           hjust   = -0.1,
      #           vjust   = -1)
    
    print(p)
    
  }
    dev.off()
  }
  
  ################################################
  # up/down in all celltypes in each cancer
  cancers <- NULL
  cancers <- sort(unique(x$Cancer))
  for(k in 1:length(cancers)){
  subx <- NULL
  subx <- x[which(x$Cancer == cancers[k]),]
  ccombis <- NULL
  ccombis <- unique(paste(subx$Manual_Annot, sep = "|"))
  #  | (x$AD < 0 & x$log2FC_Expression < 1)
  csum <- NULL
  csum <- data.frame(sort(table(subx[which(subx$AD > 0 & subx$AD > 0.2 & subx$log2FC_Expression > 1),"Regulon"])))
  csum <- csum[which(csum$Freq > length(ccombis)/2),]
  # csum <- csum[order(csum$Freq, decreasing = T),]
  cregs <- NULL
  if(length(csum) > 0){
    cregs <- sort(unique(csum$Var1))
  }
  csum <- NULL
  csum <- data.frame(sort(table(subx[which(subx$AD < 0 & abs(subx$AD) > 0.2 & subx$log2FC_Expression < 1),"Regulon"]))) # quantile(subx$AD, probs = 0.75)
  csum <- csum[which(csum$Freq > length(ccombis)/2),]
  # csum <- csum[order(csum$Freq, decreasing = T),]
  if(length(csum) > 0){
    if(!is.null(cregs)){
      cregs <- unique(c(cregs, sort(unique(csum$Var1))))
    }else{
      cregs <- unique(sort(unique(csum$Var1)))
    }
  }
  
  if(length(cregs) > 0){
    pdf(paste(cdir,"Figure3B_Top_Common_Regulons/CancerSpecific_Top_Common/Myeloid_Figure3B_",cancers[k],"_Specific_",groups[i],"_Top_Common_Regulons_Across_Celltypes.pdf", sep = ""), width = 12, height = 6, pointsize = 12)
    for(j in 1:length(cregs)){
      plotx <- NULL
      plotx <- subx[which(subx$Regulon == cregs[j]),]
      plotx <- melt(plotx, id.vars = c("ID","Comparison","Cancer","Manual_Annot"), measure.vars = c("Denominator_Activity","Numerator_Activity","Denominator_Expression","Numerator_Expression"))
      plotx$Group <- gsub(".*(Expression|Activity)","\\1",plotx$variable, ignore.case = T)
      plotx$Category <- as.character(plotx$variable)
      plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Pre"
      plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post"
      plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Post(NR)"
      plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post(R)"
      plotx$Category <- factor(plotx$Category, levels = c("Pre","Post","Post(NR)","Post(R)"))
      
      p <- NULL
      p <- ggplot(plotx, aes(Category, value)) +
        # geom_boxplot(aes(x = rep(c(-3, 3), each = 10), group = group), fill = 'steelblue') +
        geom_point(aes(fill = Manual_Annot),color = "black", shape = 21, size = 5, alpha = 0.8) +
        geom_line(aes(group = ID, color = Manual_Annot))+
        stat_compare_means(method = "t.test", paired = T)+
        facet_wrap(~Group, scales = "free")+
        # scale_shape_manual(values = 1:length(unique(plotx$Cancer)))+
        theme_classic(base_size = 15)+scale_color_manual(values = ctcols)+
        scale_fill_manual(values = ctcols)+
        theme(plot.margin = unit(c(2,2,2,2), "cm"),
              axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
              axis.text.y = element_text(size = 20),
              axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
              legend.title = element_text(size = 25),
              legend.text = element_text(size = 10),
              legend.key.size = unit(0.3, "cm"),
              strip.text.x = element_text(size = 15),
              plot.title = element_text(size = 20, face=2, hjust = 0.5),
              strip.background = element_rect(colour="white", fill="white"),
              legend.position = "right")+
        xlab("")+ylab("Expression/Activity")+
        guides(fill=guide_legend(title="Cell Type"), color=guide_legend(title="Cell Type"))+
        ggtitle(paste(cancers[k],"(",unique(plotx$Comparison),"):",cregs[j], sep = "")) # +
      # geom_text(data = ctest,  mapping = aes(x = -Inf, y = -Inf, label = Label),
      #           hjust   = -0.1,
      #           vjust   = -1)
      
      print(p)
      
    }
    dev.off()
  }
  }
  
  ################################################
  # up/down in one celltype across all cancers
  celltypes <- NULL
  celltypes <- sort(unique(x$Manual_Annot))
  for(k in 1:length(celltypes)){
    subx <- NULL
    subx <- x[which(x$Manual_Annot == celltypes[k]),]
    ccombis <- NULL
    ccombis <- unique(paste(subx$Cancer, sep = "|"))
    #  | (x$AD < 0 & x$log2FC_Expression < 1)
    csum <- NULL
    csum <- data.frame(sort(table(subx[which(subx$AD > 0 & subx$AD > 0.2 & subx$log2FC_Expression > 1),"Regulon"])))
    csum <- csum[which(csum$Freq > length(ccombis)/4),]
    # csum <- csum[order(csum$Freq, decreasing = T),]
    cregs <- NULL
    if(length(csum) > 0){
      cregs <- sort(unique(csum$Var1))
    }
    csum <- NULL
    csum <- data.frame(sort(table(subx[which(subx$AD < 0 & abs(subx$AD) > 0.2 & subx$log2FC_Expression < 1),"Regulon"])))
    csum <- csum[which(csum$Freq > length(ccombis)/4),]
    # csum <- csum[order(csum$Freq, decreasing = T),]
    if(length(csum) > 0){
      if(!is.null(cregs)){
        cregs <- unique(c(cregs, sort(unique(csum$Var1))))
      }else{
        cregs <- unique(sort(unique(csum$Var1)))
      }
    }
    
    if(length(cregs) > 0){
      pdf(paste(cdir,"Figure3B_Top_Common_Regulons/CelltypeSpecific_Top_Common/Myeloid_Figure3B_",celltypes[k],"_Specific_",groups[i],"_Top_Common_Regulons_Across_Celltypes.pdf", sep = ""), width = 12, height = 6, pointsize = 12)
      for(j in 1:length(cregs)){
        plotx <- NULL
        plotx <- subx[which(subx$Regulon == cregs[j]),]
        plotx <- melt(plotx, id.vars = c("ID","Comparison","Cancer","Manual_Annot"), measure.vars = c("Denominator_Activity","Numerator_Activity","Denominator_Expression","Numerator_Expression"))
        plotx$Group <- gsub(".*(Expression|Activity)","\\1",plotx$variable, ignore.case = T)
        plotx$Category <- as.character(plotx$variable)
        plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Pre"
        plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post"
        plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Post(NR)"
        plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post(R)"
        plotx$Category <- factor(plotx$Category, levels = c("Pre","Post","Post(NR)","Post(R)"))
        
        p <- NULL
        p <- ggplot(plotx, aes(Category, value)) +
          # geom_boxplot(aes(x = rep(c(-3, 3), each = 10), group = group), fill = 'steelblue') +
          geom_point(aes(fill = Cancer),color = "black", shape = 21, size = 5, alpha = 0.8) +
          geom_line(aes(group = ID, color = Cancer))+
          stat_compare_means(method = "t.test", paired = T)+
          facet_wrap(~Group, scales = "free")+
          # scale_shape_manual(values = 1:length(unique(plotx$Cancer)))+
          theme_classic(base_size = 15)+scale_color_manual(values = cancercols)+
          scale_fill_manual(values = cancercols)+
          theme(plot.margin = unit(c(2,2,2,2), "cm"),
                axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
                axis.text.y = element_text(size = 20),
                axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
                axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
                legend.title = element_text(size = 25),
                legend.text = element_text(size = 10),
                legend.key.size = unit(0.3, "cm"),
                strip.text.x = element_text(size = 15),
                plot.title = element_text(size = 20, face=2, hjust = 0.5),
                strip.background = element_rect(colour="white", fill="white"),
                legend.position = "right")+
          xlab("")+ylab("Expression/Activity")+
          guides(fill=guide_legend(title="Cell Type"), color=guide_legend(title="Cell Type"))+
          ggtitle(paste(celltypes[k],"(",unique(plotx$Comparison),"):",cregs[j], sep = "")) # +
        # geom_text(data = ctest,  mapping = aes(x = -Inf, y = -Inf, label = Label),
        #           hjust   = -0.1,
        #           vjust   = -1)
        
        print(p)
        
      }
      dev.off()
    }
  }
}

################################################################################
################################################################################
# Figure 2D: Common regulon (Within Each Lineage)

library(ggpubr)

ctcols <- readRDS("Settings/Myeloid_Celltype_Cols_Pan_Cancer.RDS")

allsummary <- readRDS("Myeloid_Figure3_All_Groups_Regulon_Activity_Difference_Celltypes.RDS")
allsummary$ID <- paste(allsummary$Cancer, allsummary$Regulon, allsummary$Manual_Annot, allsummary$Comparison, sep = "_")
allsummary$Manual_Annot_Large <- as.character(allsummary$Manual_Annot)
table(allsummary$Manual_Annot_Large)

cts <- unique(allsummary$Manual_Annot_Large)
groups <- unique(allsummary$Comparison)
# allstats <- NULL

for(i in 1:length(groups)){
  for(cct in 1:length(cts)){
    x <- NULL
    x <- allsummary[which(allsummary$Comparison == groups[i] & allsummary$Manual_Annot_Large == cts[cct]),]
    
    ################################################
    # up/down in each celltype categories across all cancers
    ccombis <- NULL
    ccombis <- unique(paste(x$Cancer, x$Manual_Annot, sep = "|"))
    #  | (x$AD < 0 & x$log2FC_Expression < 1)
    csum <- NULL
    csum <- data.frame(sort(table(x[which(x$AD > 0),"Regulon"]))) # & x$log2FC_Expression > 0
    print(length(ccombis))
    csum <- csum[which(csum$Freq > length(ccombis)*0.75),]
    # csum <- csum[order(csum$Freq, decreasing = T),]
    cregs <- NULL
    if(length(csum) > 0){
      cregs <- sort(unique(csum$Var1))
    }
    csum <- NULL
    csum <- data.frame(sort(table(x[which(x$AD < 0),"Regulon"]))) #  & x$log2FC_Expression < 0
    csum <- csum[which(csum$Freq > length(ccombis)*0.75),]
    # csum <- csum[order(csum$Freq, decreasing = T),]
    if(length(csum) > 0){
      if(!is.null(cregs)){
        cregs <- unique(c(cregs, sort(unique(csum$Var1))))
      }else{
        cregs <- unique(sort(unique(csum$Var1)))
      }
    }
    
    print(paste(i,cct))
    print(cregs)
    
    if(length(cregs) > 0){
      pdf(paste(cdir,"Figure3B_Top_Common_Regulons/PanCancer_Lineage_Top_Common/Myeloid_Figure3B_",groups[i],"_",cts[cct],"_Lineage_Top_Common_Regulons_Across_Cancers.pdf", sep = ""), width = 18, height = 16, pointsize = 12)
      for(j in 1:length(cregs)){
        plotx <- NULL
        plotx <- x[which(x$Regulon == cregs[j]),]
        plotx <- melt(plotx, id.vars = c("ID","Comparison","Cancer","Manual_Annot"), measure.vars = c("Denominator_Activity","Numerator_Activity","Denominator_Expression","Numerator_Expression"))
        plotx$Group <- gsub(".*(Expression|Activity)","\\1",plotx$variable, ignore.case = T)
        plotx$Category <- as.character(plotx$variable)
        plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Pre"
        plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post"
        plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Post(NR)"
        plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post(R)"
        plotx$Category <- factor(plotx$Category, levels = c("Pre","Post","Post(NR)","Post(R)"))
        plotx[which(plotx$Group == "Expression"),"value"] <- log2(plotx[which(plotx$Group == "Expression"),"value"])
        p <- NULL
        p <- ggplot(plotx, aes(Category, value)) +
          # geom_boxplot(aes(x = rep(c(-3, 3), each = 10), group = group), fill = 'steelblue') +
          geom_point(aes(fill = Manual_Annot),color = "black", shape = 21, size = 5, alpha = 0.8) +
          geom_line(aes(group = ID, color = Manual_Annot))+
          stat_compare_means(method = "t.test", paired = T)+
          facet_wrap(~Cancer+Group, scales = "free")+
          # scale_shape_manual(values = 1:length(unique(plotx$Cancer)))+
          theme_classic(base_size = 15)+scale_color_manual(values = ctcols)+
          scale_fill_manual(values = ctcols)+
          theme(plot.margin = unit(c(2,2,2,2), "cm"),
                axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
                axis.text.y = element_text(size = 20),
                axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
                axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
                legend.title = element_text(size = 25),
                legend.text = element_text(size = 10),
                legend.key.size = unit(0.3, "cm"),
                strip.text.x = element_text(size = 15),
                plot.title = element_text(size = 20, face=2, hjust = 0.5),
                strip.background = element_rect(colour="white", fill="white"),
                legend.position = "right")+
          xlab("")+ylab("log2Expression/Activity")+
          guides(fill=guide_legend(title="Cell Type"), color=guide_legend(title="Cell Type"))+
          ggtitle(paste(cts[cct],"\n",unique(plotx$Comparison),":",cregs[j], sep = "")) # +
        # geom_text(data = ctest,  mapping = aes(x = -Inf, y = -Inf, label = Label),
        #           hjust   = -0.1,
        #           vjust   = -1)
        
        print(p)
        
      }
      dev.off()
    }
    
  }
}

selected <- c("Macro_FOLR2-APOE+","Macro_FOLR2+APOE-","Macro_FOLR2+APOE+","cDC3_LAMP3","Mono_INHBA","Macro_IFI27")

allsummary <- allsummary[which(allsummary$Manual_Annot %in% selected),]

for(i in 1:length(groups)){
  for(cct in 1:length(selected)){
    x <- NULL
    x <- allsummary[which(allsummary$Comparison == groups[i] & allsummary$Manual_Annot == selected[cct]),]
    
    ################################################
    # up/down in each celltype categories across all cancers
    ccombis <- NULL
    ccombis <- unique(paste(x$Cancer, x$Manual_Annot, sep = "|"))
    #  | (x$AD < 0 & x$log2FC_Expression < 1)
    csum <- NULL
    csum <- data.frame(sort(table(x[which(x$AD > 0),"Regulon"]))) #  & x$log2FC_Expression > 0
    csum <- csum[which(csum$Freq == length(ccombis)),]
    # csum <- csum[order(csum$Freq, decreasing = T),]
    cregs <- NULL
    if(length(csum) > 0){
      cregs <- sort(unique(csum$Var1))
    }
    csum <- NULL
    csum <- data.frame(sort(table(x[which(x$AD < 0),"Regulon"]))) #  & x$log2FC_Expression < 0
    csum <- csum[which(csum$Freq == length(ccombis)),]
    # csum <- csum[order(csum$Freq, decreasing = T),]
    if(length(csum) > 0){
      if(!is.null(cregs)){
        cregs <- unique(c(cregs, sort(unique(csum$Var1))))
      }else{
        cregs <- unique(sort(unique(csum$Var1)))
      }
    }
    
    print(paste(i,cct))
    print(cregs)
    
    if(length(cregs) > 0){
      pdf(paste(cdir,"Figure3B_Top_Common_Regulons/PanCancer_Lineage_Top_Common/Myeloid_Figure3B_",groups[i],"_",selected[cct],"_Selected_Celltypes_Top_Common_Regulons_Across_Cancers.pdf", sep = ""), width = 12, height = 6, pointsize = 12)
      for(j in 1:length(cregs)){
        plotx <- NULL
        plotx <- x[which(x$Regulon == cregs[j]),]
        plotx <- melt(plotx, id.vars = c("ID","Comparison","Cancer","Manual_Annot"), measure.vars = c("Denominator_Activity","Numerator_Activity","Denominator_Expression","Numerator_Expression"))
        plotx$Group <- gsub(".*(Expression|Activity)","\\1",plotx$variable, ignore.case = T)
        plotx$Category <- as.character(plotx$variable)
        plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Pre"
        plotx$Category[which(plotx$Comparison == "Post_vs_Pre" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post"
        plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Denominator_Activity","Denominator_Expression"))] <- "Post(NR)"
        plotx$Category[which(plotx$Comparison == "Post(R)_vs_Post(NR)" & plotx$Category %in% c("Numerator_Activity","Numerator_Expression"))] <- "Post(R)"
        plotx$Category <- factor(plotx$Category, levels = c("Pre","Post","Post(NR)","Post(R)"))
        plotx[which(plotx$Group == "Expression"),"value"] <- log2(plotx[which(plotx$Group == "Expression"),"value"]+1)
        p <- NULL
        p <- ggplot(plotx, aes(Category, value)) +
          # geom_boxplot(aes(x = rep(c(-3, 3), each = 10), group = group), fill = 'steelblue') +
          geom_point(aes(fill = Cancer),color = "black", shape = 21, size = 5, alpha = 0.8) +
          geom_line(aes(group = ID, color = Cancer))+
          facet_wrap(~Group, scales = "free")+
          stat_compare_means(method = "t.test", paired = T)+
          # scale_shape_manual(values = 1:length(unique(plotx$Cancer)))+
          theme_classic(base_size = 15)+scale_color_manual(values = cancercols)+
          scale_fill_manual(values = cancercols)+
          theme(plot.margin = unit(c(2,2,2,2), "cm"),
                axis.text.x = element_text(angle = 45,hjust=1,vjust = 1,size = 20),
                axis.text.y = element_text(size = 20),
                axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
                axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
                legend.title = element_text(size = 25),
                legend.text = element_text(size = 10),
                legend.key.size = unit(0.3, "cm"),
                strip.text.x = element_text(size = 15),
                plot.title = element_text(size = 20, face=2, hjust = 0.5),
                strip.background = element_rect(colour="white", fill="white"),
                legend.position = "right")+
          xlab("")+ylab("log2Expression/Activity")+
          guides(fill=guide_legend(title="Cell Type"), color=guide_legend(title="Cell Type"))+
          ggtitle(paste(selected[cct],"\n",unique(plotx$Comparison),":",cregs[j], sep = "")) # +
        # geom_text(data = ctest,  mapping = aes(x = -Inf, y = -Inf, label = Label),
        #           hjust   = -0.1,
        #           vjust   = -1)
        
        print(p)
        
      }
      dev.off()
    }
    
  }
}
################################################################################
################################################################################
# Cell-Cell (Myeloid vs Myeloid)

library(liana)
files <- NULL
files <- list.files(path = "CellCell/Myeloid/Myeloid_Myeloid/", full.names = T)
for(i in 1:length(files)){
  x <- NULL
  x <- readRDS(files[i])
  x <- liana_aggregate(x)
  cname <- NULL
  cname <- gsub(".*\\/(.*)_Myeloid_(Pre|Post)_Liana_Output_(.*).RDS","\\1_\\2_\\3",files[i])
  cname <- gsub(".*\\/(.*)_Myeloid_(Pre|Post)_Liana_Output.RDS","\\1_\\2",cname, ignore.case = T)
  cname <- gsub("Post_Post","Post_",cname, ignore.case = T)
  if(length(grep("RvsNR", cname, ignore.case = T)) > 0){
    x <- x[(grepl("_Responder", x$source, ignore.case = T) & grepl("_Responder", x$target, ignore.case = T)) | (grepl("_Non-responder", x$source, ignore.case = T) & grepl("_Non-responder", x$target, ignore.case = T)),]
    x$Label <- gsub(".*(Responder|Non-Responder).*","Post(\\1)",x$source, ignore.case = T)
    x$Label <- paste(gsub("^(.*?)_.*","\\1",cname, ignore.case = T), ifelse(x$Label == "Post(Responder)", "Post(R)", "Post(NR)"), sep = "_")
  }else{
    x$Label <- cname
  }
  x$source <- gsub("^.*?_(.*)","\\1",x$source, ignore.case = T)
  x$target <- gsub("^.*?_(.*)","\\1",x$target, ignore.case = T)
  x$source <- gsub("(.*)_(Pre|Post)$","\\1",x$source, ignore.case = T)
  x$target <- gsub("(.*)_(Pre|Post)$","\\1",x$target, ignore.case = T)
  x$source <- gsub("(.*)_(Responder|Non-responder)$","\\1",x$source, ignore.case = T)
  x$target <- gsub("(.*)_(Responder|Non-responder)$","\\1",x$target, ignore.case = T)
  
  # x <- filter(x, aggregate_rank < 0.01)
  saveRDS(x, paste("Myeloid_Figure3A_CellCell_",cname,".RDS", sep = ""))
  
  labels <- NULL
  labels <- unique(x$Label)
  
  for(j in 1:length(labels)){
    somePDFPath = paste(cdir, "Myeloid_Figure3A_Myeloid_Myeloid_CellCell/CellCell_Top_Interactions/Myeloid_Figure3A_CellCell_",labels[j],"_",cname,"_Top_20_Interactions.pdf", sep = "")
    pdf(file=somePDFPath, width=20, height=8,pointsize=12)
    plotx <- NULL
    plotx <- filter(x[which(x$Label == labels[j]),], aggregate_rank < 0.01)
    cts <- NULL
    cts <- sort(unique(plotx$source))
    for(k in 1:length(cts)){
      if(nrow(plotx[which(plotx$source == cts[k]),]) > 0){
        p <- NULL
        p <- liana_dotplot(plotx, source_groups = cts[k],target_groups = unique(c(plotx$source,plotx$target)),ntop=20)+RotatedAxis()+theme(axis.text.x=element_text(colour="black", size = 18, face = "plain"), plot.title = element_text(hjust = 0.5, size = 30))+scale_color_gradientn(colours = color_conditions$gradient)+ggtitle(gsub("_",":",labels[j], ignore.case = T))
        print(p)
      }
    }
    dev.off()
  }
}

################################################################################
################################################################################
# CellCell (Myeloid vs Myeloid Network Frequency for Each Cell Type)
library(liana)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

ctcols <- NULL
ctcols <- readRDS("Settings/Myeloid_Celltype_Cols_Pan_Cancer.RDS")
# names(ctcols) <- gsub(".*?_(.*)","\\1",names(ctcols), ignore.case = T)

dir.create(paste(cdir,"Myeloid_Figure3A_Myeloid_Myeloid_CellCell/CellCell_Quantile75_Interactions_Network_Test/", sep = ""))

allstats <- readRDS("Myeloid_CellCell_Quantile75_Interactions_Network_Freq.RDS")
allstats <- split(allstats, allstats$Group)
allstats <- lapply(allstats, function(x){
  x$Proportion <- x$Freq/sum(x$Freq)
  return(x)
})
allstats <- do.call(rbind.data.frame, allstats)

cname <- NULL
cname <- gsub("NonSpecific_","",allstats$Group)
cname <- gsub("Post_R|Post_NR","Post_R_vs_NR",cname)
cname <- gsub("Post$|Pre","Post_vs_Pre",cname)

groups <- NULL
groups <- split(allstats$Group, cname)
freqsummary <- NULL
for(i in 1:length(groups)){
  x <- NULL
  x <- allstats[which(allstats$Group %in% groups[[i]]),]
  ccts <- NULL
  ccts <- sort(unique(c(x$source, x$target)))
  for(j in 1:length(ccts)){
    y <- NULL
    y <- x[which(x$source == ccts[j] | x$target == ccts[j]),]
    crows <- NULL
    crows <- which(y$target == ccts[j])
    y[crows,"target"] <- y[crows,"source"]
    y[crows,"source"] <- ccts[j]
    y$ID <- paste(y$source, y$target, sep = "_and_")
    # y <- plyr::count(y, vars = c("Group","source","target"), wt_var = "Freq")
    
    cids <- NULL
    cids <- unique(y$ID)
    for(k in 1:length(cids)){
      suby <- NULL
      suby <- y[which(y$ID == cids[k]),]
      suby <- data.frame(Comparison = names(groups)[i],
                         source = unique(suby$source),
                         target = unique(suby$target),
                         FC_Numbers = (sum(suby[which(suby$Group %in% c("NonSpecific_Post_R","NonSpecific_Post")),"Freq"])+1)/(sum(suby[which(suby$Group %in% c("NonSpecific_Post_NR","NonSpecific_Pre")),"Freq"])+1),
                         FC_Prop = (sum(suby[which(suby$Group %in% c("NonSpecific_Post_R","NonSpecific_Post")),"Proportion"])+0.0001)/(sum(suby[which(suby$Group %in% c("NonSpecific_Post_NR","NonSpecific_Pre")),"Proportion"])+0.0001),
                         Freq_Denominator = (sum(suby[which(suby$Group %in% c("NonSpecific_Post_NR","NonSpecific_Pre")),"Freq"])),
                         Freq_Numerator = (sum(suby[which(suby$Group %in% c("NonSpecific_Post_R","NonSpecific_Post")),"Freq"])),
                         Prop_Denominator = (sum(suby[which(suby$Group %in% c("NonSpecific_Post_NR","NonSpecific_Pre")),"Proportion"])),
                         Prop_Numerator = (sum(suby[which(suby$Group %in% c("NonSpecific_Post_R","NonSpecific_Post")),"Proportion"])))
      freqsummary <- rbind(freqsummary, suby)
    }
  }
}

saveRDS(freqsummary,"Myeloid_Myeloid_CellCell_Celltype_Specific_Freq_FC.RDS")

freqsummary <- readRDS("Myeloid_Myeloid_CellCell_Celltype_Specific_Freq_FC.RDS")

plotx <- NULL
plotx <- freqsummary[which(freqsummary$Comparison == "Post_R_vs_NR"),]
plotx$log10FC <- log10(plotx$FC_Numbers)
p <- NULL
p <- ggplot(plotx, aes(log10FC, source, size = abs(log10FC), color = target))+
  geom_point(alpha = 0.8)+
  theme_linedraw(base_size = 15)+
  scale_color_manual(values = ctcols)+
  scale_size_continuous(range = c(1,6))+
  geom_vline(xintercept = 0.5, size = 1, linetype = "dotted")+
  geom_vline(xintercept = -0.5, size = 1, linetype = "dotted")+
  ggtitle("Myeloid vs Myeloid\nCell-Cell Cell Type Specific Interactions Fold-Change")+
  guides(color = guide_legend(title = "Cell Type"))

somePDFPath = paste(cdir, "Myeloid_Figure3A_Myeloid_Myeloid_CellCell/CellCell_Quantile75_Interactions_Network_Test/Myeloid_Myeloid_Cell_Cell_Cell_Type_Specific_Frequency_ScatterPlot.pdf", sep = "")
pdf(file=somePDFPath, width=12, height=6,pointsize=12)
print(p)
dev.off()

################################################################################
################################################################################
# CellCell (Myeloid vs T Network Frequency for Each Cell Type)

library(liana)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

ctcols <- NULL
ctcols <- readRDS("Settings/Myeloid_Celltype_Cols_Pan_Cancer.RDS")
ctcols <- c(ctcols, readRDS("Settings/Final_Celltype_Cols_Pan_Cancer.RDS"))
# names(ctcols) <- gsub(".*?_(.*)","\\1",names(ctcols), ignore.case = T)

dir.create(paste(cdir,"Myeloid_Figure3A_Myeloid_T_CellCell/CellCell_Quantile75_Interactions_Network_Test/", sep = ""))

allstats <- readRDS("Myeloid_T_CellCell_Quantile75_Interactions_Network_Freq.RDS")
allstats <- split(allstats, allstats$Group)
allstats <- lapply(allstats, function(x){
  x$Proportion <- x$Freq/sum(x$Freq)
  return(x)
})
allstats <- do.call(rbind.data.frame, allstats)

cname <- NULL
cname <- gsub("NonSpecific_","",allstats$Group)
cname <- gsub("Post_R|Post_NR","Post_R_vs_NR",cname)
cname <- gsub("Post$|Pre","Post_vs_Pre",cname)

groups <- NULL
groups <- split(allstats$Group, cname)
freqsummary <- NULL
for(i in 1:length(groups)){
  x <- NULL
  x <- allstats[which(allstats$Group %in% groups[[i]]),]
  ccts <- NULL
  ccts <- sort(unique(c(x$source, x$target)))
  for(j in 1:length(ccts)){
    y <- NULL
    y <- x[which(x$source == ccts[j] | x$target == ccts[j]),]
    crows <- NULL
    crows <- which(y$target == ccts[j])
    y[crows,"target"] <- y[crows,"source"]
    y[crows,"source"] <- ccts[j]
    y$ID <- paste(y$source, y$target, sep = "_and_")
    # y <- plyr::count(y, vars = c("Group","source","target"), wt_var = "Freq")
    
    cids <- NULL
    cids <- unique(y$ID)
    for(k in 1:length(cids)){
      suby <- NULL
      suby <- y[which(y$ID == cids[k]),]
      suby <- data.frame(Comparison = names(groups)[i],
                         source = unique(suby$source),
                         target = unique(suby$target),
                         FC_Numbers = (sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_R","NonSpecific_Myeloid_vs_CD8_Post_R","NonSpecific_Myeloid_vs_CD4_Post","NonSpecific_Myeloid_vs_CD8_Post")),"Freq"])+1)/(sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_NR","NonSpecific_Myeloid_vs_CD8_Post_NR","NonSpecific_Myeloid_vs_CD4_Pre","NonSpecific_Myeloid_vs_CD8_Pre")),"Freq"])+1),
                         FC_Prop = (sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_R","NonSpecific_Myeloid_vs_CD8_Post_R","NonSpecific_Myeloid_vs_CD4_Post","NonSpecific_Myeloid_vs_CD8_Post")),"Proportion"])+0.0001)/(sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_NR","NonSpecific_Myeloid_vs_CD8_Post_NR","NonSpecific_Myeloid_vs_CD4_Pre","NonSpecific_Myeloid_vs_CD8_Pre")),"Proportion"])+0.0001),
                         Freq_Denominator = (sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_NR","NonSpecific_Myeloid_vs_CD8_Post_NR","NonSpecific_Myeloid_vs_CD4_Pre","NonSpecific_Myeloid_vs_CD8_Pre")),"Freq"])),
                         Freq_Numerator = (sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_R","NonSpecific_Myeloid_vs_CD8_Post_R","NonSpecific_Myeloid_vs_CD4_Post","NonSpecific_Myeloid_vs_CD8_Post")),"Freq"])),
                         Prop_Denominator = (sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_NR","NonSpecific_Myeloid_vs_CD8_Post_NR","NonSpecific_Myeloid_vs_CD4_Pre","NonSpecific_Myeloid_vs_CD8_Pre")),"Proportion"])),
                         Prop_Numerator = (sum(suby[which(suby$Group %in% c("NonSpecific_Myeloid_vs_CD4_Post_R","NonSpecific_Myeloid_vs_CD8_Post_R","NonSpecific_Myeloid_vs_CD4_Post","NonSpecific_Myeloid_vs_CD8_Post")),"Proportion"])))
      freqsummary <- rbind(freqsummary, suby)
    }
  }
}

saveRDS(freqsummary,"Myeloid_T_CellCell_Celltype_Specific_Freq_FC.RDS")

freqsummary <- readRDS("Myeloid_T_CellCell_Celltype_Specific_Freq_FC.RDS")
freqsummary <- freqsummary[grep("CD4|CD8",freqsummary$source, ignore.case = T, invert = T),]
plotx <- NULL
plotx <- freqsummary[which(freqsummary$Comparison %in% c("Myeloid_vs_CD4_Post_R_vs_NR","Myeloid_vs_CD8_Post_R_vs_NR")),]
plotx$log10FC <- log10(plotx$FC_Numbers)
p <- NULL
p <- ggplot(plotx, aes(log10FC, source, size = abs(log10FC), color = target))+
  geom_point(alpha = 0.8)+
  theme_linedraw(base_size = 15)+
  scale_color_manual(values = ctcols)+
  scale_size_continuous(range = c(1,6))+
  geom_vline(xintercept = 0.5, size = 1, linetype = "dotted")+
  geom_vline(xintercept = -0.5, size = 1, linetype = "dotted")+
  facet_wrap(~Comparison)+
  ggtitle("Myeloid vs Myeloid\nCell-Cell Cell Type Specific Interactions Fold-Change")+
  guides(color = guide_legend(title = "Cell Type", ncol = 1))+
  theme(legend.key.size = unit(0.1, "cm"))

somePDFPath = paste(cdir, "Myeloid_Figure3A_Myeloid_T_CellCell/CellCell_Quantile75_Interactions_Network_Test/Myeloid_T_Cell_Cell_Cell_Type_Specific_Frequency_ScatterPlot.pdf", sep = "")
pdf(file=somePDFPath, width=14, height=7,pointsize=12)
print(p)
dev.off()

################################################################################
################################################################################
# CellCell CRC Redefine Response (Myeloid vs T Network)
library(liana)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

ctcols <- NULL
ctcols <- readRDS("Settings/Myeloid_Celltype_Cols_Pan_Cancer.RDS")
ctcols <- c(ctcols, readRDS("Settings/Final_Celltype_Cols_Pan_Cancer.RDS"))

dir.create(paste(cdir,"Myeloid_Figure3A_Myeloid_T_CellCell_Individual/CellCell_Quantile75_Interactions_Network/", sep = ""), recursive = T)

files <- NULL
files <- list.files(path = "CellCell/Myeloid/Myeloid_T_CellCell_Individual/", full.names = T)

for(i in 1:length(files)){
  print(i)
  x <- NULL
  x <- readRDS(files[i])
  x <- liana_aggregate(x)
  cname <- NULL
  cname <- gsub(".*\\/(.*)_(Myeloid_vs_.*)_(Pre|Post)_Liana_Output_(.*).RDS","\\1_\\2",files[i])
  cname <- paste(cname, "_", unique(cmeta$Response[match(gsub("CRC_(.*)_Myeloid.*","\\1",cname, ignore.case = T), cmeta$SampleID)]), sep = "")
  
  x <- x[(grepl("Mono|Macro|DC|Mast", x$source, ignore.case = T) & grepl("CD4|CD8", x$target, ignore.case = T)) | (grepl("CD4|CD8", x$source, ignore.case = T) & grepl("Mono|Macro|DC|Mast", x$target, ignore.case = T)),]
  
  saveRDS(x, paste("Myeloid_Figure3A_CellCell_Myeloid_T_Individuals_",cname,".RDS", sep = ""))
}

files <- NULL
files <- list.files(pattern = "Myeloid_Figure3A_CellCell_Myeloid_T_Individuals_", full.names = T)

for(k in 1:length(files)){
  print(k)
  cfiles <- NULL
  cfiles <- files[k]
  cgroup <- NULL
  cgroup <- gsub(".*Myeloid_Figure3A_CellCell_Myeloid_T_Individuals_(.*)\\.RDS","\\1",files[k], ignore.case = T)
  x <- NULL
  x <- readRDS(cfiles)
  x <- filter(x, aggregate_rank < 0.01)
  
  somePDFPath = paste(cdir, "Myeloid_Figure3A_Myeloid_T_CellCell_Individual/CellCell_Quantile75_Interactions_Network/Myeloid_Figure3A_",cgroup,"_Myeloid_T_CellCell_Individual_CellCell_Quantile90_Interactions_Network.pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=8,pointsize=12)
  
  mynet <- NULL
  mynet <- data.frame(table(x[,c("source","target")]))
  mynet <- mynet[which(mynet$Freq > quantile(mynet$Freq, probs = 0.50)),]
  mynet$source <- as.character(mynet$source)
  mynet$target <- as.character(mynet$target)
  
  net <- NULL
  net<- graph_from_data_frame(mynet)
  deg <- NULL
  deg <- igraph::degree(net, mode="all")
  coords <- NULL
  coords <- layout_in_circle(net) # , order = order(membership(karate_groups)))
  
  ccols <- NULL
  ccols <- ctcols[names(ctcols) %in% unique(c(mynet$source, mynet$target))]
  cmiss <- NULL
  cmiss <- unique(c(mynet$source, mynet$target))[which(!unique(c(mynet$source, mynet$target)) %in% names(ctcols))]
  if(length(cmiss) > 0){
    cmisscols <- NULL
    cmisscols <- gen_colors(n = length(cmiss))
    names(cmisscols) <- cmiss
    ccols <- c(ccols, cmisscols)
  }
  ccols <- ccols[match(names(V(net)),names(ccols))]
  
  E(net)$width  <- E(net)$Freq*0.2
  for (m in 1:length(unique(mynet$source))){
    E(net)[purrr::map(unique(mynet$source),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$source)[m],x))
    })%>% unlist()]$color <- adjustcolor(ccols[which(names(ccols) == unique(mynet$source)[m])], alpha.f = .9)
  }
  
  plot(net, edge.arrow.size=1, 
       vertex.size=20,
       edge.curved=0.3,
       vertex.color= ccols[match(names(V(net)),names(ccols))],
       vertex.shape="sphere", # "sphere", "circle"
       vertex.frame.color="#555555",
       vertex.label.color=ccols[match(names(V(net)),names(ccols))],
       layout = coords,
       vertex.label.family="Helvetica",
       vertex.label.degree=0,
       vertex.label.dist= 2, # c(8,8,5,5,-5,-5,-5,-5,-8,-5,-6,-5,5,6,5),
       vertex.label.cex=1)
  title(gsub("_",":",cgroup, ignore.case = T),cex.main=2,col.main="black")
  dev.off()
}

################################################################################
################################################################################
# CellCell Modify CRC (Myeloid vs T Network)
library(liana)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

ctcols <- NULL
ctcols <- readRDS("Settings/Myeloid_Celltype_Cols_Pan_Cancer.RDS")
ctcols <- c(ctcols, readRDS("Settings/Final_Celltype_Cols_Pan_Cancer.RDS"))

dir.create(paste(cdir,"Myeloid_Figure3A_Myeloid_T_Modify_CRC/CellCell_Quantile75_Interactions_Network/", sep = ""), recursive = T)

files <- NULL
files <- list.files(path = "CellCell/Myeloid/Myeloid_T_Modify_CRC/", full.names = T)
# files <- files[grep("NonSpecific_Myeloid_vs_CD.*_Melanoma|_Melanoma_Post_Liana_Output_RvsNR", files, ignore.case = T, invert = T)]

for(i in 1:length(files)){
  print(i)
  x <- NULL
  x <- readRDS(files[i])
  x <- liana_aggregate(x)
  cname <- NULL
  cname <- gsub(".*\\/(.*)_(Myeloid_vs_.*)_(Pre|Post)_Liana_Output_(.*).RDS","\\1_\\2_\\3_\\4",files[i])
  cname <- gsub(".*\\/(.*)_(Myeloid_vs_.*)_(Pre|Post)_Liana_Output.RDS","\\1_\\2_\\3",cname)
  
  x <- x[(grepl("Mono|Macro|DC|Mast", x$source, ignore.case = T) & grepl("CD4|CD8", x$target, ignore.case = T)) | (grepl("CD4|CD8", x$source, ignore.case = T) & grepl("Mono|Macro|DC|Mast", x$target, ignore.case = T)),]
  
  saveRDS(x, paste("Myeloid_Figure3A_CellCell_Myeloid_T_Modify_CRC_",cname,".RDS", sep = ""))
}

files <- NULL
files <- list.files(pattern = "Myeloid_Figure3A_CellCell_Myeloid_T_Modify_CRC_", full.names = T)

for(k in 1:length(files)){
  print(k)
  cfiles <- NULL
  cfiles <- files[k]
  cgroup <- NULL
  cgroup <- gsub(".*Myeloid_Figure3A_CellCell_Myeloid_T_Modify_CRC_(.*)\\.RDS","\\1",files[k], ignore.case = T)
  x <- NULL
  x <- readRDS(cfiles)
  x <- filter(x, aggregate_rank < 0.01)
  
  somePDFPath = paste(cdir, "Myeloid_Figure3A_Myeloid_T_Modify_CRC/CellCell_Quantile75_Interactions_Network/Myeloid_Figure3A_",cgroup,"_Myeloid_T_Modify_CRC_CellCell_Quantile90_Interactions_Network.pdf", sep = "")
  pdf(file=somePDFPath, width=10, height=8,pointsize=12)
  
  mynet <- NULL
  mynet <- data.frame(table(x[,c("source","target")]))
  mynet <- mynet[which(mynet$Freq > quantile(mynet$Freq, probs = 0.90)),]
  mynet$source <- as.character(mynet$source)
  mynet$target <- as.character(mynet$target)
  
  net <- NULL
  net<- graph_from_data_frame(mynet)
  deg <- NULL
  deg <- igraph::degree(net, mode="all")
  coords <- NULL
  coords <- layout_in_circle(net) # , order = order(membership(karate_groups)))
  
  ccols <- NULL
  ccols <- ctcols[names(ctcols) %in% unique(c(mynet$source, mynet$target))]
  cmiss <- NULL
  cmiss <- unique(c(mynet$source, mynet$target))[which(!unique(c(mynet$source, mynet$target)) %in% names(ctcols))]
  if(length(cmiss) > 0){
    cmisscols <- NULL
    cmisscols <- gen_colors(n = length(cmiss))
    names(cmisscols) <- cmiss
    ccols <- c(ccols, cmisscols)
  }
  ccols <- ccols[match(names(V(net)),names(ccols))]
  
  E(net)$width  <- E(net)$Freq*0.2
  for (m in 1:length(unique(mynet$source))){
    E(net)[purrr::map(unique(mynet$source),function(x) {
      get.edge.ids(net,vp = c(unique(mynet$source)[m],x))
    })%>% unlist()]$color <- adjustcolor(ccols[which(names(ccols) == unique(mynet$source)[m])], alpha.f = .9)
  }
  
  plot(net, edge.arrow.size=1, 
       vertex.size=20,
       edge.curved=0.3,
       vertex.color= ccols[match(names(V(net)),names(ccols))],
       vertex.shape="sphere", # "sphere", "circle"
       vertex.frame.color="#555555",
       vertex.label.color=ccols[match(names(V(net)),names(ccols))],
       layout = coords,
       vertex.label.family="Helvetica",
       vertex.label.degree=0,
       vertex.label.dist= 2, # c(8,8,5,5,-5,-5,-5,-5,-8,-5,-6,-5,5,6,5),
       vertex.label.cex=1)
  title(gsub("_",":",cgroup, ignore.case = T),cex.main=2,col.main="black")
  dev.off()
}

################################################################################
################################################################################
# Cell-Cell (Myeloid vs T) Replot receptor-ligand plots

library(liana)

files <- NULL
files <- list.files(path = "CellCell/Myeloid/Myeloid_T/", full.names = T)
cname <- NULL
cname <- gsub(".*\\/(.*)_(Myeloid_vs_.*)_(Pre|Post)_Liana_Output_(.*).RDS","\\1_\\2_\\3_\\4",files)
cname <- gsub(".*\\/(.*)_(Myeloid_vs_.*)_(Pre|Post)_Liana_Output.RDS","\\1_\\2_\\3",cname, ignore.case = T)
cname <- gsub("Post_Post","Post_",cname, ignore.case = T)
cname <- gsub("(.*_Myeloid_vs_CD.*_Post_)(NR|R)","\\1RvsNR",cname, ignore.case = T)
cname <- gsub("(.*)(Pre|Post)$","\\1PrevsPost",cname, ignore.case = T)
cname <- gsub("(Melanoma_Myeloid_vs_CD.*)_Melanoma(_PrevsPost)","\\1\\2",cname, ignore.case = T)
cname <- gsub("RvsNRvsNR","RvsNR",cname, ignore.case = T)
files <- split(files, cname)
dir.create(paste(cdir,"Myeloid_Figure3A_Myeloid_T_CellCell/CellCell_Top_Interactions_Replot/", sep = ""))
dir.create(paste("Myeloid_Myeloid_vs_T_CellCell_Top_Interactions_Replot/", sep = ""))

for(i in 1:length(files)){
  gc()
  print(i)
  x <- NULL
  for(j in 1:length(files[[i]])){
    subx <- NULL
    subx <- readRDS(files[[i]][[j]])
    subx <- liana_aggregate(subx)
    cname <- NULL
    cname <- gsub(".*\\/(.*)_(Myeloid_vs_.*)_(Pre|Post)_Liana_Output_(.*).RDS","\\1_\\2_\\3_\\4",files[[i]][[j]])
    cname <- gsub(".*\\/(.*)_(Myeloid_vs_.*)_(Pre|Post)_Liana_Output.RDS","\\1_\\2_\\3",cname, ignore.case = T)
    cname <- gsub("Post_Post","Post_",cname, ignore.case = T)
    # cname <- gsub("Non-responder_Post_RvsNR","Post_NR",cname, ignore.case = T)
    # cname <- gsub("Responder_Post_RvsNR","Post_R",cname, ignore.case = T)
    if(length(grep("RvsNR|PostNR|PostR|Post_R|Post_NR", cname, ignore.case = T)) > 0){
      if(length(grep("Responder", subx$source, ignore.case = T)) > 0){
        subx <- subx[(grepl("_Responder", subx$source, ignore.case = T) & grepl("_Responder", subx$target, ignore.case = T)) | (grepl("_Non-responder", subx$source, ignore.case = T) & grepl("_Non-responder", subx$target, ignore.case = T)),]
        subx$Label <- gsub(".*(Responder|Non-Responder).*","Post(\\1)",subx$source, ignore.case = T)
        subx$Label <- paste(cname, ifelse(subx$Label == "Post(Responder)", "Post(R)", "Post(NR)"), sep = "_")
      }else{
        # subx$Label <- gsub(".*Liana_Output_Post(R|NR).*","Post(\\1)",files[[i]][[j]], ignore.case = T)
        subx$Label <- gsub(".*\\/(.*_Myeloid_vs_CD.*?)_Post_Liana_Output_Post(R|NR).*","\\1_Post(\\2)",files[[i]][[j]], ignore.case = T)
      }
      
      subx$Label <- gsub("Post_RvsNR_","",subx$Label, ignore.case = T)
      subx$Label <- gsub("Post_Post","Post",subx$Label, ignore.case = T)
      # subx$Label <- gsub("(.*_Myeloid_vs_CD.*)_Post_(R|NR)$","\\1_Post(\\2)",subx$Label, ignore.case = T)
      
      subx$source <- gsub("(.*)(_Responder|_Non-Responder)(.*)","\\1",subx$source, ignore.case = T)
      subx$target <- gsub("(.*)(_Responder|_Non-Responder)(.*)","\\1",subx$target, ignore.case = T)
      subx$source <- gsub("(Responder_|Non-Responder_)(.*)","\\1",subx$source, ignore.case = T)
      subx$target <- gsub("(Responder_|Non-Responder_)(.*)","\\1",subx$target, ignore.case = T)
    }else{
      subx$Label <- cname
      if(length(grep("_Pre|_Post", subx$source, ignore.case = T)) > 0){
        subx <- subx[(grepl("_Pre", subx$source, ignore.case = T) & grepl("_Pre", subx$target, ignore.case = T)) | (grepl("_Post", subx$source, ignore.case = T) & grepl("_Post", subx$target, ignore.case = T)),]
      }
    }
    subx <- subx[(grepl("Mono|Macro|DC|Mast", subx$source, ignore.case = T) & grepl("CD4|CD8", subx$target, ignore.case = T)) | (grepl("CD4|CD8", subx$source, ignore.case = T) & grepl("Mono|Macro|DC|Mast", subx$target, ignore.case = T)),]
    
    subx$source <- gsub("(Melanoma_|Breast Cancer_|BCC_|ccRCC_|CRC_|HCC_|iCCA_|HNSCC_)(.*)","\\2",subx$source, ignore.case = T)
    subx$target <- gsub("(Melanoma_|Breast Cancer_|BCC_|ccRCC_|CRC_|HCC_|iCCA_|HNSCC_)(.*)","\\2",subx$target, ignore.case = T)
    
    subx <- filter(subx, aggregate_rank < 0.01)
    x <- rbind(x, subx)
  }
  
  print(table(x$Label))
  print(cname)
  
  if(length(files[[i]]) > 1){
    cname <- NULL
    cname <- names(files)[i]
    # cname <- paste(unique(x$Label), collapse = "_vs_")
    # cname <- gsub("(.*)_Post_NR_vs_.*_Post_R","\\1_Post_RvsNR",cname, ignore.case = T)
    # cname <- gsub("(.*)_Post_vs_.*_Pre","\\1_PrevsPost",cname, ignore.case = T)
  }
  
  cout <- NULL
  cout <- list(x)
  names(cout) <- cname
  saveRDS(cout, paste("Myeloid_Myeloid_vs_T_CellCell_Top_Interactions_Replot/Myeloid_Figure3A_CellCell_",cname,"_Myeloid_T_Top_5_Interactions.RDS", sep = ""))
  
  print(head(cout, n=2))
  print(table(x$Label))
  
}

files <- NULL
files <- list.files(path = "Myeloid_Myeloid_vs_T_CellCell_Top_Interactions_Replot", pattern = "Myeloid_Figure3A_CellCell_", full.names = T)

for(i in 1:length(files)){
  
  print(i)
  
  x <- NULL
  x <- readRDS(files[i])
  cname <- NULL
  cname <- names(x)
  x <- x[[1]]
  
  cts <- NULL
  cts <- unique(sort(c(x$source,x$target)))
  cts <- cts[grep("Macro|Mono|DC|Mast", cts, ignore.case = T)]
  
  somePDFPath = paste(cdir, "Myeloid_Figure3A_Myeloid_T_CellCell/CellCell_Top_Interactions_Replot/Myeloid_Figure3A_CellCell_",cname,"_Myeloid_T_Top_5_Interactions.pdf", sep = "")
  pdf(file=somePDFPath, width=16, height=11,pointsize=12)
  
  for(k in 1:length(cts)){
    creceptor <- NULL
    cligand <- NULL
    p1 <- NULL
    p2 <- NULL
    
    if(nrow(x[which(x$source == cts[k]),]) > 0){
      cligand <- x[which(x$source == cts[k]),]
      cligand <- split(cligand, cligand$Label)
      cligand <- lapply(cligand, function(x){
        x <- x[order(x$aggregate_rank, decreasing = F),]
        x <- split(x, x$target)
        x <- lapply(x, function(y){
          cgenes <- NULL
          cgenes <- paste(y$ligand.complex, y$receptor.complex, sep = "_")
          cgenes <- cgenes[1:ifelse(length(cgenes) > 5, 5, length(cgenes))]
          y <- y[which(paste(y$ligand.complex, y$receptor.complex, sep = "_") %in% cgenes),]
          return(y)
        })
        x <- do.call(rbind.data.frame, x)
        return(x)
      })
      cligand <- do.call(rbind.data.frame, cligand)
      cligand$source <- gsub(".*?_(.*)","\\1",cligand$Label)
      p1 <- liana_dotplot(cligand, source_groups = unique(cligand$source),target_groups = unique(cligand$target),ntop=length(unique(paste(cligand$ligand.complex, cligand$receptor.complex))))+RotatedAxis()+theme(axis.text.x=element_text(colour="black", size = 20, face = "plain"), axis.text.y=element_text(colour="black", size = 12, face = "plain"), plot.title = element_text(hjust = 0.5, size = 30))+scale_color_gradientn(colours = color_conditions$gradient)+ggtitle(paste(unique(gsub("(.*)_.*","\\1",cligand$Label)), " - ", cts[k], "\n(As Ligand)", sep = ""))
      print(p1)
    }
    
    if(nrow(x[which(x$target == cts[k]),]) > 0){
      creceptor <- x[which(x$target == cts[k]),]
      creceptor <- split(creceptor, creceptor$Label)
      creceptor <- lapply(creceptor, function(x){
        x <- x[order(x$aggregate_rank, decreasing = F),]
        x <- split(x, x$source)
        x <- lapply(x, function(y){
          cgenes <- NULL
          cgenes <- paste(y$ligand.complex, y$receptor.complex, sep = "_")
          cgenes <- cgenes[1:ifelse(length(cgenes) > 5, 5, length(cgenes))]
          y <- y[which(paste(y$ligand.complex, y$receptor.complex, sep = "_") %in% cgenes),]
          return(y)
        })
        x <- do.call(rbind.data.frame, x)
        return(x)
      })
      creceptor <- do.call(rbind.data.frame, creceptor)
      creceptor$target <- paste(creceptor$target, ":", gsub(".*?_(.*)","\\1",creceptor$Label), sep = "")
      creceptor$target <- gsub(".*?_(.*)","\\1",creceptor$Label)
      creceptor$source_actual <- creceptor$source
      creceptor$target_actual <- creceptor$target
      creceptor$source <- creceptor$target_actual
      creceptor$target <- creceptor$source_actual
      
      p2 <- liana_dotplot(creceptor, source_groups = unique(creceptor$source),target_groups = unique(creceptor$target),ntop=length(unique(paste(creceptor$ligand.complex, creceptor$receptor.complex))))+RotatedAxis()+theme(axis.text.x=element_text(colour="black", size = 20, face = "plain"), axis.text.y=element_text(colour="black", size = 12, face = "plain"), plot.title = element_text(hjust = 0.5, size = 30))+scale_color_gradientn(colours = color_conditions$gradient)+ggtitle(paste(unique(gsub("(.*)_.*","\\1",cligand$Label)), " - ", cts[k], "\n(As Receptor)", sep = ""))
      print(p2)
    }
  }
  dev.off()
}


########################################################################################################
########################################################################################################
# Figure 3B: GSEA
files <- list.files(path = "GSEA/GSEA/", pattern = "RDS", full.names = T, recursive = T)
files <- files[grep("Myeloid", files, ignore.case = T)]
allgsea <- NULL
allgo <- NULL
allkegg <- NULL

for(i in 1:length(files)){
  print(paste("i:",i,",",files[i], sep =))
  x <- NULL
  x <- readRDS(files[i])
  cname <- NULL
  cname <- gsub(".*\\/(.*)_Cluster.*Pan_Cancer_.*.RDS","\\1",files[i], ignore.case = T)
  if(!is.null(x)){
    if(length(grep("GSEA_Selected", files[i], ignore.case = T)) > 0){
      x$BH <- p.adjust(x$pval, method = "BH")
      x <- x[which(x$pval < 0.05),]
      allgsea <- rbind(allgsea, x)
    }else if(length(grep("_GO_", files[i], ignore.case = T)) > 0){
      csubgroup <- NULL
      if(length(grep("_Down_Regulated_GO_", files[i], ignore.case = T)) > 0){
        csubgroup <- "Down"
      }else if(length(grep("_Up_Regulated_GO_", files[i], ignore.case = T)) > 0){
        csubgroup <- "Up"
      }
      x$BH <- p.adjust(x$P.DE, method = "BH")
      x <- x[which(x$P.DE < 0.05),]
      allgo <- rbind(allgo, data.frame(Regulation = csubgroup, x))
    }else if(length(grep("_KEGG_", files[i], ignore.case = T)) > 0){
      csubgroup <- NULL
      if(length(grep("_Down_Regulated_KEGG_", files[i], ignore.case = T)) > 0){
        csubgroup <- "Down"
      }else if(length(grep("_Up_Regulated_KEGG_", files[i], ignore.case = T)) > 0){
        csubgroup <- "Up"
      }
      x$BH <- p.adjust(x$P.DE, method = "BH")
      x <- x[which(x$P.DE < 0.05),]
      allkegg <- rbind(allkegg, data.frame(Regulation = csubgroup, x))
    }
  }
}

saveRDS(allgsea, "Myeloid_Figure3B_GSEA_Combined_BH005.RDS")
saveRDS(allgo, "Myeloid_Figure3B_GO_Combined_BH005.RDS")
saveRDS(allkegg, "Myeloid_Figure3B_KEGG_Combined_BH005.RDS")

allgsea <- readRDS("Myeloid_Figure3B_GSEA_Combined_BH005.RDS")
allgo <- readRDS("Myeloid_Figure3B_GO_Combined_BH005.RDS")
allkegg <- readRDS("Myeloid_Figure3B_KEGG_Combined_BH005.RDS")

