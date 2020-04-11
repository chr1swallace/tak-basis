#!/usr/bin/env Rscript

## basis p values
## remotes::install_github('ollyburren/cupcake')
library(data.table)
library(magrittr)
library(ggplot2)
library(cupcake)
library(parallel)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(dendextend)
library(factoextra)

setwd("~/A") # root of https://github.com/chr1swallace/imd-basis
source("R/cw-utils.R")
source("R/cw-palette.R")
source("R/cw-colours2.R") # colored_labelled_boxes

proj <- fread("~/tak+fg+ukbb+.csv")
proj[,stars:=ifelse(newfdr<0.01,"*","")]
sigs <- proj[fdrcat=="general"][
  fdr.overall <= 0.01 |
  trait %in% c("TAK","UKBB_NEALE:SRD:ankylosing.spondylitis")]$lab  %>%
  unique()
"Sawalha TAK" %in% (sigs)

proj <- proj[lab %in% sigs]
stars <- dt2mat(proj[category!="FinnGen"], trait ~ PC, value.var="stars")[,paste0("PC",1:13)]
obsData <- dt2mat(proj[category!="FinnGen"], trait ~ PC, value.var="delta")[,paste0("PC",1:13)]
obsVars <- dt2mat(proj[category!="FinnGen"], trait ~ PC, value.var="var.proj")[,paste0("PC",1:13)]^2
fdr <- dt2mat(proj[category!="FinnGen"], trait ~ PC, value.var="newfdr")[,1,drop=FALSE]
b <- obsData

renamer <- unique(proj[,.(category,trait,trait.label,category.label,lab)])
renamer[is.na(lab)| is.na(trait)]

m <- match(rownames(b), renamer$trait)
summary(m)
rownames(stars) <- rownames(b) <- sub("_combined","combined",renamer$lab[m])

## clust plot
tb <- t(b)
zp <- (tb/apply(abs(tb),1,max))  %>% t()  %>% 
    cut(. ,breaks=seq(-1,1,by=0.1),include.lowest=TRUE)  %>%  t()  %>% as.numeric()  
cols <- matrix(grnvi(21)[zp], nrow(b),ncol(b),dimnames=dimnames(b))

D <-  get_dist(b,method="euclidean")  
cl <- hclust(D,method="ward.D2")
dd <- as.dendrogram(cl)
## this group are largely null - push them towards one end of display rather than in middle (topologically this doesn't change dendrogram)
nulls <-  c("UKBB basal cell carcinoma", "UKBB malignant melanoma", 
"UKBB allergy hypersensitivity anaphylaxis", "UKBB eczema dermatitis", 
"UKBB psoriasis", "methotrexate", "UKBB high cholesterol", "UKBB hypertension", 
"UKBB diabetes", "UKBB emphysema chronic bronchitis")
dd %<>% rotate(., order=c(setdiff(labels(dd),nulls),intersect(labels(dd),nulls)))

k <- 4
col4 <- tol5qualitative[c(1,2,3,5)]
## col4 <- c("#015501","#010155","#015501","#010155")
## col4 <- rep(c("dodgerblue","grey30"),length.out=k)
col4 <- rep(c(mygreen,"grey10"),length.out=k)
pch4 <- rep(c(19,18),length.out=k)
cuts <- cutree(dd,k=k)[labels(dd)]
## colbb <- ifelse(grepl("UKBB",labels(dd)),"grey20", "steelblue") #tol5qualitative[1])
colbb <- ifelse(grepl("UKBB",labels(dd)),"grey10", mygreen) #darken(mygreen,1.2)) #tol5qualitative[1])
pchbb <- ifelse(grepl("UKBB",labels(dd)), 18,19)



pdf("~/Projects/basis-tak/cluster.pdf",height=10,width=8,pointsize=10)

par( mar = c(1,0,0,29))
## par(mar = c(28,2,1,2))
dd  %>% 
    dendextend::set("leaves_pch",pchbb)  %>%
    dendextend::set("labels_cex",1)  %>%
    dendextend::set("leaves_col",colbb)  %>% 
  dendextend::set("labels_col",colbb)  %>%
  plot(.,axes=FALSE,horiz=TRUE)
M <- ncol(b)
colored_labelled_boxes(colors = cols[labels(dd),M:1],
                       labels=stars[labels(dd),M:1],
                       dend=dd,sort_by_labels_order=FALSE,horiz=TRUE)
rectcol <- "grey30" #tol5qualitative[1]
rectx <- 0.15 #0.45 #0.25
recth <- 1 #0.85
recthdiff=-0.3
r <- rect.dendrogram(dd, k=4,border=rectcol,horiz=TRUE,prop_k_height=0.8,lty=2,
                     lower_rect=-0.61) # -1.46 for 6in wide, -0.5 for 8in wide
text(rectx,length(r[[1]])*recth + recthdiff,labels="I",col=rectcol)
text(rectx,length(r[[1]]) + length(r[[2]])*recth + recthdiff,labels="II",col=rectcol)
text(rectx,length(unlist(r[1:2])) + length(r[[3]])*recth + recthdiff,labels="III",col=rectcol)
text(rectx,length(unlist(r[1:3])) + length(r[[4]])*recth + recthdiff,labels="IV",col=rectcol)

dev.off()    
