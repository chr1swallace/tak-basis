#!/usr/bin/env Rscript

## basis p values
## remotes::install_github('ollyburren/cupcake')

## library(randomFunctions)
library(cluster)
library(cowplot)
library(cupcake)
library(data.table)
library(dendextend)
library(factoextra) # get_dist
library(ggplot2)
library(gridGraphics) # recordPlot
library(magrittr)
library(parallel)
library(pheatmap)
source("~/A/R/cw-utils.R")

## load tak data
tak <- fread("~/MetalResults_bases.txt")
tak[,pid:=paste(CHR,POS,sep=":")]
takp1 <- fread("~/tak-projection.csv")

## tak v
btak <- with(tak, log(OR) * shrinkage[pid] - beta.centers[pid])
ptak <- btak %*% rot.pca[tak$pid, ]
v <- with(tak, SE * shrinkage[pid] * rot.pca[pid, ])
vtak <- t(v) %*% LD.int[tak$pid, tak$pid] %*% v

################################################################################

## load finngen projections
g <- fread("/home/cew54/share/Data/Tidy_GWAS-summary/Projecting_dbs/Projected_table_finnmas_20200403.tsv")
qc <- fread("/home/cew54/share/Data/Tidy_GWAS-summary/Projecting_dbs/QC_table_finnmas_20200403.tsv")
head(g)
grep("K11_IBD",g$Trait,value=TRUE) %>% unique()
grep("IBD",g$Trait,value=TRUE) %>% unique()

g0 <- fread("/home/cew54/share/Data/Tidy_GWAS-summary/Projecting_dbs/Projected_table_20200403.tsv")
grep("K11",g0$Trait,value=TRUE)  %>% unique()

qc0 <- fread("/home/cew54/share/Data/Tidy_GWAS-summary/Projecting_dbs/QC_table_20200403.tsv")
g0 <- g0[grep("^CEL_Finn|^IBD_Finn",Trait)]
table(g0$Trait)
g0 <- merge(g0, qc0[,.(Trait,overall_p)],by="Trait")
g0[,Trait:=paste0("K11_",Trait)]
g <- merge(g, qc[,.(Trait,overall_p)],by="Trait")
g <- rbind(g,g0)

setnames(g,
         c("Trait","Delta","Var.Delta"),
         c("trait","delta","var.proj"))
g <- g[grep("FinnGen",trait)]
g$category <- g$category.label <- "FinnGen"
g$trait.label <- g$trait <- sub("_FinnGen.*","",g$trait)
g$fdr.overall <- p.adjust(g$overall_p,method="BH")
## cat(trans$BP38,
## ibd <- fread("~/finngen_r2_K11_IBD.gz")
## head(ibd)
grep("K11_IBD",g$trait,value=TRUE) %>% unique()
grep("IBD",g$trait,value=TRUE) %>% unique()
## g <- g[grep("K11",trait)]

## keep disease traits
g[grep("^[D-S][0-9][0-9]?_",trait)]$trait %>%unique()  %>% sample(., 100)
g <- g[grep("^[D-S][0-9][0-9]?_",trait)]

################################################################################

## plot with other traits
takp1$category <- "TAK"
takp1$category.label <- "Sawalha"
takp1$trait <- takp1$trait.label <- "TAK"
takp1$fdrcat <- "general"
names(takp1)

################################################################################

## load burren data
setwd("~/A")
source("~/A/R/cw-reader.R")
source("~/A/R/cw-colours2.R")
source("~/A/R/cw-palette.R")
source("~/A/R/cw-renamer.R")
sparse=reader()
names(sparse)

proj <- copy(sparse)
g$fdrcat <- "general"
proj[grepl("addison",trait),trait:="addisons disease"]# too many characters
proj[grepl("addison",trait),trait.label:="addisons disease"]# too many characters
proj[,trait.label:=sub("UKBB ","",trait.label)]
proj <- rbind(proj,takp1,g, fill=TRUE)
proj[,lab:=paste(category.label,trait.label) %>% sub("^ ","",.)]
proj[,newfdr:=p.adjust(p.value,method="BH"),by=c("category","PC")]
table(proj$category)
tail(proj)
proj[trait.label=="TAK"]
fwrite(proj,file="~/tak+fg+ukbb+.csv")

################################################################################

