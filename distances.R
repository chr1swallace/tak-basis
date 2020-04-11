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
# A is root of https://github.com/chr1swallace/imd-basis
source("~/A/R/cw-utils.R")

x=fread("~/tak+fg+ukbb+.csv")
tak <- fread("~/MetalResults_bases.txt")
tak[,pid:=paste(CHR,POS,sep=":")]

## calculate distances

todo <- setdiff(x$trait,c("TAK",NA))

ttest <- function(x,y,vx,vy) {
  vp <- (vx + vy)/2
  cbind(d=(x-y), v=vp)
}

delta <- dt2mat(x[!is.na(trait)], PC ~ trait, fun=mean, value.var="delta")
V <- dt2mat(x[!is.na(trait)], PC ~ trait, fun=mean, value.var="var.proj")

pids <- tak$pid
v <- shrinkage[pids] * rot.pca[pids, ]
W <- t(v) %*% LD[pids, pids] %*% v

## mahalanobis
delta <- dt2mat(x[!is.na(trait)], PC ~ trait, fun=mean, value.var="delta")
D <- lapply(todo, function(tr) {
  data.table(trait=tr,
             DW=sqrt(t(delta[,"TAK"] - delta[,tr]) %*% solve(W) %*% (delta[,"TAK"] - delta[,tr]))[1,1],
             EUC=sqrt(sum( (delta[,"TAK"] - delta[,tr])^2 ))
             )
}) %>% rbindlist()
D[,cl:=ifelse(trait %in% x[fdr.overall < 0.01]$trait, "<0.01",">=0.01")]
D$cl %<>% factor(., levels=c(">=0.01","<0.01"))
D <- merge(D, unique(x[,.(trait,trait.label,category.label)]),by="trait")
head(D[category.label=="FinnGen"][order(DW)],20)
## tail(D[order(DW)])

delta <- dt2mat(x[!is.na(trait)], PC ~ trait, fun=mean, value.var="delta")
with(D, plot(EUC,DW)) ## euclidean and mahalanobis highly correlated

D <- D[category.label %in% c("UKBB","FinnGen","Geneatlas_ICD")]
D[,category.label:=as.character(category.label)]
D <- D[!grepl("KELA",trait)]
D[,trait2:=sub("  "," \n ",trait)][,trait2:=sub("UKBB_NEALE:SRD:","",trait2)]
D[,trait2:=sub("..not.elsewhere.classified",".other",trait2)]
D[DW<1]
D[grep("elsewhere",trait2)][cl=="<0.01"]

## show finngen, and geneatlas
D <- D[category.label!="UKBB"]
table(D$category.label)
D[category.label=="Geneatlas_ICD", category.label:="UK Biobank"]

pdf("~/Projects/basis-tak/mahalanobis.pdf",height=6,width=9)

theme_set(theme_cowplot())
ggplot(D, aes(x=DW,fill=cl)) +
  geom_histogram(binwidth=0.01,col="grey")+
  scale_fill_manual("IMD basis significance, FDR",values=c("<0.01"="grey20",">=0.01"="grey90")) +
  ## geom_label_repel(aes(label=trait),data=D[DW< 1.1], y=1,fill="white",nudge_y=1)
  geom_text_repel(aes(label=gsub("."," ",
                                 sub("UKBB_NEALE:SRD:","",trait2),fixed=TRUE)),
                  data=D[grepl("K11_IBD",trait) | (DW< 1.2 & cl=="<0.01")], y=1.5,
                  angle=90,
                  nudge_y=2,
                  hjust=0,vjust=0.5) +
  labs(x="Mahalanobis Distance",y="Count") +
  scale_y_continuous(expand=c(0,0)) +
  background_grid(major="y") +
  theme(legend.position="bottom",
        strip.background=element_blank(),
        strip.text=element_text(face="bold"),
        legend.justification=c("right","center"),
        axis.line.y=element_blank()) +
  facet_wrap(~category.label,ncol=3,scales="free")

dev.off()

