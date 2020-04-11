#!/usr/bin/env Rscript

## basis p values
## remotes::install_github('ollyburren/cupcake')
library(data.table)
library(magrittr)
library(cupcake)
source("~/A/R/cw-utils.R")

## load tak data
tak <- fread("~/MetalResults_bases.txt")
head(tak)
tak[,pid:=paste(CHR,POS,sep=":")]
takp0 <- project_sparse(beta=log(tak$OR), seb=tak$SE, pids=tak$pid)
takp0

## create LD
s.DT <- split(SNP.manifest, SNP.manifest$chr)
R <- fread("~/R_SNP_bases_ld.txt")
R[,pidA:=paste(CHR_A,BP_A,sep=":")]
R[,pidB:=paste(CHR_B,BP_B,sep=":")]
dim(R)
hist(R$BP_A-R$BP_B) # only up to 1MB. Assume r=0 outside this
plot(abs(R$BP_A-R$BP_B),R$R) # not unreasonable given this plot

LD.int <- lapply(names(s.DT), function(chr) {
  r <- matrix(0,nrow(s.DT[[chr]]),nrow(s.DT[[chr]]),
              dimnames=list(s.DT[[chr]]$pid,s.DT[[chr]]$pid))
  diag(r) <- 1
  Rsub <- R[ pidA %in% rownames(r) & pidB %in% rownames(r)]
  nrow(Rsub)
  dim(r)
  r[cbind(Rsub$pidA,Rsub$pidB)] <- Rsub$R
  r[cbind(Rsub$pidB,Rsub$pidA)] <- Rsub$R
  r
})  %>% bdiag_with_dimnames(.)
LD.int <- LD.int[rownames(rot.pca),rownames(rot.pca)]

## project with international LD
plot(as.vector(cupcake::LD),as.vector(LD.int))

## vary project_sparse function
project_sparse_int <- function (beta, seb, pids) {
    if (length(beta) != length(seb) || length(beta) != length(pids) || 
        !length(beta)) 
        stop("arguments must be equal length vectors > 0")
    if (!all(pids %in% SNP.manifest$pid)) 
        stop("all pids must be members of sparse basis (SNP.manifest$pid)")
    if (length(pids) < 0.95 * nrow(rot.pca)) 
      warning("more than 5% sparse basis snps missing:",
              100 * (1 - length(pids)/nrow(rot.pca)))
    b <- beta * shrinkage[pids] - beta.centers[pids]
    proj <- b %*% rot.pca[pids, ]
    v <- seb * shrinkage[pids] * rot.pca[pids, ]
    var.proj <- t(v) %*% LD.int[pids, pids] %*% v
    ctl <- (-beta.centers[pids]) %*% rot.pca[pids, ]
    delta <- (proj - ctl)[1, ]
    chi2 <- (t(delta) %*% solve(var.proj) %*% delta)[1, 1]
    ret <- data.table::data.table(PC = colnames(proj), proj = proj[1, 
        ], var.proj = Matrix::diag(var.proj), delta = delta, 
        p.overall = stats::pchisq(chi2, df = 13, lower.tail = FALSE))
    ret$z = ret$delta/sqrt(ret$var.proj)
    ret$p = stats::pnorm(abs(ret$z), lower.tail = FALSE) * 2
    copy(ret)
}

takp1 <- project_sparse_int(beta=log(tak$OR), seb=tak$SE, pids=tak$pid)
takp1
table(SNP.manifest$pid %in% tak$pid)

# check sparse results robust to LD
par(mfrow=c(1,2))
plot(takp0$z,takp1$z,main="z",xlab="LD",ylab="LD.int"); abline(0,1,col="red")
plot(takp0$proj,takp1$proj,main="proj",xlab="LD",ylab="LD.int"); abline(0,1,col="red")

fwrite(takp1[,.(PC,delta,var.delta=var.proj,z,p.value=p)],
       file="~/tak-projection.csv")
