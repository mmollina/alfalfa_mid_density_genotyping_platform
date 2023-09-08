require(updog)
require(mappoly)
require(tidyverse)

#### Genotype call: I195 x J432 ####
## Loading allele counts
setwd("~/repos/collaborations/alfalfa_map/")
full.dat <- read.table("combined_geno_F1_BC1F1_Alfalfa.txt", header = TRUE)
full.dat[1:10, 1:8]
## Gathering BC individuals
bc.ind <- grep("AphBC.", colnames(full.dat), value = TRUE)
x <- full.dat[,bc.ind, drop = FALSE]
BC <- data.frame(AlleleID = full.dat[,"AlleleID"],
                 MarkerName = full.dat[,"CloneID"],## Using 'MarkerName' to be compatible with probe info
                 AlleleSequence = full.dat[,"AlleleSequence"],
                 I195 = apply(full.dat[,c("I195.x","I195.y")], 
                              1, mean),
                 F1_85_209 = full.dat[,"X85.209"], x)
## Loading chromosome info and snp position
probe.info <- read_csv("~/repos/collaborations/alfalfa_map/20201030-BI-Alfalfa_SNPs_DArTag-probe-design.csv")
## Merge info
BC <- left_join(probe.info, BC, by = "MarkerName") %>%
  arrange(Chrom, Pos)
head(BC)[,1:17]
## Filtering alleles with less than 10 reads in both parents
BC <- BC[BC$I195 + BC$F1_85_209 > 10, ]
head(BC)[,1:4]
## Filtering alleles with a mean depth ranging from 10 to 600
BC <- BC[apply(BC[,bc.ind], 1, mean) < 600 & apply(BC[,bc.ind], 1, mean) > 10, ]
hist(apply(BC[,bc.ind], 1, mean), breaks = 50, main = "Read depth", xlab = "Number of reads")
###Get unique haplotype IDs
hap.id<-unique(BC$MarkerName)
## Organizing count information of SNPs within haplotypes and
## formatting read counts as 'updog' input file
cte <- 1
refmat <- sizemat <- matrix(NA, nrow = nrow(BC), ncol = length(bc.ind) + 2)
nms <- character(nrow(BC))
pos <- chrom <- ref <- alt <- NULL 
for(i in 1:length(hap.id)){
  cur.hap.id <- hap.id[i]
  cat(cur.hap.id, "-->", 100*round(i/length(hap.id), 2) ,"%\n   ")
  cur.alleles <- BC$AlleleSequence[BC$MarkerName==cur.hap.id]
  if(length(cur.alleles) == 1) next()
  cur.alleles <- sapply(cur.alleles, function(x) unlist(strsplit(x, split="")), USE.NAMES = FALSE)
  cur.var.pos <- which(apply(cur.alleles, 1, function(x) length(unique(x))) > 1)
  cur.alleles <- cur.alleles[cur.var.pos, , drop = F]
  Z <- apply(cur.alleles, 1, unique)
  ch <- unique(BC$Chrom[BC$MarkerName==cur.hap.id])
  ps <- unique(BC$Pos[BC$MarkerName==cur.hap.id]) - 50 + cur.var.pos
  for(j in 1:nrow(cur.alleles)){
    cat(j, " ")
    nms[cte] <- paste0(cur.hap.id, "_SNP", j)
    parents<-rbind(by(apply(BC[BC$MarkerName==cur.hap.id, grep(pattern = "I195", colnames(BC), value = TRUE), drop = FALSE], 1, sum), cur.alleles[j,], sum),
                   by(apply(BC[BC$MarkerName==cur.hap.id, grep(pattern = "F1_85_209", colnames(BC), value = TRUE), drop = FALSE], 1, sum), cur.alleles[j,], sum))
    alt <- c(alt, colnames(parents)[1])
    ref <- c(ref, colnames(parents)[2])
    chrom <- c(chrom, ch)
    pos <- c(pos, ps[j])
    offspring<-t(apply(BC[BC$MarkerName==cur.hap.id, bc.ind], 2, by, cur.alleles[j,], sum))
    X <- rbind(parents, offspring)
    refmat[cte,] <- X[,1]
    sizemat[cte,] <- apply(X, 1, sum)
    cte <- cte + 1
  }
  cat("\n")
}
## Shortening SNP names and attributing names to matrices
nms <- str_replace_all(nms, "alfalfaRep2vsXJDY1_shared", "Alf")
dimnames(sizemat) <- dimnames(refmat) <- list(nms,c("I195", "F1_85_209", rownames(offspring)))
sizemat <- sizemat[!apply(sizemat, 1, function(x) all(is.na(x))),]
refmat <- refmat[!apply(refmat, 1, function(x) all(is.na(x))),]
names(chrom) <- names(pos) <- names(ref) <- names(alt) <- nms[which(!apply(refmat, 1, function(x) all(is.na(x))))]
## Genotype calling using updog
mout = multidog(refmat = refmat, 
                sizemat = sizemat, 
                ploidy = 4, 
                model = "f1",
                p1_id = "I195",
                p2_id = "F1_85_209",
                nc = 8)
## Importing updog results to MAPpoly
ud.dat.bc <- import_from_updog(mout, prob.thres = 0.8)
ud.dat.bc$seq.ref <- ref[ud.dat.bc$mrk.names]
ud.dat.bc$seq.alt <- alt[ud.dat.bc$mrk.names]
ud.dat.bc$chrom <- sapply(strsplit(chrom[ud.dat.bc$mrk.names], "chr|.1"), function(x) as.numeric(x[2]))
ud.dat.bc$genome.pos <- pos[ud.dat.bc$mrk.names]
plot(ud.dat.bc)
print(ud.dat.bc, detailed = T)
s <- make_seq_mappoly(ud.dat.bc, arg = "all")
a <- get_genomic_order(s)
plot(a$ord$seq.pos)
plot(a)
save(ud.dat.bc, file = "bc_updog_data.rda")
