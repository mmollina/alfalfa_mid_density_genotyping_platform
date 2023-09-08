require(mappoly)
require(tidyverse)
require(ggsci)
if(!require(highprecHMM)){
  devtools::install_github("mmollina/highprecHMM")
  require(highprecHMM)
}
load(file = "~/repos/collaborations/alfalfa_map/bc_map_construction/bc_map_result.rda")
load(file = "~/repos/collaborations/alfalfa_map/f1_map_construction/f1_map_result.rda")
setwd("~/repos/collaborations/alfalfa_map/joint_map")
#### Functions ####
get_map_dosage <- function(map){
  p <- abs(abs(map$info$seq.dose.p1 - map$info$ploidy/2) - map$info$ploidy/2)
  q <- abs(abs(map$info$seq.dose.p2 - map$info$ploidy/2) - map$info$ploidy/2)
  s.p <- names(which(p  ==  1 & q  ==  0))
  s.q <- names(which(p  ==  0 & q  ==  1))
  ds <- names(which(p  ==  1 & q  ==  1))
  mplx <- setdiff(map$info$mrk.names, c(s.p, s.q, ds))
  freq = c(length(s.p) + length(s.q), length(ds), length(mplx))
  names(freq) <- c("simplex", "double.simplex", "multiplex")
  list(simplex.p = s.p,
       simplex.q = s.q,
       double.simplex = ds,
       multiplex = mplx,
       freq = freq)
}
homolog_correspondence <- function(map1, map2, parent1, parent2){
  ph1 <- ph_list_to_matrix(map1$maps[[1]]$seq.ph[[parent1]], map1$info$ploidy)
  ph2 <- ph_list_to_matrix(map2$maps[[1]]$seq.ph[[parent2]], map2$info$ploidy)
  rownames(ph1) <- map1$info$mrk.names
  rownames(ph2) <- map2$info$mrk.names
  id1 <- intersect(rownames(ph1), rownames(ph2))
  id2 <- names(which(apply(ph1[id1,], 1, sum) - apply(ph2[id1,], 1, sum) == 0))
  diff.dose <- setdiff(id1, id2)
  ph1.temp <- ph1[id2,]
  ph2.temp <- ph2[id2,]
  id <- permute::allPerms(1:ncol(ph2.temp))
  phase.similarity <- apply(id, 1, function(x) sum(abs(diff(ph1.temp == ph2.temp[,x]))))
  best.conf <- which.min(phase.similarity)
  ph2 <- ph2[,id[best.conf,]]
  list(ph1 = ph1, ph2 = ph2)
}
merge_data <- function(dat, P){
  id <- expand.grid(c("dosage.p1", "dosage.p2"),
                    1:length(dat), stringsAsFactors = FALSE)[,2:1]
  snp <- unique(Reduce(union, lapply(dat, function(x) x$mrk.names)))
  chrom <- unlist(lapply(dat, function(x) x$chrom))
  pos <-  unlist(lapply(dat, function(x) x$genome.pos))
  seq.ref <- unlist(lapply(dat, function(x) x$seq.ref))
  seq.alt <- unlist(lapply(dat, function(x) x$seq.alt))
  chisq.pval <- unlist(lapply(dat, function(x) x$chisq.pval))
  up <- unique(P)
  X <- data.frame(snp = snp,
                  chrom = chrom[snp],
                  pos = pos[snp],
                  seq.ref = seq.ref[snp],
                  seq.alt = seq.alt[snp],
                  chisq.pval = chisq.pval[snp],
                  matrix(NA, nrow = length(snp), ncol = length(up), dimnames = list(snp, up)))
  for(i in up){
    v <- which(P%in%i)
    Y <- matrix(NA, nrow = length(snp), ncol = length(v), dimnames = list(snp, NULL))
    for(j in 1:length(v)){
      y <- dat[[id[v[j],1]]][[id[v[j],2]]]
      Y[names(y),j] <- y
    }
    snp.temp <- unlist(apply(Y, 1, function(x) {
      a <- na.omit(unique(x))
      if(length(a) == 1) return(a)
      else return(NULL)
    }))
    X[names(snp.temp),i] <- snp.temp
  }
  for(i in 1:length(dat)){
    d <- rownames_to_column(as.data.frame(dat[[i]]$geno.dose), "snp")

    X <- left_join(X, d, by = "snp")
  }
  return(X)
}
build_hmm_joint_map <- function(map1, map2, ped, ploidy, err = 0.05, tol1 = 1e-2, tol2= 1e-3){
  ## I195
  H1 <- homolog_correspondence(map1 = map1,
                               map2 = map2,
                               parent1 = 1,
                               parent2 = 1)
  ### Merging maps based on genome and MDS
  ch.f1 <- names(sort(table(dat.f1$chrom[rownames(H1$ph1)]), decreasing = TRUE))[1]
  ch.bc <- names(sort(table(dat.bc$chrom[rownames(H1$ph2)]), decreasing = TRUE))[1]
  if(ch.f1 != ch.bc) stop()
  dont.belong.f1 <- names(which(dat.f1$chrom[rownames(H1$ph1)] != ch.f1))
  dont.belong.bc <- names(which(dat.bc$chrom[rownames(H1$ph2)] != ch.bc))
  belong.f1 <- names(which(dat.f1$chrom[rownames(H1$ph1)] == ch.f1))
  belong.bc <- names(which(dat.bc$chrom[rownames(H1$ph2)] == ch.bc))
  mrk.pos <- sort(c(dat.f1$genome.pos[belong.f1],
                    dat.bc$genome.pos[belong.bc]))
  mrk.ord <- names(mrk.pos[!duplicated(mrk.pos)])
  a1 <- mappoly::extract_map(map1)
  b1 <- abs(kronecker(a1[mrk.ord], t(a1[dont.belong.f1]), FUN = "-"))
  dimnames(b1) <- list(mrk.ord, dont.belong.f1)
  u <- apply(b1, 2, which.min)
  for (i in seq_along(u))
    mrk.ord <- append(mrk.ord, names(u)[i], after = u[i] - 1)
  if(length(setdiff(dont.belong.bc, mrk.ord))!=0){
    a2 <- mappoly::extract_map(map2)
    b2 <- abs(kronecker(a2[mrk.ord], t(a2[dont.belong.bc]), FUN = "-"))
    dimnames(b2) <- list(mrk.ord, dont.belong.bc)
    u <- apply(b2, 2, which.min)
    for (i in seq_along(u))
      mrk.ord <- append(mrk.ord, names(u)[i], after = u[i] - 1)
  }
  H.res <- NULL
  for(i in mrk.ord){
    x1 <- try(H1$ph1[i,,drop=FALSE], TRUE)
    x2 <- try(H1$ph2[i,,drop=FALSE], TRUE)
    if(!is.matrix(x1) & !is.matrix(x2))next()
    if(is.matrix(x1) & is.matrix(x2)){
      if(any(x1 != x2)) next()
      H.res <- rbind(H.res, x1)
    } else if(is.matrix(x1) & !is.matrix(x2))
      H.res <- rbind(H.res, x1)
    else if(!is.matrix(x1) & is.matrix(x2))
      H.res <- rbind(H.res, x2)
    else stop()
  }
  dim(H.res)
  I195 <- ph_matrix_to_list(H.res)
  names(I195) <- rownames(H.res)
  ## J432
  J432 <- map1$map[[1]]$seq.ph[[2]]
  names(J432) <- map1$info$mrk.names
  ## BC85.209
  BC85.209 <- map2$map[[1]]$seq.ph[[2]]
  names(BC85.209) <-  map2$info$mrk.names

  Ph <- list(I195 = I195,
             J432 = J432,
             BC85.209 = BC85.209)

  #### Joint map ####
  ############################
  input.data <- column_to_rownames(dose.dat, var = "snp")
  input.data[1:6,1:10]
  input.data[input.data == (ploidy+1)] <- NA
  ind.vec <- rownames(ped)
  ngam <- choose(ploidy, ploidy/2)
  A<-as.matrix(expand.grid(0:(ngam-1),
                           0:(ngam-1))[,2:1])
  M <- vector("list", length(ind.vec))
  names(M) <- ind.vec
  dim(A)
  length(M)
  ######################
  # create progress bar
  pb <- txtProgressBar(min = 0, max = length(ind.vec), style = 3)
  cte <- 1
  cat("\n Organizing phases to visit\n")
  for(j in ind.vec){
    setTxtProgressBar(pb, cte)
    cte<- cte +1
    #cat(j, "\n")
    Parents <- ped[j, ]
    if(any(is.na(Parents)))
      next()
    P1 <- Ph[[Parents$P1]]
    P2 <- Ph[[Parents$P2]]
    P1M <- ph_list_to_matrix(P1, ploidy)
    P2M <- ph_list_to_matrix(P2, ploidy)
    rownames(P1M) <- names(P1)
    rownames(P2M) <- names(P2)
    d <- input.data[mrk.ord,j,drop = FALSE]
    d[setdiff(rownames(d), rownames(P1M)),1] <- NA
    d[setdiff(rownames(d), rownames(P2M)),1] <- NA
    I <- vector("list", nrow(d))
    names(I) <- rownames(d)
    for(i in rownames(d)){
      qq <- pp <- rep(0, ploidy)
      if(i%in%rownames(P1M))
        pp <- P1M[i,]
      if(i%in%rownames(P2M))
        qq <- P2M[i,]
      a <- kronecker(apply(combn(pp, ploidy/2), 2, sum),
                     apply(combn(qq, ploidy/2), 2, sum), "+")
      if(is.na(d[i,1]) | !d[i,1]%in%0:ploidy)
        I[[i]] <- A
      else
        I[[i]] <- A[d[i,1] == a, , drop = FALSE]
    }
    for(i in which(sapply(I, is.null))){
      I[[i]] <- A
    }
    M[[j]] <- I
  }
  close(pb)
  A <- vector("list", length(ind.vec))
  names(A) <- ind.vec
  H <- vector("list", nrow(d))
  names(H) <- rownames(d)
  for(i in names(H))
    H[[i]] <- A
  for(i in names(H))
    for(j in names(A))
      H[[i]][[j]] <- M[[j]][[i]]
  H.out <- H
  #### Joint mapping - no error modeling ####
  if(err < 10e-5){
    map <- mappoly:::est_haplo_hmm(ploidy = ploidy,
                                   n.mrk = length(H),
                                   n.ind = length(H[[1]]),
                                   haplo = H,
                                   rf_vec = rep(0.01, length(H)-1),
                                   verbose = FALSE,
                                   use_H0 = FALSE,
                                   tol = tol1)
  } else {
    #### Joint mapping - error modeling ####
    A<-as.matrix(expand.grid(0:(ngam-1),
                             0:(ngam-1))[,2:1])
    B<-apply(A, 1, paste0, collapse = "-")
    log_E <- numeric(length(H) * length(H[[1]]) * nrow(A))
    cte <- count <- 1
    pb <- txtProgressBar(min = 0, max = length(ind.vec), style = 3)
    cat("\n Organizing phases to visit (with error)\n")
    for(j in 1:length(H[[1]])){
      setTxtProgressBar(pb, cte)
      cte<- cte +1
      for(i in 1:length(H)){
        id <- match(apply(H[[i]][[j]], 1, paste0, collapse = "-"), B)
        z <- rep(NA, ngam^2)
        if(length(id) == ngam^2){
          z[id] <- 1/(ngam^2)
        } else {
          z[id] <- (1 - err)/length(id)
          z[is.na(z)] <- err/(ngam^2 - length(id))
        }
        log_E[count:(count+nrow(A)-1)] <- log(z)
        count <- count + nrow(A)
      }
    }
    close(pb)
    n.mrk <- length(H)
    n.ind <- length(H[[1]])
    rf.vec <- rep(0.01, n.mrk-1)
    system.time(map <- est_map_R(m = ploidy,
                                 n.mrk = n.mrk,
                                 n.ind = n.ind,
                                 emit = log_E,
                                 rf_vec = rf.vec,
                                 tol = tol2))
  }
  mrk.names <- names(H)
  seq.num <- match(mrk.names, dose.dat$snp)
  seq.dose.p1 = dose.dat$I195[seq.num]
  seq.dose.p2 = dose.dat$J432[seq.num]
  seq.dose.p3 = dose.dat$BC85.209[seq.num]
  chrom  = dose.dat$chrom[seq.num]
  genome.pos = dose.dat$pos[seq.num]
  seq.ref = dose.dat$seq.ref[seq.num]
  seq.alt = dose.dat$seq.alt[seq.num]
  chisq.pval = dose.dat$chisq.pval[seq.num]
  names(seq.num) <- names(seq.ref) <- names(seq.alt) <- names(chisq.pval) <-
    names(chrom) <- names(genome.pos) <- names(seq.dose.p1) <-
    names(seq.dose.p2) <- names(seq.dose.p3) <- mrk.names


  info = list(ploidy = ploidy,
              n.mrk = length(mrk.names),
              seq.num = seq.num,
              mrk.names = mrk.names,
              seq.dose.p1 = seq.dose.p1,
              seq.dose.p2 = seq.dose.p2,
              seq.dose.p3 = seq.dose.p3,
              chrom = chrom,
              genome.pos = genome.pos,
              seq.ref = seq.ref,
              seq.alt = seq.alt,
              chisq.pval = chisq.pval,
              data.name = "dose.dat",
              ph.thresh = NULL)

  PH1 <-   ph_list_to_matrix(Ph$I195, ploidy)
  rownames(PH1) <- names(Ph$I195)
  PH2 <- ph_list_to_matrix(Ph$J432, ploidy)
  rownames(PH2) <- names(Ph$J432)
  PH3 <- ph_list_to_matrix(Ph$BC85.209, ploidy)
  rownames(PH3) <- names(Ph$BC85.209)
  ph <- vector("list", 3)
  names(ph) <- c("I195","J432","BC85.209")
  for(i in mrk.names){
    ph$I195 <- rbind(ph$I195, tryCatch(PH1[i,],
                                       error = function(e)
                                         return(rep(0, ploidy))))
    ph$J432<- rbind(ph$J432, tryCatch(PH2[i,],
                                      error = function(e)
                                        return(rep(0, ploidy))))
    ph$BC85.209 <- rbind(ph$BC85.209, tryCatch(PH3[i,],
                                               error = function(e)
                                                 return(rep(0, ploidy))))
  }
  for(i in 1:3){
    ph[[i]] <- ph_matrix_to_list(ph[[i]])
    names(ph[[i]]) <- mrk.names
  }
  maps <- list(list(seq.num = seq.num,
                    seq.rf = map[[2]],
                    seq.ph = ph,
                    loglike = map[[1]]))
  structure(list(info = info,
                 maps = maps,
                 H.out = H.out),
            class = "mappoly.map")
}
get_submap_hmm <- function(ploidy, H, mrks, tol= 10e-3){
  Htemp <- H[mrks]
  map <- mappoly:::est_haplo_hmm(ploidy = ploidy,
                                 n.mrk = length(Htemp),
                                 n.ind = length(Htemp[[1]]),
                                 haplo = Htemp,
                                 rf_vec = rep(0.01, length(Htemp)-1),
                                 verbose = FALSE,
                                 use_H0 = FALSE,
                                 tol = tol)
  map.h0 <- mappoly:::est_haplo_hmm(ploidy = ploidy,
                                 n.mrk = length(Htemp),
                                 n.ind = length(Htemp[[1]]),
                                 haplo = Htemp,
                                 rf_vec = 0.5,
                                 verbose = FALSE,
                                 use_H0 = TRUE,
                                 tol = tol)
  return(list(map = map, map.h0 = map.h0))
}
phasing_and_hmm_rf <- function(ch){
  map <- build_hmm_joint_map(map1 = MAPs.f1[[ch]],
                             map2 = MAPs.bc[[ch]],
                             ped = ped,
                             ploidy = ploidy,
                             err = 0.05)
  return(map)
}
#### Merging data ####
dat <- list(dat.f1, dat.bc)
P <- c("I195", "J432", "I195", "BC85.209")
dose.dat <- merge_data(dat, P)
dim(dose.dat)
#### Pedigree ####
F1 <- grep("F1", colnames(dose.dat), value = TRUE)
BC <- grep("AphBC", colnames(dose.dat), value = TRUE)
ped <- data.frame(P1 = c(rep(NA,3), rep("I195", length(c(F1,BC)))),
                  P2 = c(rep(NA,3), rep("J432", length(F1)), rep("BC85.209", length(BC))),
                  row.names = c("I195","J432","BC85.209", F1, BC))
#### Build joint map ####
ploidy = 4
ch <- 1:8
cl <- parallel::makeCluster(8)
parallel::clusterEvalQ(cl, require(mappoly))
parallel::clusterEvalQ(cl, require(tidyverse))
parallel::clusterEvalQ(cl, require(highprecHMM))
parallel::clusterExport(cl, c("dat.f1", "dat.bc", "MAPs.f1", "MAPs.bc", "ped", "ploidy", "dose.dat", "homolog_correspondence", "build_hmm_joint_map"))
MAPs <- parallel::parLapply(cl,ch, phasing_and_hmm_rf)
parallel::stopCluster(cl)
plot_map_list(MAPs, col = inlmisc::GetColors(n = 8, alpha = .8, reverse = TRUE))
#####Removing marker Alf_1045205_SNP1, which is causing a big gap in chromosome 3####
MAPs.bc[[3]] <- drop_marker(MAPs.bc[[3]], mrk = "Alf_1045205_SNP1")
MAPs.f1[[3]] <- drop_marker(MAPs.f1[[3]], mrk = "Alf_1045205_SNP1")
MAPs[[3]] <- phasing_and_hmm_rf(3)
######
plot_map_list(MAPs, col = inlmisc::GetColors(n = 8, alpha = .8, reverse = TRUE))
plot_genome_vs_map(MAPs, alpha = .5, size = 2)
sum.map <- summary_maps(MAPs)
x <- t(sapply(MAPs, function(x) get_map_dosage(x)$freq))
sum.map[, 5:7] <- rbind(x, apply(x, 2, sum))
save(dat, MAPs, sum.map, file = "BC_F1_alfalfa_map.rda")
