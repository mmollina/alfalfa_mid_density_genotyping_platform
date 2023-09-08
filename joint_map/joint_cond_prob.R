##Reading arguments
options(echo=TRUE)
arg <- commandArgs(trailingOnly = TRUE)
print(arg)
for(j in 1:length(arg)){
  eval(parse(text=arg[[j]]))
}
print(arg)
rm(arg)
.libPaths( c( "/home3/mmollin/singularity_R_packages" , .libPaths()))
require(highprecHMM)
require(mappoly)
setwd("~/repos/collaborations/bi-alfalfa-map/joint_map")
load("BC_F1_alfalfa_map.rda")
err <- 0.05
ploidy <- 4
genoprob <- vector("list", 8)
for(ch in 1:8){
  H <- MAPs[[ch]]$H.out
  #### preparing log_error emission ####
  ngam <- choose(ploidy, ploidy/2)
  A<-as.matrix(expand.grid(0:(ngam-1),
                           0:(ngam-1))[,2:1])
  B<-apply(A, 1, paste0, collapse = "-")
  log_E <- numeric(length(H) * length(H[[1]]) * nrow(A))
  count <- 1
  for(j in 1:length(H[[1]])){
    cat(j, "\n")
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
  mrknames <- names(H)
  indnames <- names(H[[1]])
  n.mrk <- length(mrknames)
  n.ind <- length(indnames)
  rf.vec <- MAPs[[ch]]$maps[[1]]$seq.rf
  genoprob[[ch]] <- calc_genoprob_R(m = ploidy,
                                    mrknames = mrknames,
                                    indnames = indnames,
                                    n.mrk = n.mrk,
                                    n.ind = n.ind,
                                    emit = log_E,
                                    rf_vec = rf.vec)
}
save(genoprob, ped, file = paste0("~/repos/collaborations/bi-alfalfa-map/joint_map/genoprobs.rda"),
     compress = "xz", compression_level = 9)
w <- calc_homologprob(genoprob)
plot(w, lg = "all", ind = "F1.85.47", use.plotly = FALSE)
plot(w, lg = "all", ind = "AphBC.43", use.plotly = FALSE)


