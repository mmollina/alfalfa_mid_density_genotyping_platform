setwd("/backup_mm/repos/collaborations/bi-alfalfa-map/")
load("bc_updog_data.rda")
require(mappoly)
require(reshape2)

#### Filtering ####
dat.bc <- filter_missing(ud.dat.bc, type = 'marker', filter.thres = 0.1, inter = F)
dat.bc <- filter_missing(dat.bc, type = 'individual', filter.thres = 0.05, inter = F)
plot(dat.bc)
#dat1 <- filter_individuals(dat1)
s.f.bc <- filter_segregation(dat.bc, chisq.pval.thres = 0.05/dat.bc$n.mrk, inter = F)  
s.bc <- make_seq_mappoly(s.f.bc)
b.bc <- elim_redundant(s.bc)
s.bc <- make_seq_mappoly(b.bc)
## Two-points
tpt.bc <- est_pairwise_rf(s.bc, ncpus = 16)
m.bc <- rf_list_to_matrix(tpt.bc)
## It was not possible to use UPGMA
#grs <- group_mappoly(m.bc, expected.groups = 50, comp.mat = TRUE)
plot(m.bc, ord = make_seq_mappoly(get_genomic_order(s.bc)), fact = 3)

## selecting markers ##
LGs.bc <- vector("list", 8)
names(LGs.bc) <- 1:8
load("f1_filtered_organized_data_.rda")
pdf("alfalfa_BC_rec_mats.pdf", width = 8.5, height = 11)
layout(matrix(1:8, ncol = 2))
for(i in c(1:8)){
  ## Using F1 map to assemble groups
  d <- abs(kronecker(LGs.f1[[i]]$s$genome.pos, t(s.bc$genome.pos), FUN = "-"))
  s.temp.bc <- make_seq_mappoly(dat.bc, s.bc$seq.mrk.names[unlist(apply(d, 1, function(x) which(x==0)))])
  tpt.temp.bc <- make_pairs_mappoly(tpt.bc, s.temp.bc)
  m.temp.bc <- rf_list_to_matrix(tpt.temp.bc)
  plot(m.temp.bc, ord = s.temp.bc, main.text = i)
  LGs.bc[[i]] <- list(s = s.temp.bc, tpt = tpt.temp.bc, m = m.temp.bc)
}
dev.off()
save(LGs.bc, dat.bc, file = "bc_filtered_organized_data_.rda")