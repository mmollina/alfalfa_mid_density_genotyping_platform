load("f1_updog_data.rda")
require(mappoly)
require(reshape2)
#### Filtering ####
dat.f1 <- filter_missing(ud.dat.f1, type = 'marker', filter.thres = 0.1, inter = F)
dat.f1 <- filter_missing(dat.f1, type = 'individual', filter.thres = 0.05, inter = F)
s.f <- filter_segregation(dat.f1, chisq.pval.thres = 0.05/dat.f1$n.mrk, inter = F)  
s <- make_seq_mappoly(s.f)
b <- elim_redundant(s)
s <- make_seq_mappoly(b)
## Two-points
tpt <- est_pairwise_rf(s, ncpus = 8)
m.f1 <- rf_list_to_matrix(tpt)
plot(m.f1, ord = make_seq_mappoly(get_genomic_order(s)), fact = 3)
#### Grouping ####
gr <- group_mappoly(input.mat = m.f1, 
                    expected.groups = 9, comp.mat = TRUE, inter = FALSE)
plot(gr)
id <- c(1,3,2,4,9,5,6,7)
## selecting markers ##
LGs.f1 <- vector("list", 8)
names(LGs.f1) <- 1:8
pdf("alfalfa_F1_rec_mats.pdf", width = 8.5, height = 11)
layout(matrix(1:8, ncol = 2))
for(i in c(1:8)){
  s.temp <- make_seq_mappoly(gr, id[i])
  tpt.temp <- make_pairs_mappoly(tpt, s.temp)
  m.temp <- rf_list_to_matrix(tpt.temp)
  o.mds <- mds_mappoly(m.temp)
  s.mds <- make_seq_mappoly(o.mds) 
  plot(m.temp, ord = s.mds, main.text = i)
  LGs.f1[[i]] <- list(s = s.mds, tpt = tpt.temp, m = m.temp)
}
dev.off()
save(LGs.f1, dat.f1, file = "f1_filtered_organized_data_.rda")
