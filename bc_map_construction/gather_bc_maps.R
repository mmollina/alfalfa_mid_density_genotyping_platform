require(mappoly)
load(file = "~/repos/collaborations/bi-alfalfa-map/bc_filtered_organized_data_.rda")
fl <- list.files(pattern = ".rda", path = "~/repos/collaborations/bi-alfalfa-map/bc_map_construction", full.names = TRUE)
MAPs.bc<- vector("list", length(fl))
cte<-1
for(i in 1:length(MAPs.bc)){
   load(fl[i])
  if(cte == 2 | cte == 3 | cte == 5 | cte == 6 | cte == 7 | cte == 8)
    MAPs.bc[[cte]] <- rev_map(map.bc$mds.err) 
  else
    MAPs.bc[[cte]] <- map.bc$mds.err
  cte<-cte+1
}

plot_map_list(MAPs.bc, col = "ggstyle")
plot_genome_vs_map(MAPs.bc)
summary_maps(MAPs.bc)
save(MAPs.bc, dat.bc, file = "~/repos/collaborations/bi-alfalfa-map/bc_map_construction/bc_map_result.rda")

