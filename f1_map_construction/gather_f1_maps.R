require(mappoly)
load(file = "~/repos/collaborations/bi-alfalfa-map/f1_filtered_organized_data_.rda")
fl <- list.files(pattern = ".rda", path = "~/repos/collaborations/bi-alfalfa-map/f1_map_construction", full.names = TRUE)
MAPs.f1<- vector("list", length(fl))
cte<-1
for(i in 1:length(MAPs.f1)){
   load(fl[i])
  if(cte == 2 | cte == 3 | cte == 5 | cte == 6 | cte == 7 | cte == 8)
    MAPs.f1[[cte]] <- rev_map(map.f1$mds.err) 
  else
    MAPs.f1[[cte]] <- map.f1$mds.err
  cte<-cte+1
}

plot_map_list(MAPs.f1, col = "ggstyle")
plot_genome_vs_map(MAPs.f1)
summary_maps(MAPs.f1)
save(MAPs.f1, dat.f1, file = "~/repos/collaborations/bi-alfalfa-map/f1_map_construction/f1_map_result.rda")

