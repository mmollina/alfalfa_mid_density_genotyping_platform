#R CMD BATCH --no-save --no-restore '--args ex=10, st=5, thmm=20' map_construction.R  map_construction.Rout &
print(t1 <- Sys.time())
load(file = "~/repos/collaborations/bi-alfalfa-map/bc_filtered_organized_data_.rda")
##Reading arguments
options(echo=TRUE)
arg <- commandArgs(trailingOnly = TRUE)
print(arg)
for(j in 1:length(arg)){
  eval(parse(text=arg[[j]]))
}
print(arg)
rm(arg)
require(mappoly)

## Function
map.construction <- function(LG, start.set = 5, expansion = 10, thres.hmm = 20, phlim = 20){  
  ## MDS order ##
  map <- est_rf_hmm_sequential(input.seq = LG$s,
                                      start.set = start.set,
                                      thres.twopt = 5, 
                                      thres.hmm = thres.hmm,
                                      twopt = LG$tpt,
                                      sub.map.size.diff.limit = expansion, 
                                      phase.number.limit = phlim,
                                      reestimate.single.ph.configuration = TRUE,
                                      tol = 10e-3,
                                      tol.final = 10e-5, detailed.verbose = F)
  map.unique <- filter_map_at_hmm_thres(map, thres.hmm = 0.0000001)
  map.inique.reest <- est_full_hmm_with_global_error(map.unique, 
                                                     error = 0.05, 
                                                     tol = 10e-5)
  map.new <- split_and_rephase(map.inique.reest, twopt = LG$tpt, gap.threshold = 5)
  map.reest <- est_full_hmm_with_global_error(map.new, 
                                                     error = 0.05, 
                                                     tol = 10e-5)
  list(map = map.unique, 
       mds.err = map.reest)
}
map.bc <- map.construction(LG = LGs.bc[[ch]], 
                           expansion = ex, 
                           start.set = st, 
                           thres.hmm = thmm, 
                           phlim = phlim)
save(map.bc, dat.bc, ch,  file = paste0("~/repos/collaborations/bi-alfalfa-map/bc_map_construction/expansion_", ex, "_map_ch_", ch, ".rda"))
print(t2 <- Sys.time())
print(t2-t1)