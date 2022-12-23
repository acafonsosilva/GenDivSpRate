library(RPANDA)
library(coda)
library(stringr)

folder_path <- '/data/biodiv/asilva/5.speciation_rate/outputs/'

pars <- data.frame()
ClaDS_spec_ratedt <- data.frame()
tree_MAPS <- list()

dobj <- list.files(folder_path, pattern = "upham_tree", full.names = T) ### these files are a bit heavy so stayed in the cluster
for (j in dobj){
  load(j)
  treeN = gsub(paste0(folder_path,"/fromNathan_"), "",str_extract(j, regex(".*(?=.Rdata)")))
  
  ntips = tree$Nnode + 1
  nedges = 2 * tree$Nnode
  names = c("sigma","alpha","epsilon",paste0("lambda_",0:nedges))
  names(MAPS) = names
  
  rates = MAPS[5:npar]
  
  branchRate = MAPS[(nedges+5):(nedges+4+ntips)]
  tipRate = rates[tree$edge[, 2] <= (tree$Nnode + 1)]
  
  temp <- data.frame(species = tree$tip.label, branchRate = branchRate, tipRate = tipRate, treeN = treeN)
  
  ClaDS_spec_ratedt <- rbind(ClaDS_spec_ratedt, temp) #### keeps only rates
  pars <- rbind(pars, data.frame(t(MAPS[1:4]), treeN = treeN)) #### keeps hyperparameters
  tree_MAPS[[treeN]] <- list(tree = tree, rates = MAPS)  #### keeps all parameters

  rm(coda_chain)
  print(paste(treeN, 'is done'))
}

saveRDS(tree_MAPS, paste0(folder_path, 'MCCposterior100_treeMAPS.rds'))
write.table(ClaDS_spec_ratedt, paste0(folder_path, 'MCCposterior100_tipRate.txt'), sep = "\t",quote = FALSE, row.names = FALSE)
write.table(pars, paste0(folder_path,'MCCposterior100_pars.txt'), sep = "\t",quote = FALSE, row.names = FALSE)
