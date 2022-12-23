library(tidyverse)
library(caper)

folder_path <- '/data/biodiv/asilva/'

### global ####

pglsGlobal <- list.files(paste0(folder_path, '7.pgls/outputs/'), full.names = T)

pglssum100 <- data.frame()
for(i in pglsGlobal){
  
  pgls <- readRDS(i)
  pglssum <- data.frame()                      
  for (j in names(pgls)){
    if(!is.null(pgls[[j]])){
      pglssum <- bind_rows(pglssum, 
                           data.frame(set = str_replace(word(i,-1, sep = "_"),".rds",""),
                                      analysis = j,
                                      df = summary(pgls[[j]])$df[2],
                                      term = rownames(summary(pgls[[j]])$coefficients),
                                      summary(pgls[[j]])$coefficients,
                                      sigma = summary(pgls[[j]])$sigma,
                                      lambda = pgls[[j]]$param.CI$lambda$opt,
                                      lambda.CIlow = pgls[[j]]$param.CI$lambda$ci.val[1],
                                      lambda.CIup = pgls[[j]]$param.CI$lambda$ci.val[2], 
                                      stringsAsFactors = FALSE))
    }
  }           
  pglssum100 <- rbind(pglssum100,pglssum)
}
saveRDS(pglssum100, paste0(folder_path, '7.pgls/outputs/global_gendivSpRate_PGLSresults.rds'))


### clade ####

pglsClade <- list.files(paste0(folder_path, '7.pgls/outputs/'), full.names = T)

pglssum100 <- data.frame()
for(i in pglsClade){
  pgls <- readRDS(i)
  clade <- str_replace(word(i,-1, sep = "_"),".rds","")
  for(g in names(pgls)){
    for (j in names(pgls[[g]])){
      pglssum100 <- bind_rows(pglssum100, 
                              data.frame(clade = g,
                                         set = str_replace(word(i,-1, sep = "_"),".rds",""),
                                         analysis = j,
                                         df = summary(pgls[[g]][[j]])$df[2],
                                         term = rownames(summary(pgls[[g]][[j]])$coefficients),
                                         summary(pgls[[g]][[j]])$coefficients,
                                         sigma = summary(pgls[[g]][[j]])$sigma,
                                         lambda = pgls[[g]][[j]]$param.CI$lambda$opt,
                                         lambda.CIlow = pgls[[g]][[j]]$param.CI$lambda$ci.val[1],
                                         lambda.CIup = pgls[[g]][[j]]$param.CI$lambda$ci.val[2], 
                                         stringsAsFactors = FALSE))
    }
  }
}
saveRDS(pglssum100, paste0(folder_path, '7.pgls/outputs/clade_gendivSpRate_PGLSresults.rds'))

### global - traits ####

pglsGlobal <- list.files(paste0(folder_path, '7.pgls/outputs/'), full.names = T)

pglssum100 <- data.frame()
for(i in pglsGlobal){
  
  pgls <- readRDS(i)
  pglssum <- data.frame()                      
  for (j in names(pgls)){
    if(!is.null(pgls[[j]])){
      pglssum <- bind_rows(pglssum, 
                           data.frame(set = str_replace(word(i,-1, sep = "_"),".rds",""),
                                      analysis = j,
                                      df = summary(pgls[[j]])$df[2],
                                      term = rownames(summary(pgls[[j]])$coefficients),
                                      summary(pgls[[j]])$coefficients,
                                      sigma = summary(pgls[[j]])$sigma,
                                      lambda = pgls[[j]]$param.CI$lambda$opt,
                                      lambda.CIlow = pgls[[j]]$param.CI$lambda$ci.val[1],
                                      lambda.CIup = pgls[[j]]$param.CI$lambda$ci.val[2], 
                                      stringsAsFactors = FALSE))
    }
  }           
  pglssum100 <- rbind(pglssum100,pglssum)
}
saveRDS(pglssum100, paste0(folder_path, '7.pgls/outputs/global_gendivSpRateTraits_PGLSresults.rds'))

### global - mutRate ####

pglsGlobal <- list.files(paste0(folder_path, '7.pgls/outputs/'), full.names = T)

pglssum100 <- data.frame()
for(i in pglsGlobal){
  
  pgls <- readRDS(i)
  pglssum <- data.frame()                      
  for (j in names(pgls)){
    if(!is.null(pgls[[j]])){
      pglssum <- bind_rows(pglssum, 
                           data.frame(set = str_replace(word(i,-1, sep = "_"),".rds",""),
                                      analysis = j,
                                      df = summary(pgls[[j]])$df[2],
                                      term = rownames(summary(pgls[[j]])$coefficients),
                                      summary(pgls[[j]])$coefficients,
                                      sigma = summary(pgls[[j]])$sigma,
                                      lambda = pgls[[j]]$param.CI$lambda$opt,
                                      lambda.CIlow = pgls[[j]]$param.CI$lambda$ci.val[1],
                                      lambda.CIup = pgls[[j]]$param.CI$lambda$ci.val[2], 
                                      stringsAsFactors = FALSE))
    }
  }           
  pglssum100 <- rbind(pglssum100,pglssum)
}

saveRDS(pglssum100, paste0(folder_path, '7.pgls/outputs/global_gendivSpRateMutRate_PGLSresults.rds'))


### global - perSite ####

pglsGlobal <- list.files(paste0(folder_path, '7.pgls/output/'), full.names = T)

pglssum100 <- data.frame()
for(i in pglsGlobal){
  
  pgls <- readRDS(i)
  pglssum <- data.frame()                      
  for (j in names(pgls)){
    if(!is.null(pgls[[j]])){
      pglssum <- bind_rows(pglssum, 
                           data.frame(set = str_replace(word(i,-1, sep = "_"),".rds",""),
                                      analysis = j,
                                      df = summary(pgls[[j]])$df[2],
                                      term = rownames(summary(pgls[[j]])$coefficients),
                                      summary(pgls[[j]])$coefficients,
                                      sigma = summary(pgls[[j]])$sigma,
                                      lambda = pgls[[j]]$param.CI$lambda$opt,
                                      lambda.CIlow = pgls[[j]]$param.CI$lambda$ci.val[1],
                                      lambda.CIup = pgls[[j]]$param.CI$lambda$ci.val[2], 
                                      stringsAsFactors = FALSE))
    }
  }           
  pglssum100 <- rbind(pglssum100,pglssum)
}

saveRDS(pglssum100, paste0(folder_path, '7.pgls/outputs/global_gendivSpRatePerSite_PGLSresults.rds'))