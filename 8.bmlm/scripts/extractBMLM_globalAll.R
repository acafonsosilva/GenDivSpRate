library(tidyverse)
library(tidybayes)
library(brms)

folder_path <- '/kingdoms/biodiv/workspace5/asilva/SuperCrunchClean/'
wd <- paste0(folder_path, 'bmlm/global/output/')
# wd <- paste0(folder_path, 'bmlm/GenDiv_subsample6Ind/global/output/')   ## temporarily to extract the subsampled dataset that had previously run 

mods <- list.dirs(wd, recursive = F, full.names = F)

allPosterior <- list()
for (j in mods){
#for (j in mods[3:4]){
  trees <- list.files(paste0(wd,j), recursive = F, full.names = F)
  
  for(tr in trees){
    mod <- readRDS(file = paste0(wd, j, "/",tr))
    
    p1 <- mod %>%
      spread_draws(b_logtipRate, r_clades[condition,term], n = 2000) %>%
      filter(term %in% 'logtipRate',
             condition %in% names(table(mod$data$clades)[table(mod$data$clades) > 20])) %>%
      ungroup() %>%
      mutate(.variable = condition,
             .value = b_logtipRate + r_clades,
             model = j,
             set = str_replace(word(tr,-1, sep = "_"),".rds",""),
             rhat = summary(mod)$fixed['logtipRate','Rhat']) %>%
      dplyr::select( -condition, -term)
    
    allPosterior[[j]][[tr]] <- p1
  }
}
# saveRDS(allPosterior, paste0(wd, 'global_allPosterior_subsample6Ind.rds'))

saveRDS(allPosterior, paste0(wd, 'global_allPosterior.rds'))


#### Extract traits ####

wd <- paste0(folder_path, 'bmlm/global_traits/output/')
mods <- list.dirs(wd, recursive = F, full.names = F)

allPosterior <- list()
for (j in mods){
  trees <- list.files(paste0(wd,j), recursive = F, full.names = F)
  
  for(tr in trees){
    mod <- readRDS(file = paste0(wd, j, "/",tr))
    
    rhats <- data.frame(summary(mod)$fixed[-1,'Rhat']) %>% 
      rownames_to_column() %>% 
      mutate(model = j, 
             set = str_replace(word(tr,-1, sep = "_"),".rds",""),
             .variable = paste0("b_",rowname)) %>%
      dplyr::rename(rhat = summary.mod..fixed..1...Rhat..) %>%
      select(-rowname)
   if(j %in% "piSpRateTraits"){
     p1 <- mod %>%
       gather_draws(b_logtipRate,
                    b_logBodyMassKg_notInputed, 
                    b_loggeoArea_km2, 
                    b_logmean_temp,
                    b_loglitter_or_clutch_size_n,
                    b_logGenerationLength_d, n = 2000) %>%
       ungroup() %>%
       left_join(., rhats, by = '.variable')
   } 
    
    if(j %in% c("SpRateTraits","piTraits")){
      p1 <- mod %>%
        gather_draws(b_logBodyMassKg_notInputed, 
                     b_loggeoArea_km2, 
                     b_logmean_temp,
                     b_loglitter_or_clutch_size_n,
                     b_logGenerationLength_d, n = 2000) %>%
        ungroup() %>%
        left_join(., rhats, by = '.variable')
    } 

    
    allPosterior[[j]][[tr]] <- p1
  }
}
saveRDS(allPosterior, paste0(wd, 'global_traits_allPosterior.rds'))


#### Extract mutation rates ####

wd <- paste0(folder_path, 'bmlm/global_mutRate/output/')
mods <- list.dirs(wd, recursive = F, full.names = F)[!list.dirs(wd, recursive = F, full.names = F) %in% "PiSpRatemutRateRand"]

allPosterior <- list()
for (j in mods){
  trees <- list.files(paste0(wd, j), recursive = F, full.names = F)
  
  for(tr in trees){
    mod <- readRDS(file = paste0(wd, j, "/",tr))
    
    rhats <- data.frame(summary(mod)$fixed) %>% 
      rownames_to_column() %>% 
      mutate(model = j, 
             set = str_replace(word(tr,-1, sep = "_"),".rds",""),
             .variable = paste0("b_",rowname))
   
      p1 <- mod %>%
        gather_draws(b_logtipRate, n = 2000) %>%
        ungroup() %>%
        left_join(., rhats[-1, c('model','set','.variable','Rhat')], by = '.variable')
  
    allPosterior[[j]][[tr]] <- p1
  }
}
saveRDS(allPosterior, paste0(wd, 'global_mutRate_allPosterior.rds'))
