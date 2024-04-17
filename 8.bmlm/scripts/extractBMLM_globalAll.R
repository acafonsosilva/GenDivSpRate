library(tidyverse)
library(tidybayes)
library(brms)

folder_path <- '/data/biodiv/asilva/rerunAnalyses/'
folder_path <- '/kingdoms/biodiv/workspace5/asilva/rerunAnalyses2023/'
wd <- paste0(folder_path, 'bmlm/global_mutRate/output/')

mods <- list.dirs(wd, recursive = F, full.names = F)

#allPosterior <- list()
allPosterior <- data.frame()
for (j in mods){
#for (j in mods[3:4]){
  trees <- list.files(paste0(wd,j), recursive = F, full.names = F)
  
  for(tr in trees){
    mod <- readRDS(file = paste0(wd, j, "/",tr))
    
    p1 <- mod %>%
      spread_draws(b_logtipRate, r_clades[condition,term], ndraws = 2000) %>%
      filter(term %in% 'logtipRate',
             condition %in% names(table(mod$data$clades)[table(mod$data$clades) > 20])) %>%
      ungroup() %>%
      mutate(.variable = condition,
             .value = b_logtipRate + r_clades,
             model = j,
             set = str_replace(word(tr,-1, sep = "_"),".rds",""),
             rhat = summary(mod)$fixed['logtipRate','Rhat']) %>%
      dplyr::select( -condition, -term)
    
    allPosterior <- bind_rows(allPosterior, p1)
    #allPosterior[[j]][[tr]] <- p1
    print(tr)
  }
  print(j)
}

saveRDS(allPosterior, paste0(wd, 'global_allPosterior.rds'))


#### Extract traits ####

wd <- paste0(folder_path, 'bmlm/global_traits/output/')
mods <- list.dirs(wd, recursive = F, full.names = F)

#allPosterior <- list()
allPosterior <- data.frame()
for (j in mods){
  trees <- list.files(paste0(wd,j), recursive = F, full.names = F)
  
  for(tr in trees){
    mod <- readRDS(file = paste0(wd, j, "/",tr))
    
    rhats <- data.frame(summary(mod)$fixed[-1,]) %>% 
      rownames_to_column() %>% 
      select(rowname, Rhat) %>%
      mutate(model = j, 
             set = str_replace(word(tr,-1, sep = "_"),".rds",""),
             .variable = paste0("b_",rowname)) %>%
      select(-rowname)
    
   if(j %in% "piSpRateTraits"){
     p1 <- mod %>%
       gather_draws(b_logtipRate,
                    b_logBodyMassKg_notInputed, 
                    b_loggeoArea_km2, 
                    b_logmean_temp,
                    b_loglatitude_mean,
                    b_loglitter_or_clutch_size_n,
                    b_logGenerationLength_d, ndraws = 2000) %>%
       ungroup() %>%
       left_join(., rhats, by = '.variable')
   } 
    
    if(j %in% c("SpRateTraits","piTraits")){
      p1 <- mod %>%
        gather_draws(b_logBodyMassKg_notInputed, 
                     b_loggeoArea_km2, 
                     b_logmean_temp,
                     b_loglatitude_mean,
                     b_loglitter_or_clutch_size_n,
                     b_logGenerationLength_d, ndraws = 2000) %>%
        ungroup() %>%
        left_join(., rhats, by = '.variable')
    } 

    
    allPosterior <- bind_rows(allPosterior, p1)
    print(tr)
    
    #allPosterior[[j]][[tr]] <- p1
  }
  print(j)
}
saveRDS(allPosterior, paste0(wd, 'global_traits_allPosterior.rds'))


#### Extract mutation rates ####

wd <- paste0(folder_path, 'bmlm/global_mutRate/output/')
mods <- list.dirs(wd, recursive = F, full.names = F)

#allPosterior <- list()
allPosterior <- data.frame()
for (j in mods){
  trees <- list.files(paste0(wd, j), recursive = F, full.names = F)
  
  for(tr in trees[90:101]){
    mod <- readRDS(file = paste0(wd, j, "/",tr))
    
    rhats <- data.frame(summary(mod)$fixed) %>% 
      rownames_to_column() %>% 
      mutate(model = j, 
             set = str_replace(word(tr,-1, sep = "_"),".rds",""),
             .variable = paste0("b_",rowname))
   
      p1 <- mod %>%
        gather_draws(b_logtipRate, ndraws = 2000) %>%
        ungroup() %>%
        left_join(., rhats[-1, c('model','set','.variable','Rhat')], by = '.variable')
  
      allPosterior <- bind_rows(allPosterior, p1)
      print(tr)
      
      #allPosterior[[j]][[tr]] <- p1
  }
  print(j)
  
}
saveRDS(allPosterior, paste0(wd, 'global_mutRate_allPosterior_year.rds'))

# folder_path <- '/data/biodiv/asilva/rerunAnalyses/'
# #folder_path <- '/kingdoms/biodiv/workspace5/asilva/rerunAnalyses2023/'
# wd <- paste0(folder_path, 'bmlm/global_mutRate/output/')
# saveRDS(allPosterior, paste0(wd, 'global_mutRate_allPosterior_year.rds'))
