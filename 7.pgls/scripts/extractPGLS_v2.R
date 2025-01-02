library(tidyverse)
library(caper)
library(nlme)

folder_path <- '/data/biodiv/asilva/rerunAnalyses/'

### all from nmle####

pglsGlobal <- list.files(paste0(folder_path, 'pgls/output/'), full.names = T, pattern = "global")

pglssum100 <- data.frame()
for(i in pglsGlobal){
  
  pgls <- readRDS(i)
  pglssum <- data.frame()                      
  for (j in names(pgls)){
    if(!is.null(pgls[[j]])){
      pglssum <- bind_rows(pglssum, 
                           data.frame(set = str_replace(word(i,-1, sep = "_"),".rds",""),
                                      analysis = j,
                                      df = pgls[[j]]$dims$N,
                                      #term = rownames(coef(summary(pgls[[j]]))),
                                      t(coef(summary(pgls[[j]]))[,1]),
                                      t(coef(summary(pgls[[j]]))[,2]),
                                      t(coef(summary(pgls[[j]]))[,4]),
                                      lambda = attributes(pgls[[j]]$apVar)$Pars['corStruct'], ####
                                      stringsAsFactors = FALSE))
    }
  }           
  pglssum100 <- rbind(pglssum100,pglssum)
}

colnames(pglssum100)[4:28] <- c('Estimate_intercept','Estimate_tipRate','SE_intercept', 'SE_tipRate','pvalue_intercept','pvalue_tipRate','pgls_lambda',
                                   'Estimate_BodyMass','Estimate_geoArea','Estimate_meanTemp','Estimate_meanLat','Estimate_litterSize','Estimate_genLength',
                                'SE_BodyMass','SE_geoArea','SE_meanTemp','SE_meanLat','SE_litterSize','SE_genLength',
                                   'pvalue_BodyMass','pvalue_geoArea','pvalue_meanTemp','pvalue_meanLat','pvalue_litterSize','pvalue_genLength')

PGLSglobalAllL <- pglssum100 %>%  
  pivot_longer(cols = -c(set, analysis,df, pgls_lambda), 
               names_to = c('Term','variable'), 
               names_pattern = "(.*)_(.*)") %>% 
  drop_na(value)

saveRDS(PGLSglobalAllL, paste0(folder_path, 'pgls/gendivSpRate_PGLSresultsAll.rds'))

#### to estimate R2 of a pgls
library(rr2)
pglsListR2 <- list()
for(i in pglsGlobal){
  pgls <- readRDS(i)
  tree = str_replace(word(i,-1, sep = "_"),".rds","")
  r2 <- R2(mod = pgls[['EstPiSyn']], pred = FALSE) ## pred takes a lot of time
  pglsListR2[[tree]] <- list(pgls[['EstPiSyn']], r2)
}  
saveRDS(pglsListR2, paste0(folder_path, 'pgls/gendivSpRate_PGLS_R2EstPiSyn.rds'))


### clade from caper ####

pglsClade <- list.files(paste0(folder_path, 'pgls/output/'), full.names = T, pattern = "clade")

pglssum100 <- data.frame()
for(i in pglsClade){
  pgls <- readRDS(i)
  #clade <- str_replace(word(i,-1, sep = "_"),".rds","")
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
  print(i)
}
saveRDS(pglssum100, paste0(folder_path, 'pgls/clade/clade_gendivSpRate_PGLSresults.rds'))


