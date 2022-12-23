library(tidyverse)

### make tree labels table where each matches the tsn of the found species name
## Just for species labels in the DNAonly phylogeny
folder_path <- "/data/biodiv/asilva/"

syn <- readRDS(paste0(folder_path,'nameMatching/outputs/synonymsClean_uphamUpdatedNCBIupdated.rds'))

syn0 <- syn %>% 
  filter(!nameUsage %in% 'extinct') 

uphamPC <- read.delim(paste0(folder_path,'nameMatching/inputs/MamPhy_5911sp_tipGenFamOrdCladeGenesSampPC.txt')) %>%
  mutate(species = word(tiplabel,1,2,sep = "_"),
         Family = str_to_sentence(fam),
         Order = str_to_sentence(ord)) %>%
  filter(samp %in% 'sampled', extinct. %in% 0) %>%
  select(species, Family, Order, PC)

treelabels <- uphamPC %>% 
  mutate(syn_accepted = 
           ifelse(species %in% str_replace(syn0$Accepted," ","_"),
                  'Accepted','NotAccepted'),
         syn_syn = ifelse(species %in% 
                            str_replace(syn0$Synonym," ","_"),
                          'Synonym','NotSynonym'),
         syn_syn_synSpp = ifelse(species %in% 
                                   str_replace(syn0$SynSpp," ","_"),
                                 'SynSpp','NotSynSpp'))

for (i in treelabels$species){
  dt <- treelabels[treelabels$species %in% i,] 
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('Accepted','Synonym','NotSynSpp'))){
    id <- syn0[(syn0$Synonym %in% str_replace(i, "_", " ") |
                  syn0$Accepted %in% str_replace(i, "_", " ")),] 
    if(length(id$ID) == 1){
      treelabels[treelabels$species %in% i,'syn_ID'] <- id$ID
      treelabels[treelabels$species %in% i,'nameUsage'] <- id$nameUsage
    }
    if(length(id$ID) != 1){
      id2 <- syn0[syn0$Synonym %in% str_replace(i, "_", " "),] 
      if(length(id2$ID) == 1){
        treelabels[treelabels$species %in% i,'syn_ID'] <- id2$ID
        treelabels[treelabels$species %in% i,'nameUsage'] <- id2$nameUsage
      }
      if(length(id2$ID) != 1){
        if(length(unique(id2$nameUsage)) == 1){
          treelabels[treelabels$species %in% i,'nameUsage'] <- unique(id2$nameUsage)
        }
        if(length(unique(id2$nameUsage)) != 1 & any(id2$nameUsage %in% 'valid')){
          treelabels[treelabels$species %in% i,'nameUsage'] <- 'valid'
        }
      }
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in%
         c('Accepted','Synonym','SynSpp'))){
    # nU <- syn0[!syn0$nameUsage %in% 'valid' &
    #              !syn0$tsn %in% 'iucn' &
    #              (syn0$Accepted %in% str_replace(i, "_", " ") |
    #                 syn0$Synonym %in% str_replace(i, "_", " ") |
    #                 syn0$SynSpp %in% str_replace(i, "_", " ")),]
    # if(length(unique(nU$nameUsage)) == 1){
    #   treelabels[treelabels$species %in% i,'nameUsage'] <- unique(syn0[syn0$Synonym %in% str_replace(i, "_", " "),'nameUsage'])
    # }
    # if(length(unique(nU)) != 1){
    #   print(i)
    # }
  } ##fixed but here just for logic
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('Accepted','NotSynonym','NotSynSpp'))){ 
    id <- syn0[syn0$Accepted %in% str_replace(i, "_", " "),]  
    if(length(id$ID) == 1){
      treelabels[treelabels$species %in% i,'syn_ID'] <- id$ID
      if(syn0[syn0$Accepted %in% str_replace(i, "_", " "),'tsn'] %in% 'iucn'){
        treelabels[treelabels$species %in% i,'nameUsage'] <- 
          ifelse(syn0[syn0$Accepted %in% 
                        str_replace(i, "_", " "),'nameUsage'] %in% 
                   'iucn_problSyn', 'iucn_problSyn', 
                 syn0[syn0$Accepted %in% str_replace(i, "_", " "),'syn_tsn']) ## synonym is also iucn
      }
    }
    if(length(id$ID) != 1){
      print(i)
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in%
         c('Accepted','NotSynonym','SynSpp'))){
    # id <- syn0[syn0$Accepted %in% str_replace(i, "_", " ") |
    #              syn0$SynSpp %in% str_replace(i, "_", " "),]
    # if(length(id) == 1){
    #   treelabels[treelabels$species %in% i,'syn_ID'] <- id
    #   if(syn0[syn0$Accepted %in% str_replace(i, "_", " "),'tsn'] %in% 'iucn'){
    #     treelabels[treelabels$species %in% i,'nameUsage'] <-
    #       ifelse(syn0[syn0$Accepted %in%
    #                     str_replace(i, "_", " "),'nameUsage'] %in%
    #                'iucn_problSyn', 'iucn_problSyn',
    #              syn0[syn0$Accepted %in% str_replace(i, "_", " "),'syn_tsn']) ## synonym is also iucn
    #   }
    # }
    # if(length(id) != 1){
    #   print(i)
    # }
  } ## fixed but here just for logic
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','Synonym','NotSynSpp'))){   
    id <- syn0[syn0$Synonym %in% str_replace(i, "_", " "),] 
    if(length(id$ID) == 1){
      treelabels[treelabels$species %in% i,'syn_ID'] <- id$ID
      treelabels[treelabels$species %in% i,'nameUsage'] <- syn0[syn0$Synonym %in% str_replace(i, "_", " "),'nameUsage']
    }
    if(length(id$ID) != 1){
      if(length(unique(syn0[syn0$Synonym %in% str_replace(i, "_", " "),
                            'nameUsage'])) == 1){
        treelabels[treelabels$species %in% i,'nameUsage'] <- unique(syn0[syn0$Synonym %in% str_replace(i, "_", " "),'nameUsage'])
      }
      if(length(unique(syn0[syn0$Synonym %in% str_replace(i, "_", " "),
                            'nameUsage'])) != 1){
        print(i)
      }
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','Synonym','SynSpp'))){  
    # id <- syn0[syn0$Synonym %in% str_replace(i, "_", " "),'ID'] 
  }## fixed but here just for logic
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','NotSynonym','NotSynSpp'))){ 
    treelabels[treelabels$species %in% i,'nameUsage'] <- 'neverFound'
  }  ##neverfound
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','NotSynonym','SynSpp'))){
    # id <- syn0[syn0$SynSpp %in% str_replace(i, "_", " "),] 
  } ##fixed but here just for logic
}

saveRDS(treelabels, paste0(folder_path,'nameMatching/outputs/treelabels_synIDs.rds'))
