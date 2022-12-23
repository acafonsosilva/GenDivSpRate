library(tidyverse)
select <- dplyr::select

#### Load data ####
folder_path <- '/data/biodiv/asilva/'

# genetic diversity species - to discard species without genetic diversity and species which either mean subsampled pi or theta are outside the range of 1000 subsamples to 5 individuals
gen.divCIc_sc <- read.delim(paste0(folder_path, '4.genetic_diversity/GenDiv_subsample5Ind.txt')) %>%
  filter(EstPi > 0, outPi %in% 'in', outTheta %in% 'in') 

## Body size, range size and vagility from Upham et al. 2020
uphamTraits <- read.delim(paste0(folder_path,'otherData/traits/input/Upham2020_vagilityTraits.txt'), stringsAsFactors = FALSE) %>%
  filter(species %in% gen.divCIc_sc$species)

### Other traits that need to be matched first
## Mean temperature (Rolland et al. 2014) - temperature * 10
rolland <- read.delim(paste0(folder_path,'otherData/traits/input/EnvData_Rolland2014.txt'), stringsAsFactors = FALSE) %>%
  dplyr::select(species, mean_temp) 

## Fecundity - litter_or_clutch_size_n (Myhrvold et al. 2015).
# myhrvold
myhrvold <- read.csv(paste0(folder_path,'otherData/traits/input/Amniote_Database_Myhrvold2015.csv'), stringsAsFactors = FALSE) %>% filter(class == 'Mammalia') %>% 
  mutate(species = paste(genus, species, sep= '_'),
         litter_or_clutch_size_n = na_if(litter_or_clutch_size_n,-999.000)) %>%
  select(species, litter_or_clutch_size_n) 

###Longevity - GenerationLength_d (Pacifici et al. 2013)
pacifici <- read.table(paste0(folder_path,'otherData/traits/input/GenerationLenghtMammals_Pacifici2013.txt'),
                       header = TRUE, sep = "\t") %>% 
  mutate(species = str_replace(Scientific_name, " ", "_")) %>%
  select(species, GenerationLength_d) %>%
  filter(!duplicated(species))

## to match trait species names

syn <- readRDS(paste0(folder_path,'2.nameMatching/outputs/synonymsClean_uphamUpdatedNCBIupdated.rds'))

treelabels <- readRDS(paste0(folder_path,'2.nameMatching/outputs/treelabels_synIDs.rds'))

#### Match trait species names to same genetic diversity species ####
traitlabels <- data.frame(traitSp = 
                             c(rolland$species,
                               myhrvold$species, 
                               pacifici$species)) %>%
  filter(!duplicated(traitSp), str_replace(traitSp,"_"," ") %in% 
           c(syn$Accepted, syn$Synonym, syn$SynSpp)) %>%
  mutate(traitSp = str_replace(traitSp,"_"," "),
         syn_accepted = ifelse(traitSp %in% syn$Accepted,
                               'Accepted','NotAccepted'),
         syn_syn = ifelse(traitSp %in% syn$Synonym,
                          'Synonym','NotSynonym'),
         syn_syn_synSpp = ifelse(traitSp %in% syn$SynSpp,
                                 'SynSpp','NotSynSpp'))

## Find synonyms ID for each trait species (not very efficient)
for (i in traitlabels$traitSp){
  dt <- traitlabels[traitlabels$traitSp %in% i,] 
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('Accepted','Synonym','NotSynSpp'))){
    id <- syn[(syn$Synonym %in% i |
                 syn$Accepted %in% i),] 
    if(length(id$ID) == 1){
      traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id$ID
      traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- id$nameUsage
    }
    if(length(id$ID) != 1){
      id2 <- syn[syn$Synonym %in% i,] 
      if(length(id2$ID) == 1){
        traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id2$ID
        traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- id2$nameUsage
      }
      if(length(id2$ID) != 1){
        if(length(unique(id2$nameUsage)) == 1){
          traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- unique(id2$nameUsage)
        }
        if(length(unique(id2$nameUsage)) != 1 & any(id2$nameUsage %in% 'valid')){
          traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- 'valid'
        }
      }
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in%
         c('Accepted','Synonym','SynSpp'))){
    nU <- syn[!syn$nameUsage %in% 'valid' &
                 !syn$tsn %in% 'iucn' &
                 (syn$Accepted %in% i |
                    syn$Synonym %in% i |
                    syn$SynSpp %in% i),]
    if(length(unique(nU$nameUsage)) == 1){
      traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- unique(syn[syn$Synonym %in% i,'nameUsage'])
    }
    if(length(unique(nU)) != 1){
      print(i)
    }
  } 
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('Accepted','NotSynonym','NotSynSpp'))){ 
    id <- syn[syn$Accepted %in% i,]  
    if(length(id$ID) == 1){
      traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id$ID
      if(syn[syn$Accepted %in% i,'tsn'] %in% 'iucn'){
        traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- 
          ifelse(syn[syn$Accepted %in% 
                       i,'nameUsage'] %in% 
                   'iucn_problSyn', 'iucn_problSyn', 
                 syn[syn$Accepted %in% i,'syn_tsn']) ## synonym is also iucn
      }
      if(syn[syn$Accepted %in% i,'tsn'] %in% 'otl'){ 
        traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- 'otl'
      }
      
    }
    if(length(id$ID) != 1){
      print(i)
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in%
         c('Accepted','NotSynonym','SynSpp'))){
    id <- syn[syn$Accepted %in% i |
                syn$SynSpp %in% i,]
    if(length(id$ID) == 1){
      traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id$ID
      if(syn[syn$Accepted %in% i,'tsn'] %in% 'iucn'){
        traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <-
          ifelse(syn[syn$Accepted %in%
                       i,'nameUsage'] %in%
                   'iucn_problSyn', 'iucn_problSyn',
                 syn[syn$Accepted %in% i,'syn_tsn'])
      }
    }
    if(length(id$ID) != 1){
      if(any(!id[id$tsn %in% 'iucn','Synonym'] %in% traitlabels$traitSp)){
        id1 <- id %>% filter(!tsn %in% 'iucn')
        if(length(id1$ID) == 1){
          traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id1$ID
          traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- id1$nameUsage
        }
        if(length(id1$ID) != 1){
          print(i)
        }
      }
      if(any(id[id$tsn %in% c('iucn','otl'),'Accepted'] %in% traitlabels$traitSp)){
        id1 <- syn[syn$SynSpp %in% i,] 
        if(length(id1$ID) == 1){
          traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id1$ID
          traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- id1$synSpp_tsn
        }
        if(length(id1$ID) != 1){
          print(i)
        }
      } 
    }
  } 
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','Synonym','NotSynSpp'))){   
    id <- syn[syn$Synonym %in% i,] 
    if(length(id$ID) == 1){
      traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id$ID
      traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- syn[syn$Synonym %in% i,'nameUsage']
    }
    if(length(id$ID) != 1){
      if(length(unique(syn[syn$Synonym %in% i,
                           'nameUsage'])) == 1){
        traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- syn[syn$Synonym %in% i,'ID'][1] ##doesn't matter which ID because there are multiple IDs due to multiple SynSpp, so the first is fine
        
        traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- unique(syn[syn$Synonym %in% i,'nameUsage'])
      }
      if(length(unique(syn[syn$Synonym %in% i,
                           'nameUsage'])) != 1){
        print(i)
      }
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','Synonym','SynSpp'))){  
    id <- syn[syn$Synonym %in% i,] 
    if(length(id$ID) == 1){
      traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id$ID
      traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- id$nameUsage
    }
    if(length(id$ID) != 1){
      if(all(unique(id$nameUsage) %in% 'problSyn')){
        traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- 
          unique(id$nameUsage)
      }
      if(!all(unique(id$nameUsage) %in% 'problSyn')){
        print(i)
      }
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','NotSynonym','NotSynSpp'))){ 
    traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- 'neverFound'
  }  ##neverfound
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','NotSynonym','SynSpp'))){
    id <- syn[syn$SynSpp %in% i,] 
    if(length(id$ID) == 1){
      traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id$ID
      traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- id$nameUsage
    }
    
    if(length(id$ID) != 1){
      if(length(unique(id$Accepted)) == 1){     
        id1 <- syn[syn$Synonym %in% id[id$SynSpp %in% i,'Accepted'],] %>%
          filter(nameUsage %in% 'valid') #make sure
        if(length(id1$ID) == 1){
          traitlabels[traitlabels$traitSp %in% i,'syn_ID'] <- id1$ID
          traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- id1$nameUsage
        }
        if(length(id1$ID) != 1){
          if(any(str_replace(c(id$Accepted,id$Synonym,id$SynSpp)," ","_")
                 %in% treelabels$species) & 
             is.na(traitlabels[traitlabels$traitSp %in% i,'nameUsage'])){
            print(i)
          }
        }
      }
      if(length(unique(id$Accepted)) != 1){
        if(any(str_replace(c(id$Accepted,id$Synonym,id$SynSpp)," ","_")
               %in% treelabels$species) & any(id$tsn %in% 'iucn')){
          traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- 'iucn_problSyn'
        }
        if(any(str_replace(c(id$Accepted,id$Synonym,id$SynSpp)," ","_")
               %in% treelabels$species) & !any(id$tsn %in% 'iucn')){
          print(i)
        }
      }
    }
  } 
}

## Loop through problems 
for( i in traitlabels$traitSp){
  ss <- syn[syn$Accepted %in% i | syn$Synonym %in% i | syn$SynSpp %in% i,] ##to be sure
  #traitlabels[traitlabels$traitSp %in% i,]
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'valid'){ # - if valid - ok
    if(dim(ss)[1] == 1){
      if(ss$ID %in% treelabels$syn_ID){
        if(treelabels[treelabels$syn_ID %in% ss$ID,'nameUsage'] %in% 'valid'){ ## double check that the same synID is valid in both
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in%
                                                                              traitlabels[traitlabels$traitSp %in%
                                                                                            i,'syn_ID'],
                                                                            'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in%
                                                                           traitlabels[traitlabels$traitSp %in%
                                                                                         i,'syn_ID'],
                                                                         'Family']
        }
        if(!treelabels[treelabels$syn_ID %in% ss$ID,'nameUsage'] %in% 'valid'){
          print(i)
        }
      }
      if(!ss$ID %in% treelabels$syn_ID){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
      }
    }
    if(dim(ss)[1] != 1){
      ss1 <- filter(ss, ID %in% treelabels$syn_ID) ## only one match in tree species
      
      if(dim(ss1)[1] == 1){
        if(treelabels[treelabels$syn_ID %in% ss1$ID,'nameUsage'] %in% 'valid'){ ## double check that the same synID is valid in both
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in%
                                                                              traitlabels[traitlabels$traitSp %in%
                                                                                            i,'syn_ID'],
                                                                            'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in%
                                                                           traitlabels[traitlabels$traitSp %in%
                                                                                         i,'syn_ID'],
                                                                         'Family']
        }
        if(!treelabels[treelabels$syn_ID %in% ss1$ID,'nameUsage'] %in% 'valid'){
          ##if more than one that is valid it should have lumped, etc, then find all that match the valid and check there is only one ID matching tree IDs
          ss_syn <- syn %>%
            filter(Accepted %in% ss$Accepted,
                   ID %in% treelabels$syn_ID)
          
          if(dim(ss_syn)[1] == 1){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                              %in% ss_syn$ID,'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                           %in% ss_syn$ID,'Family']
          }
          if(dim(ss_syn)[1] != 1){
            print(i)
          }
        }
      }
      if(dim(ss1)[1] == 0){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
      }
      if(dim(ss1)[1] > 1){ ## matches with 2 species names, but only one in itis and both in tree
        ### keep ID of respective tip name
        
        if(!all(ss1$ID %in% treelabels$syn_ID)){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in%
                                                                              traitlabels[traitlabels$traitSp %in%
                                                                                            i,'syn_ID'],
                                                                            'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in%
                                                                           traitlabels[traitlabels$traitSp %in%
                                                                                         i,'syn_ID'],
                                                                         'Family']
        }
        if(all(ss1$ID %in% treelabels$syn_ID)){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'subspecies'){ # - if subspecies name - ok
    if(dim(ss)[1] == 1){
      if(syn[syn$Accepted %in% ss$Accepted & syn$nameUsage %in% 'valid','ID'] %in% treelabels$syn_ID){
        id <- syn[syn$Accepted %in% ss$Accepted & syn$nameUsage %in% 'valid','ID'] 
        
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                          %in% id,'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                       %in% id,'Family']
      }
      if(any(!syn[syn$Accepted %in% ss$Accepted,'ID'] %in% treelabels$syn_ID)){
        ss_syn <- syn %>% 
          filter(Accepted %in% ss$Accepted,
                 ID %in% treelabels$syn_ID)
        
        if(dim(ss_syn)[1] == 1){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% ss_syn$ID,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% ss_syn$ID,'Family']
          
        }
        if(dim(ss_syn)[1] == 0){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
        if(dim(ss_syn)[1] > 1){
          print(i)
        }
      } 
    }
    if(dim(ss)[1] != 1){
      s <- str_replace(unique(na.omit(c(ss$Accepted,ss$Synonym,ss$SynSpp)))," ","_")
      s <- s[s %in% treelabels$species]
      
      if(length(s) == 1){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species 
                                                                          %in% s,'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                       %in% s,'Family']
      }
      if(length(s) > 1 & str_replace(i," ","_") %in% s){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species 
                                                                          %in% str_replace(i," ","_"),'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                       %in% str_replace(i," ","_"),'Family']
        
      }
      if(length(s) == 0){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
      }
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'synonymised'){ # - if synonymised - ok
    if(dim(ss)[1] == 1){
      ss_syn <- syn[syn$ID %in% syn[syn$Accepted %in% ss$Accepted ,'ID'] & 
                      syn$ID %in% treelabels$syn_ID, ]
      if(dim(ss_syn)[1] == 1){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                          %in% ss_syn$ID,'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                       %in% ss_syn$ID,'Family']
      }
      if(dim(ss_syn)[1] != 1){
        ss_syn1 <- syn %>%
          filter(ID %in% treelabels$ID)
        
        if(dim(ss_syn1)[1] == 1){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% ss_syn1$ID,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% ss_syn1$ID,'Family']
        }
        if(dim(ss_syn1)[1] == 0){
          if(ss$ID %in% treelabels$syn_ID){
            ###situation with 2 tip names where one is in itis but the other is not
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                              %in% ss$ID,'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                           %in% ss$ID,'Family']
          }
          if(!ss$ID %in% treelabels$syn_ID){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
          }
        }
        if(dim(ss_syn1)[1] > 1){
          print(i)
        }
      }
    }
    if(dim(ss)[1] != 1){
      print(i)
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'lumped'){ # - if lumped - ok
    if(dim(ss)[1] == 1){
      ss_syn <- syn[syn$ID %in% syn[syn$Accepted %in% ss$Accepted,'ID'] & 
                      syn$ID %in% treelabels$syn_ID,]
      if(dim(ss_syn)[1] == 1){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                          %in% ss_syn$ID,'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                       %in% ss_syn$ID,'Family']
      }
      
      if(dim(ss_syn)[1] == 0){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
      }
      if(dim(ss_syn)[1] > 1){ ##both syn and valid are in the tree - keep match regardless
        if(ss$ID %in% treelabels$syn_ID){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% ss$ID,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% ss$ID,'Family']
        }
        if(!ss$ID %in% treelabels$syn_ID){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
    if(dim(ss)[1] != 1){
      
      if(length(unique(ss$Accepted)) == 1){
        if(unique(ss$Accepted) %in% str_replace(treelabels$species,"_"," ")){ 
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in% syn[syn$Synonym %in% unique(ss$Accepted),'ID'],'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in% syn[syn$Synonym %in% unique(ss$Accepted),'ID'],'Family']
        }
        if(!unique(ss$Accepted) %in% str_replace(treelabels$species,"_"," ")){
          print(i)
        }
      }
      if(length(unique(ss$Accepted)) != 1){
        if(any(unique(na.omit(c(ss$Accepted,ss$Synonym, ss$SynSpp))) %in% str_replace(treelabels$species,"_"," "))){
          s <- treelabels[treelabels$species %in% str_replace(unique(na.omit(c(ss$Accepted,ss$Synonym, ss$SynSpp))), " ","_"),'species']
          if(length(s) == 1){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species %in% s,'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species %in% s,'Family']
          }
          if(length(s) != 1){
            print(i)
          }
        }
        
        if(!any(unique(na.omit(c(ss$Accepted,ss$Synonym, ss$SynSpp))) %in% str_replace(treelabels$species,"_"," "))){
          print(i)
        }
      }
    }
  }
  
  ## to check if to keep:
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'invalid'){ #   - invalid - maybe not ok
    if(dim(ss)[1] == 1){
      if(!str_replace(ss$Accepted," ","_") %in% treelabels$species){ ##if accepted otl not in tree
        if(str_replace(word(ss$Synonym,1,2)," ","_") %in% treelabels$species){ ## if species name match tree tip
          id <- treelabels[treelabels$species %in% 
                             str_replace(word(ss$Synonym,1,2)," ","_") ,'syn_ID']
          
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% id,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% id,'Family']
        }
        
        if(!str_replace(word(ss$Synonym,1,2)," ","_") %in% treelabels$species){
          if(ss$ID %in% treelabels$syn_ID){
            print(i)
          }
          if(!ss$ID %in% treelabels$syn_ID){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
          }
        }
      }
      
      if(str_replace(ss$Accepted," ","_") %in% treelabels$species){
        if(ss$ID %in% treelabels$syn_ID){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% ss$ID,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% ss$ID,'Family']
        }
        if(!ss$ID %in% treelabels$syn_ID){
          print(i)
        }
      }
    }
    if(dim(ss)[1] != 1){
      
      if(any(str_replace(ss$Accepted," ","_") %in% treelabels$species)){
        ss1 <- syn %>% filter(Accepted %in% ss$Accepted, 
                              ID %in% treelabels$syn_ID)
        if(dim(ss1)[1] == 1){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% ss1$ID,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% ss1$ID,'Family']
          
        }
        if(dim(ss1)[1] != 1){
          print(i)
        }
      }
      
      if(!any(str_replace(ss$Accepted," ","_") %in% treelabels$species)){
        if(str_replace(word(i,1,2)," ","_") %in% treelabels$species){ ##some subspecies as invalid
          ss2 <- syn[syn$Accepted %in% word(i,1,2) | syn$Synonym %in% word(i,1,2) | syn$SynSpp %in% word(i,1,2),]
          if(dim(ss2)[1] == 1){
            print(i)
          }
          if(dim(ss2)[1] != 1){
            syn_ss2 <- ss2 %>% filter(ID %in% treelabels$syn_ID)
            if(dim(syn_ss2)[1] == 1){
              traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                                %in% syn_ss2$ID,'species']
              
              traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                             %in% syn_ss2$ID,'Family']
            }
            if(dim(syn_ss2)[1] != 1){
              print(i)
            }
          }
        }
        
        if(!str_replace(word(i,1,2)," ","_") %in% treelabels$species){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'invalidspp'){ # - invalidspp - maybe not ok
    if(dim(ss)[1] == 1){
      if(any(syn[syn$Accepted %in% ss$Accepted,'ID'] %in% treelabels$syn_ID)){
        ss_syn <- syn %>% filter(Accepted %in% ss$Accepted,
                                 ID %in% treelabels$syn_ID)
        if(dim(ss_syn)[1] == 1){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% ss_syn$ID,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% ss_syn$ID,'Family']
        }
        
        if(dim(ss_syn)[1] != 1){
          if(str_replace(word(i,1,2)," ","_") %in% treelabels$species & !paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){
            
            ss_syn2 <- syn[syn$Accepted %in% word(i,1,2) & syn$nameUsage %in% 'valid',]
            if(dim(ss_syn2)[1] == 1){
              traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                                %in% ss_syn2$ID,'species']
              
              traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                             %in% ss_syn2$ID,'Family']
            }
            if(dim(ss_syn2)[1] != 1){
              print(i)
            }
          }
          
          if(!str_replace(word(i,1,2)," ","_") %in% treelabels$species){
            if(paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){
              print(i)
            } 
            if(!paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){
              traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
            } 
          }
        }
      }
      if(!any(syn[syn$Accepted %in% ss$Accepted,'ID'] %in% treelabels$syn_ID)){
        if(str_replace(word(i,1,2)," ","_") %in% treelabels$species){
          
          if(!paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){
            ss_syn2 <- treelabels[treelabels$species %in% str_replace(word(i,1,2)," ","_"),]
            
            if(dim(ss_syn2)[1] == 1){ ##some might not have an ID
              traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species 
                                                                                %in% ss_syn2$species,'species']
              
              traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                             %in% ss_syn2$species,'Family']
            }
            if(dim(ss_syn2)[1] != 1){
              print(i)
            }
          }
          if(paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){ #rearranged spp name is also in tree then don't include
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
          }
        }
        
        if(!str_replace(word(i,1,2)," ","_") %in% treelabels$species){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
        
      }
    }
    
    if(dim(ss)[1] != 1){
      print(i)
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'iucn'){ # - iucn - maybe not ok
    if(dim(ss)[1] == 1){
      if(str_replace(ss$Accepted, " ","_") %in% treelabels$species){
        if(!str_replace(ss$Synonym, " ","_") %in% treelabels$species & is.na(ss$SynSpp)){
          if(dim(syn[syn$Accepted %in% i,])[1] == 1){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                              %in% ss$ID,'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                           %in% ss$ID,'Family']
          }
          if(dim(syn[syn$Accepted %in% i,])[1] != 1){
            print(i)
          }
        }
        
        if(str_replace(ss$Synonym, " ","_") %in% treelabels$species | !is.na(ss$SynSpp)){
          if(str_replace(ss$Accepted, " ","_") %in% treelabels$species){
            if(dim(syn[syn$Accepted %in% i,])[1] == 1){
              traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in% ss$ID,'species']
              
              traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                             %in% ss$ID,'Family']
            }
            if(dim(syn[syn$Accepted %in% i,])[1] != 1){
              print(i)
            }   
          }
          if(!str_replace(ss$Accepted, " ","_") %in% treelabels$species){
            print(i)
          }
        }
      }
      
      if(!str_replace(ss$Accepted, " ","_") %in% treelabels$species){
        if(!ss$ID %in% treelabels$syn_ID){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
        
        if(ss$ID %in% treelabels$syn_ID){
          print(i) 
        }
      }
    }
    
    if(dim(ss)[1] != 1){
      if(any(ss$Accepted %in% ss$SynSpp)){
        ss1 <- filter(ss, SynSpp %in% i)
        if(dim(ss1)[1] == 1){
          traitlabels[traitlabels$traitSp %in% i, 'syn_accepted'] <- 'NotAccepted'
          
          if(str_replace(ss1$Accepted, " ","_") %in% treelabels$species){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species %in% str_replace(ss1$Accepted, " ","_"),'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species %in% str_replace(ss1$Accepted, " ","_"),'Family']
          }
          
          if(str_replace(ss1$Synonym, " ","_") %in% treelabels$species & !ss1$nameUsage %in% 'valid'){
            print(i)
            # if(!ss$ID %in% treelabels$syn_ID){
            #   traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
            # }
          }
          
          if(str_replace(ss1$SynSpp, " ","_") %in% treelabels$species){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species %in% str_replace(ss1$SynSpp, " ","_"),'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species %in% str_replace(ss1$SynSpp, " ","_"),'Family']
          }
        }
        
        if(dim(ss1)[1] != 1){
          print(i)
        }
      }
    }
  }
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'otl'){ # - otl - maybe not ok
    if(dim(ss)[1] == 1){
      if(!any(na.omit(c(ss$Accepted,ss$Synonym,ss$SynSpp))  %in% str_replace(treelabels$species,"_"," "))){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
      }
      
      if(any(na.omit(c(ss$Accepted,ss$Synonym,ss$SynSpp))  %in% str_replace(treelabels$species,"_"," "))){
        if(str_replace(ss$Accepted, " ","_") %in% treelabels$species){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species %in% str_replace(ss$Accepted, " ","_"),'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species %in% str_replace(ss$Accepted, " ","_"),'Family']
        }
        
        if(str_replace(ss$Synonym, " ","_") %in% treelabels$species & !ss$nameUsage %in% 'valid'){
          if(ss$ID %in% treelabels$syn_ID){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in% ss$ID,'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in% ss$ID,'Family']
          }
          
           if(!ss$ID %in% treelabels$syn_ID){
             if(ss$nameUsage %in% 'otl_problSyn'){
               ## this trait species name can be matched in phylogeny using synonym but in this case synonym is also otl species so should be verified in the otl_problSyn section
               traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- 'otl_problSyn'
             }
             if(!ss$nameUsage %in% 'otl_problSyn'){
          print(i)
          #   traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
             }
        }
        }
        if(str_replace(ss$SynSpp, " ","_") %in% treelabels$species){
          print(i)
          # traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species %in% str_replace(ss$SynSpp, " ","_"),'species']
          # 
          # traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species %in% str_replace(ss$SynSpp, " ","_"),'Family']
        }
        
      }
    }
    
    if(dim(ss)[1] != 1){
      print(i)
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'iucn_valid'){ # - iucn_valid - maybe not ok
    if(dim(ss)[1] == 1){
      if(!str_replace(ss$Accepted, " ","_") %in% treelabels$species){
        if(str_replace(ss$Synonym, " ","_") %in% treelabels$species & is.na(ss$SynSpp)){
          
          if(dim(syn[syn$Synonym %in% i,])[1] == 1){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                              %in% ss$ID,'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                           %in% ss$ID,'Family']
          }
          if(dim(syn[syn$Synonym %in% i,])[1] != 1){
            print(i)
          }
        }
        if(!str_replace(ss$Synonym, " ","_") %in% treelabels$species | !is.na(ss$SynSpp)){
          if(!ss$ID %in% treelabels$syn_ID){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
          }
          if(ss$ID %in% treelabels$syn_ID){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                              %in% ss$ID,'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                           %in% ss$ID,'Family']
            
          }
        }
      }
      if(str_replace(ss$Accepted, " ","_") %in% treelabels$species){
        print(i)
      }
    }
    if(dim(ss)[1] != 1){
      print(i)
    }
  }
  
  ### to likely exclude 
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'extinct'){ # - extinct
    traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'iucn_problSyn'){ # - iucn_problSyn - maybe some ok?
    if(dim(ss)[1] != 1){
      if(any(ss$ID %in% treelabels$syn_ID)){
        nss <- ss %>% filter(!ID %in% treelabels$syn_ID) ## the names not in tree
        nss0 <- c(nss$Accepted,nss$Synonym,nss$SynSpp) ## the names not in tree
        if(any(!str_replace(nss0," ","_") %in% treelabels$species)){ #if no other names in tree we can keep name
          id <- ss[ss$ID %in% treelabels$syn_ID,'ID']
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% id,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% id,'Family']
        }
        if(any(str_replace(nss0," ","_") %in% treelabels$species)){
          print(i)
        }
      }
      
      if(!any(ss$ID %in% treelabels$syn_ID)){
        print(i)
      }
    }
    if(dim(ss)[1] == 1){
      if(ss$ID %in% treelabels$syn_ID){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                          %in% ss$ID,'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                       %in% ss$ID,'Family']
      }
      
      if(!ss$ID %in% treelabels$syn_ID){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
      }
      
    }
  }    
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'otl_problSyn'){ # - otl_problSyn - maybe some ok?
    nss <- as.character(na.omit(unique(c(ss$Accepted,ss$Synonym,ss$SynSpp))))
    
    if(any(str_replace(nss," ","_") %in% treelabels$species)){
      if(any(syn[syn$Accepted %in% nss | 
                 syn$Synonym %in% nss | 
                 syn$SynSpp %in% nss,'ID'] %in% treelabels$syn_ID)){
        syn_nss <- syn[(syn$Accepted %in% nss | syn$Synonym %in% nss | syn$SynSpp %in% nss) & syn$ID %in% treelabels$syn_ID, ]
        if(dim(syn_nss)[1] == 1){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% syn_nss$ID,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% syn_nss$ID,'Family']
        }
        
        if(dim(syn_nss)[1] != 1){
          print(i)
        }
        
      }
      if(!any(syn[syn$Accepted %in% nss | 
                  syn$Synonym %in% nss | 
                  syn$SynSpp %in% nss,'ID'] %in% treelabels$syn_ID)){
        
        if(length(nss[nss %in% str_replace(treelabels$species,"_"," ")]) == 1){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species %in% str_replace(ss$Synonym, " ","_"),'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species %in% str_replace(ss$Synonym, " ","_"),'Family']
          
        }
        if(length(nss[nss %in% str_replace(treelabels$species,"_"," ")] == 1)){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
    if(!any(str_replace(nss," ","_") %in% treelabels$species)){
      traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'problSyn'){ # - problSyn
    if(dim(ss)[1] == 1){
      ## get all possible matchable species names and check if more than one is in tree
      ss1 <- syn[syn$ID %in% syn[syn$Accepted %in% 
                                   syn[syn$Synonym %in% i,'Accepted'],'ID'] & 
                   syn$ID %in% treelabels$syn_ID,]
      if(dim(ss1)[1] == 1){
        
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                          %in% ss1$ID,'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                       %in% ss1$ID,'Family']
      }
      if(dim(ss1)[1] == 0){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
      }
      
      if(dim(ss1)[1] > 1){
        print(i)
      }
    }
    if(dim(ss)[1] != 1){
      if(str_replace(i," ","_") %in% treelabels$species){
        traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                          %in% str_replace(i," ","_") ,'species']
        
        traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                       %in% str_replace(i," ","_"),'Family']
        traitlabels[traitlabels$traitSp %in% i,'nameUsage'] <- ifelse(any(ss[ss$Synonym %in% i, 'syn_tsn'] %in% 'iucn'), 'iucn', traitlabels[traitlabels$traitSp %in% i,'nameUsage'])
      }
      
      if(!str_replace(i," ","_") %in% treelabels$species){
        print(i)
      }
    }
  }
  
  if(traitlabels[traitlabels$traitSp %in% i, 'nameUsage'] %in% 'neverFound'){
    if(str_detect(i, "-") | str_detect(i, "\\.") | str_detect(i, " x ")){
      traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
    }
    
    if(!(str_detect(i, "-") | str_detect(i, "\\.") | str_detect(i, " x "))){
      if(!str_replace(word(i,1,2)," ","_") %in% treelabels$species){
        if(!paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
        if(paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                            %in% paste(word(i,1),word(i,3),
                                                                                       sep = "_") ,'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                         %in% paste(word(i,1),word(i,3),
                                                                                    sep = "_"),'Family']
        }
      }
      
      if(str_replace(word(i,1,2)," ","_") %in% treelabels$species){
        if(str_count(i, "\\S+") == 3){
          if(!paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                              %in% str_replace(word(i,1,2)," ","_"),'species']
            
            traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                           %in% str_replace(word(i,1,2)," ","_"),'Family']
          }
          if(paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
            traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip' ##don't include spp where rearranged names are both tips
          }
        }
        if(str_count(i, "\\S+") == 2){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                            %in% str_replace(i," ","_"),'species']
          
          traitlabels[traitlabels$traitSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                         %in% str_replace(i," ","_"),'Family']
        }
        if(str_count(i, "\\S+") > 3){
          traitlabels[traitlabels$traitSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
  }
}

#### Filter to genetic diversity dataset and merge datasets ####
traitlabels1 <- traitlabels %>%
  filter(uphamTips %in% gen.divCIc_sc$species) %>%
  mutate(species = str_replace(traitSp," ","_")) %>%
  left_join(., rolland, by = 'species') %>%
  left_join(., myhrvold, by = 'species') %>%
  left_join(., pacifici, by = 'species') 

## collapse NAs and exclude synonyms trait data when 1 duplicated trait species name matches tree species name
traitData <- tibble()
for (i in unique(traitlabels1$uphamTips)){
  dt <- filter(traitlabels1, uphamTips %in% i) %>%
    mutate_all(as.character)
  
  if(dim(dt)[1] == 1){
    traitData <- bind_rows(traitData, 
                           select(dt,uphamTips, Family, 
                                  mean_temp:GenerationLength_d))
  }
  
  if(dim(dt)[1] > 1){
    if(length(na.omit(c(dt$mean_temp, dt$litter_or_clutch_size_n, dt$GenerationLength_d))) < 6){ ##collapses NAs
      dt <- dt %>%
        select(.,uphamTips, Family, mean_temp:GenerationLength_d)
      dt2 <- matrix(ncol = length(colnames(dt)))
      colnames(dt2) <- colnames(dt)
      
      for(j in 1:length(colnames(dt))){
        dt2[,j] <- ifelse(!all(is.na(dt[,j])), na.omit(unique(dt[,j])), NA)
      }
    }
    
    if(length(na.omit(c(dt$mean_temp,dt$litter_or_clutch_size_n, dt$GenerationLength_d))) > 5){  ##keeps only rows that trait data matches tree (not worth merging based on synonyms)
      
      dt2 <- filter(dt, species %in% i) %>%
        select(uphamTips, Family, mean_temp:GenerationLength_d)
    }
    traitData <- bind_rows(traitData, data.frame(dt2, stringsAsFactors = FALSE))
  }
}

traitData <- traitData %>%
  rename(species = uphamTips) %>%
  right_join(., uphamTraits, by = 'species') 

write.table(traitData, paste0(folder_path,'otherData/traits/matchedTraits.txt'), quote = F, col.names = T, row.names = F, sep = '\t')
