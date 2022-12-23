library(tidyverse)

### Make seqLabels list
folder_path <- "/data/biodiv/asilva/"

inputSeq <- readRDS(paste0(folder_path,'nameMatching/inputs/genBank_metadata_complete.rds'))
syn <- readRDS(paste0(folder_path,'nameMatching/outputs/synonymsClean_uphamUpdatedNCBIupdated.rds'))
treelabels <- readRDS(paste0(folder_path,'nameMatching/outputs/treelabels_synIDs.rds'))

seqlabels <- inputSeq %>% 
  select(seqSp) %>%
  distinct() %>%
  mutate(syn_accepted = ifelse(seqSp %in% syn$Accepted,
                               'Accepted','NotAccepted'),
         syn_syn = ifelse(seqSp %in% syn$Synonym,
                          'Synonym','NotSynonym'),
         syn_syn_synSpp = ifelse(seqSp %in% syn$SynSpp,
                                 'SynSpp','NotSynSpp'))

for (i in seqlabels$seqSp){
  dt <- seqlabels[seqlabels$seqSp %in% i,] 
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('Accepted','Synonym','NotSynSpp'))){
    id <- syn[(syn$Synonym %in% i |
                 syn$Accepted %in% i),] 
    if(length(id$ID) == 1){
      seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id$ID
      seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- id$nameUsage
    }
    if(length(id$ID) != 1){
      id2 <- syn[syn$Synonym %in% i,] 
      if(length(id2$ID) == 1){
        seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id2$ID
        seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- id2$nameUsage
      }
      if(length(id2$ID) != 1){
        if(length(unique(id2$nameUsage)) == 1){
          seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- unique(id2$nameUsage)
        }
        if(length(unique(id2$nameUsage)) != 1 & any(id2$nameUsage %in% 'valid')){
          seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- 'valid'
        }
      }
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in%
         c('Accepted','Synonym','SynSpp'))){
    # nU <- syn[!syn$nameUsage %in% 'valid' &
    #              !syn$tsn %in% 'iucn' &
    #              (syn$Accepted %in% i |
    #                 syn$Synonym %in% i |
    #                 syn$SynSpp %in% i),]
    # if(length(unique(nU$nameUsage)) == 1){
    #   seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- unique(syn[syn$Synonym %in% i,'nameUsage'])
    # }
    # if(length(unique(nU)) != 1){
    #   print(i)
    # }
  } ##fixed but here just for logic
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('Accepted','NotSynonym','NotSynSpp'))){ 
    id <- syn[syn$Accepted %in% i,]  
    if(length(id$ID) == 1){
      seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id$ID
      if(syn[syn$Accepted %in% i,'tsn'] %in% 'iucn'){
        seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- 
          ifelse(syn[syn$Accepted %in% 
                       i,'nameUsage'] %in% 
                   'iucn_problSyn', 'iucn_problSyn', 
                 syn[syn$Accepted %in% i,'syn_tsn']) ## synonym is also iucn
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
      seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id$ID
      if(syn[syn$Accepted %in% i,'tsn'] %in% 'iucn'){
        seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <-
          ifelse(syn[syn$Accepted %in%
                       i,'nameUsage'] %in%
                   'iucn_problSyn', 'iucn_problSyn',
                 syn[syn$Accepted %in% i,'syn_tsn'])
      }
    }
    if(length(id$ID) != 1){
      if(any(!id[id$tsn %in% 'iucn','Synonym'] %in% seqlabels$seqSp)){
        id1 <- id %>% filter(!tsn %in% 'iucn')
        if(length(id1$ID) == 1){
          seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id1$ID
          seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- id1$nameUsage
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
      seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id$ID
      seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- syn[syn$Synonym %in% i,'nameUsage']
    }
    if(length(id$ID) != 1){
      if(length(unique(syn[syn$Synonym %in% i,
                           'nameUsage'])) == 1){
        seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- unique(syn[syn$Synonym %in% i,'nameUsage'])
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
      seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id$ID
      seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- id$nameUsage
    }
    if(length(id$ID) != 1){
      if(all(unique(id$nameUsage) %in% 'problSyn')){
        seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- 
          unique(id$nameUsage)
      }
      if(!all(unique(id$nameUsage) %in% 'problSyn')){
        print(i)
      }
    }
  }
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','NotSynonym','NotSynSpp'))){ 
    seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- 'neverFound'
  }  ##neverfound
  if(all(select(dt,syn_accepted:syn_syn_synSpp) %in% 
         c('NotAccepted','NotSynonym','SynSpp'))){
    id <- syn[syn$SynSpp %in% i,] 
    if(length(id$ID) == 1){
      seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id$ID
      seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- id$nameUsage
    }
    
    if(length(id$ID) != 1){
      if(length(unique(id$Accepted)) == 1){     
        id1 <- syn[syn$Synonym %in% id[id$SynSpp %in% i,'Accepted'],] %>%
          filter(nameUsage %in% 'valid') #make sure
        if(length(id1$ID) == 1){
          seqlabels[seqlabels$seqSp %in% i,'syn_ID'] <- id1$ID
          seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- id1$nameUsage
        }
        if(length(id1$ID) != 1){
          if(any(str_replace(c(id$Accepted,id$Synonym,id$SynSpp)," ","_")
                 %in% treelabels$species) & 
             is.na(seqlabels[seqlabels$seqSp %in% i,'nameUsage'])){
            print(i)
          }
        }
      }
      if(length(unique(id$Accepted)) != 1){
        if(any(str_replace(c(id$Accepted,id$Synonym,id$SynSpp)," ","_")
               %in% treelabels$species) & any(id$tsn %in% 'iucn')){
          seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- 'iucn_problSyn'
        }
        if(any(str_replace(c(id$Accepted,id$Synonym,id$SynSpp)," ","_")
               %in% treelabels$species) & !any(id$tsn %in% 'iucn')){
          print(i)
        }
      }
    }
  } 
}


#### match sequence species name to tree species names
seqlabels <- seqlabels %>% mutate(species = word(seqSp,1,2))
  
# match directly sequence species names that are found in tree but never in synonyms list
## risking to include subspecies names
naSynID <- treelabels %>% 
  filter(species %in% str_replace(seqlabels[
    word(seqlabels[is.na(seqlabels$syn_ID),'seqSp'],1,2) 
                     %in% str_replace(treelabels$species,"_"," "),'seqSp']," ","_"),
         nameUsage %in% 'neverFound') %>%
  select(species, Family) %>%
  mutate(uphamTips = species,
         species = str_replace(species, "_", " "))

seqlabels <- left_join(seqlabels, naSynID, by = 'species') %>% drop_na(seqSp)

#match tree to seqlabels through ID - loop through potential problems 
for( i in seqlabels[is.na(seqlabels$uphamTips), 'seqSp']){
  ss <- syn[syn$Accepted %in% i | syn$Synonym %in% i | syn$SynSpp %in% i,] ##to be sure
  #seqlabels[seqlabels$seqSp %in% i,]
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'valid'){ # - if valid - ok
    if(dim(ss)[1] == 1){
      if(ss$ID %in% treelabels$syn_ID){
        if(treelabels[treelabels$syn_ID %in% ss$ID,'nameUsage'] %in% 'valid'){ ## double check that the same synID is valid in both
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in%
                                                                        seqlabels[seqlabels$seqSp %in%
                                                                                    i,'syn_ID'],
                                                                      'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in%
                                                                     seqlabels[seqlabels$seqSp %in%
                                                                                 i,'syn_ID'],
                                                                   'Family']
        }
        if(!treelabels[treelabels$syn_ID %in% ss$ID,'nameUsage'] %in% 'valid'){
          print(i)
        }
      }
      if(!ss$ID %in% treelabels$syn_ID){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
      }
    }
    if(dim(ss)[1] != 1){
      ss1 <- filter(ss, ID %in% treelabels$syn_ID) ## only one match in tree species
      
      if(dim(ss1)[1] == 1){
          if(treelabels[treelabels$syn_ID %in% ss1$ID,'nameUsage'] %in% 'valid'){ ## double check that the same synID is valid in both
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in%
                                                                          seqlabels[seqlabels$seqSp %in%
                                                                                      i,'syn_ID'],
                                                                        'species']
            
            seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in%
                                                                       seqlabels[seqlabels$seqSp %in%
                                                                                   i,'syn_ID'],
                                                                     'Family']
          }
          if(!treelabels[treelabels$syn_ID %in% ss1$ID,'nameUsage'] %in% 'valid'){
            ##if more than one that is valid it should have lumped, etc, then find all that match the valid and check there is only one ID matching tree IDs
            ss_syn <- syn %>%
              filter(Accepted %in% ss$Accepted,
                     ID %in% treelabels$syn_ID)
            
              if(dim(ss_syn)[1] == 1){
                seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                            %in% ss_syn$ID,'species']
                
                seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                         %in% ss_syn$ID,'Family']
              }
            if(dim(ss_syn)[1] != 1){
              print(i)
            }
          }
      }
      if(dim(ss1)[1] == 0){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
      }
      if(dim(ss1)[1] > 1){ ## matches with 2 species names, but only one in itis and both in tree
        ### keep ID of respective tip name
        
        if(!all(ss1$ID %in% treelabels$syn_ID)){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID %in%
                                                                      seqlabels[seqlabels$seqSp %in%
                                                                                  i,'syn_ID'],
                                                                    'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID %in%
                                                                   seqlabels[seqlabels$seqSp %in%
                                                                               i,'syn_ID'],
                                                                 'Family']
        }
        if(all(ss1$ID %in% treelabels$syn_ID)){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
  }
}
    
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'subspecies'){ # - if subspecies name - ok
    if(dim(ss)[1] == 1){
      if(syn[syn$Accepted %in% ss$Accepted & syn$nameUsage %in% 'valid','ID'] %in% treelabels$syn_ID){
        id <- syn[syn$Accepted %in% ss$Accepted & syn$nameUsage %in% 'valid','ID'] 
        
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                    %in% id,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                 %in% id,'Family']
      }
      if(any(!syn[syn$Accepted %in% ss$Accepted,'ID'] %in% treelabels$syn_ID)){
        ss_syn <- syn %>% 
          filter(Accepted %in% ss$Accepted,
                 ID %in% treelabels$syn_ID)
        
        if(dim(ss_syn)[1] == 1){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% ss_syn$ID,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                   %in% ss_syn$ID,'Family']
          
        }
        if(dim(ss_syn)[1] == 0){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
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
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$species 
                                                                    %in% s,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                 %in% s,'Family']
      }
      if(length(s) > 1 & str_replace(i," ","_") %in% s){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$species 
                                                                    %in% str_replace(i," ","_"),'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                 %in% str_replace(i," ","_"),'Family']
        
      }
      if(length(s) == 0){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
      }
    }
  }
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'synonymised'){ # - if synonymised - ok
    if(dim(ss)[1] == 1){
      ss_syn <- syn[syn$ID %in% syn[syn$Accepted %in% ss$Accepted ,'ID'] & 
                      syn$ID %in% treelabels$syn_ID, ]
      if(dim(ss_syn)[1] == 1){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                    %in% ss_syn$ID,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                 %in% ss_syn$ID,'Family']
      }
      if(dim(ss_syn)[1] != 1){
        ss_syn1 <- syn %>%
          filter(ID %in% treelabels$ID)
        
        if(dim(ss_syn1)[1] == 1){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% ss_syn1$ID,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                   %in% ss_syn1$ID,'Family']
        }
        if(dim(ss_syn1)[1] == 0){
          if(ss$ID %in% treelabels$syn_ID){
            ###situation with 2 tip names where one is in itis but the other is not
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                        %in% ss$ID,'species']
            
            seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                     %in% ss$ID,'Family']
          }
          if(!ss$ID %in% treelabels$syn_ID){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
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
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'lumped'){ # - if lumped - ok
    if(dim(ss)[1] == 1){
      ss_syn <- syn[syn$ID %in% syn[syn$Accepted %in% ss$Accepted,'ID'] & 
                      syn$ID %in% treelabels$syn_ID,]
      if(dim(ss_syn)[1] == 1){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                    %in% ss_syn$ID,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                 %in% ss_syn$ID,'Family']
      }
      
      if(dim(ss_syn)[1] == 0){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
      }
      if(dim(ss_syn)[1] > 1){ ##both syn and valid are in the tree - keep match regardless
        if(ss$ID %in% treelabels$syn_ID){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                    %in% ss$ID,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                 %in% ss$ID,'Family']
        }
        if(!ss$ID %in% treelabels$syn_ID){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
    if(dim(ss)[1] != 1){
      print(i)
    }
  }
  
  ## to check if to keep:
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'invalid'){ #   - invalid - maybe not ok
    if(dim(ss)[1] == 1){
      if(!str_replace(ss$Accepted," ","_") %in% treelabels$species){ ##if accepted otl not in tree
        if(str_replace(word(ss$Synonym,1,2)," ","_") %in% treelabels$species){ ## if species name match tree tip
          id <- treelabels[treelabels$species %in% 
                             str_replace(word(ss$Synonym,1,2)," ","_") ,'syn_ID']
          
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% id,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                   %in% id,'Family']
        }
        
        if(!str_replace(word(ss$Synonym,1,2)," ","_") %in% treelabels$species){
          if(ss$ID %in% treelabels$syn_ID){
            print(i)
          }
          if(!ss$ID %in% treelabels$syn_ID){
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
          }
        }
      }
      
      if(str_replace(ss$Accepted," ","_") %in% treelabels$species){
        if(ss$ID %in% treelabels$syn_ID){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% ss$ID,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
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
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% ss1$ID,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
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
              seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                          %in% syn_ss2$ID,'species']
              
              seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                       %in% syn_ss2$ID,'Family']
            }
            if(dim(syn_ss2)[1] != 1){
              print(i)
            }
          }
        }
        
        if(!str_replace(word(i,1,2)," ","_") %in% treelabels$species){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
  }
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'invalidspp'){ # - invalidspp - maybe not ok
    if(dim(ss)[1] == 1){
      if(any(syn[syn$Accepted %in% ss$Accepted,'ID'] %in% treelabels$syn_ID)){
        ss_syn <- syn %>% filter(Accepted %in% ss$Accepted,
                                 ID %in% treelabels$syn_ID)
        if(dim(ss_syn)[1] == 1){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% ss_syn$ID,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                   %in% ss_syn$ID,'Family']
        }
        
        if(dim(ss_syn)[1] != 1){
          if(str_replace(word(i,1,2)," ","_") %in% treelabels$species & !paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){
            
            ss_syn2 <- syn[syn$Accepted %in% word(i,1,2) & syn$nameUsage %in% 'valid',]
            if(dim(ss_syn2)[1] == 1){
              seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                          %in% ss_syn2$ID,'species']
              
              seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
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
              seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
            } 
          }
        }
      }
      if(!any(syn[syn$Accepted %in% ss$Accepted,'ID'] %in% treelabels$syn_ID)){
        if(str_replace(word(i,1,2)," ","_") %in% treelabels$species){
          
          if(!paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){
            ss_syn2 <- treelabels[treelabels$species %in% str_replace(word(i,1,2)," ","_"),]
            
            if(dim(ss_syn2)[1] == 1){ ##some might not have an ID
              seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$species 
                                                                          %in% ss_syn2$species,'species']
              
              seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                       %in% ss_syn2$species,'Family']
            }
            if(dim(ss_syn2)[1] != 1){
              print(i)
            }
        }
          if(paste(word(i,1),word(i,3), sep = "_") %in% treelabels$species){ #rearranged spp name is also in tree then don't include
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
          }
        }
        
        if(!str_replace(word(i,1,2)," ","_") %in% treelabels$species){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
        }
        
      }
    }
    
    if(dim(ss)[1] != 1){
      print(i)
    }
  }
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'iucn'){ # - iucn - maybe not ok
    if(dim(ss)[1] == 1){
      if(str_replace(ss$Accepted, " ","_") %in% treelabels$species){
        if(!str_replace(ss$Synonym, " ","_") %in% treelabels$species & is.na(ss$SynSpp)){
          if(dim(syn[syn$Accepted %in% i,])[1] == 1){
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                        %in% ss$ID,'species']
            
            seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                     %in% ss$ID,'Family']
          }
          if(dim(syn[syn$Accepted %in% i,])[1] != 1){
            print(i)
          }
        }
        
        if(str_replace(ss$Synonym, " ","_") %in% treelabels$species | !is.na(ss$SynSpp)){
          print(i)
        }
      }
      
      if(!str_replace(ss$Accepted, " ","_") %in% treelabels$species){
        if(!ss$ID %in% treelabels$syn_ID){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
        }
        
        if(ss$ID %in% treelabels$syn_ID){
          print(i) 
        }
      }
    }
    
    if(dim(ss)[1] != 1){
      print(i)
    }
  }
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'iucn_valid'){ # - iucn_valid - maybe not ok
    if(dim(ss)[1] == 1){
      if(!str_replace(ss$Accepted, " ","_") %in% treelabels$species){
        if(str_replace(ss$Synonym, " ","_") %in% treelabels$species & is.na(ss$SynSpp)){
          
          if(dim(syn[syn$Synonym %in% i,])[1] == 1){
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                        %in% ss$ID,'species']
            
            seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                     %in% ss$ID,'Family']
          }
          if(dim(syn[syn$Synonym %in% i,])[1] != 1){
            print(i)
          }
        }
        if(!str_replace(ss$Synonym, " ","_") %in% treelabels$species | !is.na(ss$SynSpp)){
          if(!ss$ID %in% treelabels$syn_ID){
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
          }
          if(ss$ID %in% treelabels$syn_ID){
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                        %in% ss$ID,'species']
            
            seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
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
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'extinct'){ # - extinct
    seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
  }
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'iucn_problSyn'){ # - iucn_problSyn - maybe some ok?
    if(dim(ss)[1] != 1){
      if(any(ss$ID %in% treelabels$syn_ID)){
        nss <- ss %>% filter(!ID %in% treelabels$syn_ID) ## the names not in tree
        nss0 <- c(nss$Accepted,nss$Synonym,nss$SynSpp) ## the names not in tree
        if(any(!str_replace(nss0," ","_") %in% treelabels$species)){ #if no other names in tree we can keep name
          id <- ss[ss$ID %in% treelabels$syn_ID,'ID']
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% id,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
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
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                    %in% ss$ID,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                 %in% ss$ID,'Family']
      }
      
      if(!ss$ID %in% treelabels$syn_ID){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
      }
      
    }
  }    
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'otl_problSyn'){ # - otl_problSyn - maybe some ok?
    nss <- as.character(na.omit(unique(c(ss$Accepted,ss$Synonym,ss$SynSpp))))
    
    if(any(str_replace(nss," ","_") %in% treelabels$species)){
      if(any(syn[syn$Accepted %in% nss | 
                 syn$Synonym %in% nss | 
                 syn$SynSpp %in% nss,'ID'] %in% treelabels$syn_ID)){
        syn_nss <- syn[(syn$Accepted %in% nss | syn$Synonym %in% nss | syn$SynSpp %in% nss) & syn$ID %in% treelabels$syn_ID, ]
        if(dim(syn_nss)[1] == 1){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                      %in% syn_nss$ID,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                   %in% syn_nss$ID,'Family']
        }
        
        if(dim(syn_nss)[1] != 1){
          print(i)
        }
        
      }
      if(!any(syn[syn$Accepted %in% nss | 
                 syn$Synonym %in% nss | 
                 syn$SynSpp %in% nss,'ID'] %in% treelabels$syn_ID)){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
      }
    }
    if(!any(str_replace(nss," ","_") %in% treelabels$species)){
      seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
    }
  }
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'problSyn'){ # - problSyn
    if(dim(ss)[1] == 1){
      ## get all possible matchable species names and check if more than one is in tree
      ss1 <- syn[syn$ID %in% syn[syn$Accepted %in% 
                                   syn[syn$Synonym %in% i,'Accepted'],'ID'] & 
                   syn$ID %in% treelabels$syn_ID,]
      if(dim(ss1)[1] == 1){
        
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$syn_ID 
                                                                    %in% ss1$ID,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$syn_ID 
                                                                 %in% ss1$ID,'Family']
      }
      if(dim(ss1)[1] == 0){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
      }
      
      if(dim(ss1)[1] > 1){
        print(i)
      }
    }
    if(dim(ss)[1] != 1){
      if(str_replace(i," ","_") %in% treelabels$species){
        seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                    %in% str_replace(i," ","_") ,'species']
        
        seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                 %in% str_replace(i," ","_"),'Family']
        seqlabels[seqlabels$seqSp %in% i,'nameUsage'] <- ifelse(any(ss[ss$Synonym %in% i, 'syn_tsn'] %in% 'iucn'), 'iucn', seqlabels[seqlabels$seqSp %in% i,'nameUsage'])
      }
      
      if(!str_replace(i," ","_") %in% treelabels$species){
        print(i)
      }
    }
  }
  
  if(seqlabels[seqlabels$seqSp %in% i, 'nameUsage'] %in% 'neverFound'){
    if(str_detect(i, "-") | str_detect(i, "\\.") | str_detect(i, " x ")){
      seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
    }
    
    if(!(str_detect(i, "-") | str_detect(i, "\\.") | str_detect(i, " x "))){
      if(!str_replace(word(i,1,2)," ","_") %in% treelabels$species){
        if(!paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
        }
        if(paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                      %in% paste(word(i,1),word(i,3),
                                                                                 sep = "_") ,'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                   %in% paste(word(i,1),word(i,3),
                                                                              sep = "_"),'Family']
        }
      }
      
      if(str_replace(word(i,1,2)," ","_") %in% treelabels$species){
        if(str_count(i, "\\S+") == 3){
          if(!paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                        %in% str_replace(word(i,1,2)," ","_"),'species']
            
            seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                     %in% str_replace(word(i,1,2)," ","_"),'Family']
          }
          if(paste(word(i,1),word(i,3),sep = "_") %in% treelabels$species){
            seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip' ##don't include spp where rearranged names are both tips
          }
        }
        if(str_count(i, "\\S+") == 2){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- treelabels[treelabels$species
                                                                      %in% str_replace(i," ","_"),'species']
          
          seqlabels[seqlabels$seqSp %in% i,'Family'] <- treelabels[treelabels$species 
                                                                   %in% str_replace(i," ","_"),'Family']
        }
        if(str_count(i, "\\S+") > 3){
          seqlabels[seqlabels$seqSp %in% i,'uphamTips'] <- 'noTip'
        }
      }
    }
  }
}

write.table(seqlabels,paste0(folder_path,'nameMatching/outputs/seqlabels_matchedUpham.txt'), quote = F, col.names = T, row.names = F, sep = '\t')

#### assigning all identified sequence species names to sequences

inputSeq1 <- left_join(inputSeq, seqlabels[,c('seqSp','Family','uphamTips')],
                       by ='seqSp') %>% 
  filter(!duplicated(accession)) %>%
  rename(genBank_seqSp = seqSp,
         uphamFamily = Family)

write.table(inputSeq1,paste0(folder_path,'nameMatching/outputs/genBank_metadata_complete_matchedUpham0.txt'), quote = F, col.names = T, row.names = F, sep = '\t')

write.table(inputSeq1,paste0(folder_path,'supercrunch_family/genBank_metadata_complete_matchedUpham.txt'), quote = F, col.names = T, row.names = F, sep = '\t')

save.image(paste0(folder_path,'nameMatching/outputs/AllinputSeq_matchedUpham.RData'))
