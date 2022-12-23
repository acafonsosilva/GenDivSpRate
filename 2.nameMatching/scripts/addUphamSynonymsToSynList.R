library(knitr)
library(tidyverse)
library(taxize)
library(rentrez)
keyiucn <- 'd374fbdef6525e72f41962487622c73f4225bc48e7154c7f405d15ddd87ec59a'
set_entrez_key("78915277c9890b626a8b082c889ceca86a08")

folder_path <- "/data/biodiv/asilva/"

upham <- 
  read.table(paste0(folder_path,'nameMatching/inputs/upham_NCBI_matchup.txt'), header = T, sep = "\t")
  
syn <- readRDS(paste0(folder_path,'nameMatching/outputs/synonymsClean.rds')) %>% 
  distinct() %>%
  rename(Accepted = scientificName, Synonym = synonym)

#### Find missing species in phylogeny ####

treelabels <- upham %>% 
  transmute(treeSp = MasterTax_SciName) %>%
  distinct %>%
  mutate(syn_accepted = ifelse(treelabels$treeSp %in% 
                          str_replace(str_to_sentence(syn$Accepted)," ","_"),
                        'accepted','notAccepted'),
  syn_syn = ifelse(treelabels$treeSp %in% 
                     str_replace(str_to_sentence(syn$Synonym)," ","_"),
                   'synonym','notSynonym'),
  syn_syn_synSpp = ifelse(treelabels$treeSp %in% 
                            str_replace(str_to_sentence(syn$SynSpp)," ","_"),
                          'SynSpp','notsynSpp'))

treeNotFound <- treelabels %>% 
  filter(syn_accepted %in% 'notAccepted', 
         syn_syn %in% 'notSynonym',
         syn_syn_synSpp %in% 'notsynSpp') 


### Customized functions to make repetitive steps easier
### get whatever species it can find

get_tsn_binded <- function(sp,db, accepted){
  id <- c()
  for(i in unique(sp)){
    tryCatch({
      id <- c(id,suppressMessages(get_tsn_(i, db = db, accepted = accepted))) 
    }, error=function(e){cat("Error in",conditionMessage(e), "\n")})
  }
  id <- id %>%
    lapply(., function(x) data.frame(as.list(x), stringsAsFactors=F)) %>%
    bind_rows()
  
  return(id)
}

id2name_binded <- function(sp,db){
  id <- c()
  for(i in unique(sp)){
    tryCatch({
      id <- c(id,suppressMessages(id2name(i, db = db))) 
    }, error=function(e){cat("Error in",conditionMessage(e), "\n")})
  }
  id <- id %>%
    lapply(., function(x) data.frame(as.list(x), stringsAsFactors=F)) %>%
    bind_rows()
  
  return(id)
}

synonyms_binded <- function(sp, db){

  syns <- c()
  for(i in unique(sp)){
    tryCatch({
      syns <- c(syns,taxize::synonyms(i, db = db)) 
    }, error=function(e){cat("Error in",conditionMessage(e), "\n")})
  }
  syns <- syns %>%
    lapply(., function(x) data.frame(as.list(x), stringsAsFactors=F)) %>%
    bind_rows() %>%
    distinct()
  
  return(syns)
}


### link itis synonyms to treeNotFound by SynSpp
syntemp0 <- synonyms_binded(
  str_replace(unique(treeNotFound$treeSp),"_"," "), db = 'itis') %>% 
  filter(syn_name %in% str_replace(treeNotFound$treeSp,"_"," ")) %>%
  rename(SynSpp = syn_name, id = syn_tsn,
         synonym = acc_name, syn_tsn = acc_tsn) %>%
  select(-matches("sub_tsn|parent_tsn|acc_author|syn_author|NA_character_.")) %>%
  drop_na() 

syntemp <- id2name_binded(unique(syntemp0$id), db = 'itis') %>%
  right_join(., syntemp0, by = 'id') %>%                  
  select(-matches("sub_tsn|name|parent_tsn|acc_author|syn_author|NA_character_.")) %>%
  rename(synSpp_tsn = id, synSpp_rank = rank, synSpp_status = status)

syntemp <- id2name_binded(syntemp$syn_tsn, db = 'itis') %>%
  rename(syn_tsn = id, id = parent_tsn) %>%
  left_join(., id2name_binded(.$id, db = 'itis'), by = 'id') %>%
  rename(tsn = id, scientificName = name.y, rank = rank.x, status = status.x) %>%
  select(-matches("name.x|rank.y|status.y|parent_tsn|NA_character_.")) %>%
  right_join(.,syntemp, by = 'syn_tsn') %>%
  mutate(treeSp = str_replace(SynSpp," ","_")) 

### link iucn synonyms to treeNotFound with iucn_name
synIUCN0 <- synonyms_binded(
  str_replace(unique(treeNotFound$treeSp),"_"," "), db = 'iucn') %>% 
  filter(!synonym %in% str_replace(treeNotFound$treeSp,"_"," ")) %>%
  rename(iucn_syn = synonym, 
         iucn_id = accepted_id,
         iucn_name = accepted_name) %>% 
  select(-matches("authority|syn_authority|NA_character_.")) %>%
  drop_na() %>%
  distinct()

synIUCN1 <- get_tsn_binded(unique(synIUCN0$iucn_syn), 
                          db = 'itis', accepted = TRUE) 

synIUCN <- data.frame()
for(i in unique(synIUCN1$scientificName)){
  if(i %in% synIUCN0$iucn_name){
    dt0 <- synIUCN0[synIUCN0$iucn_name %in% i,]
    dt1 <- synIUCN1[synIUCN1$scientificName %in% 
                      dt0[dt0$iucn_name %in% i, 'iucn_name'],] %>%
      transmute(iucn_syn_tsn = tsn, iucn_name = scientificName, 
                iucnSyn_nameUsage = nameUsage) %>%
      right_join(., dt0, by = 'iucn_name') %>%
      mutate(treeSp = str_replace(iucn_name," ","_"))
    synIUCN <- bind_rows(synIUCN,dt1)
  }
  
  if(i %in% synIUCN0$iucn_syn){
    dt0 <- synIUCN0[synIUCN0$iucn_syn %in% i,]
    dt1 <- synIUCN1[synIUCN1$scientificName %in% 
                      dt0[dt0$iucn_syn %in% i, 'iucn_syn'],] %>%
      transmute(iucn_syn_tsn = tsn, iucn_syn = scientificName, 
                iucnSyn_nameUsage = nameUsage) %>%
      right_join(., dt0, by = 'iucn_syn') %>%
      mutate(treeSp = str_replace(iucn_name," ","_"))
    synIUCN <- bind_rows(synIUCN,dt1)
  }
}


### previous bit was not effective so this was just to fix invalids
for(i in unique(synIUCN[synIUCN$iucnSyn_nameUsage %in% 'invalid', 'iucn_syn'])){

  dt <- synIUCN[synIUCN$iucn_syn %in% i,]
  synIUCN1 <- get_tsn_binded(i,db = 'itis', accepted = TRUE) %>%
    filter(scientificName %in% dt[dt$iucn_syn %in% i,'iucn_name'])

  if(!is_empty(synIUCN1) & dim(synIUCN1)[1] == 1){
    synIUCN[synIUCN$iucn_name %in% dt[dt$iucn_syn %in% i,'iucn_name'], 
            c('iucn_syn_tsn', 'iucnSyn_nameUsage')] <-
      c(synIUCN1$tsn, synIUCN1$nameUsage)
  }
}


### link tol synonyms to treeNotFound with treeSp

synOTL <- data.frame()
for(i in str_replace(unique(treeNotFound$treeSp),"_"," ")){
  tolid <- suppressWarnings(get_tolid_(i, messages = FALSE)[[1]])
  tolid <- tolid[tolid$unique_name %in% i | tolid$matched_name %in% i, ]
  
  syn_otl <- NULL
  if(!is_empty(tolid[[1]])){
    if(!is.null(tolid[[1]])){
      syn_otl <- unlist(
        rotl::taxonomy_taxon_info(ott_ids = tolid$ott_id)[[1]]$synonyms)
      
      if(any(!tolid$unique_name %in% tolid$matched_name)){
        syn_otl <- unique(c(syn_otl, tolid$unique_name)
                          [!c(syn_otl, tolid$unique_name) %in% i])
        
        if(dim(tolid)[1] > 1){
          tolid <- tolid %>% filter(!unique_name %in% matched_name)
        }
      }
    }
  }

  if(!is.null(syn_otl)){
    if(!all(tolid$unique_name %in% syn_otl)){
      synOTL <- bind_rows(synOTL, 
                          data.frame(otl_id = tolid$ott_id,
                                     otl_name = tolid$unique_name,
                                     otl_syn = syn_otl[!syn_otl %in% 
                                                         tolid$unique_name],
                                     treeSp = str_replace(i," ","_"), 
                                     stringsAsFactors = FALSE))
    }
  }
}

 
##### match different synonyms to treeNotFound
# otl_id only meaningful when there is an iucn id
treeNotFound2 <- treeNotFound %>%
  left_join(. ,syntemp, by = 'treeSp') %>%
  left_join(., synIUCN, by = 'treeSp') %>%
  left_join(., synOTL, by = 'treeSp') %>%
  mutate(neverFound = case_when(is.na(syn_tsn) & is.na(tsn) & 
                                is.na(synSpp_tsn) & is.na(iucn_syn_tsn) & 
                                is.na(iucn_id) ~ "neverFound",
                                TRUE ~ 'somethingFound')) %>%
  rename(Accepted = scientificName,
         Synonym = synonym) %>% distinct()
  
## some fuckup before when keeping the species names only kept the genus
l <- unique(treeNotFound2[which(str_count(treeNotFound2$Accepted, "\\S+") == 1),'treeSp'])
for(i in l){
  id <- id2name_binded(unique(treeNotFound2[treeNotFound2$treeSp %in% i, 'syn_tsn']), db = 'itis')
  treeNotFound2[treeNotFound2$treeSp %in% i, 'tsn'] <- id$id
  treeNotFound2[treeNotFound2$treeSp %in% i, 'Accepted'] <- id$name
  treeNotFound2[treeNotFound2$treeSp %in% i, 'syn_tsn'] <- treeNotFound2[treeNotFound2$treeSp %in% i, 'synSpp_tsn']
  treeNotFound2[treeNotFound2$treeSp %in% i, 'Synonym'] <- treeNotFound2[treeNotFound2$treeSp %in% i, 'SynSpp']
  treeNotFound2[treeNotFound2$treeSp %in% i, 'synSpp_tsn'] <- NA
  treeNotFound2[treeNotFound2$treeSp %in% i, 'SynSpp'] <- NA
}


### not all synonyms match the valid subspecies synonyms
#to try to improve check if rearranging the iucn name with synonym gives an accepted itis subspecies name that may or may not be listed in the synonyms names
for (i in unique(treeNotFound2[!is.na(treeNotFound2$iucn_name),'iucn_name'])){
  dt <- treeNotFound2[treeNotFound2$iucn_name %in% i,]
  
  if(any(word(dt$iucn_syn,1) == word(dt$iucn_name,1))){
    
    spp <- paste(dt$iucn_syn,word(dt$iucn_name,2)) 
    spp <- spp[word(spp,1) %in% unique(word(dt$iucn_name,1))]
    
    t <- get_tsn_binded(spp, db = 'itis', accepted = TRUE)

    if(!is_empty(t)){
      if(dim(t)[1] != 0){
      if(t$scientificName %in% syn$Synonym & t$nameUsage %in% 'valid'){
        syn[syn$Synonym %in% t$scientificName,c('SynSpp', 'synSpp_tsn')] <-
          c(dt$iucn_name, 'iucn')

        treeNotFound2[treeNotFound2$iucn_name %in% i,
                      c('tsn','Accepted','Synonym',
                        'syn_tsn', 'SynSpp','synSpp_tsn')] <-
          syn[syn$Synonym %in% t$scientificName,
              c('tsn','Accepted','Synonym',
                'syn_tsn', 'SynSpp','synSpp_tsn')]
      }
      }
    }
  }
}

### some random fuckup with iucn_syn_tsn %in% 'valid'
## one should be retrieved the tsn and the other excluded
## treeNotFound2[treeNotFound2$treeSp %in% treeNotFound2[treeNotFound2$iucn_syn_tsn %in% 'valid','treeSp'],]
treeNotFound2 <- treeNotFound2[!treeNotFound2$iucn_syn_tsn %in% 'valid',]
treeNotFound2[!treeNotFound2$iucnSyn_nameUsage %in% c('invalid','valid') & !is.na(treeNotFound2$iucnSyn_nameUsage ),'iucnSyn_nameUsage'] <- 'valid'

toAdd <- data.frame()
c <- as.character(colnames(syn))[!as.character(colnames(syn)) %in% c('commonNames', 'Family','nameUsage')]
for(i in unique(treeNotFound2$treeSp)){
  dt <- treeNotFound2[treeNotFound2$treeSp %in% i,] %>% distinct()

  ## add tol synonyms as SynSpp to synonyms already there, ignore the rest
  if(all(dt$neverFound %in% 'neverFound')){
    dt <- dt[!is.na(dt$otl_id),]
    
    if(dim(dt)[1] != 0){ ## the only labels that nothing can be done...
    if(any(dt$otl_syn %in% syn$Synonym)){
      if(any(is.na(syn[syn$Synonym %in% dt$otl_syn, 'SynSpp']))){
        syn[syn$Synonym %in% dt$otl_syn, 'SynSpp'] <- dt[dt$otl_syn %in% syn$Synonym,'otl_name']
        syn[syn$Synonym %in% dt$otl_syn, 'synSpp_tsn'] <- 'otl'
      }
    }
    
    if(all(dt$otl_syn %in% syn$Accepted)){
      dt1 <- dt %>% filter(otl_syn %in% syn$Accepted) %>%
        transmute(tsn = unique(syn[syn$Accepted %in% dt$otl_syn,'tsn']),
                  Accepted = otl_syn,
                  syn_tsn = 'otl',
                  Synonym = otl_name,
                  nameUsage = 'invalid') %>%
        distinct()
      toAdd <- bind_rows(toAdd, dt1)
    }
    if(all(!dt$otl_syn %in% syn$Accepted & !dt$otl_syn %in% syn$Synonym)){
      dt1 <- dt %>%
        transmute(tsn = 'otl',
                  Accepted = otl_syn,
                  syn_tsn = 'otl',
                  Synonym = otl_name,
                  nameUsage = 'invalid') %>%
        distinct()
      if(!str_detect(dt1$Accepted, ' sp. ')){
        toAdd <- bind_rows(toAdd, dt1)
      }
    }
    }
  }

  
  if(any(dt$neverFound %in% 'somethingFound')){
    ### if valid tsn needs to be included in syn
    
    if(all(is.na(dt$Synonym))){

      if(any(dt$iucn_syn %in% syn$Synonym)){
        dtu <- dt[dt$iucn_syn %in% syn$Synonym,] %>% 
          select(iucn_syn, iucn_name) %>% distinct()
        if(dim(dtu)[1] == 1){
          if(!dtu$iucn_syn %in% syn$Accepted){  ## make sure that they match with synonyms
          syn[syn$Synonym %in% dtu$iucn_syn,c('SynSpp','synSpp_tsn')] <- c(dtu$iucn_name, 'iucn') ### when iucn name matches syn synonyms that already exists
          }
        }
        if(dim(dtu)[1] != 1){
          if(any(dtu$iucn_syn %in% syn$Accepted)){
            dt <- dt %>% filter(iucn_syn_tsn %in% syn$tsn)
            if(all(dt$iucn_syn %in% syn$Accepted)){
              dt1 <- dt %>%
                transmute(tsn = iucn_syn_tsn,
                          Accepted = iucn_syn,
                          syn_tsn = 'iucn',
                          Synonym = iucn_name,
                          nameUsage = 'problSyn')
              toAdd <- bind_rows(toAdd, dt1)
            }
            
          }
        }
      }
      
      if(any(!dt$iucn_syn %in% syn$Accepted & 
             dt$iucnSyn_nameUsage %in% 'valid')){
        dt1 <- dt %>% 
          transmute(tsn = iucn_syn_tsn,
                    Accepted = iucn_syn,
                    syn_tsn = 'iucn',
                    Synonym = iucn_name,
                    nameUsage = iucnSyn_nameUsage) %>%
          distinct()
        
        if(dim(dt1)[1] == 1){
          toAdd <- bind_rows(toAdd, dt1)
        }
        if(dim(dt1)[1] != 1){
          dt1 <- dt %>% filter(iucnSyn_nameUsage %in% 'valid') %>%
            transmute(tsn = iucn_syn_tsn,
                      Accepted = iucn_syn,
                      syn_tsn = 'iucn',
                      Synonym = iucn_name,
                      nameUsage = iucnSyn_nameUsage) %>%
            distinct()
          
          if(dim(dt1)[1] == 1){
            toAdd <- bind_rows(toAdd, dt1)
          }
          if(dim(dt1)[1] != 1){
            dt1 <- dt %>% filter(iucnSyn_nameUsage %in% 'valid') %>%
              transmute(tsn = iucn_syn_tsn,
                        Accepted = iucn_syn,
                        syn_tsn = 'iucn',
                        Synonym = iucn_name,
                        nameUsage = 'problSyn') %>%
              distinct()
            toAdd <- bind_rows(toAdd, dt1)
          }
        }
      }
    }

    if(any(dt$status %in% 'valid')){
      dt1 <- dt %>% select(c) %>% distinct()
      if(!all(is.na(dt1$Synonym))){
        toAdd <- bind_rows(toAdd, dt1)
      }
      
    }
    
    if(any(dt$iucnSyn_nameUsage %in% 'invalid' & is.na(dt$tsn))){
      dt1 <- dt %>% 
        transmute(Accepted = iucn_name,
                  tsn = "iucn",
                  syn_tsn = "iucn",
                  Synonym = iucn_syn,
                  nameUsage = 'invalid') %>%
        distinct()
      toAdd <- bind_rows(toAdd, dt1)
      
    }
  }
}

syn2 <- syn
toAdd <- unique(toAdd) %>%
  filter(!str_detect(Accepted, " sp. "), 
         !str_detect(Accepted, "\\?"), 
         !str_detect(Accepted, "-"),
         !str_detect(Accepted, "\\("),
         !str_detect(Synonym, "\\["),
         !str_detect(Synonym, "\\<i>"),
         !str_detect(Synonym, "\\."),
         str_count(Accepted, "\\S+") <= 2) %>%
  group_by(Accepted) %>% 
  summarise_all(coalesce_by_column) %>%
  data.frame()

syn <- bind_rows(syn, toAdd)

### define extinct species, at least based on Upham list
syn[syn$Synonym %in% 
      str_replace(unique(c( 
        as.character(upham[upham$Source %in% 'newSpExt', 'MasterTax_SciName']), 
        as.character(upham[upham$Source %in% 'newSpExt', 'NCBI_SciName']))),
        "_"," "), 'nameUsage'] <- 'extinct'

### treat duplicates in synonyms (only problems should be dups)
c <- syn[syn$Synonym %in% syn[duplicated(syn$Synonym),'Synonym'] & 
           !syn$nameUsage %in% c('problSyn','extinct'),] %>%
  arrange(Synonym) 
c <- data.table::setDT(c)[, lapply(.SD, na.omit), by = 'Accepted'] %>%
  distinct() %>%
  data.frame()
syn <- bind_rows(syn[!syn$Synonym %in% c$Synonym,], c)
  
### merge iucn 'accepted' that match with available synonyms
### change nameUsage
for (i in unique(syn[syn$tsn %in% 'iucn','Synonym'])){
  dt1 <- NULL
  dt0 <- syn[syn$Synonym %in% i,]
  if(dim(dt0)[1] > 1){
    dt <- dt0 %>% filter(!Accepted %in% SynSpp)
    if(dim(dt)[1] > 1 & any(!dt$nameUsage %in% 'invalid')){
      dt1 <- dt %>% filter(!tsn %in% 'iucn') %>%
        bind_rows(., tibble(
          dt[!dt$tsn %in% 'iucn',1:6], 
          SynSpp = dt[dt$tsn %in% 'iucn','Accepted'], 
          nameUsage = dt[!dt$tsn %in% 'iucn','nameUsage'],
          synSpp_tsn = 'iucn'))
      
      syn <- bind_rows(syn[!syn$Synonym %in% i,],dt1)
    }
    
    if(dim(dt)[1] > 1 & all(!dt$nameUsage %in% 'valid') & is.null(dt1)){
      syn[syn$Synonym %in% i,'nameUsage'] <- 'iucn_problSyn'
    }
    
    dt1 <- NULL
    if(dim(dt)[1] == 1 & all(dt$nameUsage %in% 'valid')){
      syn <- syn[!syn$Accepted %in% 
                   dt0[dt0$Accepted %in% dt0$SynSpp,'Accepted'], ]
    }
  }
}

c <- syn[syn$Synonym %in% syn[duplicated(syn$Synonym),'Synonym'] & 
           !syn$nameUsage %in% c('problSyn','extinct'),] %>%
  arrange(Synonym) 
c <- data.table::setDT(c)[, lapply(.SD, na.omit), by = 'Accepted'] %>%
  distinct() %>%
  data.frame() %>%
  filter(nameUsage %in% 'invalid')

syn[syn$Accepted %in% c$Accepted, 'nameUsage'] <- paste0(c$tsn,'_problSyn')

### there are still NAs in nameUsage for valid species

for(i in syn[is.na(syn$nameUsage),'Accepted']){
  dt <- syn[syn$Accepted %in% i,]
  
  if(any(str_count(dt$Synonym, "\\S+") == 3 & is.na(dt$nameUsage))){
    syn[syn$Accepted %in% i &
          str_count(syn$Synonym, "\\S+") == 3 &
          is.na(syn$nameUsage), 'nameUsage'] <- 'subspecies'
  }
  
  if(any(str_count(dt$Synonym, "\\S+") == 2 & is.na(dt$nameUsage))){
    if(dim(syn[syn$Synonym %in% dt$Synonym,])[1] == 1){
      syn[syn$Accepted %in% i &
            str_count(syn$Synonym, "\\S+") == 2 &
            is.na(syn$nameUsage), 'nameUsage'] <- 'synonymised'
    }
  }
}

### Add family to missing rows
for(i in syn[is.na(syn$Family) & 
             !syn$nameUsage %in% 'extinct' &
             !syn$tsn %in% 'iucn','Accepted']){
  dt <- syn[syn$Accepted %in% i,]
  
  if(dim(dt)[1] > 1 & any(!is.na(dt$Family))){
    syn[syn$Accepted %in% i & 
          is.na(syn$Family), 'Family'] <- unique(dt$Family)[!is.na(unique(dt$Family))]
    
  }
}

### there are Accepted names that should be in Synonym as valid names as well
d <- syn[!syn$Accepted %in% syn$Synonym & !syn$tsn %in% c('iucn','otl'),] %>%
  mutate(syn_tsn = tsn,
         Synonym = Accepted,
         nameUsage = 'valid')

syn[!syn$Accepted %in% syn$Synonym & 
      !syn$tsn %in% c('iucn','otl') & 
      syn$nameUsage %in% 'valid', 'nameUsage'] <- 'iucn_valid'

syn <- bind_rows(syn,d) 

#### still issues with the iucn species 
## remove any iucn species that can be added as SynSpp
for (i in syn[syn$tsn %in% 'iucn','Accepted']){
  s <- syn[syn$tsn %in% 'iucn' & syn$Accepted %in% i,'Synonym']
  if(any(is.na(syn[syn$Synonym %in% i & 
                   !syn$nameUsage %in% 'valid','SynSpp'])) &
     length(s) == 1){
    
    syn[syn$Synonym %in% i & 
          !syn$nameUsage %in% 'valid' & 
          is.na(syn$SynSpp),'SynSpp'] <- s
    
    if(syn[syn$tsn %in% 'iucn' & syn$Accepted %in% i,'syn_tsn'] %in% 'iucn'){
      syn[syn$Synonym %in% i & 
            !syn$nameUsage %in% 'valid' & 
            is.na(syn$SynSpp),'synSpp_tsn'] <-  'iucn'
    }
    syn <- syn[!(syn$tsn %in% 'iucn' & syn$Accepted %in% i),]
  }
}

## add an ID to syn to make it easier to match datasets
syn <- syn %>% 
  arrange(Accepted, syn_tsn) %>% 
  mutate(ID = seq(1,length(syn$tsn)))

## a couple of more issues...
##syn[!is.na(syn$synSpp_tsn) & !syn$synSpp_tsn %in% 'iucn' & !syn$synSpp_tsn %in% 'otl', ]
syn[c(2844,9774), 'synSpp_tsn'] <- NA

saveRDS(syn, paste0(folder_path,'nameMatching/outputs/synonymsClean_uphamUpdated.rds'))
