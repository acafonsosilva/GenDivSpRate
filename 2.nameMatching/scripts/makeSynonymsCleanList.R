library(tidyverse)
library(taxize)

### Meyer et al. (2015) synonyms list is a bit problematic ####
##Accepted species are not updated: remove species that are not valid, include all subspecies in synonyms and not include valid species (that are not subspecies) in the synonym column 

folder_path <- "/data/biodiv/asilva/"

syn <- read.delim(paste0(folder_path,'nameMatching/inputs/df_Syns_mamm_update2016.txt')
                  , sep = "\t", stringsAsFactors = FALSE)

#####
synSp <- unique(syn$Accepted)

temp <- data.frame()
for(i in synSp){
  t0 <- get_tsn_(i, db = 'itis', 
                 accepted = TRUE, 
                 messages = FALSE)[[i]]
  
  t <- data.frame(Accepted = i,tsn = NA,scientificName = NA, commonNames = NA, nameUsage = 'invalid', stringsAsFactors = FALSE)
  
  if(!is.null(t0)){
    if(dim(t0)[1] != 0){
      t1 <- data.frame(Accepted = i, t0)
      t1 <- t1[str_to_lower(t1$scientificName) %in% i & 
                 !t1$nameUsage %in% 'accepted',]
      if(dim(t1)[1] != 0){
        t <- t1
      }
    }
  }
  temp <- bind_rows(temp, t) 
  print(paste(which(synSp %in% i),'of',length(synSp)))
}

synClean <- temp[!temp$nameUsage %in% 'accepted',] %>%
  mutate(rank = case_when(
    str_count(scientificName, "\\S+") == 2 ~ 'species', 
    str_count(scientificName, "\\S+") == 3 ~ 'subspecies',
    TRUE ~ 'species'))

### loop through each accepted name and remove the rows that are subspecies as well the ones that the scientific name doesn't match the accepted

for (i in unique(synClean$Accepted)){
  r <- as.numeric(row.names(synClean[synClean$Accepted %in% i,]))
  if(any(synClean[r,'nameUsage'] %in% 'valid')){
    a <- synClean[r,] %>% 
      filter(rank %in% 'species', str_to_lower(scientificName) %in% Accepted)
    synClean[r[!r %in% which(synClean$tsn %in% a$tsn)],] <- NA
  }
  print(which(unique(synClean$Accepted) %in% i))
}

synClean <- synClean %>% drop_na(Accepted) %>% distinct()

#table(synClean$nameUsage) ##441 supposedly accepted species are invalid...

### add family information ####

for(i in synClean[synClean$nameUsage %in% 'valid','scientificName']){
  
  a <- tax_name(query = i, 
                get = "family", 
                db = "itis", 
                accepted = TRUE, 
                messages = FALSE)['family']
  
  synClean[synClean$scientificName %in% i,'Family'] <- a
  
  print(paste(which(synClean$scientificName %in% i),'of',length(synClean$scientificName)))
}

### manage subspecies information: ####

temp <- data.frame()
for(i in synClean[!is.na(synClean$tsn),'tsn']){
  t0 <- downstream(i, 'itis', downto = 'subspecies')[[i]]
  
  if(any(!t0$tsn %in% 'No data')){
    temp <- bind_rows(temp, t0) 
  }
  print(paste(which(synClean$tsn %in% i),'of',length(synClean$tsn)))
}
subspp <- temp[-6]
colnames(subspp) <- c('syn_tsn','scientificName','tsn','rank','synonym')

### manage synonyms information: ####
## needs to search for subspecies synonyms as well

temp1 <- data.frame()
for(i in synClean[!is.na(synClean$tsn),'tsn']){
  t0 <- synonyms(i, db = 'itis', message = FALSE)[[i]]
  
  if(!is_empty(t0)){
    temp1 <- bind_rows(temp1, t0) 
  }
  print(paste(which(synClean$tsn %in% i),'of',length(synClean$tsn)))
}
syno <- temp1[,-c(1,3)]
colnames(syno) <- c('tsn','synonym','syn_tsn')

### find accepted names of invalid species ####
temp2 <- data.frame()
for(i in synClean[synClean$nameUsage %in% 'invalid','Accepted']){
  t0 <- synonyms(i, db = 'itis', message = FALSE)[[i]]
  
  if(!is.na(t0)){
    temp2 <- bind_rows(temp2, data.frame(Accepted = i,t0, stringsAsFactors = FALSE )) 
  }
  print(paste(which(synClean$Accepted %in% i),'of',length(synClean$Accepted)))
}
invalids <- temp2[,c(3,4,7,8)]
colnames(invalids) <- c('scientificName','tsn','synonym','syn_tsn')

# check acc_name corresponds to valid species
v <- id2name(unique(temp2$acc_tsn), db = 'itis') %>%
  bind_rows()

#### the non valid names ####
## bind subspp, syno and invalids ## 
nonvalid <- bind_rows(subspp, invalids,syno) %>%
  rename(nameUsage = rank) %>%
  distinct()

nonvalid1 <- left_join(nonvalid,synClean[,c('tsn','scientificName','commonNames','Family')], by = 'tsn') 

nonvalid1$scientificName <- NA
for (i in 1:length(nonvalid1$scientificName.x)){
  nonvalid1[i,'scientificName'] <- ifelse(
    is.na(nonvalid1[i,'scientificName.x']), nonvalid1[i,'scientificName.y'],
    nonvalid1[i,'scientificName.x'])
}

nonvalid1 <- nonvalid1 %>% 
  select(tsn,scientificName,commonNames,Family,syn_tsn,synonym,nameUsage) %>%
  distinct()

#### join synClean with nonvalids ####
## synClean needs to have a synonym species in which valid means scientificName = synonym

synClean <- synClean %>% 
  filter(nameUsage %in% 'valid') %>%
  mutate(syn_tsn = tsn,
         synonym = scientificName) %>%
  select(tsn,scientificName,commonNames,Family,syn_tsn,synonym,nameUsage)

synCleanComp <- bind_rows(synClean,nonvalid1) %>%
  distinct()

#### complete family information missing: ####
for(i in synCleanComp[is.na(synCleanComp$Family),'scientificName']){
  a <- tax_name(query = i, get = "family", db = "itis", 
                accepted = TRUE, messages = FALSE)['family']
  
  synCleanComp[synCleanComp$scientificName %in% i,'Family'] <- a
}

### exclude the rows that just differ in tsn_syn, just because of different sources for the same synonym name 
synCleanComp <- synCleanComp %>% 
  distinct_at(vars(tsn, scientificName, commonNames, Family, synonym), .keep_all = TRUE)


#### set nameUsage according with rules ####

## synonymised - if matching scientificName with synonym and synonym with scientificName gives only one hit
## but a couple of synonyms might be just due to more than one change in name

tosynonymise <- synCleanComp %>% 
  filter(scientificName %in% synonym, !nameUsage %in% c('valid','subspecies')) %>%
  select(scientificName) %>%
  distinct()


tosynonymise1 <- data.frame()
for(i in tosynonymise$scientificName){
  a <- synCleanComp[synCleanComp$scientificName %in% 
                      i,]
  if(dim(a)[1] > 1){
    b <- synCleanComp[synCleanComp$synonym %in% a$synonym, ]
    if(any(is.na(b$nameUsage))){
      if(dim(b[is.na(b$nameUsage),])[1] == 1){
        tosynonymise1 <- bind_rows(tosynonymise1,
                                   data.frame(
                                     b[is.na(b$nameUsage),
                                       c('scientificName','syn_tsn')], 
                                     nameUsage1 = 'synonymised', stringsAsFactors = FALSE))
      }
      if(dim(b[is.na(b$nameUsage),])[1] != 1){
        tosynonymise1 <- bind_rows(tosynonymise1,
                                   data.frame(
                                     b[is.na(b$nameUsage),
                                       c('scientificName','syn_tsn')], 
                                     nameUsage1 = 'lumped', stringsAsFactors = FALSE))
      }
    }
  }
}

tosynonymise <- left_join(tosynonymise,tosynonymise1, by = 'scientificName')

#### remove subspecies names in the scientificName column ####
## the subspecies that match synonym names that are not valid
## make list of subspecies synonym species names of accepted synonyms

rmScinam <- id2name(unique(synCleanComp[synCleanComp$scientificName %in% 
                                          tosynonymise[is.na(tosynonymise$nameUsage1),'scientificName'],'tsn']), 
                    db = 'itis') %>% bind_rows()  ### all subspecies are valid

### check synonyms names of subspecies names are invalid species
rmScinam0 <- id2name(synCleanComp[synCleanComp$scientificName %in% rmScinam$name, 'syn_tsn'],
                     db = 'itis') %>% bind_rows() %>% 
  rename(SynSpp = name, synSpp_tsn = id, synSpp_rank = rank, synSpp_status = status) %>%
  mutate(name = rmScinam$name) %>% 
  select(-parent_tsn)

rmScinam <- left_join(rmScinam, rmScinam0, by = 'name') %>% 
  mutate(scientificName = word(name,1,2)) %>%
  rename(synonym = name, syn_tsn = id, tsn = parent_tsn)

## remove subspecies from scientificName 
## later include flag for these synonym names
synCleanComp <- synCleanComp[!synCleanComp$scientificName %in% rmScinam$synonym,]
tosynonymise <- tosynonymise[!is.na(tosynonymise$nameUsage1),]

## for 27 species 

### there are subspecies names not flagged as subspecies because they have invalid status. For some of these the species don't have valid subspecies names, so should keep these but flag them as invalidspp

synCleanComp[is.na(synCleanComp$nameUsage) & str_count(synCleanComp$synonym, "\\S+") == 3,'nameUsage'] <- 'invalidspp'  ### don't forget these are likely real synonyms but possibly matching with other valid species name... need to check again rmScinam


### some synonym names have a different syn_tsn - seems to be the case where these synonym names match different valid species. Need to match all synonyms, valid and invalid, to these synonyms, so the valid species with this synonym can be flagged. 

## pull problematic synonyms  - 17 species names
problSyn <- synCleanComp[synCleanComp$synonym %in% synCleanComp[duplicated(synCleanComp$synonym),'synonym'],'synonym'] %>% 
  unique()

temp <- data.frame()
for (i in problSyn){
  ## extract all tsn numbers matching that synonym
  a <- get_tsn_(i, db = 'itis', message = FALSE, accepted = FALSE)[[i]]
  ## extract all synonyms matching these tsn ids and flag the problematic syns
  b <- synonyms(a$tsn, db = 'itis') %>% bind_rows() %>% mutate(nameUsage = ifelse(syn_tsn %in% a$tsn, 'problSyn',NA))
  
  ## back to main list, flag synonyms. 
  synCleanComp[synCleanComp$syn_tsn %in% b[!is.na(b$nameUsage),'syn_tsn'],'nameUsage'] <- 'problSyn'
  
  ##If not in the list yet, make synonyms list of other synonym names (both with flag and no flag).
  
  subb <- b[!b$syn_tsn %in% synCleanComp$syn_tsn,]
  
  if(dim(subb)[1] != 0){
    if(any(is.na(subb$acc_name))){
      asubb <- id2name(subb[is.na(subb$acc_name), 'acc_tsn'], db = 'itis') %>% 
        bind_rows() %>%
        rename(acc_tsn = id, acc_name = name)
      subb1 <- left_join(subb,asubb[1:2], by = "acc_tsn") %>% 
        mutate(acc_name = coalesce(acc_name.x, acc_name.y)) %>% 
        select(-acc_name.x, -acc_name.y)
    }
    subb1 <- id2name(subb[, 'acc_tsn'], db = 'itis') %>% 
      bind_rows() %>%
      rename(acc_tsn = id, acc_name = name) %>%
      select(acc_tsn,acc_name) %>%
      right_join(., subb, by = "acc_tsn") %>% 
      mutate(acc_name = coalesce(acc_name.x, acc_name.y)) %>% 
      select(-acc_name.x, -acc_name.y)
    
    f <- suppressWarnings(
      tax_name(query = unique(subb1$acc_name), get = "family", db = "itis",
               accepted = TRUE, messages = FALSE)) %>%
      rename(acc_name = query)
    
    if(any(f$family %in% synCleanComp$Family)){
      subb1 <- subb1 %>% 
        left_join(., f[-1], by = 'acc_name') %>%
        filter(family %in% synCleanComp$Family)
      temp <- bind_rows(temp,subb1)
    }
  }
}

### only 6 possible new problems remain
temp <- temp %>% select(-sub_tsn) %>%
  distinct() %>%
  rename(id = syn_tsn) %>%
  left_join(., bind_rows(id2name(.$id, db = 'itis')), by = 'id') %>%
  select(-name, -parent_tsn, -acc_author,-syn_author) 

## a couple of subspecies name should be added to rmScinam
rmScinam0 <- temp %>% 
  filter(str_count(acc_name, "\\S+") == 3) %>%
  rename(synonym = acc_name, syn_tsn = acc_tsn,SynSpp = syn_name,synSpp_tsn = id, synSpp_rank = rank, synSpp_status = status) 

rmScinam0 <- id2name(rmScinam0$syn_tsn, db = 'itis') %>% 
  bind_rows() %>%
  rename(syn_tsn = id, id = parent_tsn) %>%
  left_join(., bind_rows(id2name(.$id, db = 'itis')), by = 'id') %>%
  rename(tsn = id, scientificName = name.y, rank = rank.x, status = status.x) %>%
  select(-name.x, -rank.y, -status.y, -parent_tsn) %>%
  right_join(.,rmScinam0, by = 'syn_tsn') 

rmScinam <- bind_rows(rmScinam, rmScinam0) %>% select(-nameUsage,-family)

## add other synonyms names as invalid names to main list
## 1 name is invalid, the others are from a new species not in main list
tempinvsp <- temp %>% 
  filter(str_count(acc_name, "\\S+") != 3) %>%
  rename(synonym = acc_name, syn_tsn = acc_tsn,SynSpp = syn_name,synSpp_tsn = id, synSpp_rank = rank, synSpp_status = status) 

synCleanComp[synCleanComp$synonym %in% tempinvsp$SynSpp,'nameUsage'] <- 'synonymised' ## after checking that there was no other accepted name with the same synonym

tempinvsp <- tempinvsp %>%
  filter(!synSpp_tsn %in% tempinvsp[tempinvsp$SynSpp %in% synCleanComp$synonym,'synSpp_tsn']) %>%
  rename(tsn = syn_tsn, scientificName = synonym, Family = family, syn_tsn = synSpp_tsn, synonym = SynSpp) %>%
  mutate(nameUsage = c('synonymised','invalidspp','invalidspp')) %>%
  select(-synSpp_rank, -synSpp_status) %>%
  bind_rows(., data.frame(tsn = '930330', synonym = 'Prosciurillus topapuensis', 
                          syn_tsn = '930330', nameUsage = 'valid', 
                          scientificName = 'Prosciurillus topapuensis',
                          Family = 'Sciuridae', stringsAsFactors = FALSE))

synCleanComp <- bind_rows(synCleanComp, tempinvsp)

## all of these problems sorted - proceed to check remaining nameUsage = NA

##### check synonyms %in% scientificName ####

tosynonymiseRev <- synCleanComp %>% 
  filter(is.na(nameUsage)) %>%
  select(synonym) %>%
  distinct()

tosynonymiseRev1 <- data.frame()   
for (i in tosynonymiseRev$synonym){
  a <- synCleanComp[synCleanComp$synonym %in% 
                      i,]
  b <- synCleanComp[synCleanComp$synonym %in% a$synonym, ]
  if(any(is.na(b$nameUsage))){
    if(dim(b[is.na(b$nameUsage),])[1] == 1){
      tosynonymiseRev1 <- bind_rows(tosynonymiseRev1,
                                    data.frame(
                                      b[is.na(b$nameUsage),
                                        c('scientificName','syn_tsn')], 
                                      nameUsage2 = 'synonymised', stringsAsFactors = FALSE))
    }
  }
}
length(tosynonymiseRev$synonym) == length(tosynonymiseRev1$scientificName) ## all are synonymised

### there are still valid subspecies names in accepted list - remove those and add whichever can be find

x <- synCleanComp %>% filter(is.na(nameUsage), str_count(scientificName, "\\S+") == 3 ) %>% pull(scientificName) %>% unique()

synCleanComp <- synCleanComp[!synCleanComp$scientificName %in% x, ]

toAdd <- unique(synCleanComp[synCleanComp$scientificName %in% x & is.na(synCleanComp$nameUsage),'scientificName']) %>% 
  get_tsn_(.,db = 'itis', accepted = TRUE) %>% bind_rows() %>%
  filter(!scientificName %in% synCleanComp$synonym, nameUsage %in% 'valid') %>%
  mutate(synonym = scientificName, scientificName = word(scientificName,1,2))

toAdd1 <- get_tsn_(toAdd1$scientificName, db = 'itis', accepted = TRUE) %>% bind_rows() %>%
  filter(!tsn %in% toAdd$tsn, !scientificName %in% synCleanComp$synonym) %>%
  mutate(synonym = scientificName, tsn = syn_tsn, scientificName = word(scientificName,1,2)) %>%
  distinct()

for (i in unique(toAdd1$scientificName)){
  t <- toAdd1[toAdd1$scientificName %in% i,'tsn']$tsn[1]
  toAdd1[toAdd1$scientificName %in% i,'tsn'] <- toAdd[toAdd$scientificName %in% i,'tsn']$tsn[1]
  
  toAdd1[toAdd1$synonym %in% i,'syn_tsn'] <- toAdd[toAdd$scientificName %in% i,'tsn']$tsn[1]
}

for(i in unique(toAdd1$scientificName)){
  a <- tax_name(query = i, get = "family",  db = "itis", 
                accepted = TRUE, messages = FALSE)['family']
  b <- toAdd1[toAdd1$scientificName %in% i,]
  toAdd1[toAdd1$scientificName %in% i,'Family'] <- a
  toAdd1[toAdd1$synonym %in% i,'nameUsage'] <- 'valid'
  toAdd1[toAdd1$syn_tsn %in% b$syn_tsn & is.na(toAdd1$nameUsage),'nameUsage'] <- 'subspecies'
}

synCleanComp <- bind_rows(synCleanComp, toAdd1)

### somehow there are still scientificName not in synonym with a valid as nameUsage

val <- synCleanComp[!synCleanComp$scientificName %in% synCleanComp$synonym,] %>% 
  distinct_at(vars(scientificName),.keep_all = TRUE) %>%
  mutate(syn_tsn = tsn, synonym = scientificName, nameUsage = 'valid')

synCleanComp <- bind_rows(synCleanComp, val)

## remove the species that had a different syn_tsn and now are not in the main list anymore. This step has been moved to be done after extracting family name, to avoid this problem now. 

d <- full_join(tosynonymise, tosynonymiseRev1, by = c('scientificName','syn_tsn')) %>%
  filter(syn_tsn %in% synCleanComp[is.na(synCleanComp$nameUsage), 'syn_tsn']) %>%
  mutate(status = paste(nameUsage1,nameUsage2,sep = "_"))

for (i in unique(d$scientificName)){
  a <- d[d$scientificName %in% i,]
  
  if(all(a$status %in% 'lumped_synonymised')){
    synCleanComp[synCleanComp$syn_tsn %in% a$syn_tsn,'nameUsage'] <- 'lumped'
    d[d$syn_tsn %in% a$syn_tsn, 'status'] <- 'lumped'
  }
  
  if(all(a$status %in% 'synonymised_synonymised')){
    if(length(synCleanComp[synCleanComp$scientificName %in% a$scientificName,'nameUsage']) == 2) {
      synCleanComp[synCleanComp$syn_tsn %in% a$syn_tsn,'nameUsage'] <- 'synonymised'
      d[d$syn_tsn %in% a$syn_tsn, 'status'] <- 'synonymised'
    }
    
    n <- synCleanComp[synCleanComp$scientificName %in% a$scientificName,] %>%
      filter(is.na(nameUsage))
    if(dim(n)[1] == 1){
      synCleanComp[synCleanComp$syn_tsn %in% a$syn_tsn,'nameUsage'] <- 'synonymised'
      d[d$syn_tsn %in% a$syn_tsn, 'status'] <- 'synonymised'
    }
  }
  
  if(all(a$status %in% 'NA_synonymised')){
    ### double check these synonyms don't match more than one species
    b <- synCleanComp[synCleanComp$scientificName %in% i,] %>%
      filter(!nameUsage %in% 'valid')
    
    if(dim(b)[1] > 1){
      for(j in b$synonym){
        c <- synCleanComp[synCleanComp$scientificName %in% j,] #to be sure
        if(dim(c)[1] == 0 & is.na(b[b$synonym %in% j,'nameUsage'])){
          synCleanComp[synCleanComp$syn_tsn %in% b[b$synonym %in% j,'syn_tsn'],'nameUsage'] <- 'lumped'
          d[d$syn_tsn %in% b[b$synonym %in% j,'syn_tsn'], 'status'] <- 'lumped'
        }
        if(dim(c)[1] == 0 & !is.na(b[b$synonym %in% j,'nameUsage'])){
          d[d$syn_tsn %in% b[b$synonym %in% j,'syn_tsn'], 'status'] <- synCleanComp[synCleanComp$syn_tsn %in% b[b$synonym %in% j,'syn_tsn'],'nameUsage']
        }
      }
    }
    if(dim(b)[1] == 1){ 
      synCleanComp[synCleanComp$syn_tsn %in% a$syn_tsn,'nameUsage'] <- 'synonymised'
      d[d$syn_tsn %in% a$syn_tsn, 'status'] <- 'synonymised'
    }
  }
}


#### Clean issue with the tsn of valid species ####

valTSN <- data.frame()
for (i in unique(synCleanComp[synCleanComp$nameUsage %in% 'valid','scientificName'])){
  
  valTSN <- bind_rows(valTSN, get_tsn_(i, db = 'itis', accepted = TRUE))
}

valTSN <- valTSN %>%
  filter(scientificName %in% unique(synCleanComp[synCleanComp$nameUsage %in% 'valid','scientificName'])) %>%
  distinct()

for(i in unique(synCleanComp[synCleanComp$nameUsage %in% 'valid','scientificName'])) {
  synCleanComp[synCleanComp$scientificName %in% i, 'tsn'] <- valTSN[valTSN$scientificName %in% i,'tsn']
  synCleanComp[synCleanComp$synonym %in% i, 'syn_tsn'] <- valTSN[valTSN$scientificName %in% i,'tsn']
}

if(all(rmScinam$scientificName %in% synCleanComp$scientificName) &
   all(!rmScinam$synonym %in% synCleanComp$Synonym) & 
   all(rmScinam$tsn %in% synCleanComp$tsn) & 
   all(rmScinam$status %in% 'valid' & rmScinam$rank %in% 'Subspecies')){  
  
  toAdd <- rmScinam %>% select(syn_tsn, SynSpp, synSpp_tsn) 
  synCleanComp <- left_join(synCleanComp, toAdd, by =  'syn_tsn')
}

saveRDS(synCleanComp, paste0(folder_path,'nameMatching/outputs/synonymsClean.rds'))
