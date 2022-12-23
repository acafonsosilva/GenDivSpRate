library(ape)
library(tidyverse)

folder_path <- '/data/biodiv/asilva/'
folder_path <- '/Users/acas/Dropbox/Post-docs/Morlon_Lab/manuscript_projects_info/mammals_genDiv_SpRate/manuscript_scripts_data/'
wd <- paste0(folder_path, 'mutationRate/clades/')
mutClades <- list.dirs(wd, recursive = FALSE, full.names = FALSE)

MutRateAll <- data.frame()
for ( i in (1:length(mutClades))){
 set <- list.dirs(paste0(wd, mutClades[i]), recursive = FALSE, full.names = FALSE)[-102]
  # if(i  == 6){ #Cricetidae
  #   set <- names(TreeSet)[-c(19)]
  # }
  # 
  # if(i == 9){ #Eulipotyphla
  #   set <- names(TreeSet)[-c(1,6,21)]
  # }
  # 
  # if(i == 10){ #GuineaPigRelated
  #   set <- names(TreeSet)[-c(5)]
  # }
  # 
  # if(all(i != c(6, 9, 10))){
  #   set <- names(TreeSet)
  # }
  
  for (j in set){
    tree <- read.tree(paste0(wd, mutClades[i],"/",j,"/", mutClades[i],".tree"))
    
    tipLengths <- data.frame(mutclade = mutClades[i], set = j, edge = tree$edge, time = tree$edge.length) %>%
      filter(edge.2 <= Ntip(tree)) %>%
      mutate(species = tree$tip.label,
             nodes_sister = ifelse(edge.1 %in% edge.1[duplicated(edge.1)], "sister","nonSister"))
    
    #extract trees from results
    read_in <- readLines(paste0(wd, mutClades[i],"/",j,"/", mutClades[i],"_res.txt"))
    muttree_pos <- grep(pattern = "tree length",x = read_in)
    muttree <- read.tree(text=read_in[muttree_pos+4])

    # extract tip branch lengths ~ expected number of substitutions per site at 3rd codon positions
    #which tips
    tips <- which(muttree$edge[,2]<=Ntip(muttree))
    #tip branch lengths
    MutRate <- data.frame(species = muttree$tip.label, expNsub = muttree$edge.length[tips], stringsAsFactors = FALSE) %>% 
      left_join(., tipLengths, by = 'species') %>%
      mutate(mutrate = expNsub/time)
    
    MutRateAll <- bind_rows(MutRateAll, MutRate)
  }
}

MutRateAll <- filter(MutRateAll, !expNsub %in% 4e-06) #these correspond to species that got a near zero expected number of substitutions
saveRDS(MutRateAll, paste0(folder_path, 'mutationRate/output/mutationRate_cytb3rdcodonPAML.rds'))

