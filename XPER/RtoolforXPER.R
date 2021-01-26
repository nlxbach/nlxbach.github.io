##### READ SDD FILE ####"
# This function to read sdd format exported from XPER
# And return class xpersdd

read_sdd_xper <- function(sddfile){
  xml_data <- XML::xmlParse(file = sddfile) %>% XML::xmlToList()
  class(xml_data) <- c(class(xml_data), "xpersdd")
  return(xml_data)
}

##### TAXON ####
# function to get table of taxon ---
get_df_Taxon <- function(xml_data){
  if (!"xpersdd" %in% class(xml_data)) stop("xml_data must be created by function read_sdd_xper()")

  TaxonNames <- xml_data[["Dataset"]][["TaxonNames"]]
  TaxonNames_id = sapply(TaxonNames, FUN = function(x) x$.attrs)
  TaxonNames_label = sapply(TaxonNames, FUN = function(x) x$Representation$Label)
  
  df_TaxonNames <- data.frame(taxon_id = TaxonNames_id, 
                                taxon_label = TaxonNames_label)
    
  return(df_TaxonNames)

}


####### CHARACTER #########
# function to get state of one character and return to a vector
#Characters <- xml_data[["Dataset"]][["Characters"]]
#Character <- Characters[[4]]
get_state <- function(Character){
  if(any(names(Character)=="States")){
    list_states = Character$States
    state_id <- sapply(list_states, FUN = function(x) x$.attrs)
    state_label <- sapply(list_states, FUN = function(x) x$Representation$Label)
    state <- paste0("{",state_id,":",state_label,"}")
    res <- paste(state, collapse = " ; ")
  }else{
    res <- NA
  }
  return(res)
}
# function to get table of character
get_df_Characters <- function(xml_data){
  if (!"xpersdd" %in% class(xml_data)) stop("xml_data must be created by function read_sdd_xper()")
  
  Characters <- xml_data[["Dataset"]][["Characters"]]
  Characters_Label = sapply(Characters, FUN = function(x) x$Representation$Label)
  Characters_id = sapply(Characters, FUN = function(x) x$.attrs)
  Characters_type = names(Characters)
  Characters_states = sapply(Characters, get_state)
  
  df_Characters <- data.frame(ch_id = Characters_id,
                              ch_label = Characters_Label,
                              ch_type = Characters_type,
                              ch_state = Characters_states)
  
  return(df_Characters)
}

####### STATE #########
# function to get state of one character and return to a df
get_df_state <- function(Character){
  if(any(names(Character) == "States")){
    list_states = Character$States
    state_id <- sapply(list_states, FUN = function(x) x$.attrs)
    state_label <- sapply(list_states, FUN = function(x) x$Representation$Label)
    ch_id <- rep(Character$.attrs, length(state_id))
    res <- data.frame(state_id, state_label, ch_id )
  } else {
    res <- NULL
  }
  return(res)
}
# function to get table of states
get_df_States <- function(xml_data){
  if (!"xpersdd" %in% class(xml_data)) stop("xml_data must be created by function read_sdd_xper()")
  
  Characters <- xml_data[["Dataset"]][["Characters"]]
  df_States <- do.call("rbind", sapply(Characters, get_df_state)) 
  row.names(df_States) = c(1:nrow(df_States))
  return(df_States)
}



######## DESCRIPTION ########
# function to get description of one character of one taxon
#CodedDescriptions <- xml_data[["Dataset"]][["CodedDescriptions"]]
#Character = CodedDescriptions[[1]]$SummaryData[[2]]
get_descr_character <- function(Character){
  
  if(any(names(Character) == "State")){
    index <- which(names(Character) == "State")
    states_checked <- sapply(index, FUN = function(x) Character[x])
    res <- paste(unlist(states_checked), collapse = " & ")
    
    
  } else if(any(names(Character) == "Measure")){
    index <- which(names(Character) == "Measure")
    states_checked <- sapply(index, FUN = function(x) Character[x])
    type <- sapply(states_checked, FUN = function(x) x[1])
    value <- sapply(states_checked, FUN = function(x) x[2])
    min <- value[type == "Min"]
    max <-value[type == "Max"]
    res <- paste0("[",min,";",max,"]")
  } else {
    res <- NA
  }
  
  return(res)
}

# function to get description of all character of one taxon
#Descr_taxon <- CodedDescriptions[[1]]
get_descr_taxon <- function(Descr_taxon){
  SummaryData <- Descr_taxon$SummaryData
  
  ch_id = sapply(SummaryData, FUN = function(x) x$.attrs)
  descr = sapply(SummaryData, get_descr_character)
  taxon_id = rep(Descr_taxon$Scope$TaxonName, length(ch_id))
  
  df_descr <- data.frame(taxon_id, ch_id, descr)
  return(df_descr)
}




# function to get matrix of character and taxon
get_df_Descriptions <- function(xml_data){
  if (!"xpersdd" %in% class(xml_data)) stop("xml_data must be created by function read_sdd_xper()")
  
  CodedDescriptions <- xml_data[["Dataset"]][["CodedDescriptions"]]
  df_Descriptions <- do.call("rbind", lapply(CodedDescriptions, get_descr_taxon))
  df_Descriptions <- tidyr::pivot_wider(df_Descriptions, names_from = "ch_id", values_from = "descr")

  
  df_rule <- get_df_rule(xml_data)
  
  for(i in 1:nrow(df_rule)){
    ch_fa <- df_rule$ch_id_fa[i]
    ch_child <- df_rule$ch_id_child[i]
    state <- df_rule$state_inapp[i]
    index <- sapply(df_Descriptions[ch_fa], FUN = function(x) str_detect(x, pattern = state))
    df_Descriptions[ch_child][index] <- "nonapp"
  }
  
  return(df_Descriptions)
}




####### RULE #########
# Function to get rule of one character
#CharNode <- Nodes[[2]]
get_rule_char <- function(CharNode){
  if(any(names(CharNode) == "DependencyRules")){
    state_inapp <- sapply(CharNode$DependencyRules$InapplicableIf, FUN = function(x) x)
    ch_id_child <- rep(CharNode$Character, length(state_inapp))
    res <- data.frame(ch_id_child, state_inapp)
  }else{
    res <- NULL
  }
  return(res)
}

# Function to get table of all child character
get_df_rule <- function(xml_data){
  if (!"xpersdd" %in% class(xml_data)) stop("xml_data must be created by function read_sdd_xper()")
  
  df_States <- get_df_States(xml_data) %>% select(state_id,"ch_id_fa"=ch_id)
  Nodes <- xml_data$Dataset$CharacterTrees$CharacterTree$Nodes
  
  res <- do.call("rbind", sapply(Nodes, get_rule_char)) %>% 
    left_join(df_States, by = c("state_inapp"="state_id"))
  return(res)
}

############### MERGEMOD ############### 
#' This function is rewrite from Algorithm MERGEMOD copyright by:
#' Lebbe, J., and R. Vignes. 2003. 
#' Utilitaires XPER: bibliothèque de programmes C et leur documentation pour l’analyse des bases de connaissances taxonomiques.


#' Funtion to get mergemod for one character
#' trait = c("A","B&C","D&E","B&D&F","A&B&E","C&E&F")
#' list_state <- c("A", "B", "C", "D", "E", "F")
#' b.mergemod.one.trait(trait,list_state)

mergemod_one_trait <- function(trait, list_state, sep = "&"){
  #get number of taxon (n) and state (m)
  n <- length(trait) #number of taxon
  m <- length(list_state) #number of state
  
  # create vector result
  list_pair_merge <- c()
  
  #Run loop 
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      #test if i and j can be merged
      pair_state <- c(list_state[i], list_state[j])
      res = TRUE
      # for each couple of taxa a and b
      for(a in 1:(n-1)){
        for (b in (a+1):n){
          #creat a vector state of each taxa a and b
          sp_a <- unlist(strsplit(trait[a],sep))
          sp_b <- unlist(strsplit(trait[b],sep))
          # check length of vector intersect
          len <- length(intersect(sp_a,sp_b))
          
          # If len = 0 => the two taxa are distinguished
          if (!is.na(sp_a) && !is.na(sp_b) && len == 0){
            if(length(intersect(sp_a, pair_state)) > 0 & length(intersect(sp_b, pair_state)) > 0){
              res = FALSE
              break
            }
          }
        }
        if(res == FALSE){break}
      }
      
      if(res == TRUE){
        list_pair_merge <- append(list_pair_merge, paste(pair_state, collapse = "-"))
      }
      
    }
  }
  return(list_pair_merge)
}


#' Function to get results from mergemod analysis for 
#' all categories characters in database

mergemod <- function(xml_data){
  if (!"xpersdd" %in% class(xml_data)) stop("xml_data must be created by function read_sdd_xper()")
  
  df_Descriptions <- get_df_Descriptions(xml_data) %>% as.data.frame()
  df_Characters <- get_df_Characters(xml_data)
  df_States <- get_df_States(xml_data)
  
  ch_id_list <- c()
  state_merge <- c()
  
  for(i in 1:nrow(df_Characters)){
    ch_id <- df_Characters$ch_id[i]
    trait = df_Descriptions[,ch_id]
    ch_id_list[i] <- ch_id
    
    if(df_Characters$ch_type[i] == "CategoricalCharacter"){
      list_state <- df_States$state_id[df_States$ch_id == ch_id]
      state_merge[i] <-  paste0(mergemod_one_trait(trait, list_state, sep = " & "), 
                                collapse = " ; ")
    }else{
      state_merge[i] <-"Non Application"
    }
  }
  
  res_df_Mergemod <- data.frame(ch_id = ch_id_list, state_merge) 
  df_Characters <- left_join(df_Characters, res_df_Mergemod, by="ch_id")
  
  return(df_Characters)
  
}
