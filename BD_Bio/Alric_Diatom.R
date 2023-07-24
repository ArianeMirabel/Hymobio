
Alr_Names_Diatom <- fread("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/3_DIATOMS/AMOBIO_DIAT_transcodage.csv")
Alr_Names_Diatom <- Alr_Names_Diatom[,lapply(.SD, function(x) return(iconv(x, "ISO-8859-15", "UTF-8")))]

Alr_Abd_Diatom <- fread("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/3_DIATOMS/AMOBIO_DIAT_RLQtable.csv",
                       encoding="Latin-1")

Alr_Traits_Diatom <- fread("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/3_DIATOMS/traitf_diat.csv",
                           encoding="Latin-1")

Correspondance_Diatom <- function(Alric_Traits, Alric_Abd, Naiades_Abd, Naiades_Station, Abs_Nai){
  
  return(lapply(unique(Abs_Nai$tax.GENUS), function(gen){
    
    print(gen)
  
if(nrow(Alric_traits[GENUS == gen]) == 1){
  return(c(Nai.GENUS = gen, Alr.SPECIES = Alric_traits[GENUS == gen,SCIENTIFIC_NAME]))} 
  
if(nrow(Alric_traits[GENUS == gen]) > 1){
      return(tryCatch({Find_close_match_DIA(Alr_rep = gen,
                                  Alric_traits = Alric_traits, Alric_Abd = Alric_Abd, 
                                  Naiades_Abd = Naiades_Abd, Naiades_Station = Naiades_Station, Abs_Nai= Abs_Nai)},
                      error = function(e) print(paste0("!!!!!!!!!!!!!!!!!!!!!!!!!!! Error with ",gen)) 
                      ))
  }
  
  if(nrow(Alric_traits[GENUS == gen]) == 0){
    
    Sim <- sort(unique(agrep(gen, Alric_traits$GENUS, max.distance = 0.2, ignore.case = FALSE, value = TRUE)))
    
    if(!identical(Sim, character(0))){
      
      print(paste0("Similar to ", gen , " from Alric : ", paste(paste0(Sim ," (",1:length(Sim),")"), collapse=", ", " ")))
      resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep: "))
      
      while(is.na(resp) | resp < 0 | resp > length(Sim)){ 
        resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep :"))}
      
      if(resp == 0){ return(character(0))} else {
        
          if(length(unique(Alric_traits[GENUS == Sim[resp], SCIENTIFIC_NAME])) == 1){
                    return(c(Nai.GENUS = gen, Alr.SPECIES = unique(Alric_traits[GENUS == Sim[resp]][,SCIENTIFIC_NAME])))}
        
              else {return(tryCatch({Find_close_match_DIA(Alr_rep = Sim[resp], same = FALSE, Nai_init = gen,
                                                Alric_traits = Alric_traits, Alric_Abd = Alric_Abd, 
                                                Naiades_Abd = Naiades_Abd, Naiades_Station = Naiades_Station, Abs_Nai= Abs_Nai)},
                error = function(e) print(paste0("!!!!!!!!!!!!!!!!!!!!!!!!!!! Error second round with ",gen)) 
              ))
                           }
        }
    }
      
    else {return(character(0))}
    } 
  }))}
