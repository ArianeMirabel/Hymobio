
Alr_traits_Fish <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/2_FISHES/AMOBIO_FISH_Traits.csv",
                      sep=";"))[,.(Nom.vernaculaire, English.name, sciname, Code_TAXON, 
                                   tax.CLASS, tax.ORDER, tax.FAMILY, tax.GENUS, tax.SPECIES)]
Alr_traits_Fish <- Alr_traits_Fish[,lapply(.SD, function(x) return(iconv(x, "ISO-8859-15", "UTF-8")))]
Alr_traits_Fish[sciname == "Gymnocephalus cernua", sciname := "Gymnocephalus cernuus"]

Alr_Abd_Fish <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/2_FISHES/AMOBIO_FISH_Abundance.csv",
                                  sep=";"))
Alr_Abd_Fish <- Alr_Abd_Fish[!is.na(X_L93)]
Alr_Abd_Fish <- Alr_Abd_Fish[!is.na(CODE_STATION),]

Alr_Abd_Fish[, year := format(as.Date(Alr_Abd_Fish[,DATE], format =  "%d/%m/%Y"), "%Y")]

Correspondance_Fish <- function(Alric_traits , Alric_Abd, Naiades_Abd, Naiades_Station, Abs_Nai){
  
  return(lapply(unique(Abs_Nai$tax.GENUS), function(gen){
  
if(nrow(Alric_traits[tax.GENUS == gen]) == 1){return(c(Nai.GENUS = gen,
         Alr.SPECIES = Alric_traits[tax.GENUS == gen][,paste(tax.GENUS,tax.SPECIES)]))} 
  
if(nrow(Alric_traits[tax.GENUS == gen]) > 1){
    return(Find_close_match_F(Alr_rep = gen,
                              Alric_traits = Alric_traits, Alric_Abd = Alric_Abd, 
                              Naiades_Abd = Naiades_Abd, Naiades_Station = Naiades_Station, Abs_Nai= Abs_Nai))
  }
  
  if(nrow(Alric_traits[tax.GENUS == gen]) == 0){
    
    Sim <- unique(agrep(gen, Alric_traits$tax.GENUS, max.distance = 0.4, ignore.case = FALSE, value = TRUE))
    
    if(!identical(Sim, character(0))){
      
      print(paste0("Similar to ", gen , " from Alric : ", paste(paste0(Sim ," (",1:length(Sim),")"), collapse=", ", " ")))
      resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep: "))
      
      while(is.na(resp) | resp<0 | resp>length(Sim)){ resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep :"))}
      
      if(resp == 0){ return(character(0))} else {
        
          if(length(Sim) == 1){return(c(Nai.GENUS = gen, Alr.SPECIES = Alric_traits[tax.GENUS == Sim[resp]][,paste(tax.GENUS,tax.SPECIES)]))}
        
              else {return(Find_close_match_F(Alr_rep = Sim[resp], same = FALSE, Nai_init = gen,
                                              Alric_traits = Alric_traits, Alric_Abd = Alric_Abd, 
                                              Naiades_Abd = Naiades_Abd, Naiades_Station = Naiades_Station, Abs_Nai= Abs_Nai))}
        }
    }
    
    SimVern <- unique(agrep(gen, Alric_traits$Nom.vernaculaire, max.distance = 0.5, ignore.case = FALSE, value = TRUE),
                      agrep(gen, Alric_traits$English.name, max.distance = 0.5, ignore.case = FALSE, value = TRUE))
    SimVern <- unique(Alric_traits[Nom.vernaculaire %in% SimVern, tax.GENUS], Alric_traits[English.name %in% SimVern, tax.GENUS])
    
    if(!identical(SimVern, character(0))){
      
      print(paste0("Vernacular name similar to ", gen , " from Alric : ", paste(paste0(SimVern ," (",1:length(SimVern),")"), collapse=", ", " ")))
      resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep: "))
      
      while(is.na(resp) | resp<0 | resp>length(SimVern)){ resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep :"))}
      
      if(resp == 0){ return(character(0))} else {
        
        if(length(SimVern) == 1){return(c(Nai.GENUS = gen, Alr.SPECIES = Alric_traits[tax.GENUS == SimVern[resp]][,paste(tax.GENUS,tax.SPECIES)]))}
        
        else {return(Find_close_match_F(Alr_rep = SimVern[resp], same = FALSE, Nai_init = gen,
                                        Alric_traits = Alric_traits, Alric_Abd = Alric_Abd, 
                                        Naiades_Abd = Naiades_Abd, Naiades_Station = Naiades_Station, Abs_Nai= Abs_Nai))}
      }
    }
      
    else {return(character(0))}
    } 
}))}
