
Alr_Names_Macroinv <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/1_MACROINVERTEBRATES/2_BASE_AMOBIO/AMOBIO_INV_transcodage.csv",
                                       sep=";"))

Alr_Names_Macroinv <- Alr_Names_Macroinv[,lapply(.SD, function(x) return(iconv(x, "ISO-8859-15", "UTF-8")))]
Alr_Names_Macroinv <- Alr_Names_Macroinv[ ,tax.GENUS := sub(" .*", "", TAXON)]

Alr_Traits_Macroinv <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/1_MACROINVERTEBRATES/2_BASE_AMOBIO/AMOBIO_INV_traits.csv",
                      sep=";"))
Alr_Abd_Macroinv <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/1_MACROINVERTEBRATES/2_BASE_AMOBIO/AMOBIO_INV_RLQtable.csv",
                                  sep=";"))
Alr_Abd_Macroinv <- Alr_Abd_Macroinv[!is.na(x)]
Alr_Abd_Macroinv <- Alr_Abd_Macroinv[,1:(grep("m",colnames(Alr_Abd_Macroinv))[1]-1), with = FALSE]

setnafill(Alr_Abd_Macroinv[, .SD, .SDcols = is.integer], fill = 0)
setnafill(Alr_Abd_Macroinv[, .SD, .SDcols = is.numeric], fill = 0) 

Alr_Names_Macroinv[grep("Helodes", TAXON), TAXON := gsub("Helodes", "Elodes", TAXON)]

Correspondance_Macroinv <- function(Alric_traits , Alric_Abd, Naiades_Abd, Naiades_Station, Abs_Nai){
  
  return(lapply(unique(Abs_Nai$tax.GENUS), function(gen){
if(nrow(Alric_traits[tax.GENUS == gen]) == 1){return(c(Nai.GENUS = gen,
                                                   Alr.SPECIES = Alric_traits[tax.GENUS == gen,TAXON]))} 
  
if(nrow(Alric_traits[tax.GENUS == gen]) > 1){
    return(tryCatch({Find_close_match_MI(Alr_rep = gen,
                                         Alric_traits = Alric_traits, Alric_Abd = Alric_Abd, 
                                         Naiades_Abd = Naiades_Abd, Naiades_Station = Naiades_Station, Abs_Nai= Abs_Nai)},
                    
                    error = function(e) print(gen) 
    ))
  }
  
  if(nrow(Alric_traits[tax.GENUS == gen]) == 0){
    
    Sim <- sort(unique(agrep(gen, Alric_traits$tax.GENUS, max.distance = 0.3, ignore.case = FALSE, value = TRUE)))
    
    if(!identical(Sim, character(0))){
      
      print(paste0("Similar to ", gen , " from Alric : ", paste(paste0(Sim ," (",1:length(Sim),")"), collapse=", ", " ")))
      resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep: "))
      
      while(is.na(resp) | resp<0 | resp>length(Sim)){ resp <- as.integer(readline(prompt = "Enter 0 or n° of genus to keep :"))}
      
      if(resp == 0){ return(character(0))} else {
        
          if(nrow(Alric_traits[tax.GENUS == Sim[resp]]) == 1){return(c(Nai.GENUS = gen, Alr.SPECIES = Alric_traits[tax.GENUS == Sim[resp]][,TAXON]))}
        
              else {return(Find_close_match_MI(Alr_rep = Sim[resp], same = FALSE, Nai_init = gen,
                                               Alric_traits = Alric_traits, Alric_Abd = Alric_Abd, 
                                               Naiades_Abd = Naiades_Abd, Naiades_Station = Naiades_Station, Abs_Nai= Abs_Nai))}
        }
    }
      
    else {return(character(0))}
    } 
}))
}






