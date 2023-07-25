invisible(lapply(c("data.table", "httr", "jsonlite","stringr", "sp", "parallel", "foreach", "doParallel", "snow", "ggplot2", "dplyr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Hydro")

# Get Stations
#####
Fields <- paste0("&fields=code_site,code_station,libelle_station,longitude_station,latitude_station,libelle_projection,",
                "libelle_departement,qualification_donnees_station,en_service,date_ouverture_station,date_fermeture_station,",
                "code_finalite_station,date_maj_station,code_projection")
Ref_url <- fromJSON(paste0("https://hubeau.eaufrance.fr/api/v1/hydrometrie/referentiel/stations?",Fields))
Hydro_Station <- data.table()
while(nrow(Hydro_Station) != Ref_url$count){ Hydro_Station <- rbind(Hydro_Station,Ref_url$data)
  if(!is.null(Ref_url$`next`)) Ref_url <- fromJSON(Ref_url$`next`)}
rm(Ref_url)
Hydro_Station <- Hydro_Station[!libelle_departement %in% c("GUADELOUPE", "MARTINIQUE", "MAYOTTE","LA REUNION", "GUYANE")]
Hydro_Station <- Hydro_Station[longitude_station <= 40]
#####

extract_function <- function(adress){
 ret <- tryCatch(
   { fromJSON(adress) 
     }, error = function(e) { e
      }, warning= function(w) w )

 if(!is.null(ret$message) && !any(grep("500 Internal Server Error", ret$message))){
   ret <-list(count = "Other error")
   } else { 
     while(!is.null(ret$message) && any(grep("500 Internal Server Error", ret$message))){
   print("Server error")
   Sys.sleep(time = 2)
   ret <- fromJSON(adress)} 
   }
 return(ret)
}


# Debit Mensuel
#####
Hydro_mensuel <- data.table()
No.data <- character()
for(station in unique(Hydro_Station$code_station)){
  ret <- data.table()
  
  extract <- extract_function(paste0("https://hubeau.eaufrance.fr/api/v1/hydrometrie/obs_elab?code_entite=", station,
  "&grandeur_hydro_elab=QmJ", sep = ""))

  if(extract$count != 0){
    while(nrow(ret) != extract$count){
      ret <- rbind(ret,extract$data)
      if(!is.null(extract$`next`)) {extract <- extract_function(extract$`next`)}
    }
    
    Hydro_mensuel <- rbind(Hydro_mensuel,ret[,.(code_site, code_station, longitude,latitude, 
                                                date_obs_elab, grandeur_hydro_elab, resultat_obs_elab)])
  }
  if(extract$count == 0) {No.data <- c(No.data, station)}
  if(extract$count == "Other error") {No.data <- c(No.data, paste("Other pb with", station))}
}

Hydro_mensuel[, year := format(as.Date(date_obs_elab), "%Y")][, month := format(as.Date(date_obs_elab), "%m")]
save(Hydro_mensuel, file = "Extract_HydrologieMensuelle")
#####

# Debit Journalier
#####
Station.run <- split(Hydro_Station$code_station, ceiling(seq_along(Hydro_Station$code_station) / 500))

rm(list = setdiff(ls(),c("Station.run", "extract_function")))

clust <- makeCluster(detectCores()-2)
registerDoParallel(clust)

Extract_slices <- mclapply(1:length(Station.run), function(i){
  Station.run.slice <- Station.run[[i]]
  writeLines(paste("Chunck ", i, "/", length(Station.run)))
  
Hydro_journalier <- data.table()
No.data <- character()

pb <- txtProgressBar(min = 0, max = length(Station.run.slice), style = 3, char = "=")

  for(station in Station.run.slice){
  
  ret <- data.table()
    
  extract <- extract_function(paste0("https://hubeau.eaufrance.fr/api/v1/hydrometrie/obs_elab?code_entite=", station,
                                      "&grandeur_hydro_elab=QmJ", sep = ""))
  
  if(extract$count == "Other error") {
    print(paste("Error station ", station))
    No.data <- c(No.data, paste("Other error with", station))
  } else {
    if(extract$count != 0){
    while(nrow(ret) != extract$count){
      ret <- rbind(ret,extract$data)
      if(!is.null(extract$`next`)) {extract <- extract_function(extract$`next`)}
    }
    
    Hydro_journalier <- rbind(Hydro_journalier,ret[,.(code_site, code_station, longitude,latitude, 
                                                date_obs_elab, grandeur_hydro_elab, resultat_obs_elab)])
    rm(ret)
  } 
  if(extract$count == 0) {No.data <- c(No.data, station)}   
  }
  
    setTxtProgressBar(pb, which(Station.run.slice == station))
}
close(pb)

save(Hydro_journalier, file = paste0("Extract_HydrologieJournaliere_", i))
return(No.data)
})
  
stopCluster(clust)

Hydro_journaliere <- do.call(rbind,lapply(grep("Extract_HydrologieJournaliere_",list.files(), value = TRUE), function(f) {
  load(f); return(Hydro_journalier)
}))

Hydro_journaliere[, Date := as.Date(date_obs_elab)][, year := format(as.Date(date_obs_elab, format =  "%Y-%m-%d"), "%Y")][
  , month := format(as.Date(date_obs_elab, format =  "%Y-%m-%d"), "%m")][, day := format(as.Date(date_obs_elab, format =  "%Y-%m-%d"), "%d")]

save(Hydro_journaliere, file = "Hydro_journaliere")
#####

# Complétude des données
#####
load("Hydro_journaliere")

nrow(Hydro_journaliere[resultat_obs_elab < 0])/nrow(Hydro_journaliere)*100
uniqueN(Hydro_journaliere[resultat_obs_elab < 0, code_site])/uniqueN(Hydro_journaliere[, code_site])*100


DayLim <- 1
clust <- makeCluster(detectCores()-2)
registerDoParallel(clust)

Run_periods <- do.call(rbind,mclapply(unique(Hydro_journaliere$code_station), function(st){
  x <- Hydro_journaliere[code_station == st, .(Date)]
  x <- x[order(Date)]
  run <- 1
  x[, Run_period := run][, Date_check := seq.Date(from = first(Date), length.out = length(Date), by = "days")-Date]
  x[Date_check > -DayLim, Date_check := 0]
  
  while(any(x[, Date_check] != 0)){
    x[Date_check != 0, c("Run_period", "Date_check") := list(run + 1, seq.Date(from = first(Date), length.out = length(Date), by = "days")-Date)]
    run <- run + 1
    x[Date_check > -DayLim, Date_check := 0]
  }
  return(x[, code_station := st][,.(code_station, Date, Run_period)])}))

stopCluster(clust)
gc()

Hydro_journaliere[, Run_period := 1]
Hydro_journaliere[Run_periods, on = c("code_station", "Date"), Run_period := i.Run_period]

Hydro_journaliere[, Length_period := .N, by = c("code_station", "Run_period")]
Hydro_journaliere[, Cover := "All_periods"]
Hydro_journaliere[Length_period < 1826 & Length_period >= 365, Cover := "One_year"]
Hydro_journaliere[Length_period < 365 & Length_period >= 182, Cover := "Six_months"]
Hydro_journaliere[Length_period < 182 & Length_period >= 91, Cover := "Three_months"]
Hydro_journaliere[Length_period < 91 & Length_period >= 30, Cover := "One_month"]
Hydro_journaliere[Length_period < 30 , Cover := "Too_short"]

PiePeriod <- unique(Hydro_journaliere[, Eff_periods := .N, by = Cover][,.(Eff_periods, Cover)])[, 
                    Cover := factor(Cover, levels = c("All_periods", "One_year","Six_months", "Three_months", "One_month", "Too_short"))]

ggplot(PiePeriod,aes(x = "", y = Eff_periods, fill = Cover)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Accent") +
  theme_void() + guides(fill = guide_legend(title="")) + ggtitle(paste("Authorized gap", DayLim, "days"))  

Hydro_journaliere[, Eff_periods := NULL]; Hydro_journaliere <- Hydro_journaliere[resultat_obs_elab > 0]
save(Hydro_journaliere, file = "Hydro_journaliere_raw")

#####

# Link to Bio Stations

Hydro_station_DataPresence <- Hydro_Station[, Presence_data := "Yes"][!code_station %in% unique(Hydro_journaliere$code_station), Presence_data := "No"]
save(Hydro_station_DataPresence, file = "HydroStations_DataPresence")

load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")
Correspondance_station <- as.data.table(list_station_filter5_clean)[COMPARTIMENT_START %in% c("DIATOM", "MACROINVERTEBRATE", "FISH") & 
                          COMPARTIMENT_TARGET == "HYDRO"][,.(COMPARTIMENT_START, ID_AMOBIO_START, COMPARTIMENT_TARGET, ID_TARGET, to_keep)][
                          , Nstation := .N, by = "ID_AMOBIO_START"]
load("Hydro_journaliere")
Correspondance_station <- Correspondance_station[ID_TARGET %in% unique(Hydro_journaliere$code_station)]
rm(list = c("list_station_filter5_clean", "Hydro_journaliere"))

load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/Fish_Inventories")
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Inventories")
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/Diatom_Inventories")

BioInv <- rbind(unique(AllInv_Fish[, ID_AMOBIO_START := paste0("FISH_", CdStation)][,.(CdStation, ID_AMOBIO_START, Date_PrelBio)]),
                unique(AllInv_Macroinvertebrate[, ID_AMOBIO_START := paste0("MACROINVERTEBRATE_", CdStation)][,.(CdStation, ID_AMOBIO_START, Date_PrelBio)]),
                unique(AllInv_Diatom[, ID_AMOBIO_START := paste0("DIATOM_", CdStation)][,.(CdStation, ID_AMOBIO_START, Date_PrelBio)]))
rm(list=c("AllInv_Fish", "AllInv_Macroinvertebrate", "AllInv_Diatom"))

Correspondance_station <- Correspondance_station[ID_AMOBIO_START %in% BioInv$ID_AMOBIO_START]
Correspondance_station <- left_join(Correspondance_station, BioInv, by = "ID_AMOBIO_START", relationship = "many-to-many")

load("CheckSingleSpeciesStations")
SingleSpeciesStations <- paste0("FISH_", unique(SingleSpeciesStations))

Correspondance_station <- Correspondance_station[!ID_AMOBIO_START %in% SingleSpeciesStations]
save(Correspondance_station, file = "Hydro_BioInventories_correspondance")








