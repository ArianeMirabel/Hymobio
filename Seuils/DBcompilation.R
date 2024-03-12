invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel", "dplyr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

load("AMOBIO_WP3_2_FISH_INV_DIA_ENTROPIE_20231123.Rdata")
RFD_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN)[, H_Ncrue_S1_ALL := NULL]
RFD_all[,grep("B_INV|B_FISH|B_DIA", colnames(RFD_all), value = T) := NULL]
RFD_all[, Year := as.character(format(DATE_OPERATION, format = "%Y"))]

directory <- "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/"
lapply(paste0(directory,
    grep("DiversityFunctional_2401|DiversityTaxonomic_2401",list.files(directory), value =T)),
    load, .GlobalEnv)

for(compartment in c("Fish", "Macroinvertebrate", "Diatom")) {
  RFD_all <- left_join(RFD_all,  get(paste0("DiversityTaxonomic_", compartment))[, ID_AMOBIO_START := NULL][
    , Year := as.character(Year)], by = c("COMPARTMENT" = "COMPARTMENT", "ID" = "CdStation", "Year" = "Year"))
  RFD_all <- left_join(RFD_all,  get(paste0("DiversityFunctional_", compartment))[, ID_AMOBIO_START := NULL][
    , Year := as.character(Year)],by = c("COMPARTMENT" = "COMPARTMENT", "ID" = "CdStation", "Year" = "Year"))
}


load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Hydro/HydroIndex_36125all_AM_20232911")
RFD_all <- left_join(RFD_all, Hydrolaps[, ID := sub('.*_', '', ID_AMOBIO_START)][, COMPARTMENT := sub('_.*', '', ID_AMOBIO_START)][
  , c("ID", "ID_AMOBIO_START", "COMPARTMENT", "Samp_date", grep("Ncrues|Netiage", colnames(Hydrolaps), value = T)), with = F],
  by = c("ID" = "ID", "DATE_OPERATION" = "Samp_date", "COMPARTMENT" = "COMPARTMENT"))
colnames(RFD_all) <- sub(" ", "_", colnames(RFD_all))
colnames(RFD_all) <- sub("[^ -~]", "", colnames(RFD_all))
colnames(RFD_all) <- sub("(", "", colnames(RFD_all), fixed = T)
colnames(RFD_all) <- sub(")", "", colnames(RFD_all), fixed = T)

save(RFD_all, file = "AMOBIO_WP3_2_FISH_INV_DIA_ENTROPIE_20240212.RData")

Catalog <- fread("../MetricsCatalogue.csv")[which(TokeepThreshold)]

RFD_all <- RFD_all[, c("ID","Year",grep(paste0("B_FISH|B_INV|B_DIA|",paste(Catalog$NameOrigin, collapse = "|")), colnames(RFD_all), value = T)), with = F]
RFD_all <- RFD_all[, lapply(.SD, function(X) {X[X == "-Inf"] <- NA; return (X)})]
save(RFD_all, file = "../HYMOBIO_FULLDATA_202403.RData")
rm(list= c("MATRIX_AMOBIO_WP3_CLEAN", "Hydrolaps"))


load("AMOBIO_WP3_2_FISH_INV_DIA_ENTROPIE_20240212.RData")
RFD_all[COMPARTMENT == "FISH" & !is.na(B_FISH_Richness),grep("B_FISH|ID|Year", colnames(RFD_all), value = T), with = F ]

Count <- do.call(rbind,lapply(c("FISH", "INV", "DIA"), function(comp){
return(data.frame(Nsites = RFD_all[grep(comp,COMPARTMENT)][!is.na(get(paste0("B_",comp,"_Richness"))), uniqueN(ID)],
           Nopes = RFD_all[grep(comp,COMPARTMENT)][!is.na(get(paste0("B_",comp,"_Richness"))), uniqueN(paste0(ID, Year))],
           Dstart = RFD_all[grep(comp,COMPARTMENT)][!is.na(get(paste0("B_",comp,"_Richness"))), min(Year)],
           Dend = RFD_all[grep(comp,COMPARTMENT)][!is.na(get(paste0("B_",comp,"_Richness"))), max(Year)],
           Compartment = comp))
  }))


source("C:/Users/armirabel/Documents/INRAE/Hymobio/Hymobio_GitHub/BD_Bio/Indices/FunctionalDiversity_Functions.R")

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

directory <- paste0(getwd(), "/FinalMAJ_Naiades.Alric.Traits/")

Count <- do.call(rbind,lapply(c("Fish", "Macroinvertebrate", "Diatom"), function(comp){
  
  invisible(load_CompartmentFiles(Directory = directory, Compartment = comp))
  
  return(data.frame(Ntaxon = get(paste0("Traits_", comp))[, uniqueN(CdAppelTaxon)],
                    Ntraits = uniqueN(get(paste0("Code_", comp))[,Trait]),
                    Compartment = comp))
}))



setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

directory <- "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/"
lapply(paste0(directory,
              grep("DiversityFunctional_2401|DiversityTaxonomic_2401",list.files(directory), value =T)),
       load, .GlobalEnv)

Count <- do.call(rbind,lapply(c("Fish", "Macroinvertebrate", "Diatom"), function(comp){
  code <- toupper(comp)
  if(comp == "Macroinvertebrate"){code <- "INV"}
  if(comp == "Diatom"){code <- "DIA"}
  
  return(data.frame(Ntaxo = length(grep(paste0("B_", code, "_"),colnames(get(paste0("DiversityTaxonomic_", comp))))),
                    Nfunct = length(grep(paste0("B_", code, "_"),colnames(get(paste0("DiversityFunctional_", comp))))),
                    Compartment = comp))
}))


