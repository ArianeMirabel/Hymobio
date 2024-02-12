#!/usr/bin/Rscript 
invisible(lapply(c("data.table", "ggplot2", "entropart"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

load("SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

files <- list.files("FinalMAJ_Naiades.Alric.Traits/")
lapply(paste0("FinalMAJ_Naiades.Alric.Traits/",files[setdiff(grep("_Inventories", files),grep(".csv", files))]), load, .GlobalEnv)


Inds <- c("Richness", "Shannon", "Simpson")

Indexes_measure <- function(X){
  ret <- tryCatch({lapply(0:2, function(Q){round(expq(bcTsallis(as.AbdVector(X),q=Q,Correction = "Marcon"),q=Q),2)})},
                  
 error = function(e) {
   ret <- tryCatch({lapply(0:2, function(Q){round(expq(bcTsallis(as.AbdVector(X),q=Q,Correction = "ChaoShen"),q=Q),2)})},
                               
                               error = function(e) {print(e)})
   })
  ret[[4]] <- names(ret[[1]])
  return(ret)
}

Inventory_taxo <- function(Inv){
  setnames(Inv, names(Inv), gsub("NomLatin_Taxo|NomLatinAppelTaxon", "NomTaxo", names(Inv)))
  
  Inv <- unique(Inv[,.(Group, CdStation, Year, Date_PrelBio, NomTaxo, Abundance)][
    , Abundance := sum(Abundance), by = c("CdStation", "NomTaxo", "Date_PrelBio")])
  Inv <- unique(Inv[, Abundance := ceiling(mean(Abundance)), by =c("CdStation", "NomTaxo", "Year")][
    , .(Group, CdStation, NomTaxo, Year, Date_PrelBio, Abundance)])
  
  Inv <- unique(Inv[
    , c("Nspecies", Inds, "Bias_Estimator") := c(uniqueN(NomTaxo), Indexes_measure(Abundance)), 
    by = c("CdStation","Year")][ ,.(Group, CdStation, Year, Date_PrelBio, Nspecies, Richness, Shannon, Simpson)])
  Inv[Richness == Inf, Richness := Nspecies]
  
  Inv <- merge(Inv, 
            unique(as.data.table(list_station_filter5_clean)[COMPARTIMENT_START == toupper(unique(Inv$Group)), 
                   .(ID_AMOBIO_START, ID_START)]), by.x = "CdStation", by.y = "ID_START", all.x = T)[
          ,c("Group", "CdStation", "ID_AMOBIO_START", "Year", "Nspecies", "Richness", "Shannon", "Simpson"), with = F]
  
  return(Inv)
}

Index_Taxo_Diatom <- Inventory_taxo(Inv = AllInv_Diatom)
save(Index_Taxo_Diatom, file = "TaxonomicBiodiv_Janv24_Diatom")

Index_Taxo_Macroinvertebrate <- Inventory_taxo(AllInv_Macroinvertebrate)
save(Index_Taxo_Macroinvertebrate, file = "TaxonomicBiodiv_Janv24_Macroinvertebrate")

Index_Taxo_Fish <- Inventory_taxo(AllInv_Fish)
save(Index_Taxo_Fish, file = "TaxonomicBiodiv_Janv24_Fish")
