invisible(lapply(c("data.table", "ggplot2", "entropart"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

dir <- "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/"
files <- list.files(dir)
lapply(paste0(dir,files[setdiff(grep("_Inventories", files),grep(".csv", files))]), load, .GlobalEnv)
load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

#Old, lapply version
#####
DiatomTaxo <- lapply(unique(AllInv_Diatom$CdStation), function(stat){
  Stat <- AllInv_Diatom[CdStation == stat]
  
  Ret <- do.call(rbind,lapply(unique(Stat$Year), function(yr){
    Stat_yr <- Stat[Year == yr]
    ret <- tryCatch({unlist(lapply(0:2, function(Q){round(expq(bcTsallis(as.AbdVector(Stat_yr$Abundance),q=Q,Correction = "UnveilC"),q=Q),2)}))},
                    
                    error = function(e) tryCatch({unlist(lapply(0:2, function(Q){round(expq(bcTsallis(as.AbdVector(Stat_yr$Abundance),q=Q,Correction = "ChaoShen"),q=Q),2)}))},
                                                 error = function(e) {print(paste("Other problem", stat))}))
    names(ret) <- c("Richness", "Shannon", "Simpson")
    return(ret)
  }))
  Ret$Cdstation <- stat
  return(Ret)})
save(DiatomTaxo, file = "TaxonomicBiodiv_Diatom")
#####

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
    , .(Group, CdStation, NomTaxo, Year, Abundance)])
  
  Inv <- unique(Inv[
    , c("Nspecies", Inds, "Bias_Estimator") := c(uniqueN(NomTaxo), Indexes_measure(Abundance)), 
    by = c("CdStation","Year")][ ,.(Group, CdStation, Year, Nspecies, Richness, Shannon, Simpson)])
  Inv[Richness == Inf, Richness := Nspecies]
  
  Inv <- merge(Inv, 
            unique(as.data.table(list_station_filter5_clean)[COMPARTIMENT_START == toupper(unique(Inv$Group)), 
                   .(ID_AMOBIO_START, ID_START)]), by.x = "CdStation", by.y = "ID_START", all.x = T)[
          ,c("Group", "CdStation", "ID_AMOBIO_START", "Year", "Nspecies", "Richness", "Shannon", "Simpson"), with = F]
  
  return(Inv)
}

Index_Taxo <- lapply(c("Fish", "Macroinvertebrate", "Diatom"), function(compartment){
  assign(paste0("DiversityTaxo_", compartment), Inventory_taxo(Inv = get(paste0("AllInv_", compartment))))
  save(list = paste0("DiversityTaxo_", compartment), 
       file = paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/DiversityTaxonomic_", compartment))
  })

Index_Taxo_Diatom <- Inventory_taxo(Inv = AllInv_Diatom)
save(Index_Taxo_Diatom, file = "TaxonomicBiodiv_Diatom")

Index_Taxo_Macroinvertebrate <- Inventory_taxo(AllInv_Macroinvertebrate)
save(Index_Taxo_Macroinvertebrate, file = "TaxonomicBiodiv_Macroinvertebrate")

Index_Taxo_Fish <- Inventory_taxo(AllInv_Fish)
save(Index_Taxo_Fish, file = "TaxonomicBiodiv_Fish")

##
# Diatoms supplementary taxo
dir <- "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/"
load(paste0(dir,"Diatom_Inventories"))
load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

Inv <- unique(AllInv_Diatom[,.(Group, CdStation, Year, Date_PrelBio, NomLatin_Simple, Genus, Abundance)][
  , Abundance := sum(Abundance), by = c("CdStation", "NomLatin_Simple", "Genus", "Date_PrelBio")])
Inv <- unique(Inv[, Abundance := ceiling(mean(Abundance)), by =c("CdStation", "NomLatin_Simple", "Year")][
  , .(Group, CdStation, NomLatin_Simple, Genus, Date_PrelBio, Year, Abundance)])

Inv <- unique(Inv[, Ngenus := sum(Abundance), by = c("CdStation", "Year", "Genus")][
    , Abdgenus := Ngenus/sum(Abundance), by = c("CdStation", "Year")][
    , Rgenus := as.numeric(uniqueN(NomLatin_Simple)), by = c("CdStation", "Year", "Genus")][
    , Rgenus := as.numeric(Rgenus/uniqueN(NomLatin_Simple)), by = c("CdStation", "Year")][
    , .(Group, CdStation, Date_PrelBio, Year, Genus, Abdgenus, Rgenus)])

Inv <- dcast(Inv, Group + CdStation + Date_PrelBio + Year ~ Genus, value.var = c("Abdgenus", "Rgenus"))
Inv[is.na(Inv),] <-0

setnames(Inv, grep("Abdgenus|Rgenus", colnames(Inv), value = T), 
         paste0("B_DIA_", grep("Abdgenus|Rgenus", colnames(Inv), value = T)))
Inv[, grep("Abnormal", colnames(Inv), value = T) := NULL]
cols <- setdiff(colnames(Inv), c("Group", "CdStation", "Date_PrelBio", "Year"))
Inv <- Inv[,  lapply(.SD, mean), .SDcols = cols, by = c("Group", "CdStation","Year")]

SuppTaxo_Diatom <- merge(Inv, 
             unique(as.data.table(list_station_filter5_clean)[COMPARTIMENT_START == toupper(unique(Inv$Group)), 
                    .(ID_AMOBIO_START, ID_START)]), by.x = "CdStation", by.y = "ID_START", all.x = T)[
                     ,c("Group", "CdStation", "ID_AMOBIO_START", "Year", grep("B_DIA", colnames(Inv), value = T)), with = F]

save(SuppTaxo_Diatom, file = 
       "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/TaxonomicBiodiv_Diatom_SuppTaxo")


##
# Data plots
Toplot <- rbind(Index_Taxo_Fish, Index_Taxo_Macroinvertebrate, Index_Taxo_Diatom)[ , as.list(unlist(lapply(.SD, 
      function(x) list(mean = mean(x, na.rm = T), 
                       up = mean(x, na.rm = T)+sd(x, na.rm = T)/2, 
                       low = mean(x, na.rm = T)-sd(x, na.rm = T)/2)))),
      .SDcols = Inds, by = c("Group", "Year")]
Toplot <- melt(Toplot, measure = list(paste(Inds, "mean", sep = "."), paste(Inds, "up", sep = "."), paste(Inds, "low", sep = ".")), 
      value.name = c("mean", "up", "low"), variable.name = "Index", variable.factor = F)
Toplot <- Toplot[data.table(IndName = Inds, Index = as.character(1:3)), on = "Index", Index := i.IndName]
Toplot[is.finite(mean), quant075 := quantile(mean, .75), by = c("Group", "Index")]
Toplot <- Toplot[is.finite(mean) & mean < quant075]


ggplot(Toplot, aes(x = Year, y = mean, group = Group, color = Group, fill = Group)) +
  geom_line() + scale_x_discrete(breaks = seq(min(Toplot$Year), max(Toplot$Year), by = 5)) +
  geom_ribbon(aes(ymax = up, ymin = low), color = "grey", alpha = 0.3) + facet_wrap( ~ Index, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("Stations Mean")
  

lapply(grep("TaxonomicBiodiv_", list.files(), value = T), load, .GlobalEnv)














