invisible(lapply(c("data.table", "ggplot2", "entropart"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

dir <- "C:/Users/armirabel/Documents/DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/FinalMAJ_Naiades.Alric.Traits/"
files <- list.files(dir)
lapply(paste0(dir,files[setdiff(grep("_Inventories", files),grep(".csv", files))]), load, .GlobalEnv)

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
  
  Inv <- unique(Inv[,.(Group, CdStation, Year, NomTaxo, Abundance)][
    , Abundance := ceiling(mean(Abundance)), by = c("CdStation", "NomTaxo", "Year")])
  
  Inv <- unique(Inv[
    , c("Nspecies", Inds, "Bias_Estimator") := c(uniqueN(NomTaxo), Indexes_measure(Abundance)), 
    by = c("CdStation","Year")][ ,.(Group, CdStation, Year, Nspecies, Richness, Shannon, Simpson)])
  Inv[Richness == Inf, Richness := Nspecies]
  return(Inv)
}

Index_Taxo_Diatom <- Inventory_taxo(AllInv_Diatom)
save(Index_Taxo_Diatom, file = "TaxonomicBiodiv_Diatom")

Index_Taxo_Macroinvertebrate <- Inventory_taxo(AllInv_Macroinvertebrate)
save(Index_Taxo_Macroinvertebrate, file = "TaxonomicBiodiv_Macroinvertebrate")

Index_Taxo_Fish <- Inventory_taxo(AllInv_Fish)
save(Index_Taxo_Fish, file = "TaxonomicBiodiv_Fish")


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














