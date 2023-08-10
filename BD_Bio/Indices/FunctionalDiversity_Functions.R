invisible(lapply(c("data.table", "ggplot2", "entropart", "cluster", "dplyr", "StatMatch"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

Code_Diatom <- fread(file = paste0(getwd(),"/FinalMAJ_Naiades.Alric.Traits/Diatom_CodeTraits.csv"))
Code_Macroinvertebrate <- fread(file = paste0(getwd(),"/FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_CodeTraits.csv"))
Code_Fish <- fread(file = paste0(getwd(),"/FinalMAJ_Naiades.Alric.Traits/Fish_CodeTraits.csv"))


load_CompartmentFiles <- function(Directory, Compartment){
files <- list.files(Directory)
lapply(paste0(Directory,files[setdiff(grep(Compartment, files),grep(".csv", files))]), load, .GlobalEnv)
}

Split_Inv <- function(Inventory, max = 1e+03){
  return(lapply(split(unique(Inventory$CdStation),ceiling(seq_along(unique(Inventory$CdStation)) /max)), function(slice){
    return(Inventory[CdStation %in% slice]) }))
}

Indexes_measure_fun <- function(Abce, Names, Dissim){
  
  ret <- tryCatch({lapply(0:2, function(Q){round(expq(Hqz(as.AbdVector(setNames(Abce, Names)), q=Q, 
                                                          Dissim, Correction = "None"), q=Q),2)})},
                  
                  error = function(e) {print(e)})
  return(ret)
}

Inventory_Functional_long <- function(Code, Inventory, Functional_traits, level) {
  
  Functional_tree <- unique(Functional_traits)[, c("CdAppelTaxon", intersect(Code$Trait, colnames(Functional_traits))), with = F]
  
  pb <- txtProgressBar(min = 0, max = uniqueN(Code$Category), style = 3)
  
ret <- lapply(unique(Code$Category), function (trait) {
  
    traits <- unique(Code[Category == trait, Abbreviations])
    Tree <- Functional_tree[, traits, with = F][, (traits) := lapply(.SD, function(i){i[is.na(i)] <- 0; i}), .SDcols = traits][
      , lapply(.SD, as.factor)]
    Inv <- Inventory
    dissim_trait <- as.matrix(daisy(data.frame(Tree, row.names = Functional_tree[,CdAppelTaxon]), metric="gower"))
    dissim_trait <- 1 - dissim_trait/max(dissim_trait)
    
    Inv[, c("Richness_Fun", "Shannon_Fun", "Simpson_Fun") := Indexes_measure_fun(Abce = Abundance, Names = CdAppelTaxon_Join, 
          Dissim = dissim_trait), by = c("CdStation","Year")][, Category := gsub(" ","",trait)]
    
    setTxtProgressBar(pb, which(unique(Code$Category) == trait))
    
   return(unique(Inv[ ,c("Group", "CdStation", "Year", "Richness_Fun", "Shannon_Fun", "Simpson_Fun", "Category"), with = F]))
})
  return(do.call(rbind, ret))
}














