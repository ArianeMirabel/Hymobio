#load packages
invisible(lapply(c("data.table", "ggplot2", "entropart", "cluster", "dplyr", "StatMatch"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

# Create functions
#####

load_CompartmentFiles <- function(Directory, Compartment){
files <- list.files(Directory)
lapply(paste0(Directory,files[setdiff(grep(Compartment, files),grep(".csv", files))]), load, .GlobalEnv)

if(Compartment == "Fish"){
  Traits_Fish <- unique(Traits_Fish)[, c("CdAppelTaxon", intersect(Code_Fish$Trait, colnames(Traits_Fish))), with = F]
  ToCast <- names(Traits_Fish[, .SD, .SDcols = -c("CdAppelTaxon", grep("food.DET|food.PLANT|food.ANI", names(Traits_Fish), value = T))])
  
  Traits_Fish_modif <- lapply(ToCast, function(col){
    ret <- as.data.table(dcast(Traits_Fish[, list(V1=1, CdAppelTaxon, get(col))], CdAppelTaxon ~ V3,
                               fun=sum, value.var="V1", drop=c(TRUE, FALSE)))
    ret <- ret[,.SD,.SDcols = setdiff(names(ret), "NA")]
    setnames(ret, names(ret[,!"CdAppelTaxon"]), paste0(col, "_", names(ret[,!"CdAppelTaxon"])))
    return(ret)})
  Traits_Fish_modif <- Reduce(function(...) merge(..., by = "CdAppelTaxon"), Traits_Fish_modif)
  Traits_Fish <- merge(Traits_Fish[,.SD,.SDcols = c("CdAppelTaxon", grep("food.DET|food.PLANT|food.ANI", names(Traits_Fish), value = T))], 
                       Traits_Fish_modif, by = "CdAppelTaxon")
  assign("Traits_Fish", Traits_Fish, envir = .GlobalEnv)
  
  Code_Fish[, Modalite_simple := Modalite][!grep("food.DET|food.PLANT|food.ANI",Trait), Modalite := paste0(Trait, "_", Modalite)]
  assign("Code_Fish", Code_Fish, envir = .GlobalEnv)
}

if(Compartment == "Diatom"){
  setnames(Traits_Diatom, c( "guilde_h","guilde_l","guilde_m","guilde_e","acidobi","acidoph","neutro","alcaliph","alcalib",
                             "acid_indif","douce","douce_saum","mod_saum","saumatre","N_auto1","N_auto2","N_hetero1","N_hetero2",
                             "ox_high","ox_forte","ox_mod","ox_low","ox_verylow","oligosaprobe","B_meso","A_meso","A_meso_poly",
                             "polysaprobe","oligotroph","oligo_meso","meso","meso_eutro","eutrophe","hypereutro","troph_indif",
                             "centriques","pennees_araphidees","pennees_monoraphidees","pennees_biraphidees","planktonic","benthic",
                             "mob","no_mob","jamais_attach","attach_vu","solitaire","col_muc","col_zigzag","col_ruban","col_etoile",
                             "col_arbuscule","C1","C2","C3","C4","C5","C6","C7","gf1","gf2","gf3","gf4"), Code_Diatom$Modalite)}
}


Split_Inv <- function(Inventory, Slice, N_Slices){
  SplitList <- split(unique(Inventory$CdStation),ceiling(seq_along(unique(Inventory$CdStation)) /
                                                           ceiling(uniqueN(Inventory$CdStation)/N_slices)))
   return(Inventory[CdStation %in% SplitList[[Slice]]])}


Inventory_TraitsFrequencies <- function(Code, Inventory, Traits) {
  
  pb <- txtProgressBar(min = 0, max = uniqueN(Code$Trait), style = 3)
  
  # get functional tree
  Functional_tree <- unique(Traits)[, c("CdAppelTaxon", grep(paste(Code$Trait, collapse = "|"), colnames(Traits), value = T)),
                                    with = F]
  
  # lapply for each trait
  Indices_Modalite <- lapply(unique(Code$Trait), function (trait) {
    
    # lapply for each modality
    ret <- lapply(Code[Trait == trait, Modalite], function(modalite){
      
      # get the species x modality matrix en merge with inventories
      Tree <- Functional_tree[, c("CdAppelTaxon",modalite), with = F][, (modalite) := lapply(.SD, function(i){i[is.na(i)] <- 0; i}), 
                                                    .SDcols = modalite]
      
      Inv <- merge(Inventory[,.(CdStation, Year, CdAppelTaxon, Abundance)], Tree, by = "CdAppelTaxon", all.x = T, all.y = F)
      
      # For the abundance or frequency of modality in the inventory: multiply species abundance (or frequency) by the modality 
      # (if the species doesn't have the modaloty, it's null), then sum by station and year
      
      Inv[, Freq := Abundance/sum(Abundance), by = c("CdStation","Year")]
      
      Inv[, paste0("N_", modalite) := sum(Abundance * as.numeric(get(modalite))), by = c("CdStation","Year")][
          , paste0("Freq_", modalite) := sum(Freq * as.numeric(get(modalite))), by = c("CdStation","Year")]
      
      return(unique(Inv[ ,c("CdStation", "Year", paste0("N_", modalite), paste0("Freq_", modalite)), with = F]))
    })
    
    setTxtProgressBar(pb, which(unique(Code$Trait) == trait))
    
    # merge the lapply return that is a inventory x modality frequency matrix for each trait
    ret <- tryCatch({
      Reduce(function(...) merge(..., by = c("CdStation", "Year"), all =T), ret)},
      warning = function(w) {print(paste(w, "\n", "On trait", trait))},
      error = function(e) {print(paste(e, "\n", "On trait", trait))})
    
    return(ret)
    
  })
  
  # merge the lapply return that is a inventory x trait modality frequency matrix
  Indices_Modalite <- tryCatch({
    Reduce(function(...) merge(..., by = c("CdStation", "Year")), Indices_Modalite)},
    warning = function(w){ print(paste(w, "\n", "On Modalite"))},
    error = function(e){ print(paste(e, "\n", "On Modalite"))})

  
  return(Indices_Modalite)}

#####

# Process

#####

# Set direction for files. My file is a general "Bio" file with the subfile "/FinalMAJ_Naiades.Alric.Traits/" for the final inventories and traits data base
setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

# Set direction for inventories
directory <- paste0(getwd(), "/FinalMAJ_Naiades.Alric.Traits/")

# If Genomig, make one file by compartment and i.slice, and replace lapply by simple function and final save
# Compartment to study
# compartment <- "Macroinvertebrate"
# Slice number
# i.slice <- 1

# Number of database subdivision (for computation speed)
N_slices <- 5

# Load files
load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

Modality_save <- lapply(c("Fish", "Macroinvertebrate", "Diatom"), function(compartment){
  

invisible(load_CompartmentFiles(Directory = directory, Compartment = compartment))
  
  print(compartment)
  
FreqMod <- lapply(1:N_slices, function(i.slice){
  
  print(paste0("Slice ", i.slice, "/", N_slices))
  
AllInv_slice <- Split_Inv(Inventory = get(paste0("AllInv_",compartment))[!is.na(CdAppelTaxon_Join),],
                          Slice = i.slice, N_Slices = N_slices)

Index_Fun <- Inventory_TraitsFrequencies(
  Code = get(paste0("Code_",compartment))[Modalite %in% colnames(get(paste0("Traits_",compartment)))],
  Inventory = AllInv_slice, 
  Traits = get(paste0("Traits_",compartment)))[, Group := compartment]

# Add Amobio reference code
Index_Fun <- left_join(Index_Fun, 
            unique(as.data.table(list_station_filter5_clean)[COMPARTIMENT_START == toupper(compartment), 
            .(ID_AMOBIO_START, ID_START)]), by = c("CdStation" = "ID_START"))

return(Index_Fun)

})

FreqMod <- do.call(rbind,FreqMod)


save(FreqMod, file = paste0("FreqModalite_", compartment))

rm(list = grep(compartment, ls()))
})









