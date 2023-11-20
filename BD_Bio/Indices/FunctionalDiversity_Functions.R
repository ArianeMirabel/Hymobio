invisible(lapply(c("data.table", "ggplot2", "entropart", "cluster", "dplyr", "StatMatch", "vegan"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

#Code_Diatom <- fread(file = paste0(getwd(),"/FinalMAJ_Naiades.Alric.Traits/Diatom_CodeTraits.csv"))
#Code_Macroinvertebrate <- fread(file = paste0(getwd(),"/FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_CodeTraits.csv"))
#Code_Fish <- fread(file = paste0(getwd(),"/FinalMAJ_Naiades.Alric.Traits/Fish_CodeTraits.csv"))


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
  Traits_Fish <- merge(Traits_Fish[,c("CdAppelTaxon", grep("food.DET|food.PLANT|food.ANI", names(Traits_Fish), value = T)), with = F], 
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

Inventory_Functional_TSI <- function(Code, Inventory, Traits) {
  
  TSI_tree <- unique(Traits)[, c("CdAppelTaxon", Code[Type == "Ecology_Food", Modalite]), with = F]
  TSI_tree <- TSI_tree[, lapply(.SD, as.numeric), .SDcols = setdiff(names(TSI_tree), "CdAppelTaxon"), by = CdAppelTaxon]
  
  TSI_tree <- unique(TSI_tree[, Fsize := (sum(.SD^2) - (1/nrow(Code[Category == "Food_Size"])))/(1-(1/nrow(Code[Category == "Food_Size"])))
  , .SDcols = Code[Category == "Food_Size", Modalite], by = CdAppelTaxon][,.(CdAppelTaxon, Fsize)])
  
  Inventory <- merge(Inventory, TSI_tree, by.x = "CdAppelTaxon_Join", by.y = "CdAppelTaxon", all.x = T)
  
  Inventory[, TSI := sum(Fsize * log(Abundance +1)) / sum(log(Abundance + 1)) , by = c("CdStation","Year")]
  
  return(unique(Inventory[ ,.(CdStation,Year, TSI)]))
}


Inventory_Functional_NicheOverlap <- function(Code, Inventory, Traits) {
  
  NO_tree <- unique(Traits)[, c("CdAppelTaxon", Code[Type == "Ecology_Food", Modalite]), with = F]
  NO_tree <- NO_tree[, lapply(.SD, as.numeric), .SDcols = setdiff(names(NO_tree), "CdAppelTaxon"), by = CdAppelTaxon]
  NO_tree <- as.matrix(as.matrix(designdist(NO_tree[, Code[Category == "Food_Size", Modalite], with = F], 
                       method = "J/sqrt(A*B)", terms = "quadratic")))
  rownames(NO_tree) <- unique(Traits[,CdAppelTaxon]); colnames(NO_tree) <- unique(Traits[,CdAppelTaxon])
  
  return(unique(Inventory[, 
                      meanNO := mean(NO_tree[rownames(NO_tree) %in% CdAppelTaxon_Join, colnames(NO_tree) %in% CdAppelTaxon_Join])
                      , by = c("CdStation","Year")][ ,.(CdStation, Year, meanNO)]))
}

Indexes_measure_fun <- function(Abce, Names, Dissim, Hill){
  
  ret <- tryCatch({lapply(0:2, function(Q){
    ifelse(Hill,round(expq(Hqz(as.AbdVector(setNames(Abce, Names)), q=Q, 
                Dissim, Correction = "None"), q=Q),2),
           round(Hqz(as.AbdVector(setNames(Abce, Names)), q=Q, 
                          Dissim, Correction = "None"),2))
    })},
                  error = function(e) {print(e)})
  return(ret)
}


Inventory_Functional_long <- function(Code, Inventory, Traits, hill = TRUE) {
  
  Functional_tree <- unique(Traits)[, c("CdAppelTaxon", grep(paste(Code$Trait, collapse = "|"), colnames(Traits), value = T)),
                                    with = F]
  
  pb <- txtProgressBar(min = 0, max = uniqueN(Code$Trait)*2 + 1, style = 3)
  
   Indices_Modalite <- lapply(unique(Code$Trait), function (trait) {
      
      ret <- lapply(Code[Trait == trait, Modalite], function(modalite){
        
        Inv <- copy(Inventory)
        
        Tree <- Functional_tree[, modalite, with = F][, (modalite) := lapply(.SD, function(i){i[is.na(i)] <- 0; i}), 
                                                      .SDcols = modalite][, lapply(.SD, as.factor)]
        
        dissim_trait <- as.matrix(daisy(data.frame(Tree, row.names = Functional_tree[,CdAppelTaxon]), metric="gower"))
        dissim_trait <- 1 - dissim_trait/max(dissim_trait)
        
        Inv[, paste0(c("R_Mod_", "Sh_Mod_", "Si_Mod_"), modalite) := Indexes_measure_fun(Abce = Abundance, Names = CdAppelTaxon_Join, 
             Dissim = dissim_trait, Hill = hill), by = c("CdStation","Year")]
        
        return(unique(Inv[ ,c("CdStation", "Year", grep(modalite, colnames(Inv), value = T)), with = F]))
        })
        
      setTxtProgressBar(pb, which(unique(Code$Trait) == trait))
      
      ret <- tryCatch({
        Reduce(function(...) merge(..., by = c("CdStation", "Year")), ret)},
              warning = function(w) {print(paste(w, "\n", "On trait", trait))},
              error = function(e) {print(paste(e, "\n", "On trait", trait))})
      
      return(ret)
      
      })
      
   Indices_Modalite <- tryCatch({
     Reduce(function(...) merge(..., by = c("CdStation", "Year")), Indices_Modalite)},
     warning = function(w){ print(paste(w, "\n", "On Modalite"))},
     error = function(e){ print(paste(e, "\n", "On Modalite"))})
   
   
   Indices_Trait <- lapply(unique(Code$Trait), function (trait) {
     
     Inv <- copy(Inventory)
     
       Tree <- Functional_tree[, Code[Trait == trait, Modalite], with = F][
         , (Code[Trait == trait, Modalite]) := lapply(.SD, function(i){i[is.na(i)] <- 0; i}), 
         .SDcols = Code[Trait == trait, Modalite]][, lapply(.SD, as.factor)]
       
       dissim_trait <- as.matrix(daisy(data.frame(Tree, row.names = Functional_tree[,CdAppelTaxon]), metric="gower"))
       dissim_trait <- 1 - dissim_trait/max(dissim_trait)
       
       Inv[, paste0(c("R_Tr_", "Sh_Tr_", "Si_Tr_"), trait) := Indexes_measure_fun(Abce = Abundance,
                Names = CdAppelTaxon_Join, Dissim = dissim_trait), by = c("CdStation","Year")]
       
       setTxtProgressBar(pb, which(unique(Code$Trait) == trait) + uniqueN(Code$Trait))
       
       return(unique(Inv[ , c("CdStation", "Year",  grep(trait, colnames(Inv), value = T)), with = F]))
     })
     
   Indices_Trait <- tryCatch({
     Reduce(function(...) merge(..., by = c("CdStation", "Year")), Indices_Trait)},
                                warning = function(w){ print(paste(w, "\n", "On Trait"))},
                                error = function(e){ print(paste(e, "\n", "On Trait"))}) 
   
   Tree <- Functional_tree[, !"CdAppelTaxon"][, lapply(.SD, function(i){i[is.na(i)] <- 0; as.factor(i)})]
     
   dissim_trait <- as.matrix(daisy(data.frame(Tree, row.names = Functional_tree[,CdAppelTaxon]), metric="gower"))
   dissim_trait <- 1 - dissim_trait/max(dissim_trait)
     
   Indices_Total <- Inventory[, c("R_all", "Sh_all", "Si_all") := Indexes_measure_fun(Abce = Abundance,
                                Names = CdAppelTaxon_Join, Dissim = dissim_trait), by = c("CdStation","Year")]
     
   Indices_Total <- unique(Indices_Total[ ,c("CdStation", "Year", "R_all", "Sh_all", "Si_all"), with = F])


  return(Reduce(function(...) merge(..., by = c("CdStation", "Year")), 
                list(Indices_Modalite, Indices_Trait, Indices_Total)))
}




Inventory_TraitsFrequencies <- function(Code, Inventory, Traits, hill = TRUE) {
  
  Functional_tree <- unique(Traits)[, c("CdAppelTaxon", grep(paste(Code$Trait, collapse = "|"), colnames(Traits), value = T)),
                                    with = F]
  
  pb <- txtProgressBar(min = 0, max = uniqueN(Code$Trait)*2 + 1, style = 3)
  
  Indices_Modalite <- lapply(unique(Code$Trait), function (trait) {
    
    ret <- lapply(Code[Trait == trait, Modalite], function(modalite){
      
      Tree <- Functional_tree[, c("CdAppelTaxon",modalite), with = F][, (modalite) := lapply(.SD, function(i){i[is.na(i)] <- 0; i}), 
                                                    .SDcols = modalite]
      
      Inv <- merge(Inventory[,.(CdStation, Year, CdAppelTaxon, Abundance)], Tree, by = "CdAppelTaxon", all.x = T, all.y = F)
      Inv[, Freq := Abundance/sum(Abundance), by = c("CdStation","Year")]
      
      Inv[, paste0("N_", modalite) := sum(Abundance * as.numeric(get(modalite))), by = c("CdStation","Year")][
          , paste0("Freq_", modalite) := sum(Freq * as.numeric(get(modalite))), by = c("CdStation","Year")]
      
      return(unique(Inv[ ,c("CdStation", "Year", paste0("N_", modalite), paste0("Freq_", modalite)), with = F]))
    })
    
    setTxtProgressBar(pb, which(unique(Code$Trait) == trait))
    
    ret <- tryCatch({
      Reduce(function(...) merge(..., by = c("CdStation", "Year"), all =T), ret)},
      warning = function(w) {print(paste(w, "\n", "On trait", trait))},
      error = function(e) {print(paste(e, "\n", "On trait", trait))})
    
    return(ret)
    
  })
  
  Indices_Modalite <- tryCatch({
    Reduce(function(...) merge(..., by = c("CdStation", "Year")), Indices_Modalite)},
    warning = function(w){ print(paste(w, "\n", "On Modalite"))},
    error = function(e){ print(paste(e, "\n", "On Modalite"))})

  
  return(Indices_Modalite)}












