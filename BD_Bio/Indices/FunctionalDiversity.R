invisible(lapply(c("data.table", "ggplot2", "entropart", "cluster", "dplyr", "StatMatch"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

dir <- "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/"
files <- list.files(dir)
lapply(paste0(dir,files[setdiff(grep("_Inventories|_Traits", files),grep(".csv", files))]), load, .GlobalEnv)


#Inventory_Functional_wide <- function(Code, Inventory, Functional_tree) {Inv <- Inventory;for(trait in unique(Code$Category)){
  traits <- Code[Category == trait, Abbreviations]
  dissim_trait <- as.matrix(daisy(data.frame(Functional_tree[, traits, with = F][, lapply(.SD, as.factor)], 
                                             row.names = Functional_tree$code), metric="gower"))
  dissim_trait <- 1 - dissim_trait/max(dissim_trait)
  
  Inv <- Inv[, paste(c("Richness_Fun", "Shannon_Fun", "Simpson_Fun"), gsub(" ","",trait), sep = "_") := Indexes_measure_fun(Abundance, CdAppelTaxon_Join, dissim_trait),
      by = c("CdStation","Year")]
}
#return(unique(Inv[ ,c("Group", "CdStation", "Year", grep("Richness|Shannon|Simpson", colnames(Inv), value = T)), with = F]))

#Abce = Inv$Abundance; Names=Inv$CdAppelTaxon_Join; Dissim = dissim_trait
Indexes_measure_fun <- function(Abce, Names, Dissim){
  ret <- tryCatch({lapply(0:2, function(Q){round(expq(Hqz(as.AbdVector(setNames(Abce, Names)), q=Q, 
                                                          Dissim, Correction = "None"), q=Q),2)})},
                  
                  error = function(e) {print(e)})
  return(ret)
}

#trait<-"Mig_SwimMode"
Inventory_Functional_long <- function(Code, Inventory, Functional_tree) {
  
ret <- lapply(unique(Code$Category), function (trait) {
  
    traits <- unique(Code[Category == trait, Abbreviations])
    Tree <- Functional_tree[, traits, with = F][, (traits) := lapply(.SD, function(i){i[is.na(i)] <- 0; i}), .SDcols = traits][
      , lapply(.SD, as.factor)]
    Inv <- Inventory
    dissim_trait <- as.matrix(daisy(data.frame(Tree, row.names = Functional_tree[,CdAppelTaxon]), metric="gower"))
    dissim_trait <- 1 - dissim_trait/max(dissim_trait)
    
    Inv[, c("Richness_Fun", "Shannon_Fun", "Simpson_Fun") := Indexes_measure_fun(Abce = Abundance, Names = CdAppelTaxon_Join, 
          Dissim = dissim_trait), by = c("CdStation","Year")][, Category := gsub(" ","",trait)]
    
   return(unique(Inv[ ,c("Group", "CdStation", "Year", "Richness_Fun", "Shannon_Fun", "Simpson_Fun", "Category"), with = F]))
})
  return(do.call(rbind, ret))
}

# Diatom
#####
Code_diatom <- fread(file = "Traits_code_Diatom.csv")
FunctionalTree_diatom <- unique(Traits_Diatom[, Code_diatom$Abbreviations := lapply(.SD, function(x){
  as.numeric(gsub(",",".",x))}), .SDcols = Code_diatom$Abbreviations])[, c("code", Code_diatom$Abbreviations), with = F]
setnames(FunctionalTree_diatom, "code", "CdAppelTaxon")

Index_Fun_Diatom <- Inventory_Functional_long(Code_diatom, AllInv_Diatom[!is.na(CdAppelTaxon_Join),], FunctionalTree_diatom)
#####

# Macroinv
#####
Code_macroinv <- fread(file = "Traits_code_Macroinv.csv")
FunctionalTree_macroinv <- unique(Traits_Macroinvertebrate)[
  , c("CdAppelTaxon", intersect(Code_macroinv$Abbreviations, colnames(Traits_Macroinvertebrate))), with = F]

Code = Code_macroinv; Inventory = AllInv_Macroinvertebrate[!CdAppelTaxon_Join %in% c("0", "Unmatched"),][1:1000,]
Functional_tree = FunctionalTree_macroinv
Index_Fun_Macroinv <- Inventory_Functional_long(Code_macroinv[Abbreviations %in% colnames(FunctionalTree_macroinv)],
                                              AllInv_Macroinvertebrate[!CdAppelTaxon_Join %in% c("0", "Unmatched"),][1:1000,],
                                              FunctionalTree_macroinv)
#####

#Fish
#####
Code_fish <- fread(file = "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Traits_code_Fish.csv")
FunctionalTree_fish <- unique(Traits_Fish)[, c("CdAppelTaxon", intersect(Code_fish$Abbreviations, colnames(Traits_Fish))), with = F]

Code = Code_fish; Inventory = AllInv_Fish[!is.na(CdAppelTaxon_Join),][1:1000,];Functional_tree = FunctionalTree_fish

Index_Fun_Fish <- Inventory_Functional_long(Code = Code_fish[Abbreviations %in% colnames(FunctionalTree_fish)],
                                            Inventory = AllInv_Fish[!is.na(CdAppelTaxon_Join),],
                                            Functional_tree = FunctionalTree_fish)


#####

# Plot
#####
Toplot <- Index_Fun_Macroinv
Toplot[, as.vector(outer(c("Mean", "Up", "Low"), grep("_Fun", colnames(Toplot), value = T), paste0)) := 
       do.call(c,lapply(.SD, function(x) {list(mean = mean(x, na.rm = T), 
                                               up = mean(x, na.rm = T)+sd(x, na.rm = T)/2, 
                                               low = mean(x, na.rm = T)-sd(x, na.rm = T)/2)})),
     .SDcols = grep("_Fun", colnames(Toplot), value = T), by = c("Category", "Year")]
Toplot <- melt(Toplot, 
    measure = list(grep("Mean", colnames(Toplot), value = T), grep("Up", colnames(Toplot), value = T), grep("Low", colnames(Toplot), value = T)), 
     value.name = c("Mean", "Up", "Low"), variable.name = "Index", variable.factor = F)

ggplot(Toplot, aes(x = Year, y = Mean, group = Category, color = Category, fill = Category)) +
  geom_line() + scale_x_discrete(breaks = seq(min(Toplot$Year), max(Toplot$Year), by = 5)) +
  geom_ribbon(aes(ymax = Up, ymin = Low), color = "grey", alpha = 0.3) + facet_wrap( ~ Index, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("Stations Mean")



Toplot <- Toplot[data.table(IndName = Inds, Index = as.character(1:3)), on = "Index", Index := i.IndName]
Toplot[is.finite(mean), quant075 := quantile(mean, .75), by = c("Group", "Index")]
Toplot <- Toplot[is.finite(mean) & mean < quant075]


ggplot(Toplot, aes(x = Year, y = mean, group = Group, color = Group, fill = Group)) +
  geom_line() + scale_x_discrete(breaks = seq(min(Toplot$Year), max(Toplot$Year), by = 5)) +
  geom_ribbon(aes(ymax = up, ymin = low), color = "grey", alpha = 0.3) + facet_wrap( ~ Index, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("Stations Mean")
  

lapply(grep("TaxonomicBiodiv_", list.files(), value = T), load, .GlobalEnv)














