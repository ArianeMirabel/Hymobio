source("C:/Users/armirabel/Documents/INRAE/Hymobio/Hymobio_GitHub/BD_Bio/Indices/FunctionalDiversity_Functions.R")

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

directory <- paste0(getwd(), "/FinalMAJ_Naiades.Alric.Traits/")
compartment <- "Macroinvertebrate"
slice <- 1

N_slices <- 5

invisible(load_CompartmentFiles(Directory = directory, Compartment = compartment))

AllInv_slice <- Split_Inv(Inventory = get(paste0("AllInv_",compartment))[!is.na(CdAppelTaxon_Join),],
                    Slice = slice, N_Slices = N_slices)

Index_Fun <- Inventory_Functional_long(
  Code = get(paste0("Code_",compartment))[Modalite %in% colnames(get(paste0("Traits_",compartment)))],
  Inventory = AllInv_slice, 
  Traits = get(paste0("Traits_",compartment)),
  hill = FALSE)[, Group := compartment]

Index_Fun <- left_join(Index_Fun, 
                       unique(as.data.table(list_station_filter5_clean)[COMPARTIMENT_START == toupper(compartment), 
                       .(ID_AMOBIO_START, ID_START)]), by = c("CdStation" = "ID_START"))[
            ,c("Group", "CdStation", "ID_AMOBIO_START", "Date_PrelBio", "Year", grep("R_|Sh_|Si_", colnames(Index_Fun), value = T)), with = F]

if(compartment == "Fish"){pref = "B_FISH"}
if(compartment == "Macroinvertebrate"){pref = "B_INV"}
if(compartment == "Diatom"){pref = "B_DIA"}
setnames(Index_Fun, grep("R_|Sh_|Si_", colnames(Index_Fun), value = T), 
         paste0(pref, grep("R_|Sh_|Si_", colnames(Index_Fun), value = T)))

save(Index_Fun, file = paste("Dfun_2401", compartment, slice, sep = "_"))


lapply(c("Fish", "Macroinvertebrate", "Diatom"), function(compartment){
  
  if(compartment == "Fish"){pref<-"B_FISH_"}
  if(compartment == "Macroinvertebrate"){pref<-"B_INV_"}
  if(compartment == "Diatom"){pref<-"B_DIA_"}
  
  assign(paste0("TSI_", compartment),do.call(rbind,lapply(1:N_slices, function(slice){
  invisible(load_CompartmentFiles(Directory = directory, Compartment = compartment))
  
  AllInv_slice <- Split_Inv(Inventory = get(paste0("AllInv_",compartment))[!is.na(CdAppelTaxon_Join),],
                            Slice = slice, N_Slices = N_slices)
  
  return(Inventory_Functional_TSI(
    Code = get(paste0("Code_",compartment))[Modalite %in% colnames(get(paste0("Traits_",compartment)))],
    Inventory = AllInv_slice, 
    Traits = get(paste0("Traits_",compartment))))
})))
  setnames(get(paste0("TSI_", compartment)), grep("TSI", colnames(get(paste0("TSI_", compartment))), value = T),
           paste0(pref, grep("TSI", colnames(get(paste0("TSI_", compartment))), value = T)))
  save(list = paste0("TSI_", compartment), 
       file = paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/TSI_2401_", compartment))
  

assign(paste0("NicheOverlap_", compartment),do.call(rbind,lapply(1:N_slices, function(slice){
  invisible(load_CompartmentFiles(Directory = directory, Compartment = compartment))
  
  AllInv_slice <- Split_Inv(Inventory = get(paste0("AllInv_",compartment))[!is.na(CdAppelTaxon_Join),],
                            Slice = slice, N_Slices = N_slices)
  
  return(Inventory_Functional_NicheOverlap(
    Code = get(paste0("Code_",compartment))[Modalite %in% colnames(get(paste0("Traits_",compartment)))],
    Inventory = AllInv_slice, 
    Traits = get(paste0("Traits_",compartment))))
})))
setnames(get(paste0("NicheOverlap_", compartment)), grep("meanNO", colnames(get(paste0("NicheOverlap_", compartment))), value = T),
         paste0(pref, grep("meanNO", colnames(get(paste0("NicheOverlap_", compartment))), value = T)))

save(list = paste0("NicheOverlap_", compartment), 
     file = paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/NicheOverlap_2401_", compartment))


})

#Save functional
#####
lapply(c("Fish", "Macroinvertebrate", "Diatom"),function(compartment){
  assign(paste0("DiversityFunctional_", compartment),
  do.call(rbind, lapply(paste("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Dfun_2401_",
  compartment, 1:N_slices, sep = "_"), function(file){load(file); return(Index_Fun)})
))
  
  load(paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/NicheOverlap_2401_", compartment))
  load(paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/TSI_2401_", compartment))
  
  assign(paste0("DiversityFunctional_", compartment),
         merge(merge(
           get(paste0("DiversityFunctional_", compartment)), 
           get(paste0("TSI_", compartment)), by = c("CdStation", "Year")),
           get(paste0("NicheOverlap_", compartment)), by = c("CdStation", "Year")))
  
  
  if(compartment == "Fish"){get(paste0("DiversityFunctional_", compartment))[, COMPARTMENT := "FISH"]}
  if(compartment == "Macroinvertebrate"){get(paste0("DiversityFunctional_", compartment))[, COMPARTMENT := "MACROINVERTEBRATE"]}
  if(compartment == "Diatom"){get(paste0("DiversityFunctional_", compartment))[, COMPARTMENT := "DIATOM"]}
  
  bys <- c("CdStation", "Year", "Group", "COMPARTMENT", "ID_AMOBIO_START")
  get(paste0("DiversityFunctional_", compartment))[
    ,Date_PrelBio := NULL][,setdiff(names(get(paste0("DiversityFunctional_",compartment))), bys) := lapply(.SD, function(X){mean(X,na.rm = T)}), by = bys, 
                           .SDcols = setdiff(names(get(paste0("DiversityFunctional_",compartment))), bys)]

  save(list = paste0("DiversityFunctional_", compartment), 
  file = paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_2401_", compartment))

  })
#####
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_2401_Fish")
DiversityFunctional_Fish <- unique(DiversityFunctional_Fish)
save(DiversityFunctional_Fish, 
     file = "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_2401_Fish")

load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_2401_Macroinvertebrate")
DiversityFunctional_Macroinvertebrate <- unique(DiversityFunctional_Macroinvertebrate)
save(DiversityFunctional_Macroinvertebrate, 
     file = "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_2401_Macroinvertebrate")

load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_2401_Diatom")
DiversityFunctional_Diatom <- unique(DiversityFunctional_Diatom)
save(DiversityFunctional_Diatom, 
     file = "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_2401_Diatom")

# Test with entropy
#####
lapply(c("Fish", "Macroinvertebrate", "Diatom"),function(compartment){
  assign("Diversity", get(load(paste0(
           "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_"
           , compartment))))
  
  Rich <- grep("R_", colnames(Diversity)); Shan <- grep("Sh_", colnames(Diversity)); Sim <- grep("Si_", colnames(Diversity))
  Diversity[, (Rich) := lapply(.SD, function(x) {lnq(x, q = 0)}), .SDcols = Rich][
    , (Shan) := lapply(.SD, function(x) {lnq(x, q = 1)}), .SDcols = Shan][
      , (Sim) := lapply(.SD, function(x) {lnq(x, q = 2)}), .SDcols = Sim]
  
  assign(paste0("DiversityFunctional_Entropy_", compartment), Diversity)
  
  save(list = paste0("DiversityFunctional_Entropy_", compartment),
       file = paste0(
         "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_Entropy",
         compartment))
})

#####

#Save all
#####
lapply(c("Fish", "Macroinvertebrate", "Diatom"),function(compartment){
  
  load(paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/DiversityTaxonomic_", compartment))
  
  get(paste0("DiversityTaxo_", compartment))[, Richness := lnq(Richness, q = 0)][, Shannon := lnq(Shannon, q = 1)][
    , Simpson := lnq(Simpson, q = 2)]
  
  assign("Fun", get(load(paste0(
    "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/DiversityFunctional_Entropy",
    compartment))))
  
  #Fun <- do.call(rbind, lapply(paste("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Dfun",
  #                                   compartment, 1:N_slices, sep = "_"), function(file){load(file); return(Index_Fun)}))
  
  assign(paste0("Biodiversity_", compartment), merge(Fun, get(paste0("DiversityTaxo_", compartment))[, .(CdStation, Year, Nspecies, Richness, Shannon, Simpson)], 
                                                    by = c("CdStation", "Year"), all = T))
  
  save(list = paste0("Biodiversity_", compartment), 
       file = paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Biodiversity_Entropy", compartment))
  
})

  
#####

#####

#Test with gawdis
#####
source("C:/Users/armirabel/Documents/INRAE/Hymobio/Hymobio_GitHub/BD_Bio/Indices/FunctionalDiversity_Functions.R")
library("gawdis")

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

directory <- paste0(getwd(), "/FinalMAJ_Naiades.Alric.Traits/")
compartment <- "Macroinvertebrate"
slice <- 1


N_slices <- 5

invisible(load_CompartmentFiles(Directory = directory, Compartment = compartment))

AllInv_slice <- Split_Inv(Inventory = get(paste0("AllInv_",compartment))[!is.na(CdAppelTaxon_Join),],
                          Slice = slice, N_Slices = N_slices)

Code = get(paste0("Code_",compartment))[Modalite %in% colnames(get(paste0("Traits_",compartment)))]
Inventory = AllInv_slice
Traits = get(paste0("Traits_",compartment))


Functional_tree <- unique(Traits)[, c("CdAppelTaxon", grep(paste(Code$Trait, collapse = "|"), colnames(Traits), value = T)),
                                    with = F]

Tree <- Functional_tree[, !"CdAppelTaxon"][, lapply(.SD, function(i){i[is.na(i)] <- 0; as.factor(i)})]
 
Gawdis_code <- merge(data.table(Modalite = colnames(Traits)), Code, by = "Modalite")[, Group := as.numeric(factor(Trait))][
  colnames(Tree), ] 

dissim_gawdis <- gawdis(data.frame(Tree, row.names = Functional_tree[,CdAppelTaxon]), w.type = "equal",
                        groups = Gawdis_code$Group, fuzzy = unique(Gawdis_code$Group))
dissim_Gawdis <- as.matrix(dissim_gawdis)
dissim_Gawdis <- 1 - dissim_Gawdis/max(dissim_Gawdis)
Indices_Total_Gawdis <- Inventory[, c("R_all_G", "Sh_all_G", "Si_all_G") := Indexes_measure_fun(Abce = Abundance,
                Names = CdAppelTaxon_Join, Dissim = dissim_Gawdis), by = c("CdStation","Year")]
Indices_Total_Gawdis <- unique(Indices_Total_Gawdis[ ,c("CdStation", "Year", "R_all_G", "Sh_all_G", "Si_all_G"), with = F])


dissim_trait <- as.matrix(daisy(data.frame(Tree, row.names = Functional_tree[,CdAppelTaxon]), metric="gower"))
dissim_trait <- 1 - dissim_trait/max(dissim_trait)
Indices_Total <- Inventory[, c("R_all", "Sh_all", "Si_all") := Indexes_measure_fun(Abce = Abundance,
                 Names = CdAppelTaxon_Join, Dissim = dissim_trait), by = c("CdStation","Year")]
Indices_Total <- unique(Indices_Total[ ,c("CdStation", "Year", "R_all", "Sh_all", "Si_all"), with = F])

Indices <- merge(Indices_Total, Indices_Total_Gawdis, by = c("CdStation", "Year"))
Indices[, as.vector(outer(c("Mean", "Up", "Low"), grep("_all", colnames(Indices), value = T), paste0)) := 
          do.call(c,lapply(.SD, function(x) {list(mean = mean(x, na.rm = T), 
                                                  up = mean(x, na.rm = T)+sd(x, na.rm = T)/2, 
                                                  low = mean(x, na.rm = T)-sd(x, na.rm = T)/2)})),
         .SDcols = grep("_all", colnames(Indices), value = T), by = "Year"]
Indices <- unique(Indices[,c("Year", grep("Mean|Up|Low", colnames(Indices), value = T)), with = F])
  
Toplot <- melt(Indices, 
               measure = list(grep("Mean", colnames(Indices), value = T), grep("Up", colnames(Indices), value = T), 
                              grep("Low", colnames(Indices), value = T)), 
               value.name = c("Mean", "Up", "Low"), variable.name = "Index", variable.factor = F)
Toplot[Index == 1, Index := "R_daisy"][Index == 2, Index := "Sh_daisy"][Index == 3, Index := "Si_daisy"][
  Index == 4, Index := "R_gaw"][Index == 5, Index := "Sh_gaw"][Index == 6, Index := "Si_gaw"]

ggplot(Toplot, aes(x = Year, y = Mean, group = Index, color = Index, fill = Index)) +
  geom_line() + scale_x_discrete(breaks = seq(min(Toplot$Year), max(Toplot$Year), by = 5)) +
  geom_ribbon(aes(ymax = Up, ymin = Low), color = "grey", alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("Stations Mean")
#####

# Plot
#####

Toplot <- lapply(c("Fish", "Macroinvertebrate", "Diatom"),function(compartment){
  return(assign(paste0("FunctionalDiversity_", compartment),
         do.call(rbind, lapply(paste("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Biodiversity",
                                     compartment, 1:N_slices, sep = "_"), function(file){load(file); return(Index_Fun)})
         )))
})

invisible(lapply(c("data.table", "ggplot2","grid", "gridExtra", "cowplot"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))


Toplot <- lapply(c("Fish", "Macroinvertebrate", "Diatom"),function(compartment){
  load(paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Biodiversity_Entropy",
              compartment))
  return(get(paste0("Biodiversity_", compartment))[,.(CdStation, Year, Group, ID_AMOBIO_START, 
                                                      R_all, Sh_all, Si_all, Nspecies, Richness, Shannon, Simpson)])})
names(Toplot) <- c("Fish", "Macroinvertebrate", "Diatom")
Toplot[["Fish"]][,Year := as.numeric(Year)]

plot_Biodiv <- function(Compartment, Inds, type){
  toplot <- Toplot[[Compartment]][, c("Group", "CdStation", "ID_AMOBIO_START", "Year", Inds), with = F]
  
  toplot[, as.vector(outer(c("Low", "Median", "Up"), Inds, paste0)) := 
           do.call(c, lapply(.SD, function(x) {as.list(quantile(x, probs = c(0.4, 0.5, 0.6), na.rm = T))})),
         .SDcols = Inds, by =  "Year"]
  toplot <- melt(toplot, 
                 measure = list(grep("Median", colnames(toplot), value = T), grep("Up", colnames(toplot), value = T), grep("Low", colnames(toplot), value = T)), 
                 value.name = c("Median", "Up", "Low"), variable.name = "Index", variable.factor = F)
  toplot[, Index := gsub("1", "Richness", gsub("2", "Shannon", gsub("3", "Simpson", Index)))]
  
  lgd <- (type == "Functional" & Compartment == "Macroinvertebrate")
  plot <- ggplot(toplot, aes(x = Year, y = Median, group = Index, color = Index, fill = Index)) +
    geom_line(show.legend = lgd, linewidth = 1) + scale_x_continuous(breaks = seq(2005, 2020, by = 5)) +
    geom_ribbon(aes(ymax = Up, ymin = Low), color = "grey", alpha = 0.3, show.legend = lgd) + ylab("Stations level\nIndices") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none") + 
    facet_wrap(~Index, ncol = 1, scales = "free_y") + 
    theme(strip.background = element_blank(), strip.text.x = element_blank(), panel.spacing.y = unit(0.7, "lines")) +
    scale_y_continuous(limits = ~ c(0-min(.x)*0.1, ceiling(max(.x))), breaks = ~ seq(0, .x[2], length.out = 3), expand = c(0, 0))
    
  if(type == "Functional"){plot <- plot + theme(axis.title.y = element_blank())}
  if(Compartment == "Fish"){ plot <- plot + ggtitle(paste(type, "Diversity"))}
  if(Compartment != "Diatom"){ plot <- plot + theme(axis.title.x = element_blank())}
  
  
  return(plot)
}

plot <- lapply(c("Fish", "Macroinvertebrate", "Diatom"),function(compartment){
  
  fun <- plot_Biodiv(Compartment = compartment, Inds = c("R_all", "Sh_all", "Si_all"), type = "Functional")
  
  taxo <- plot_Biodiv(Compartment = compartment, Inds = c("Richness", "Shannon", "Simpson"), type = "Taxonomic")
  
return(list(textGrob(compartment,gp=gpar(fontsize=20,font=3)), taxo, fun))
})

g <- ggplotGrob(plot[[2]][[3]] + theme(legend.position = "right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lwidth <- sum(legend$width)


grid.draw(arrangeGrob(arrangeGrob(grobs = do.call(c,plot), layout_matrix = matrix(c(1,1,2,3,4,4,5,6,7,7,8,9),ncol=2, byrow=T),
             heights = c(1,5,1,5,1,5)), legend, ncol = 2,
            widths = unit.c(unit(1, "npc") - lwidth, lwidth)))



lapply(grep("TaxonomicBiodiv_", list.files(), value = T), load, .GlobalEnv)
#####
