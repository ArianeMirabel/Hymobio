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
  Traits = get(paste0("Traits_",compartment)))[, Group := compartment]

Index_Fun <- left_join(Index_Fun, 
                       unique(as.data.table(list_station_filter5_clean)[COMPARTIMENT_START == toupper(compartment), 
                       .(ID_AMOBIO_START, ID_START)]), by = c("CdStation" = "ID_START"))[
            ,c("Group", "CdStation", "ID_AMOBIO_START", "Year", grep("R_|Sh_|Si_", colnames(Index_Fun), value = T)), with = F]

save(Index_Fun, file = paste("Dfun", compartment, slice, sep = "_"))


lapply(c("Fish", "Macroinvertebrate", "Diatom"),function(compartment){
  assign(paste0("FunctionalDiversity_", compartment),
  do.call(rbind, lapply(paste("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Dfun",
  compartment, 1:N_slices, sep = "_"), function(file){load(file); return(Index_Fun)})
))
  save(list = paste0("FunctionalDiversity_", compartment), 
  file = paste0("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Dfun_tot_", compartment))

  })


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
         do.call(rbind, lapply(paste("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/Dfun",
                                     compartment, 1:N_slices, sep = "_"), function(file){load(file); return(Index_Fun)})
         )))
})

toplot <- Toplot[[1]][, c("Group", "CdStation", "ID_AMOBIO_START", "Year", "R_all", "Sh_all", "Si_all"), with = F]

toplot[, as.vector(outer(c("Mean", "Up", "Low"), grep("_al", colnames(toplot), value = T), paste0)) := 
         do.call(c,lapply(.SD, function(x) {list(mean = mean(x, na.rm = T), 
                                                 up = mean(x, na.rm = T)+sd(x, na.rm = T)/2, 
                                                 low = mean(x, na.rm = T)-sd(x, na.rm = T)/2)})),
       .SDcols = grep("_all", colnames(toplot), value = T), by =  "Year"]
toplot <- melt(toplot, 
               measure = list(grep("Mean", colnames(toplot), value = T), grep("Up", colnames(toplot), value = T), grep("Low", colnames(toplot), value = T)), 
               value.name = c("Mean", "Up", "Low"), variable.name = "Index", variable.factor = F)
toplot[, Index := gsub("1", "Richness", gsub("2", "Shannon", gsub("3", "Simpson", Index)))]

ggplot(toplot, aes(x = Year, y = Mean, group = Index, color = Index, fill = Index)) +
  geom_line() + scale_x_discrete(breaks = seq(min(toplot$Year), max(toplot$Year), by = 5)) +
  geom_ribbon(aes(ymax = Up, ymin = Low), color = "grey", alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("Stations Mean")



Toplot <- Toplot[data.table(IndName = Inds, Index = as.character(1:3)), on = "Index", Index := i.IndName]
Toplot[is.finite(mean), quant075 := quantile(mean, .75), by = c("Group", "Index")]
Toplot <- Toplot[is.finite(mean) & mean < quant075]


ggplot(Toplot, aes(x = Year, y = mean, group = Group, color = Group, fill = Group)) +
  geom_line() + scale_x_discrete(breaks = seq(min(Toplot$Year), max(Toplot$Year), by = 5)) +
  geom_ribbon(aes(ymax = up, ymin = low), color = "grey", alpha = 0.3) + facet_wrap( ~ Index, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylab("Stations Mean")


lapply(grep("TaxonomicBiodiv_", list.files(), value = T), load, .GlobalEnv)
#####
