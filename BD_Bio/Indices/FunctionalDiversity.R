source("C:/Users/armirabel/Documents/INRAE/Hymobio/Hymobio_GitHub/BD_Bio/Indices/FunctionalDiversity_Functions.R")

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

directory <- paste0(getwd(), "/FinalMAJ_Naiades.Alric.Traits/")
compartment <- "Diatom"
slice <- 1


invisible(load_CompartmentFiles(Directory = directory, Compartment = compartment))

AllInv_slice <- Split_Inv(Inventory = get(paste0("AllInv_",compartment))[!is.na(CdAppelTaxon_Join),],
                    Slice = slice, max = 100)

Index_Fun <- Inventory_Functional_long(
  Code = get(paste0("Code_",compartment))[Modalite %in% colnames(get(paste0("Traits_",compartment)))],
  Inventory = AllInv_slice, 
  Traits = get(paste0("Traits_",compartment)))[, Group := compartment]

Index_Fun <- left_join(Index_Fun, 
                       unique(as.data.table(list_station_filter5_clean)[COMPARTIMENT_START == toupper(compartment), 
                       .(ID_AMOBIO_START, ID_START)]), by = c("CdStation" = "ID_START"))[
            ,c("Group", "CdStation", "ID_AMOBIO_START", "Year", grep("R_|Sh_|Si_", colnames(Index_Fun), value = T)), with = F]

save(Index_Fun, file = paste("Dfun", compartment, slice, sep = "_"))






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
#####
