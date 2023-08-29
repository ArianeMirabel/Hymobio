#!/usr/bin/Rscript 
source("FunctionalDiversity_Functions_Bash.R")

load("SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")

directory <- paste0(getwd(), "/FinalMAJ_Naiades.Alric.Traits/")
compartment <- "Diatom"
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
