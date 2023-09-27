invisible(lapply(c("data.table", "ggplot2", "stringr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

load("AMOBIO_WP3_2_FISH_INV_DIA_20230817.Rdata")
Carhyce_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN); rm("MATRIX_AMOBIO_WP3_CLEAN")

metrics <- data.table(Name = colnames(Carhyce_all))[grep("NC_|POE_|P_|H_|T_|HM_",substr(Name, 1,5))][, NameOrigin := Name][
  , Category := sub("_.*", "",Name)][, Nscales := str_count(Name, "_")][, Scale1 := sub(".*\\_", "", Name)][
  Scale1 == "ALL", Scale1 := 99][!Scale1 %in% c("NF","WD")][
  Scale1 == "F", Name := sub("_[^_]*$", "", Name)][Scale1 == "F", Scale1 := sub(".*\\_", "", Name)][
  , Numeric1 := as.numeric(gsub("\\D", "", Scale1))][
  is.na(Numeric1), Type := Scale1][!is.na(Type), Name := sub("_[^_]*$", "", Name)][!is.na(Type), Scale1 := sub(".*\\_", "", Name)][
  , Numeric1 := as.numeric(gsub("\\D", "", Scale1))][!is.na(Numeric1)][,Name := sub("_[^_]*$", "", Name)][
  , Nscales := str_count(Name, "_")][Nscales > 1, Scale2 := sub(".*\\_", "", Name)][, Numeric2 := as.numeric(gsub("\\D", "", Scale2))][
  !is.na(Numeric2), Name := sub('_[^_]*$', '', Name)][!is.na(Type), Name := paste(Name, Type, sep = "_")][,Type := NULL][,
  c("Min1", "Min2") := list(min(Numeric1), min(Numeric2)), by = Name]
metrics[is.na(metrics), ] <- 0 

Finale <- metrics[Numeric1 == Min1 & Numeric2 == Min2]

Catalog <- fread("ThresholdMetricCatalogue.csv", sep= ";")[,Name]

#Carhyce_all <- Carhyce_all[,grep("NC_|POE_|P_|H_|T_|HM_",substr(colnames(Carhyce_all), 1,5)), with = F]
Carhyce_all <- Carhyce_all[,grep(paste0(Catalog, collapse = "|"),colnames(Carhyce_all)), with = F]



Metric <- "HM_HABBED_Sh"
Data <- Carhyce_all[, grep(Metric, colnames(Carhyce_all)), with = F]
colnames(Data) <- "metric"
ggplot(data = Data, aes(metric)) +
  geom_density() + theme_classic() + xlab(Metric) + ylab("Density")
