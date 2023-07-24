invisible(lapply(c("data.table", "rgeos", "sp", "dplyr", "gridExtra", "lubridate", "bit64"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("../../DataBase_treatment/BD_Bio")

### Naiades DB
Nai_Station <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/stations.csv", 
               sep=";"))[,.(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, CodeRegion, LbRegion)]
Stations_metropole <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")), CdStationMesureEauxSurface]
Nai_Station <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")),]

Nai_Abd <- fread(paste0("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/fauneflore.csv"))[
  MnTypTaxRep == "NbrTax"][MnemoRqNbrTaxRep == "Domaine de validité"][LbSupport %in% c("Poissons", "Macroinvertébrés aquatiques","Diatomées benthiques")]

Nai_Abd[, Correction := as.numeric(Nai_Abd$CdStationMesureEauxSurface)][, Correction := as.character(Correction)]
Nai_Abd[is.na(Correction), Correction := CdStationMesureEauxSurface]
Nai_Abd[,CdStationMesureEauxSurface := Correction][,Correction := NULL]

Nai_Abd <- Nai_Abd[CdStationMesureEauxSurface %in% Stations_metropole]

Nai_Abd[, year := as.numeric(format(DateDebutOperationPrelBio, "%Y"))]
Nai_Abd <- Nai_Abd[year >= 2000]

Nai_Abd[, NomLatinAppelTaxon.origin := NomLatinAppelTaxon]

Nai_Abd[grep("Hybrides", NomLatinAppelTaxon), NomLatinAppelTaxon := gsub("Hybrides de ", "", NomLatinAppelTaxon)]
Nai_Abd[grep("Hybride", NomLatinAppelTaxon), NomLatinAppelTaxon := gsub("Hybride ", "", NomLatinAppelTaxon)]
Nai_Abd[grep("Appellation", NomLatinAppelTaxon), NomLatinAppelTaxon := "Inconnu"]

Nai_Abd_MI <- Nai_Abd[LbSupport == "Macroinvertébrés aquatiques"][,CdStationMesureEauxSurface]

Nai_Abd[grepl("Procambarus|Orconectes|Pacifastacus", NomLatinAppelTaxon) & CdStationMesureEauxSurface %in% Nai_Abd_MI, LbSupport := "Macroinvertébrés aquatiques"]

Nai_Abd_Fish <- Nai_Abd[LbSupport == "Poissons"][,.(CdStationMesureEauxSurface, RefOperationPrelBio, DateDebutOperationPrelBio, year, CdAppelTaxon, NomLatinAppelTaxon, RsTaxRep)]
Nai_Abd_Macroinv <- Nai_Abd[LbSupport == "Macroinvertébrés aquatiques"][,.(CdStationMesureEauxSurface, RefOperationPrelBio, DateDebutOperationPrelBio, year, CdAppelTaxon, NomLatinAppelTaxon, RsTaxRep)]
Nai_Abd_Diatom <- Nai_Abd[LbSupport == "Diatomées benthiques"][,.(CdStationMesureEauxSurface, RefOperationPrelBio, DateDebutOperationPrelBio, year, CdAppelTaxon, NomLatinAppelTaxon, RsTaxRep)]

save(Nai_Abd_Fish, file = "Naiades_Fish"); save(Nai_Abd_Macroinv, file = "Naiades_Macroinvertebrate"); save(Nai_Abd_Diatom, file = "Naiades_Diatom")








