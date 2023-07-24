invisible(lapply(c("data.table", "rgeos", "sp", "dplyr", "gridExtra", "lubridate", "sf"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))


source("Fish_function.R"); source("Macroinv_function.R"); source("Diatom_function.R")

AllStations <- unique(merge(
  Nai_Abd[,.(CdStationMesureEauxSurface, DateDebutOperationPrelBio, LbSupport)], 
  Nai_Station[,.(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface)]))
setnames(AllStations, c("LbSupport", "DateDebutOperationPrelBio", "CdStationMesureEauxSurface", "CoordXStationMesureEauxSurface", "CoordYStationMesureEauxSurface"),
         c("Group", "Date","CdStation", "CoordX_L93", "CoordY_L93"))
AllStations[Group == "Poissons", Group := "Fish"][Group == "Macroinvertébrés aquatiques", Group := "Macroinvertebrates"][
  Group == "Diatomées benthiques", Group := "Diatoms"]
AllStations <- AllStations[, IDstation := paste(Group, CdStation, sep = "_")][, Date := as.character(Date)][,
               Dateop := format(as.Date(Date, format = "%Y-%m-%d"), "%Y.%m.%d")][,IDoperation := paste(IDstation, Dateop, sep="_")][,
               BDsource := "Naiades_FranceEntiere"][, MAJ := "Ariane_March2023"][, Dateop := NULL]
AllStations <- AllStations[,.(IDoperation, IDstation, Group, CdStation, CoordX_L93, CoordY_L93, BDsource, Date, MAJ)][order(IDoperation)]

Supp_FishAlr <- unique(Alr_Abd_Fish[CODE_STATION %in% setdiff(Alr_Abd_Fish$CODE_STATION, AllStations[Group == "Fish",CdStation]),
                                    .(DATE, CODE_STATION, X_L93, Y_L93)])
setnames(Supp_FishAlr, c("CODE_STATION", "DATE","X_L93", "Y_L93"), c("CdStation", "Date","CoordX_L93", "CoordY_L93"))
Supp_FishAlr[, Date := format(as.Date(Date, format = "%d/%m/%Y"), "%Y-%m-%d")][,
                                  Dateop := as.character(format(as.Date(Date, format = "%Y-%m-%d"), "%Y.%m.%d"))]
Supp_FishAlr[, Group := "Fish"][, BDsource := "Alric"][, IDstation := paste(Group, CdStation, sep = "_")][, 
                                  IDoperation := paste(IDstation, Dateop, sep="_")][,MAJ := "Ariane_March2023"][, Dateop := NULL]
Supp_FishAlr <- Supp_FishAlr[,.(IDoperation, IDstation, Group, CdStation, CoordX_L93, CoordY_L93, BDsource, Date, MAJ)][order(IDoperation)]

Supp_MacroinvertebratesAlr <- unique(Alr_Abd_MI[cd_site %in% setdiff(Alr_Abd_MI$cd_site,AllStations[Group == "Macroinvertebrates",CdStation]),
                                                .(cd_site, date_opecont, x, y)])
setnames(Supp_MacroinvertebratesAlr, c("cd_site", "date_opecont", "x", "y"), c("CdStation", "Date","CoordX_L93", "CoordY_L93"))
Supp_MacroinvertebratesAlr[, Date := format(as.Date(Date, format = "%d/%m/%Y"), "%Y-%m-%d")][,
                                  Dateop := as.character(format(as.Date(Date, format = "%Y-%m-%d"), "%Y.%m.%d"))]
Supp_MacroinvertebratesAlr[, Group := "Macroinvertebrates"][, BDsource := "Alric"][, IDstation := paste(Group, CdStation, sep = "_")][, 
                                  IDoperation := paste(IDstation, Dateop, sep="_")][,MAJ := "Ariane_March2023"][, Dateop := NULL]
Supp_MacroinvertebratesAlr <- Supp_MacroinvertebratesAlr[,.(IDoperation, IDstation, Group, CdStation, CoordX_L93, CoordY_L93, BDsource, Date, MAJ)][order(IDoperation)]

Supp_DiatomsAlr <- unique(Alr_Abd_DIA[cd_site %in% setdiff(Alr_Abd_DIA$cd_site, AllStations[Group == "Diatoms",CdStation]),
                                      .(cd_site, date_opecont, x, y)])
setnames(Supp_DiatomsAlr, c("cd_site", "date_opecont", "x", "y"), c("CdStation", "Date","CoordX_L93", "CoordY_L93"))
Supp_DiatomsAlr[, Date := format(as.Date(Date, format = "%d/%m/%Y"), "%Y-%m-%d")][,
                                  Dateop := as.character(format(as.Date(Date, format = "%Y-%m-%d"), "%Y.%m.%d"))]
Supp_DiatomsAlr[, Group := "Diatoms"][, BDsource := "Alric"][, IDstation := paste(Group, CdStation, sep = "_")][, 
                                  IDoperation := paste(IDstation, Dateop, sep="_")][,MAJ := "Ariane_March2023"][, Dateop := NULL]
Supp_DiatomsAlr <- Supp_DiatomsAlr[,.(IDoperation, IDstation, Group, CdStation, CoordX_L93, CoordY_L93, BDsource, Date, MAJ)][order(IDoperation)]


AllStations[CdStation %in% intersect(AllStations$CdStation, c(Alr_Abd_Fish$CODE_STATION, Alr_Abd_MI$cd_site, Alr_Abd_DIA$cd_site)),
            BDsource := "Naiades_Fentiere + Alric"]

AllStations <- rbind(AllStations, Supp_FishAlr, Supp_MacroinvertebratesAlr, Supp_DiatomsAlr)

save(AllStations, file = "AllStation_ArianeMarch2023")



##### Errors
Nai_Station_Fentiere <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/stations.csv", 
            sep=";"))[,.(CdStationMesureEauxSurface, LbStationMesureEauxSurface, CoordXStationMesureEauxSurface, 
            CoordYStationMesureEauxSurface, CodeRegion, LbRegion)][!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")),]

Nai_Abd_Fentiere <- fread(paste0("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/fauneflore.csv"))[
                          MnTypTaxRep == "NbrTax"][MnemoRqNbrTaxRep == "Domaine de validité"][
                          LbSupport %in% c("Poissons", "Macroinvertébrés aquatiques","Diatomées benthiques")]

Nai_Abd_Fentiere[, Correction := as.numeric(Nai_Abd_Fentiere$CdStationMesureEauxSurface)][, Correction := as.character(Correction)]
Nai_Abd_Fentiere[is.na(Correction), Correction := CdStationMesureEauxSurface]
Nai_Abd_Fentiere[,CdStationMesureEauxSurface := Correction][,Correction := NULL]


Nai_Station_Critere <- data.table()
Nai_Abd_Critere <- data.table()

for(i in 1:5){
  Nai_Station_Critere <- rbind(Nai_Station_Critere,
        as.data.table(read.csv(paste0("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/EXPORT_VIA_CRITERE_RECHERCHE/export",
        i, "/Stations.csv"), sep = ";", quote = ""))[,.(CdStationMesureEauxSurface, LbStationMesureEauxSurface, 
        CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, CodeRegion, LbRegion)])
  
  Nai_Abd_Critere <- rbind(Nai_Abd_Critere,
        as.data.table(read.csv(paste0("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/EXPORT_VIA_CRITERE_RECHERCHE/export",
        i, "/ListesFauneFlore.csv"), sep = ";", quote = "")))[MnTypTaxRep == "NbrTax"][MnemoRqNbrTaxRep == "Domaine de validité"][
        LbSupport %in% c("Poissons", "Macroinvertébrés aquatiques","Diatomées benthiques")]
}
Nai_Station_Critere <- Nai_Station_Critere[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")),]
Nai_Station_Critere[, Correction := as.numeric(Nai_Station_Critere$CdStationMesureEauxSurface)][, Correction := as.character(Correction)]
Nai_Station_Critere[is.na(Correction), Correction := CdStationMesureEauxSurface]
Nai_Station_Critere[,CdStationMesureEauxSurface := Correction][,Correction := NULL]

Nai_Abd_Critere[, Correction := as.numeric(Nai_Abd_Critere$CdStationMesureEauxSurface)][, Correction := as.character(Correction)]
Nai_Abd_Critere[is.na(Correction), Correction := CdStationMesureEauxSurface]
Nai_Abd_Critere[,CdStationMesureEauxSurface := Correction][,Correction := NULL]



Nai_Station_Critere[, geo := paste(CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, sep = "_")]
Nai_Station_Fentiere[, geo := paste(CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, sep = "_")]

unique(Nai_Station_Fentiere[,check := uniqueN(CdStationMesureEauxSurface), by = "geo"][check >= 2,])

uniqueN(setdiff(Nai_Abd_Critere$CdStationMesureEauxSurface, Nai_Abd_Fentiere$CdStationMesureEauxSurface))
uniqueN(setdiff(Nai_Abd_Fentiere$CdStationMesureEauxSurface, Nai_Abd_Critere$CdStationMesureEauxSurface))

uniqueN(setdiff(Nai_Station_Critere$geo, Nai_Station_Fentiere$geo))
uniqueN(setdiff(Nai_Station_Fentiere$geo, Nai_Station_Critere$geo))


france <- read_sf(dsn = "regions_2015_metropole_region.shp")
france <- st_geometry(st_set_crs(france, CRS("+init=epsg:2154")))

missings_Fentiere <- unique(Nai_Station_Critere[CdStationMesureEauxSurface %in% 
                     setdiff(Nai_Abd_Critere$CdStationMesureEauxSurface, Nai_Abd_Fentiere$CdStationMesureEauxSurface),
                    .(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface)])
missings_Fentiere <- st_as_sf(missings_Fentiere, coords = c('CoordXStationMesureEauxSurface', 'CoordYStationMesureEauxSurface'), crs = 2154)

plot(st_geometry(france))
plot(st_geometry(missings_Fentiere), add= T)


plot(missings_Fentiere, main = "Missings in France Entiere", col = "black", add = T)
missings_Fentiere <- st_as_sf(missings_Fentiere, CRS("+init=epsg:2154"))

ggplot() + geom_sf(data = missings_Fentiere) + theme(axis.text.x = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     rect = element_blank())


missings_Critere <- unique(Nai_Station_Fentiere[CdStationMesureEauxSurface %in% 
                    setdiff(Nai_Abd_Fentiere$CdStationMesureEauxSurface, Nai_Abd_Critere$CdStationMesureEauxSurface),
                    .(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface)])
coordinates(missings_Critere) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
plot(missings_Critere)


