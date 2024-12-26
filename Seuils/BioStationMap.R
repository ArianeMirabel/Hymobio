invisible(lapply(c("data.table", "rgdal", "sf", "ggplot2"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils")

load("../HYMOBIO_FULLDATA_202405.RData")

load("D:/Mes Donnees/Hymobio/DataBase_treatment/BD_Bio/AllStation_ArianeMarch2023")

CompartmentCols <- c("Fish" = "cornflowerblue", "Macroinvertebrates" = "tan2", "Diatoms" = "chartreuse4")

AllStations <- AllStations[CdStation %in% RFD_all$ID]

FrenchOutline <- st_read("BassinHydrographique_TOPAGE_UNION_20240301.shp")
st_as_sf(AllStations, coords = c("CoordX_L93", "CoordY_L93"))

ggplot() +
  geom_sf(data = FrenchOutline, fill="white", color="black") +
  theme_void() +
  geom_point(data = AllStations, aes(x = CoordX_L93, y = CoordY_L93, color = Group), shape = 3,
             position=position_dodge(width = 0.8))
