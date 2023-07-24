invisible(lapply(c("data.table", "httr", "jsonlite"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))



Numero_Stations <- c("I925301001","I928000101")

for (NumStat in Numero_Stations){
Station <- data.table()
extract <- fromJSON(paste0("https://hubeau.eaufrance.fr/api/v1/hydrometrie/obs_elab?code_entite=", NumStat,
                           "&grandeur_hydro_elab=QmJ", sep = ""))
while(nrow(Station) != extract$count){
      Station <- rbind(Station,extract$data)
      if(!is.null(extract$`next`)) {extract <- fromJSON(extract$`next`)}
}
Station[, date := as.Date(date_obs_elab)][, `:=`(year = format(date, "%Y"), month = format(date, "%m"), day = format(date, "%d"))]
Station <- Station[year %in% 2020:2022]

assign(x = paste0("QmJ_", NumStat), value = Station)
write.table(get(paste0("QmJ_", NumStat)), file = paste0("DebitJournalier_",NumStat, ".csv"), sep = ";", row.names = F)
}



extract <- fromJSON("https://wxs.ign.fr/mudrzvr5b9uj6pq8jkc8hk42/geoportail/wmts?SERVICE=WMTS&REQUEST=GetCapabilities")
