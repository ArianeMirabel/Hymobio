#####
invisible(lapply(c("devtools"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))
remotes::install_github("doi-usgs/EflowStats@v5.1.1")
#####
invisible(lapply(c("EflowStats", "data.table", "lmomco", "dplyr", "lubridate", "ggpubr", "ggplot2"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Temp")

# Link to Bio Stations
#####
load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")
Correspondance_station <- as.data.table(list_station_filter5_clean)[COMPARTIMENT_START %in% c("DIATOM", "MACROINVERTEBRATE", "FISH") & 
                          COMPARTIMENT_TARGET == "TEMPERATURE"][,.(COMPARTIMENT_START, ID_AMOBIO_START, COMPARTIMENT_TARGET, ID_TARGET)]

Temp_journaliere <- fread("C:/Users/armirabel/Documents/DB/TEMPERATURE/temperature_20220712.csv")
Correspondance_station <- Correspondance_station[ID_TARGET %in% unique(Temp_journaliere$stationID)]
rm(list = c("list_station_filter5_clean", "Temp_journaliere"))

load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/Fish_Inventories")
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Inventories")
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/Diatom_Inventories")

BioInv <- rbind(unique(AllInv_Fish[, ID_AMOBIO_START := paste0("FISH_", CdStation)][,.(CdStation, ID_AMOBIO_START, Date_PrelBio)]),
                unique(AllInv_Macroinvertebrate[, ID_AMOBIO_START := paste0("MACROINVERTEBRATE_", CdStation)][,.(CdStation, ID_AMOBIO_START, Date_PrelBio)]),
                unique(AllInv_Diatom[, ID_AMOBIO_START := paste0("DIATOM_", CdStation)][,.(CdStation, ID_AMOBIO_START, Date_PrelBio)]))
rm(list=c("AllInv_Fish", "AllInv_Macroinvertebrate", "AllInv_Diatom"))

Correspondance_station <- Correspondance_station[ID_AMOBIO_START %in% BioInv$ID_AMOBIO_START]
Correspondance_station <- left_join(Correspondance_station, BioInv, by = "ID_AMOBIO_START", relationship = "many-to-many")

load("CheckSingleSpeciesStations")
SingleSpeciesStations <- paste0("FISH_", unique(SingleSpeciesStations))

Correspondance_station <- Correspondance_station[!ID_AMOBIO_START %in% SingleSpeciesStations]
save(Correspondance_station, file = "Temp_BioInventories_correspondance")
#####

# mean, coefficient of variation, skewness, kurtosis, autoregressive lag-one (AR(1)) correlation coefficient, amplitude, phase of the seasonal signal
# use l-moments for the 4 first indicators
# deseasonalize for AR(1), substract the monthly mean stremfolw to each month value.
# Seasonality is measured using formula with sin, cas and tan

load("Temp_BioInventories_correspondance")
Temp_journaliere <- fread("C:/Users/armirabel/Documents/DB/TEMPERATURE/temperature_20220712.csv")

Temp_journaliere <- Temp_journaliere[!is.na(Tw_corr)]
Temp_journaliere[grep("/", Date), Samp_date := as.Date(Date, format = "%d/%m/%Y")][
  grep("-", Date), Samp_date := as.Date(Date, format = "%Y-%m-%d")][, Date := Samp_date][, Samp_date := NULL]
Temp_journaliere[, year := format(as.Date(Date, format =  "%Y-%m-%d"), "%Y")][
  , month := format(as.Date(Date, format =  "%Y-%m-%d"), "%m")][, day := format(as.Date(Date, format =  "%Y-%m-%d"), "%d")]


## Centered indices
######
# Amplitude measure, based on EflowStats functions
Amplitude <- function(x){
  
  x[, Deseas := (Tw_corr - mean(Tw_corr)), by = "month"][
    , Deseas := (Deseas - mean(Deseas)) / sd(Deseas)][
    , corrAR1 := round(ar(Deseas, aic = FALSE, order.max = 1, method = "yule-walker")$ar, 2)]
  
  jday <- lubridate::yday(x$Date)
  dec.year <- as.numeric(x$year) + (jday/365.25)
  
  scaled <- (x$Tw_corr - mean(x$Tw_corr)) / sd(x$Tw_corr)
  x_mat = cbind(1, sin(2 * pi * dec.year), cos(2 * pi * dec.year))
  mod <- .lm.fit(x_mat, scaled)
  a <- mod$coefficients[2]
  b <- mod$coefficients[3]
  
  Ampli <- round(sqrt((a^2) + (b^2)), digits = 2)
  
  Pha <- ifelse(a > 0, 365.25 * ((pi/2) - atan(b/a))/(2 * pi), 365.25 + 365.25 * ((pi/2) - pi - atan(b/a))/(2 * pi))
  if(a == 0 & b > 0){ Pha <- 365.25 }
  if(a == 0 & b < 0) { Pha <- 365.25/2 }
  if(a == 0 & b == 0) { Pha <- NA }
  
  return(data.table(period = unique(x$period),Amplitude = Ampli, Phase = Pha))
}

# Measure IHA for the sampling dates, for the given time steps
Index_timestep <- function(Tstep, Temp_Serie, station, tol_threshold){
 
   if(Tstep != 0){
  PrelDate <- data.table(Samp_date = Correspondance_station[ID_AMOBIO_START == station, Date_PrelBio])
  PrelDate[, period := 1:nrow(PrelDate)][, Start := Samp_date %m-% months(Tstep)]
  PrelDate <- left_join(PrelDate, PrelDate[, list(Date = seq.Date(from = Start, to = Samp_date, by = 'day')), by = "period"], by = "period")

Temp_serie <- merge(Temp_Serie, PrelDate, by = "Date")
if(nrow(Temp_serie) == 0) {return(NA)}
Temp_serie[, Tcover:= uniqueN(Date)/length(seq.Date(from = unique(Start), to = unique(Samp_date), by = 'day')), by = "period"]
Temp_serie <- Temp_serie[Tcover >= 0.75][,Tcover := NULL]
Temp_serie <- Temp_serie[, c("count", "unique") := list(.N, uniqueN(Tw_corr)), by = "period"][
  count >= 10 & unique > 1,][, c("count", "unique") := NULL]
if(nrow(Temp_serie) == 0) {return(NA)}
   } else {
     Temp_serie <- Temp_Serie[, c("count", "unique") := list(.N, uniqueN(Tw_corr))][
     count >= 10 & unique > 1,][, c("count", "unique") := NULL]
     if(nrow(Temp_serie) == 0) {return(NA)}
    Temp_serie = Temp_serie[, period := 0]
    PrelDate <- data.table(Samp_date = NA_Date_, period = 0, Start = NA_Date_)
    Tstep <- "all"}


Temp_serie[, c("Lmean", "Lscale", "Lskew", "Lkurt") := as.list(lmoms(Tw_corr, nmom = 4)$ratios[1:4]), by = "period"][
  , Lmean := mean(Tw_corr), by = "period"]

Temp_serie[, Deseas := (Tw_corr - mean(Tw_corr)), by = c("period", "month")][
  , Deseas := (Deseas - mean(Deseas)) / sd(Deseas), by = "period"][, corrAR1 := NA_real_]

Temp_serie[, corrAR1 := round(ar(Deseas, aic = FALSE, order.max = 1, method = "yule-walker")$ar, 2), by = "period"]

if(uniqueN(Temp_serie[,stationID])>1){print(paste("problem more sites", unique(Temp_serie[,stationID])))}

Ampli <- do.call(rbind,lapply(unique(Temp_serie$period), function(ti){
  return(Amplitude(Temp_serie[period == ti,]))}))

Temp_serie <- left_join(unique(Temp_serie[,.(stationID, period, Lmean, Lscale, Lskew, Lkurt, corrAR1)]),
               Ampli, by = "period")

Temp_serie <- right_join(unique(PrelDate[,.(Samp_date, period, Start)]), Temp_serie,
               by = "period")[
                 ,.(stationID, Samp_date, period, Lmean, Lscale, Lskew, Lkurt, corrAR1, Amplitude, Phase)]

setnames(Temp_serie, setdiff(colnames(Temp_serie), c("stationID", "Samp_date")), paste0(setdiff(colnames(Temp_serie), c("stationID", "Samp_date")), "_", Tstep))

return(Temp_serie)

}


rm(list=setdiff(ls(), c("Correspondance_station", "Temp_journaliere", "stationBio")))


Tol_threshold <- 0.75
severalDates <-vector()
Templaps <- lapply(unique(Correspondance_station[, ID_AMOBIO_START]), function(stationBio){

if(uniqueN(Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET]) == 1){
  Xtemp <- Temp_journaliere[stationID %in% unique(Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET])][
  , Date := as.Date(Date)]
} else {print(paste("Several Stations for", stationBio))}
  
Temp_laps <- tryCatch({lapply(c(3,6,12,0), function(tstep){
  return(Index_timestep(Tstep = tstep, Temp_Serie = Xtemp, station = stationBio, tol_threshold = Tol_threshold))})},
  error = function(e) print(paste("Error with", stationBio, e)))

Temp_all <- Temp_laps[[3]]

ifelse(any(unlist(lapply(Temp_laps[1:3],is.data.frame))), 
       Temp_laps <- Reduce(function(...) left_join(..., by = c("stationID", "Samp_date")), 
                           Temp_laps[1:3][which(unlist(lapply(Temp_laps[1:3],is.data.frame)))]),
       Temp_laps <- NA)


return(list(Temp_laps, Temp_all))

})
names(Templaps) <- unique(Correspondance_station[, ID_AMOBIO_START])

save(Hydrolaps, file = "HydroIndex_3612all")
#####
