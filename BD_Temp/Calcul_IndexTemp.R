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

# mean, coefficient of variation, skewness, kurtosis, autoregressive lag-one (AR(1)) correlation coefficient, amplitude, phase of the seasonal signal
# use l-moments for the 4 first indicators
# deseasonalize for AR(1), substract the monthly mean stremfolw to each month value.
# Seasonality is measured using formula with sin, cas and tan
fread("C:/Users/armirabel/Documents/DB/TEMPERATURE/temperature_20220712.csv")

##Cleaning data
#####
Hydro_journaliere[, resultat_1 := shift(resultat_obs_elab, type = "lead")][, delta := resultat_1 - resultat_obs_elab]

Xzero <- Hydro_journaliere[, .(code_station, Date, resultat_obs_elab, resultat_1, delta)][delta == 0]

Run_periods <- lapply(unique(Xzero$code_station), function(st){
  x <- Xzero[code_station == st, .(Date)]
  x <- x[order(Date)]
  run <- 1
  x[, Run_period := run][, Date_check := seq.Date(from = first(Date), length.out = length(Date), by = "days")-Date]
  
  while(any(x[, Date_check] != 0)){
    x[Date_check != 0, c("Run_period", "Date_check") := list(run + 1, seq.Date(from = first(Date), length.out = length(Date), by = "days")-Date)]
    run <- run + 1
  }
  return(x[, Nruns := .N, by = Run_period])
  }
)
names(Run_periods) <- unique(Xzero$code_station)
Run_periods <- lapply(1:length(Run_periods), function(i){
  st <- Run_periods[[i]]; st$code_station <- names(Run_periods)[i]; return(st)})
Run_periods <- do.call(rbind, Run_periods)
save(Run_periods, file = "Runperiods_ValuesCheck")

load("Runperiods_ValuesCheck")
Run_periods[, Validation_steady := "OK"][Nruns > 4, Validation_steady := "Constant"]

Hydro_journaliere[, Validation_steady := "OK"]
Hydro_journaliere[Run_periods, on = c("code_station", "Date"), Validation_steady := i.Validation_steady]
Hydro_journaliere_corrected <- Hydro_journaliere[Validation_steady == "OK"][code_station != "J571211002"]
rm("Hydro_journaliere")

# Problematic stations J571211002

DebitMoyen <- unique(Hydro_journaliere_corrected [, Qmoyen := mean(resultat_obs_elab), by = code_station][,.(code_station, Qmoyen)])
Breaks <- seq(0, ceiling(max(DebitMoyen$Qmoyen)/10)*10, length.out = 15)
DebitMoyen[, Intervals := cut(Qmoyen, breaks = Breaks, include.lowest = TRUE, labels = formatC(Breaks[-1], format = "e", digits = 0))][, EffInt := .N, by = Intervals]
Breaks2 <- seq(0, ceiling(max(DebitMoyen[Qmoyen < 1.92e+05, Qmoyen])/10)*10, length.out = 12)
Debit2 <- DebitMoyen[Qmoyen < 1.92e+05][, 
            Intervals := cut(Qmoyen, breaks = Breaks2, include.lowest = TRUE, labels = formatC(Breaks2[-1], format = "e", digits = 0))]
DebitMoyen[Debit2, on = "code_station", Intervals := i.Intervals][, EffInt := .N, by = Intervals]

plot(DebitMoyen$Qmoyen)
ggplot(data = unique(unique(DebitMoyen[,.(EffInt, Intervals)])), aes(x = Intervals, y = EffInt)) +
  geom_bar(stat= "identity")



#####

## Centered indices
######
# Amplitude measure, based on EflowStats functions
Amplitude <- function(x){
  
  x[, Deseas := (resultat_obs_elab - mean(resultat_obs_elab)), by = "month"][
    , Deseas := (Deseas - mean(Deseas)) / sd(Deseas)][
    , corrAR1 := round(ar(Deseas, aic = FALSE, order.max = 1, method = "yule-walker")$ar, 2)]
  
  jday <- lubridate::yday(x$date_obs_elab)
  dec.year <- as.numeric(x$year) + (jday/365.25)
  
  scaled <- (x$resultat_obs_elab - mean(x$resultat_obs_elab)) / sd(x$resultat_obs_elab)
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
Index_timestep <- function(Tstep, Hydro_Serie, station, tol_threshold){
  print(Tstep)
   if(Tstep != 0){
  PrelDate <- data.table(Samp_date = Correspondance_station[ID_AMOBIO_START == station, Date_PrelBio])
  PrelDate[, period := 1:nrow(PrelDate)][, Start := Samp_date %m-% months(Tstep)]
  PrelDate <- left_join(PrelDate, PrelDate[, list(Date = seq.Date(from = Start, to = Samp_date, by = 'day')), by = "period"], by = "period")

Hydro_serie <- inner_join(Hydro_Serie, PrelDate, by = "Date")
if(nrow(Hydro_serie) == 0) {return(NA)}
Hydro_serie[, Tcover:= uniqueN(Date)/length(seq.Date(from = unique(Start), to = unique(Samp_date), by = 'day')), by = "period"]
Hydro_serie <- Hydro_serie[Tcover >= 0.75][,Tcover := NULL]
Hydro_serie <- Hydro_serie[, c("count", "unique") := list(.N, uniqueN(resultat_obs_elab)), by = "period"][
  count >= 10 & unique > 1,][, c("count", "unique") := NULL]
if(nrow(Hydro_serie) == 0) {return(NA)}
   } else {
     Hydro_serie <- Hydro_Serie[, c("count", "unique") := list(.N, uniqueN(resultat_obs_elab))][
     count >= 10 & unique > 1,][, c("count", "unique") := NULL]
     if(nrow(Hydro_serie) == 0) {return(NA)}
    Hydro_serie = Hydro_serie[, period := 0]
    PrelDate <- data.table(Samp_date = NA_Date_, period = 0, Start = NA_Date_)
    Tstep <- "all"}


Hydro_serie[, c("Lmean", "Lscale", "Lskew", "Lkurt") := as.list(lmoms(resultat_obs_elab, nmom = 4)$ratios[1:4]), by = "period"][
  , Lmean := mean(resultat_obs_elab), by = "period"]

Hydro_serie[, Deseas := (resultat_obs_elab - mean(resultat_obs_elab)), by = c("period", "month")][
  , Deseas := (Deseas - mean(Deseas)) / sd(Deseas), by = "period"][, corrAR1 := NA_real_]

Hydro_serie[, corrAR1 := round(ar(Deseas, aic = FALSE, order.max = 1, method = "yule-walker")$ar, 2), by = "period"]

if(uniqueN(Hydro_serie[,code_site])>1){print(paste("problem more sites", unique(Hydro_serie[,code_site])))}

Ampli <- do.call(rbind,lapply(unique(Hydro_serie$period), function(ti){
  return(Amplitude(Hydro_serie[period == ti,]))}))

Hydro_serie <- left_join(unique(Hydro_serie[,.(code_site, period, Lmean, Lscale, Lskew, Lkurt, corrAR1)]),
               Ampli, by = "period")

Hydro_serie <- right_join(unique(PrelDate[,.(Samp_date, period, Start)]), Hydro_serie,
               by = "period")[
                 ,.(code_site, Samp_date, period, Lmean, Lscale, Lskew, Lkurt, corrAR1, Amplitude, Phase)]

setnames(Hydro_serie, setdiff(colnames(Hydro_serie), c("code_site", "Samp_date")), paste0(setdiff(colnames(Hydro_serie), c("code_site", "Samp_date")), "_", Tstep))

return(Hydro_serie)

}


rm(list=setdiff(ls(), c("Correspondance_station", "Hydro_journaliere", "stationBio")))


Tol_threshold <- 0.75
severalDates <-vector()
Hydrolaps <- lapply(unique(Correspondance_station[, ID_AMOBIO_START]), function(stationBio){

if(uniqueN(Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET]) == 1){
  Xhydro <- Hydro_journaliere[code_station %in% unique(Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET])][
  , Date := as.Date(date_obs_elab)]
} else {
    Xhydro <- Hydro_journaliere[code_station %in% unique(Correspondance_station[ID_AMOBIO_START ==  stationBio & to_keep == TRUE, ID_TARGET])][
        , Date := as.Date(date_obs_elab)]
    #Xhydro <- Xhydro[code_station == Xhydro[Ndates == max(Ndates), first(code_station)]]
    complete <- Hydro_journaliere[code_station %in% unique(Correspondance_station[ID_AMOBIO_START ==  stationBio & to_keep == FALSE, 
                                                                                  ID_TARGET])]
    Xhydro <- rbind(Xhydro, complete[!Date %in% Xhydro$Date])[order(Date)][, Ndates := .N , by = Date]
    if(any(Xhydro$Ndates>1)){severalDates <- c(severalDates,stationBio)}
    Xhydro[,Ndates := NULL]
}
  
Hydro_laps <- tryCatch({lapply(c(3,6,12,0), function(tstep){
  return(Index_timestep(Tstep = tstep, Hydro_Serie = Xhydro, station = stationBio, tol_threshold = Tol_threshold))})},
  error = function(e) print(paste("Error with", stationBio, e)))

#if(!any(unlist(lapply(Hydro_laps,is.data.frame)))){print(stationBio)}

Hydro_all <- Hydro_laps[[3]]

ifelse(any(unlist(lapply(Hydro_laps[1:3],is.data.frame))), 
       Hydro_laps <- Reduce(function(...) left_join(..., by = c("code_site", "Samp_date")), 
                            Hydro_laps[1:3][which(unlist(lapply(Hydro_laps[1:3],is.data.frame)))]),
       Hydro_laps <- NA)


return(list(Hydro_laps, Hydro_all))

})
names(Hydrolaps) <- unique(Correspondance_station[, ID_AMOBIO_START])

save(Hydrolaps, file = "HydroIndex_3612all")
#####
