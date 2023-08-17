#####
invisible(lapply(c("devtools"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))
remotes::install_github("doi-usgs/EflowStats@v5.1.1")
#####
invisible(lapply(c("EflowStats", "data.table", "lmomco", "dplyr", "lubridate", "ggpubr", "ggplot2"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Hydro")

# mean, coefficient of variation, skewness, kurtosis, autoregressive lag-one (AR(1)) correlation coefficient, amplitude, phase of the seasonal signal
# use l-moments for the 4 first indicators
# deseasonalize for AR(1), substract the monthly mean stremfolw to each month value.
# Seasonality is measured using formula with sin, cos and tan
load("Hydro_journaliere")
load("Hydro_BioInventories_correspondance")


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
Index_timestep <- function(Tstep, Hydro_Serie, Station, Sampling_Date, Tol_threshold){
  
  if (Tstep > 12){Tstep <- time_length(years(Tstep/10), unit = "months")}
  
   if(Tstep != 0){
     Sampling_Date <- data.table(Samp_date = Sampling_Date)
     Sampling_Date[, period := 1:nrow(Sampling_Date)][, Start := Samp_date %m-% months(Tstep)]
   tryCatch({ Sampling_Date <-
    left_join(Sampling_Date, Sampling_Date[, list(Date = seq.Date(from = Start, to = Samp_date, by = 'day')), by = "period"], by = "period")},
       warning = function(w){ print(paste(w, "\n", Station))}                
)
Hydro_serie <- merge(Hydro_Serie, Sampling_Date, by = "Date")

if(nrow(Hydro_serie) == 0) {return(NA)}

Hydro_serie[, Tcover:= uniqueN(Date)/length(seq.Date(from = unique(Start), to = unique(Samp_date), by = 'day')), by = "period"]
Hydro_serie <- Hydro_serie[Tcover >= Tol_threshold][,Tcover := NULL]
Hydro_serie <- Hydro_serie[, c("count", "unique") := list(.N, uniqueN(resultat_obs_elab)), by = "period"][
  count >= 10 & unique > 1,][, c("count", "unique") := NULL]

if(nrow(Hydro_serie) == 0) {return(NA)}

   } else {
     Hydro_serie <- Hydro_Serie[, c("count", "unique") := list(.N, uniqueN(resultat_obs_elab))][
     count >= 10 & unique > 1,][, c("count", "unique") := NULL]
     if(nrow(Hydro_serie) == 0) {return(NA)}
    Hydro_serie = Hydro_serie[, c("period", "Samp_date", "Start") := list(0, NA_Date_, NA_Date_)]}

Hydro_serie[, c("Lmean", "Lscale", "Lskew", "Lkurt") := as.list(lmoms(resultat_obs_elab, nmom = 4)$ratios[1:4]), by = "period"][
  , Lmean := mean(resultat_obs_elab), by = "period"]

Hydro_serie[, Deseas := (resultat_obs_elab - mean(resultat_obs_elab)), by = c("period", "month")][
  , Deseas := (Deseas - mean(Deseas)) / sd(Deseas), by = "period"][, corrAR1 := NA_real_]

Hydro_serie[, corrAR1 := round(ar(Deseas, aic = FALSE, order.max = 1, method = "yule-walker")$ar, 2), by = "period"]

Ampli <- do.call(rbind,lapply(unique(Hydro_serie$period), function(ti){
  return(Amplitude(Hydro_serie[period == ti,]))}))

Hydro_serie <- merge(Hydro_serie, Ampli, by = "period")

if(Tstep == 0){
  Tstep <- "all"
  Hydro_serie[, crue := resultat_obs_elab >= quantile(Hydro_serie$resultat_obs_elab, 0.9)*3][, Ncrue := length(which(crue))][
    , crue := NULL]
}

Hydro_serie <- unique(Hydro_serie[,intersect(colnames(Hydro_serie), 
               c("code_site", "Samp_date", "period","Lmean","Lscale","Lskew", "Lkurt","corrAR1" ,"Amplitude","Phase","Ncrue")), with = F])

setnames(Hydro_serie, setdiff(colnames(Hydro_serie), c("code_site", "Samp_date")),
         paste0(setdiff(colnames(Hydro_serie), c("code_site", "Samp_date")), "_", Tstep))

return(Hydro_serie)

}

#rm(list=setdiff(ls(), c("Correspondance_station", "Hydro_journaliere", "stationBio")))

tol_threshold <- 0.75
laps <- c(3,6,12,50,0)

Hydrolaps <- lapply(unique(Correspondance_station[, ID_AMOBIO_START]), function(stationBio){
  
  sampling_date <- unique(Correspondance_station[ID_AMOBIO_START ==  stationBio & to_keep == TRUE, Date_PrelBio])
  
  if(any(setdiff(Correspondance_station[ID_AMOBIO_START ==  stationBio & to_keep == FALSE, Date_PrelBio], sampling_date))){
    print(paste("Missing dates in", stationBio))}

  Xhydro <- Hydro_journaliere[code_station %in% Correspondance_station[ID_AMOBIO_START ==  stationBio & to_keep == TRUE, ID_TARGET]][
    , Date := as.Date(date_obs_elab)][code_station == code_station[1]]
  
  while(length(setdiff(Hydro_journaliere[code_station %in% Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET], Date],
    Xhydro$Date))>0){
    Xhydro <- rbind(Xhydro, 
    Hydro_journaliere[code_station %in% Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET]][
      !Date %in% Xhydro$Date][code_station == code_station[1]])
  }
  
  Xhydro[, Ndates := .N , by = Date]
  if(any(Xhydro$Ndates > 1)){print(paste("Several dates in", stationBio))}
  Xhydro[,Ndates := NULL]
  
Hydro_laps <- tryCatch({lapply(laps, function(tstep){
  return(Index_timestep(Tstep = tstep, Hydro_Serie = Xhydro, Station = stationBio, Sampling_Date = sampling_date, 
                        Tol_threshold = tol_threshold))})},
  error = function(e) print(paste("Error with", stationBio, e)))

names(Hydro_laps) <- laps

#if(!any(unlist(lapply(Hydro_laps,is.data.frame)))){print(stationBio)}

Hydro_all <- Hydro_laps[["0"]]

ifelse(any(unlist(lapply(Hydro_laps[as.character(head(laps, -1))],is.data.frame))), 
       tryCatch({
         Hydro_laps <- Reduce(function(...) left_join(..., by = c("code_site", "Samp_date")), 
        Hydro_laps[as.character(head(laps, -1))][which(unlist(lapply(Hydro_laps[as.character(head(laps, -1))],is.data.frame)))])
       },
                warning = function(w){ print(paste(stationBio, "\n", w))}),
       Hydro_laps <- NA)

if(is.data.table(Hydro_laps)) {Hydro_laps$ID_AMOBIO_START <- stationBio}
if(is.data.table(Hydro_all))  {Hydro_all$ID_AMOBIO_START <- stationBio}

return(list(Hydro_laps, Hydro_all))

})

Hydrolaps <- lapply(Hydrolaps, function(si){
  if(!is.data.table(si[[1]])){return(si[[2]])
  } else {
  return(merge(si[[1]], si[[2]][,c("code_site",grep("all", names(si[[2]]), value = T)), with = F], by = "code_site", all = T))}
 })

Hydrolaps <- Hydrolaps[which(unlist(lapply(Hydrolaps, is.data.table)))]
Hydrolaps <- rbindlist(Hydrolaps, fill = T)
setnames(Hydrolaps, setdiff(colnames(Hydrolaps), c("ID_AMOBIO_START", "code_site", "Samp_date")), 
         paste0("H_", setdiff(colnames(Hydrolaps), c("ID_AMOBIO_START", "code_site", "Samp_date"))))
Hydrolaps <- Hydrolaps[,c("ID_AMOBIO_START", "code_site", "Samp_date", grep("H_", colnames(Hydrolaps), value = T))]

save(Hydrolaps, file = "HydroIndex_36125all_AM_20230816")
#####

## Plot validité données
#####
Plot_comp <- function(Compartment, Threshold){
  Nstat <- sum(Pie_selection[Comp == Compartment & Equivalent != "No overlap", Eff_coverage])
  Nperc <- sum(Pie_selection[Comp == Compartment & Equivalent != "No overlap", Eff_coverage])/sum(Pie_selection[Comp == Compartment, Eff_coverage])
  assign(paste0("plot", substring(Compartment,1, 3), Threshold), 
         ggplot(Pie_selection[Comp == Compartment, .(Eff_coverage, Equivalent)],aes(x = "", y = Eff_coverage, fill = Equivalent)) +
         geom_col(color = "black") + coord_polar(theta = "y") +
         scale_fill_brewer(palette = "Spectral", direction=-1) +
         theme_void() + ggtitle(paste0("\n", substring(Compartment,1, 8))) +
         labs(subtitle = paste("N(Valid Stations) =", length(Hydrolaps_valid), "\nTime filter = ", 
                                 round(length(Hydrolaps_valid)/length(Hydrolaps)*100, 0), "%")) +
         theme(legend.position = "none", plot.subtitle = element_text(size = 8)) ,
         envir = parent.frame())
}

for(Tol_threshold in c(75, 95, 100)){
  load(paste0("HydroIndex_3612_", Tol_threshold, "thresh"))
  
  Hydrolaps_valid <- Hydrolaps[which(unlist(lapply(Hydrolaps, is.data.frame)))]
  Hydrolaps_miss <- Hydrolaps[which(!unlist(lapply(Hydrolaps, is.data.frame)))]
  
  Pie_selection <- rbind(data.table(Stations = names(Hydrolaps_valid),
                                    Coverage = unlist(lapply(Hydrolaps_valid, ncol))),
                         data.table(Stations = names(Hydrolaps_miss), 
                                    Coverage = 0))
  Cover <- data.frame(Equivalent = c("All periods", "Two periods", "Single period", "No overlap"),
                      Coverage = c(26, 18, 10, 0))
  Pie_selection[, Comp := gsub("\\_.*", "", Stations)][, Eff_coverage := .N, by = c("Coverage", "Comp")][, Eff_coverage_tot := .N, by = "Coverage"]
  Pie_selection <- unique(left_join(Pie_selection, Cover, by = "Coverage")[, .(Comp, Eff_coverage, Eff_coverage_tot, Equivalent)])
  Pie_selection[, Equivalent := factor(Equivalent, levels = c("All periods", "Two periods", "Single period", "No overlap"))]
 
  assign(paste0("Pgen", Tol_threshold), 
         ggplot(unique(Pie_selection[, .(Eff_coverage_tot, Equivalent)]),aes(x = "", y = Eff_coverage_tot, fill = Equivalent)) +
         geom_col(color = "black") + coord_polar(theta = "y") +
         scale_fill_brewer(palette = "Spectral", direction=-1) + theme_void() + guides(fill = guide_legend(title="")) + 
         ggtitle(ifelse(Tol_threshold ==100, paste(1, "Coverage threshold\nAll"), 
                        paste(paste0("0.", Tol_threshold), "Coverage threshold\nAll"))) +
         labs(subtitle = paste("N(Valid Stations) =", length(Hydrolaps_valid), "\nTime filter = ", 
                              round(length(Hydrolaps_valid)/length(Hydrolaps)*100, 0), "%")) +
         theme(plot.subtitle = element_text(size = 8))
    )
  
  Plot_comp("DIATOM", Tol_threshold); Plot_comp("MACROINVERTEBRATE", Tol_threshold); Plot_comp("FISH", Tol_threshold)
  
}

ggarrange(plotlist = c(sapply(grep("100",ls(), value = T), get), sapply(grep("95",ls(), value = T), get),
                       sapply(grep("75",ls(), value = T), get)), ncol = 4, nrow = 3, widths = c(2,1,1,1))
#####

## Indices globaux
#####

Hydro_journaliere <- Hydro_journaliere[,c("count", "unique") := list(.N,uniqueN(resultat_obs_elab)), by = "code_site"][
  count>=4 & unique > 1,][,c("count", "unique") := NULL]
Hydro_journaliere[, c("Year_Lmean", "Year_Lscale", "Year_Lskew", "Year_Lkurt") := as.list(lmoms(resultat_obs_elab, nmom = 4)$ratios[1:4]), by = "code_site"][
  , Year_Lmean := mean(resultat_obs_elab), by = "code_site"]

Hydro_journaliere[, countM := uniqueN(resultat_obs_elab), by = c("code_site", "month")][, corrAR1 := NA_real_]
Hydro_journaliere[, Deseas := (resultat_obs_elab - mean(resultat_obs_elab)), by = c("code_site","month")][
  , Deseas := (Deseas - mean(Deseas)) / sd(Deseas), by = "code_site"]

Hydro_journaliere[countM != 1,
                  Year_corrAR1 := round(ar(Deseas, aic = FALSE, order.max = 1, method = "yule-walker")$ar, 2), by = "code_site"]

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
  
  amplitude <- round(sqrt((a^2) + (b^2)), digits = 2)
  
  phase <- ifelse(a > 0, 365.25 * ((pi/2) - atan(b/a))/(2 * pi), 365.25 + 365.25 * ((pi/2) - pi - atan(b/a))/(2 * pi))
  if(a == 0 & b > 0){ phase <- 365.25 }
  if(a == 0 & b < 0) { phase <- 365.25/2 }
  if(a == 0 & b == 0) { phase <- NA }
  
  return(cbind(code_site = unique(x$code_site), amplitude, phase))
}

Ampli <- as.data.table(do.call(rbind,lapply(unique(Hydro_journaliere[countM != 1,code_site]), function(si){
  return(Amplitude(Hydro_journaliere[code_site == si,]))
})))
Hydro_journaliere <- tryCatch(left_join(Hydro_journaliere, Ampli, by = "code_site"),
                              warning(w){ print(paste("warnings with", stationBio))})

Hydro_journaliere[, c("Year_MinJ", "Year_MaxJ") := list(min(resultat_obs_elab), max(resultat_obs_elab)), by = c("code_site", "year")]
Hydro_journaliere[,  c("Month_MinJ", "Month_MaxJ") := list(min(resultat_obs_elab), max(resultat_obs_elab)), by = c("code_site", "year", "month")]
Hydro_journaliere[, MeanMonth := mean(resultat_obs_elab), by = c("code_site", "month")][, 
                                                                                        c("Year_MinM", "Year_MaxM") := list(min(MeanMonth), max(MeanMonth)), by = c("code_site", "year")][, MeanMonth := NULL]

save(Hydro_journaliere, file = "Hydro_journaliere")



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
