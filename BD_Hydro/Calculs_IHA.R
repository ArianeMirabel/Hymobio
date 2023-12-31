#####
invisible(lapply(c("devtools"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))
remotes::install_github("doi-usgs/EflowStats@v5.1.1")
#####
invisible(lapply(c("EflowStats", "data.table", "lmomco", "dplyr", "lubridate", "ggpubr", "ggplot2", "hydroTSM"),function(pk){
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

# Count the number of low water and floods events based on the QMNA qnd equivalent high water levels
CountPeriods <- function(Q, module, dates, type, threshold = 7){
  
  if(is.na(unique(module))){ return(NA) }
  
  ifelse(type == "etiage", 
                ret <- data.table(dates = dates[which(Q <= module)]), 
                ret <- data.table(dates = dates[which(Q  >= module)]))
  
  ret <- ret[order(dates)][, Run_period := 1][
    , Date_check := seq.Date(from = first(dates), length.out = length(dates), by = "days")-dates]
  
  while(any(ret[, Date_check] < -threshold )){ret[Date_check != 0 & Date_check > -threshold, dates := 
        seq.Date(from = first(dates), length.out = length(dates), by = "days")][Date_check != 0 & Date_check > -threshold, 
        Date_check := seq.Date(from = first(dates), length.out = length(dates), by = "days")-dates][
          Date_check < -threshold,  c("Run_period", "Date_check") := list(Run_period + 1, 
          seq.Date(from = first(dates), length.out = length(dates), by = "days")-dates)]
  }
  
  return(uniqueN(ret[,Run_period]))
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
  
  QM <- cbind(
    unique(Hydro_serie[, AnnualMax := max(resultat_obs_elab, na.rm = T), by = year][, Indexyear := (as.numeric(year) - min(as.numeric(year))) +1][
      ,.(AnnualMax, year, Indexyear)][, rank := dense_rank(AnnualMax)])[, Prob := (.N - rank +1)/(.N +1)][, ReturnPeriod := 1/Prob][
        order(ReturnPeriod)][, paste0("QMMax", c(2, 5, 10)) := lapply(c(2, 5, 10), function(D){AnnualMax[which(ReturnPeriod > D)[1]]})],
    
    unique(Hydro_serie[, AnnualMin := min(resultat_obs_elab, na.rm = T), by = year][, Indexyear := (as.numeric(year) - min(as.numeric(year))) + 1][
      ,.(AnnualMin, year, Indexyear)][, rank := dense_rank(-AnnualMin)])[, Prob := (.N - rank +1)/(.N +1)][, ReturnPeriod := 1/Prob][
        order(ReturnPeriod)][, paste0("QMNA", c(2, 5, 10)) := lapply(c(2, 5, 10), function(D){AnnualMin[which(ReturnPeriod > D)[1]]})]
  )
  
  Hydro_serie[, grep("QM", colnames(QM), value = T) := unique(QM[, grep("QM", colnames(QM), value = T), with = F])]
  
  Hydro_serie[, paste0("Netiage", c(2, 5, 10)) := lapply(.SD, function(Module) {return(
    CountPeriods(Q = resultat_obs_elab, module = Module, dates = Date, type = "etiage"))}), 
    .SDcols = grep("QMNA", colnames(Hydro_serie), value = T)]
  
  Hydro_serie[, paste0("Ncrues", c(2, 5, 10)) := lapply(.SD, function(Module) {return(
    CountPeriods(Q = resultat_obs_elab, module = Module, dates = Date, type = "crue"))}), 
    .SDcols = grep("QMMax", colnames(Hydro_serie), value = T)]
  
}

Hydro_serie <- unique(Hydro_serie[,intersect(colnames(Hydro_serie), c("Samp_date", "period","Lmean","Lscale","Lskew", "Lkurt",
                      "corrAR1" ,"Amplitude","Phase",grep("QM|Netiage|Ncrues", colnames(Hydro_serie), value = T))), with = F])

setnames(Hydro_serie, setdiff(colnames(Hydro_serie), c("Samp_date")),
         paste0("H_",setdiff(colnames(Hydro_serie), c("Samp_date")), "_", Tstep))

return(Hydro_serie)

}



#rm(list=setdiff(ls(), c("Correspondance_station", "Hydro_journaliere")))

tol_threshold <- 0.75
laps <- c(3,6,12,50,0)

Hydrolaps <- lapply(unique(Correspondance_station[, ID_AMOBIO_START]), function(stationBio){
  
  pb <- txtProgressBar(min = 0, max = uniqueN(Correspondance_station[, ID_AMOBIO_START]), style = 3)
  
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
  
  Xhydro[, c("code_site", "code_station") := NULL]
  
Hydro_laps <- tryCatch({lapply(laps, function(tstep){
  return(Index_timestep(Tstep = tstep, Hydro_Serie = Xhydro, Station = stationBio, Sampling_Date = sampling_date, 
                        Tol_threshold = tol_threshold))})},
  error = function(e) print(paste("Error with", stationBio, e)))

names(Hydro_laps) <- laps

Hydro_all <- Hydro_laps[["0"]]

ifelse(any(unlist(lapply(Hydro_laps[as.character(head(laps, -1))],is.data.frame))), 
       tryCatch({
         Hydro_laps <- Reduce(function(...) left_join(..., by = c("Samp_date")), 
        Hydro_laps[as.character(head(laps, -1))][which(unlist(lapply(Hydro_laps[as.character(head(laps, -1))],is.data.frame)))])
       },
                warning = function(w){ print(paste(stationBio, "\n", w))}),
       Hydro_laps <- NA)

if(is.data.table(Hydro_laps)) {Hydro_laps$ID_AMOBIO_START <- stationBio}
if(is.data.table(Hydro_all))  {Hydro_all$ID_AMOBIO_START <- stationBio}

return(list(Hydro_laps, Hydro_all))

setTxtProgressBar(pb, which(unique(Correspondance_station[, ID_AMOBIO_START])) == stationBio)

})

Hydrolaps <- lapply(Hydrolaps, function(si){
  if(!is.data.table(si[[1]])){return(si[[2]])
  } else {
  return(cbind(si[[1]], si[[2]][,setdiff(names(si[[2]]), names(si[[1]])), with = F]))}
 })

Hydrolaps <- Hydrolaps[which(unlist(lapply(Hydrolaps, is.data.table)))]
Hydrolaps <- rbindlist(Hydrolaps, fill = T)
Hydrolaps <- Hydrolaps[,c("ID_AMOBIO_START", "Samp_date", grep("H_", colnames(Hydrolaps), value = T)), with = F]

save(Hydrolaps, file = "HydroIndex_36125all_AM_20232911")




#####

## Plot validité données
#####
load("HydroIndex_36125all_AM_20232911")
laps <- c("3","6","12","60")
module <- c("2", "5", "10")

Plot_valid <- Hydrolaps[, Compartment := sub("_.*", "", ID_AMOBIO_START)][, Year := format(as.Date(Samp_date, format =  "%Y-%m-%d"), "%Y")][,
            paste0("Nstation_",laps) := lapply(laps, function(i){uniqueN(ID_AMOBIO_START)}), by = "Year"]
Plot_valid[, paste0("MeanLap_",laps) := lapply(laps, function(i){median(get(paste0("H_Lmean_",i)), na.rm = T)}), by = "Year"]
Plot_valid[, paste0("SdLap_",laps) := lapply(laps, function(i){sd(get(paste0("H_Lmean_",i)), na.rm = T)}), by = "Year"]
Plot_valid[, paste0("UpLap_",laps) := lapply(laps, function(i){quantile(get(paste0("H_Lmean_",i)), probs = 0.6, na.rm = T)}), by = "Year"]
Plot_valid[, paste0("LowLap_",laps) := lapply(laps, function(i){quantile(get(paste0("H_Lmean_",i)), probs = 0.4, na.rm = T)}), by = "Year"]

Plot_valid[, paste0("MeanEtiages_",module) := lapply(module, function(i){median(get(paste0("H_Netiage",i, "_all")), na.rm = T)}), by = "Year"]
Plot_valid[, paste0("SdEtiages_",module) := lapply(module, function(i){sd(get(paste0("H_Netiage",i, "_all")), na.rm = T)}), by = "Year"]
Plot_valid[, paste0("UpEtiages_",module) := lapply(module, function(i){quantile(get(paste0("H_Netiage",i, "_all")), probs = 0.6, na.rm = T)}), by = "Year"]
Plot_valid[, paste0("LowEtiages_",module) := lapply(module, function(i){quantile(get(paste0("H_Netiage",i, "_all")), probs = 0.4, na.rm = T)}), by = "Year"]

Plot_valid[, paste0("MeanCrues_",module) := lapply(module, function(i){median(get(paste0("H_Ncrues",i, "_all")), na.rm = T)}), by = "Year"]
Plot_valid[, paste0("SdCrues_",module) := lapply(module, function(i){sd(get(paste0("H_Ncrues",i, "_all")), na.rm = T)}), by = "Year"]
Plot_valid[, paste0("UpCrues_",module) := lapply(module, function(i){quantile(get(paste0("H_Ncrues",i, "_all")), probs = 0.6, na.rm = T)}), by = "Year"]
Plot_valid[, paste0("LowCrues_",module) := lapply(module, function(i){quantile(get(paste0("H_Ncrues",i, "_all")), probs = 0.4, na.rm = T)}), by = "Year"]


Plot_valid <- unique(Plot_valid[,c("Compartment","Year", grep("MeanLap|SdLap|UpLap|LowLap|Nstation|Crues|Etiages", names(Plot_valid), value = T)),
                                with = F])[!is.na(Year)]

Qlong <- merge(
  merge(
  merge(melt(Plot_valid[,c("Compartment", "Year", grep("MeanLap", names(Plot_valid), value = T)), with = F], 
             id.vars = c("Compartment", "Year"), variable.name = "Lap", value.name = "Mean")[, Lap := sub(".*_", "", Lap)],
        melt(Plot_valid[,c("Compartment", "Year", grep("UpLap", names(Plot_valid), value = T)), with = F], 
             id.vars = c("Compartment", "Year"), variable.name = "Lap", value.name = "Upper")[, Lap := sub(".*_", "", Lap)], 
        by = c("Compartment", "Year", "Lap")),
  melt(Plot_valid[,c("Compartment", "Year", grep("LowLap", names(Plot_valid), value = T)), with = F], 
       id.vars = c("Compartment", "Year"), variable.name = "Lap", value.name = "Lower")[, Lap := sub(".*_", "", Lap)], 
  by = c("Compartment", "Year", "Lap")),
  melt(Plot_valid[,c("Compartment", "Year", grep("Nstation", names(Plot_valid), value = T)), with = F], 
       id.vars = c("Compartment", "Year"), variable.name = "Lap", value.name = "Nstation")[, Lap := sub(".*_", "", Lap)], 
        by = c("Compartment", "Year", "Lap"))[!is.na(Mean)]
Qlong[, Lap := factor(Lap, levels = c("3","6","12","60","all"))]

Extlong <-   merge(
  merge(melt(Plot_valid[,c("Compartment", "Year", grep("MeanCrues", names(Plot_valid), value = T)), with = F], 
             id.vars = c("Compartment", "Year"), variable.name = "Module", value.name = "Mean")[, Module := sub(".*_", "", Module)],
        melt(Plot_valid[,c("Compartment", "Year", grep("UpCrues", names(Plot_valid), value = T)), with = F], 
             id.vars = c("Compartment", "Year"), variable.name = "Module", value.name = "Upper")[, Module := sub(".*_", "", Module)], 
        by = c("Compartment", "Year", "Module")),
  melt(Plot_valid[,c("Compartment", "Year", grep("LowCrues", names(Plot_valid), value = T)), with = F], 
       id.vars = c("Compartment", "Year"), variable.name = "Module", value.name = "Lower")[, Module := sub(".*_", "", Module)], 
  by = c("Compartment", "Year", "Module"))[, Type := "Crues"]

Extlong <- rbind(Extlong,
                   merge(
  merge(melt(Plot_valid[,c("Compartment", "Year", grep("MeanEtiages", names(Plot_valid), value = T)), with = F], 
             id.vars = c("Compartment", "Year"), variable.name = "Module", value.name = "Mean")[, Module := sub(".*_", "", Module)],
        melt(Plot_valid[,c("Compartment", "Year", grep("UpEtiages", names(Plot_valid), value = T)), with = F], 
             id.vars = c("Compartment", "Year"), variable.name = "Module", value.name = "Upper")[, Module := sub(".*_", "", Module)], 
        by = c("Compartment", "Year", "Module")),
  melt(Plot_valid[,c("Compartment", "Year", grep("LowEtiages", names(Plot_valid), value = T)), with = F], 
       id.vars = c("Compartment", "Year"), variable.name = "Module", value.name = "Lower")[, Module := sub(".*_", "", Module)], 
  by = c("Compartment", "Year", "Module"))[, Type := "Etiages"])
Extlong[, Module := factor(Module, levels = c("2","5","10"))]


pN <- ggplot(Qlong, aes(x = Year, y = Nstation, group = Lap, color = Lap), show.legend = F) +
  geom_line(show.legend =  F) + theme_minimal() + ylab("N\nStations") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, 600, by = 200)) + ggtitle("(a)") +
  scale_color_manual(name = "Lap", labels = labs, values =  vals)

labs <- c("3 months" , "6 months", "1 year","5 years");
vals <- c("firebrick" , "darkorange", "olivedrab3","royalblue")
pMean <- ggplot(Qlong) + 
  geom_line(aes(x = Year, y = Mean, group = Lap, color = Lap)) + 
  geom_ribbon(alpha = 0.2, aes(x = Year, y = Mean, ymin = Lower, ymax = Upper, group = Lap, fill = Lap)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Mean (l/s)") + scale_fill_manual(name = "Lap", labels = labs, values =  vals) +   
  scale_color_manual(name = "Lap", labels = labs, values =  vals) + ggtitle("(b)")

labs <- c("Limite 10 ans (décennale)" , "Limite 5 ans", "Limite 2 ans");
vals <- c( "lightskyblue", "cornflowerblue" ,"royalblue4")
pExtrem <- ggplot(Extlong) + 
  geom_line(aes(x = Year, y = Mean, group = Module, color = Module)) + 
  geom_ribbon(alpha = 0.2, aes(x = Year, y = Mean, ymin = Lower, ymax = Upper, group = Module, fill = Module)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Occurence") + scale_fill_manual(name = "Module", labels = labs, values =  vals) +   
  scale_color_manual(name = "Module", labels = labs, values =  vals) + ggtitle("(c)") + facet_wrap(~Type)
pExtrem 


ggarrange(pN, pMean, ncol = 1, common.legend = T, heights = c(1,3), 
          legend = "right")














# Old
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

Pie_selection <- rbind(data.table(Stations = names(Hydrolaps),
                                    Coverage = unlist(lapply(Hydrolaps, ncol))),
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
  

ggarrange(plotlist = c(sapply(grep("100",ls(), value = T), get), sapply(grep("95",ls(), value = T), get),
                       sapply(grep("75",ls(), value = T), get)), ncol = 4, nrow = 3, widths = c(2,1,1,1))
#####

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

#####
