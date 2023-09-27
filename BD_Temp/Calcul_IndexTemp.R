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

#min = -0.5 max = 37.908

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
Index_timestep <- function(Tstep, Temp_Serie, Station, Sampling_date, Tol_threshold){
  
  if (Tstep > 12){Tstep <- time_length(years(Tstep/10), unit = "months")}
  
   if(Tstep != 0){
     Sampling_Date <- data.table(Samp_date = Sampling_date)
     Sampling_Date[, period := 1:nrow(Sampling_Date)][, Start := Samp_date %m-% months(Tstep)]
     tryCatch({ Sampling_Date <-
       left_join(Sampling_Date, Sampling_Date[, list(Date = seq.Date(from = Start, to = Samp_date, by = 'day')), by = "period"], by = "period")},
       warning = function(w){ print(paste(w, "\n", Station))}                
     )
        
Temp_serie <- merge(Temp_Serie, Sampling_Date, by = "Date")

if(nrow(Temp_serie) == 0) {return(NA)}

Temp_serie[, Tcover:= uniqueN(Date)/length(seq.Date(from = unique(Start), to = unique(Samp_date), by = 'day')), by = "period"]
Temp_serie <- Temp_serie[Tcover >= Tol_threshold][,Tcover := NULL]
Temp_serie <- Temp_serie[, c("count", "unique") := list(.N, uniqueN(Tw_corr)), by = "period"][
  count >= 10 & unique > 1,][, c("count", "unique") := NULL]

if(nrow(Temp_serie) == 0) {return(NA)}

   } else {
     Temp_serie <- Temp_Serie[, c("count", "unique") := list(.N, uniqueN(Tw_corr))][
     count >= 10 & unique > 1,][, c("count", "unique") := NULL]
     if(nrow(Temp_serie) == 0) {return(NA)}
    Temp_serie = Temp_serie[, c("period", "Samp_date", "Start") := list(0, NA_Date_, NA_Date_)]
   }
  
Temp_serie[, c("Lmean", "Lscale", "Lskew", "Lkurt") := as.list(lmoms(Tw_corr, nmom = 4)$ratios[1:4]), by = "period"][
  , Lmean := mean(Tw_corr), by = "period"]

Temp_serie[, Deseas := (Tw_corr - mean(Tw_corr)), by = c("period", "month")][
  , Deseas := (Deseas - mean(Deseas)) / sd(Deseas), by = "period"][, corrAR1 := NA_real_]

Temp_serie[, corrAR1 := round(ar(Deseas, aic = FALSE, order.max = 1, method = "yule-walker")$ar, 2), by = "period"]

Ampli <- do.call(rbind,lapply(unique(Temp_serie$period), function(ti){
  return(Amplitude(Temp_serie[period == ti,]))}))

Temp_serie <- merge(Temp_serie, Ampli, by = "period")

if(Tstep == 0){
  Tstep <- "all"
  Temp_serie[, extreme := Tw_corr >= quantile(Temp_serie$Tw_corr, 0.9)*3][, Nextreme := length(which(extreme))][
    , extreme := NULL]
}

Temp_serie <- unique(Temp_serie[,intersect(names(Temp_serie),
    c("stationID", "Samp_date", "period", "Lmean", "Lscale", "Lskew", "Lkurt", "corrAR1", "Amplitude", "Phase", "Nextreme")), with = F])

setnames(Temp_serie, setdiff(colnames(Temp_serie), c("stationID", "Samp_date")), 
         paste0("T_",setdiff(colnames(Temp_serie), c("stationID", "Samp_date")), "_", Tstep))

return(Temp_serie)

}

tol_threshold <- 0.75
laps <- c(3,6,12,50,0)

Templaps <- lapply(unique(Correspondance_station[, ID_AMOBIO_START]), function(stationBio){
  
  pb <- txtProgressBar(min = 0, max = uniqueN(Correspondance_station[, ID_AMOBIO_START]), style = 3)
  
  sampling_date <- unique(Correspondance_station[ID_AMOBIO_START ==  stationBio , Date_PrelBio])
  
  if(uniqueN(Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET]) == 1){
    Xtemp <- Temp_journaliere[stationID %in% unique(Correspondance_station[ID_AMOBIO_START ==  stationBio, ID_TARGET])][
      , Date := as.Date(Date)]
  } else {print(paste("Several Stations for", stationBio))}
  
  Xtemp[, Ndates := .N , by = Date]
  if(any(Xtemp$Ndates > 1)){print(paste("Several dates in", stationBio))}
  Xtemp[,Ndates := NULL]
  

Temp_laps <- tryCatch({lapply(laps, function(tstep){
  return(Index_timestep(Tstep = tstep, Temp_Serie = Xtemp, Station = stationBio, Sampling_date = sampling_date,
                        Tol_threshold = tol_threshold))})},
  error = function(e) print(paste("Error with", stationBio, e)))

names(Temp_laps) <- laps

Temp_all <- Temp_laps[["0"]]

ifelse(any(unlist(lapply(Temp_laps[as.character(head(laps, -1))],is.data.frame))), 
       tryCatch({
         Temp_laps <- Reduce(function(...) left_join(..., by = c("stationID", "Samp_date")), 
         Temp_laps[as.character(head(laps, -1))][which(unlist(lapply(Temp_laps[as.character(head(laps, -1))],is.data.frame)))])
       },
       warning = function(w){ print(paste(stationBio, "\n", w))}),
       Temp_laps <- NA)

if(is.data.table(Temp_laps)) {Temp_laps$ID_AMOBIO_START <- stationBio}
if(is.data.table(Temp_all))  {Temp_all$ID_AMOBIO_START <- stationBio}

setTxtProgressBar(pb, which(unique(Correspondance_station[, ID_AMOBIO_START]) == stationBio))

return(list(Temp_laps, Temp_all))

})

Templaps <- lapply(Templaps, function(si){
  if(!is.data.table(si[[1]])){return(si[[2]])
  } else {
    return(merge(si[[1]], si[[2]][,c("stationID",grep("all", names(si[[2]]), value = T)), with = F], by = "stationID", all = T))}
})

Templaps <- Templaps[which(unlist(lapply(Templaps, is.data.table)))]
Templaps <- rbindlist(Templaps, fill = T)
Templaps <- Templaps[,c("ID_AMOBIO_START", "stationID", "Samp_date", grep("T_", colnames(Templaps), value = T)), with = F]
setnames(Templaps, "stationID", "CdStation")


save(Templaps, file = "TempIndex_36125all_AM_20230821")
#####

##Plot validation
#####
load("TempIndex_36125all_AM_20230821")

Stat_Turn <- function(yr){
  uniqueN(setdiff(Turnover[Year == as.numeric(yr), stationID],Turnover[Year == as.numeric(yr)-1, stationID])) +
    uniqueN(setdiff(Turnover[Year == as.numeric(yr)-1, stationID],Turnover[Year == as.numeric(yr), stationID])) 
}

laps <- c("3","6","12","60")
Plot_valid <- Templaps[stationID != "4044700"][, Compartment := sub("_.*", "", ID_AMOBIO_START)][, Year := format(as.Date(Samp_date, format =  "%Y-%m-%d"), "%Y")][,
             paste0("Nstation_",laps) := lapply(laps, function(i){uniqueN(stationID)}), by = "Year"]

Turnover <- unique(Plot_valid[Year > 2008,.(stationID, Year)])
Turnover <- unique(Turnover[, StatTurn := Stat_Turn(Year), by = "Year"][,.(Year, StatTurn)])[order(Year), Year := as.numeric(Year)]

Plot_valid[, paste0("MeanLap_",laps) := lapply(laps, function(i){median(get(paste0("T_Lmean_",i)), na.rm = T)}), by = "Year"]
Plot_valid[, paste0("UpLap_",laps) := lapply(laps, function(i){quantile(get(paste0("T_Lmean_",i)), probs = 0.6, na.rm = T)}), by = "Year"]
Plot_valid[, paste0("LowLap_",laps) := lapply(laps, function(i){quantile(get(paste0("T_Lmean_",i)), probs = 0.4, na.rm = T)}), by = "Year"]

Plot_valid <- unique(Plot_valid[,c("Compartment","Year", grep("MeanLap|UpLap|LowLap|Nstation", names(Plot_valid), value = T)),
                                with = F])[!is.na(Year)]

long <- merge(
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
  by = c("Compartment", "Year", "Lap"))[!is.na(Mean)][Year > 2008]
long[, Lap := factor(Lap, levels = c("3","6","12","60","all"))]

labs <- c("3 months" , "6 months", "1 year","5 years");
vals <- c("firebrick" , "darkorange", "olivedrab3","royalblue")
pMean <- ggplot(long) + 
  geom_line(aes(x = Year, y = Mean, group = Lap, color = Lap)) + 
  geom_ribbon(alpha = 0.2, aes(x = Year, y = Mean, ymin = Lower, ymax = Upper, group = Lap, fill = Lap)) + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Mean (°C)") + scale_fill_manual(name = "Lap", labels = labs, values =  vals) +   
  scale_color_manual(name = "Lap", labels = labs, values =  vals) + ggtitle("(b)")

pN <- ggplot() + theme_minimal() + ylab("N\nStations") +
  geom_col(data = Turnover[Year > 2009], aes(x = Year, y = StatTurn, fill = "Turnover")) +
  geom_line(data = long, aes(x = Year, y = Nstation, group = Lap, color = Lap), size = 2, show.legend = F) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_color_manual(name = "Lap", labels = labs, values =  vals) + ggtitle("(a)") +
  scale_fill_manual(name = "Turnover", label="N new station", values = c("Turnover" = "darkgrey")) 

ggarrange(pN, pMean, ncol = 1, common.legend = F, heights = c(1,3), 
          legend = "right")


barplot <- Turnover[, Nyrs := uniqueN(Year), by = stationID][order(Nyrs),.(Nyrs)]
ggplot(barplot, aes(Nyrs)) + geom_bar() + theme_minimal()

