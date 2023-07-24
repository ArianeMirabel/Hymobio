invisible(lapply(c("data.table", "doParallel", "ggplot2"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

## A partir d'une table "Hydro_journaliere" avec code_site/code_station/resultat_hydro/date

clust <- makeCluster(detectCores()-2)
registerDoParallel(clust)

Run_periods <- do.call(rbind,mclapply(unique(Hydro_journaliere$code_station), function(st){
  x <- Hydro_journaliere[code_station == st, .(Date)]
  x <- x[order(Date)]
  run <- 1
  x[, Run_period := run][, Date_check := seq.Date(from = first(Date), length.out = length(Date), by = "days")-Date]
  x[Date_check > -7, Date_check := 0]
  
  while(any(x[, Date_check] != 0)){
    x[Date_check != 0, c("Run_period", "Date_check") := list(run + 1, seq.Date(from = first(Date), length.out = length(Date), by = "days")-Date)]
    run <- run + 1
    x[Date_check > -7, Date_check := 0]
  }
  return(x[, code_station := st][,.(code_station, Date, Run_period)])}))

stopCluster(clust)
gc()

Hydro_journaliere[, Run_period := 1]
Hydro_journaliere[Run_periods, on = c("code_station", "Date"), Run_period := i.Run_period]

Hydro_journaliere[, Length_period := .N, by = c("code_station", "Run_period")][,Ndates := NULL]
Hydro_journaliere[, Cover := "All_periods"]
Hydro_journaliere[Length_period < 1826 & Length_period >= 365, Cover := "One_year"]
Hydro_journaliere[Length_period < 365 & Length_period >= 182, Cover := "Six_months"]
Hydro_journaliere[Length_period < 182 & Length_period >= 91, Cover := "Three_months"]
Hydro_journaliere[Length_period < 91 & Length_period >= 30, Cover := "One_months"]
Hydro_journaliere[Length_period < 30 , Cover := "Too_short"]

PiePeriod <- unique(Hydro_journaliere[, Eff_periods := .N, by = Cover][,.(Eff_periods, Cover)])

Period_pie <- ggplot(PiePeriod,aes(x = "", y = Eff_periods, fill = Cover)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("All_periods" = "olivedrab",  "One_year" = "palegreen3", "One_months" = "turquoise4", 
                               "Six_months" = "royalblue", "Three_months"= "rosybrown", "Too_short" = "darkred")) +
  theme_void() + guides(fill = guide_legend(title=""))  






