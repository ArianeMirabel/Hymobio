invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

#####
load("AMOBIO_WP3_2_FISH_INV_DIA_20230817.Rdata")
Carhyce_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN); rm("MATRIX_AMOBIO_WP3_CLEAN")

metrics <- data.table(Name = colnames(Carhyce_all))[grep("NC_|POE_|P_|H_|T_|HM_",substr(Name, 1,5))][, NameOrigin := Name][
  , Category := sub("_.*", "",Name)][, Nscales := str_count(Name, "_")][, Scale1 := sub(".*\\_", "", Name)][
  Scale1 == "ALL", Scale1 := 0][#!Scale1 %in% c("NF","WD")][
  Scale1 == "F", Name := sub("_[^_]*$", "", Name)][Scale1 == "F", Scale1 := sub(".*\\_", "", Name)][
  , Numeric1 := as.numeric(gsub("\\D", "", Scale1))][, Name := sub("_[^_]*$", "", Name)][#!is.na(Numeric1)][
  , Nscales := str_count(Name, "_")][Nscales > 1, Scale2 := sub(".*\\_", "", Name)][, Numeric2 := as.numeric(gsub("\\D", "", Scale2))][
  !is.na(Numeric2), Name := sub('_[^_]*$', '', Name)][!is.na(Scale2) & is.na(Numeric2), Scale2 := "remove"][#Scale2 != "remove"][
  ,c("Min1", "Min2") := list(min(Numeric1), min(Numeric2)), by = Name][, Name := sub("^.*?\\_", "", Name)]
metrics[is.na(metrics), ] <- 0 
write.table(Finale, file = "FinalMetrics.csv", sep=";", row.names = F)

#####

# Plot variables density distribution
#####
Catalog <- fread("MetricsCatalogue.csv")
Catalog <- Catalog[which(Tokeep)]

load("AMOBIO_WP3_2_FISH_INV_DIA_20230817.Rdata")
Carhyce_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN); rm("MATRIX_AMOBIO_WP3_CLEAN")
Carhyce_all <- Carhyce_all[, Catalog$NameOrigin, with = F]

Outliers <- 5

Categories <- data.frame(Abbrev = unique(Catalog$Category),
                       category = c("Hydrology", "Hydromorphology", "Natural Control",
                                    "Pressure", "Water Flow Obstruction", "Temperature"))

lapply(1:nrow(Categories), function(i){
  
  CatPlot <- lapply(Catalog[Category == Categories[i,"Abbrev"], NameOrigin], function(Metric){
    
  Data <- Carhyce_all[, Metric, with = F]
colnames(Data) <- "metric"

if(!is.factor(Data$metric)){
  
  pAll <- ggplot(data = Data[!is.na(metric)], aes(metric)) + #ggtitle(label = "All data") +
  geom_density(color = "firebrick", aes(y= after_stat(ndensity))) + 
  theme_classic() + ylab("") + xlab("") + scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + 
  theme(panel.background = element_rect(fill = "lightgrey"), plot.margin = unit(c(0.5, 0.5, -0.5, 0.5),"lines"))

DataNoOut <- Data[metric >= quantile(metric, Outliers/200, na.rm = T) & metric <= quantile(metric, 1-Outliers/200, na.rm = T)]

if(uniqueN(DataNoOut$metric) > 1){
  pNoOut <- ggplot(data = DataNoOut, aes(metric)) +
  geom_density(color = "firebrick", aes(y= after_stat(ndensity))) + theme_classic() + ylab("")  + 
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + #ggtitle(label = "Without\nOutliers") + 
  xlab(paste0(Catalog[NameOrigin == Metric, Name], " (",Catalog[NameOrigin == Metric, Unit],") ")) + 
  theme(plot.margin = unit(c(-0.5, 0.5, 1.5, 0.5), "lines"))
} else { pNoOut <- pAll + xlab(paste0(Catalog[NameOrigin == Metric, Name], " (",Catalog[NameOrigin == Metric, Unit],") ")) + 
                   theme(plot.margin = unit(c(-0.5, 0.5, 1.5, 0.5), "lines"))
         pAll <- NULL }

} else {
  pAll <- NULL
  pNoOut <- ggplot(data = Data[!is.na(metric)], aes(metric)) + #ggtitle(label = "All data") +
    geom_bar(fill = "firebrick", alpha = 0.65, aes(y = ..prop.., group = 1)) + 
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + theme_classic() + ylab("")  + 
    xlab(paste0(Catalog[NameOrigin == Metric, Name], " (",Catalog[NameOrigin == Metric, Unit],") ")) + 
    theme(plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "lines"))
}

if(which(Catalog[Category == Categories[i,"Abbrev"], NameOrigin] == Metric) %in% 
   seq(1, length(Catalog[Category == Categories[i,"Abbrev"], NameOrigin]), by = 3)){
  pNoOut <- pNoOut + ylab("Density")
}

return(ggarrange(pAll, pNoOut, ncol = 1, heights = c(2,5)))
})

SplitPlotCategory <- split(1:length(CatPlot), ceiling(1:length(CatPlot)/9))

lapply(1:length(SplitPlotCategory), function(j){
  PlotCategory <- ggarrange(plotlist = CatPlot[SplitPlotCategory[[j]]], ncol = 3,
                            nrow = ceiling(length(SplitPlotCategory[[j]])/3))
  
  png(paste0(sub(" ", "",Categories[i,"category"]), "_", j,".png"))
  print(annotate_figure(annotate_figure(PlotCategory,
                                        top=text_grob(paste0(j,"/",length(SplitPlotCategory)), size = 10),),
                        top=text_grob(Categories[i,"category"], size = 10, face = "bold")))
  dev.off()
})
})
#####

# Set thresholds
#####
Catalog <- fread("MetricsCatalogue.csv")
Catalog <- Catalog[which(Tokeep)]

load("../SELECTION_STATION_list_station_filter5_metadata_20230525.Rdata")
SelectedStations <- as.data.table(list_station_filter5_clean)[COMPARTIMENT_TARGET %in% c("DIATOM","FISH","MACROINVERTEBRATE")]
rm(list_station_filter5_clean)

load("AMOBIO_WP3_2_FISH_INV_DIA_20230817.Rdata")
Carhyce_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN)
Carhyce_selected <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN)[node_id %in% SelectedStations$node_id_START]; rm("MATRIX_AMOBIO_WP3_CLEAN")

Outliers <- 20
KeepOutliers <- FALSE

Categories <- data.frame(Abbrev = unique(Catalog$Category),
                         category = c("Hydrology", "Hydromorphology", "Natural Control",
                                      "Pressure", "Water Flow Obstruction", "Temperature"))
Ncut <- 9

lapply(1:nrow(Categories), function(i){
  
  CatPlot <- lapply(Catalog[Category == Categories[i,"Abbrev"], NameOrigin], function(Metric){
    
    Data <- rbind(Carhyce_all[, Metric, with = F][, Source := "All"],
                  Carhyce_selected[, Metric, with = F][, Source:= "Sel"])
    setnames(Data, Metric,"metric")
    
    if(!is.factor(Data$metric)){
      
    if(!KeepOutliers){Data <- Data[!is.na(metric) & 
                                     metric >= quantile(metric, Outliers/200, na.rm = T) &
                                     metric <= quantile(metric, 1-Outliers/200, na.rm = T)]}
     ifelse(nrow(Data)>1000, Binseq <- seq(0,1,by = 0.1), Binseq <- seq(0,1,by = 0.2))
     thresh <- unique(quantile(Data$metric, probs = Binseq))
     Data[, BinUp := cut(metric, breaks = thresh, labels = thresh[-length(thresh)], include.lowest = T)]
     
     Thresholds <- Data[,.(BinUp)][,Type := "Biology"][, BinUp := as.numeric(as.character(BinUp))]
     
    if(!is.na(Catalog[NameOrigin == Metric, LittThreshold])){
        Thresholds[BinUp == as.numeric(levels(Data$BinUp))[
          which.min(abs(as.numeric(levels(Data$BinUp)) - as.numeric(Catalog[NameOrigin == Metric, LittThreshold])))],
             c("BinUp", "Type") := list(as.numeric(Catalog[NameOrigin == Metric, LittThreshold]), "HM")]}
        
     pComp <- ggplot() + 
          stat_density(data = Data[!is.na(metric)], aes(x= metric, y = after_stat(ndensity), colour = Source),
                       linewidth = 1,  geom="line", position="identity") + 
          theme_classic() + ylab("") + scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + 
          scale_x_continuous(breaks = round(unique(Thresholds$BinUp), 2)) +
          xlab(paste(Catalog[NameOrigin == Metric, Name], Catalog[NameOrigin == Metric, Scale1],
                     "(",Catalog[NameOrigin == Metric, Unit],")")) + 
          theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"lines"), axis.text.x = element_text(angle = 70, hjust = 1, size = 8)) +
          geom_vline(data = Thresholds, aes(xintercept = BinUp, linetype = Type)) +
          scale_color_manual(name = "Stations density", breaks = c("All", "Sel"),
                             values = c("All" = "firebrick", "Sel" = "olivedrab"),
                             labels = c("All Sites", "Bio-attached sites")) +
        scale_linetype_manual(name = "Upper limit\nthresholds", breaks = c("Biology", "HM"), values = c("Biology" = 1, "HM" = 2))

      } else {
      pComp <- ggplot(data = Data[!is.na(metric)], aes(x = metric, fill = Source)) + 
        geom_bar(aes(y = ..count../sum(..count..), group = Source)) + 
        scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + theme_classic() + ylab("")  + 
        xlab(paste0(Catalog[NameOrigin == Metric, Name], " (",Catalog[NameOrigin == Metric, Unit],") ")) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))+
        scale_fill_manual(name = "", values = c("All" = "firebrick", "Sel" = "olivedrab"),labels = c("All Sites", "Bio-attached sites"))
     }
    
    if(which(Catalog[Category == Categories[i,"Abbrev"], NameOrigin] == Metric) %in% 
       seq(1, length(Catalog[Category == Categories[i,"Abbrev"], NameOrigin]), by = 3) && !is.null(pComp)){
      pComp <- pComp + ylab("Density")}
    
    return(pComp)
  })
  
  SplitPlotCategory <- split(1:length(CatPlot), ceiling(1:length(CatPlot)/9))
  
  lapply(1:length(SplitPlotCategory), function(j){
    PlotCategory <- ggarrange(plotlist = CatPlot[SplitPlotCategory[[j]]], ncol = 3,
                              nrow = ceiling(length(SplitPlotCategory[[j]])/3), common.legend = T)
    
    png(paste0(sub(" ", "",Categories[i,"category"]), "SelectComp_", j,".png"))
    print(annotate_figure(annotate_figure(PlotCategory,
                                          top=text_grob(paste0(j,"/",length(SplitPlotCategory)), size = 10),),
                          top=text_grob(Categories[i,"category"], size = 10, face = "bold")))
    dev.off()
  })
})



lapply(1:nrow(Categories), function(i){
  
  CatPlot <- lapply(Catalog[Category == Categories[i,"Abbrev"], NameOrigin], function(Metric){
    
    Data <- rbind(Carhyce_all[, Metric, with = F][, Source := "All"],
                  Carhyce_selected[, Metric, with = F][, Source:= "Sel"])
    setnames(Data, Metric,"metric")
    
    if(!is.factor(Data$metric)){
      
      if(!KeepOutliers){Data <- Data[metric >= quantile(metric, Outliers/200, na.rm = T) & metric <= quantile(metric, 1-Outliers/200, na.rm = T)]}
    
     if(uniqueN(Data[!is.na(metric),metric]) > 1){
       
       pComp <- ggplot(data = Data[!is.na(metric)], aes(x= metric, y = after_stat(ndensity), colour = Source)) + #ggtitle(label = "All data") +
         stat_density(linewidth = 1,  geom="line", position="identity") + 
         theme_classic() + ylab("") + scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + 
         xlab(paste0(Catalog[NameOrigin == Metric, Name], " (",Catalog[NameOrigin == Metric, Unit],") ")) + 
         theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"lines"))  +
         scale_color_manual(name = "", values = c("All" = "firebrick", "Sel" = "olivedrab"),labels = c("All Sites", "Bio-attached sites"))
      
       } else {pComp <- NULL}
      
    } else {
      pComp <- ggplot(data = Data[!is.na(metric)], aes(x = metric, fill = Source)) + 
        geom_bar(aes(y = after_stat(count/sum(count)), group = Source)) + 
        scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + theme_classic() + ylab("")  + 
        xlab(paste0(Catalog[NameOrigin == Metric, Name], " (",Catalog[NameOrigin == Metric, Unit],") ")) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))+
        scale_fill_manual(name = "", values = c("All" = "firebrick", "Sel" = "olivedrab"),labels = c("All Sites", "Bio-attached sites"))
      
    }
    
    if(which(Catalog[Category == Categories[i,"Abbrev"], NameOrigin] == Metric) %in% 
       seq(1, length(Catalog[Category == Categories[i,"Abbrev"], NameOrigin]), by = 3) && !is.null(pComp)){
      pComp <- pComp + ylab("Density")
    }
    
    return(pComp)
  })
  
  SplitPlotCategory <- split(1:length(CatPlot), ceiling(1:length(CatPlot)/9))
  
  lapply(1:length(SplitPlotCategory), function(j){
    PlotCategory <- ggarrange(plotlist = CatPlot[SplitPlotCategory[[j]]], ncol = 3,
                              nrow = ceiling(length(SplitPlotCategory[[j]])/3), common.legend = T)
    
    png(paste0(sub(" ", "",Categories[i,"category"]), "ThresholdsDraw_", j,".png"))
    print(annotate_figure(annotate_figure(PlotCategory,
                                          top=text_grob(paste0(j,"/",length(SplitPlotCategory)), size = 10),),
                          top=text_grob(Categories[i,"category"], size = 10, face = "bold")))
    dev.off()
  })
})
