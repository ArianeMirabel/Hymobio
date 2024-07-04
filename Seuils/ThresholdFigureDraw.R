invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid", "ggradar"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("../MetricsCatalogue.csv")[which(TokeepThreshold)]
Titles <- fread("SumUp_Seuil.csv")
Titles[, Reference_Threshold := as.numeric(gsub(",",".", Reference_Threshold))]

files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")

load("../HYMOBIO_FULLDATA_202405.RData")

toremove <- "AUC_threshold_2405_Quant95_5step_"

files <- grep(toremove, files, value = T)
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

CompartmentCols <- c("Fish" = "cornflowerblue", "Macroinvertebrate" = "tan2", "Diatom" = "chartreuse4")


# Plot difference in threshold
load("variable_thresholds")
load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
load(paste0("ThresholdPlots_list_OtherMetrics", gsub("AUC_threshold_", "", toremove)))
#load(paste0("ThresholdPlots_list_F1score", gsub("AUC_threshold_", "", toremove)))

#Threshold_plot <- F1score_plot
Threshold_plot <- Threshold_plot[which(!unlist(lapply(Threshold_plot, is.null)))]
#Threshold_max <- do.call(rbind,lapply(1:length(Threshold_plot), function(x){  if(is.data.table(Threshold_plot[x][[1]][[2]])){return(Threshold_plot[x][[1]][[2]][, Param := names(Threshold_plot)[x]])}}))

Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))

Threshold_max <- as.data.table(merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], Titles))
Threshold_max[, Category := factor(Category, levels = c("Natural Environment", "Physico-Chemistry",  "Temperature","Hydrology","Hydromorphology", 
                                                        "Flow Barriers", "Anthropic Environment"))]
Threshold_max <- Threshold_max[order(Category)]
Thresholds[, minTHR := apply(Thresholds[, !"Name_Origine", with = F],1,min)][
  , maxTHR := apply(Thresholds[, !"Name_Origine", with = F],1,max)]
Titles <- fread("SumUp_Seuil.csv")

Threshold_max <- merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], 
                       Thresholds[,.(Name_Origine, minTHR, maxTHR)])

#setnames(Threshold_max, c("THRintersect", "AccuIntersect"), c("THRmaxTest", "AUCmaxTest"))

Threshold_max[, Abs_THRmax := max(THRmaxTest), by = Param][, Abs_THRmin := min(THRmaxTest), by = Param][
  , scaleTHRmaxTest := (THRmaxTest-minTHR)/(maxTHR-minTHR), by = Param][
  #, scaleTHRmaxTrain := (THRmaxTrain-minTHR)/(maxTHR-minTHR), by = Param][
  is.na(scaleTHRmaxTest), c("scaleTHRmaxTest","scaleTHRmaxTrain") := list(0,0)][
  , Signif := 1][AUCmaxTest < 0.7, Signif := 0][, Signif := factor(Signif, levels = c(0,1))][
  , Nmax := length(which(THRmaxTest == Abs_THRmax)), by = Param][
  , Nmin := length(which(THRmaxTest == Abs_THRmin)), by = Param]

Threshold_max[, isTmax := 0][THRmaxTest == Abs_THRmax & Nmax == 1 & Signif == 1, isTmax := 1][, OccMax := sum(isTmax), by = CompartmentName][
  , isTmin := 0][THRmaxTest == Abs_THRmin & Nmin == 1 & Signif == 1, isTmin := 1][, OccMin := sum(isTmin), by = CompartmentName]

#####
ListBest08 <- unique(Threshold_max[AUCmaxTest >= 0.8, Param])
save(ListBest08, file = "Best_AUC08")
ListBest07 <- unique(Threshold_max[AUCmaxTest >= 0.7, Param])
save(ListBest07, file = "Best_AUC07")

Table <- Threshold_max[,"AUCvalid" := max(AUCmaxTest), by = "Name_Origine"][AUCmaxTest == AUCvalid]
Table <- Table[,"THRvalid" := min(THRmaxTest), by = "Name_Origine"][THRmaxTest == THRvalid]

Othermetrics <- as.data.table(do.call(rbind,lapply(OtherMetrics_plot, function(X) return(X[[2]]))))
setnames(Othermetrics, c("Param", "Threshold"),c("Name_Origine","THRvalid"))

Table <- merge(Table[,.(AUCvalid, THRvalid, Name_Origine, Description, Category, CompartmentName)],
               Othermetrics[,THRvalid := as.double(THRvalid)], by = c("Name_Origine","THRvalid", "CompartmentName"),
               all.x = T)

## Pies of occurence as first/last
#####
PpieTmin_t <- unique(Threshold_max[,.(OccMin, CompartmentName)][order(OccMin)])
PpieTmin_t[, CompartmentName := factor(CompartmentName, levels = PpieTmin_t[order(OccMin, decreasing = T),CompartmentName])][, prop := OccMin/sum(OccMin)][
  ,Ypos := cumsum(prop)- 0.5*prop]
PpieTmin <- ggplot(PpieTmin_t, aes(x="", y=prop, fill=CompartmentName)) +
  geom_bar(stat="identity", width=1) + coord_polar("y") + 
  geom_text(aes(y = Ypos, label = OccMin), color = "white") + theme_void() + 
  scale_fill_manual(values = CompartmentCols) +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank(), axis.title = element_blank(), plot.title = element_text(size =9, hjust = 1)) +
  ggtitle("Lowest threshold")

PpieTmax_t <- unique(Threshold_max[,.(OccMax, CompartmentName)][order(OccMax)])
PpieTmax_t[,CompartmentName := factor(CompartmentName, levels = PpieTmax_t[order(OccMax, decreasing = T),CompartmentName])][, prop := OccMax/sum(OccMax)][
  ,Ypos := cumsum(prop)- 0.5*prop]
PpieTmax <- ggplot(PpieTmax_t, aes(x="", y=prop, fill=CompartmentName)) +
  geom_bar(stat="identity", width=1) + coord_polar("y") + 
  geom_text(aes(y = Ypos, label = OccMax), color = "white") + theme_void() + 
  scale_fill_manual(values = CompartmentCols) +
  theme(legend.position = "none",axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank(), axis.title = element_blank(),  plot.title = element_text(size =9, hjust = 1)) +
  ggtitle("Higest threshold") 
#####


## Barplot of mean AUC and Thresholds
#####
AUCmeansCat <- merge(Threshold_max, Titles[,.(Name_Origine, Category)], by = "Name_Origine")

AUCmeansCat <- unique(Threshold_max[, c("MTest","MthreshTest") := list(round(mean(AUCmaxTest, na.rm = T),2), round(mean(scaleTHRmaxTest), 2)),
                                 by = c("CompartmentName", "Category")][,.(CompartmentName,Category, MTest, MthreshTest)])
PcatMean <- dcast(AUCmeansCat[,.(CompartmentName, Category, MTest)], CompartmentName ~ Category, value.var = "MTest")
PcatThresh <- dcast(AUCmeansCat[,.(CompartmentName, Category, MthreshTest)], CompartmentName ~ Category, value.var = "MthreshTest")

PcatMean <- ggplot(AUCmeansCat, aes(x = Category, y = MTest, fill = CompartmentName)) +
  geom_col(position = position_dodge(0.6), width = 0.5) + scale_fill_manual(values = CompartmentCols) +
  theme_classic() + labs(x = "", y = "AUC")  + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"), axis.text.x=element_blank(), legend.position = "none") +
  geom_hline(yintercept = 0.7, color = "darkgrey", lty = 1) 
PcatThresh <- ggplot(AUCmeansCat, aes(x = Category, y = MthreshTest, fill = CompartmentName)) +
  geom_col(position = position_dodge(0.6), width = 0.5) +  theme_classic() + 
  theme(plot.margin = unit(c(0, 0, 1, 0), "lines"), legend.position = "none") +
  scale_fill_manual(values = CompartmentCols) +
  labs(x = "", y = "Threshold")  

legd_Group <- ggplot_gtable(ggplot_build(PcatThresh + theme(legend.position = "right"))) 
legd_Group <- legd_Group$grobs[[which(sapply(legd_Group$grobs, function(x) x$name) == "guide-box")]]
#####

## Radar graphic
#####
Threshold_max[, NcompCat := length(which(AUCmaxTest >= 0.6))/uniqueN(Param), by = c("Category", "CompartmentName")]

Pradar <- dcast(unique(Threshold_max[,.(Category, CompartmentName, NcompCat)]), CompartmentName ~ Category, value.var = "NcompCat")

apply(Pradar[,!"CompartmentName"], 2, mean)

Pradar <- ggradar(background.circle.colour = "white", Pradar, legend.position = "none", plot.extent.x.sf = 2,
        axis.label.size = 3.5, grid.label.size = 3.5, group.line.width = 1, group.point.size = 3, fill = T,
        fill.alpha = 0.1, group.colours = CompartmentCols, label.gridline.min = F, label.gridline.mid = F,
        label.gridline.max = F) +
  ggtitle("Proportion of significant models") + 
  theme(legend.position = "none", plot.title = element_text(size = 10))
#####


## Summary table
#####
load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
Titles <- fread("SumUp_Seuil.csv")
files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
files <- grep(toremove, files, value = T)
load("variable_thresholds")
Thresholds[, minTHR := apply(Thresholds[, !"Name_Origine", with = F],1,min)][
  ,maxTHR := apply(Thresholds[, !"Name_Origine", with = F],1,max)]

Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))
Threshold_max <- as.data.table(merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], Titles))

Nsignif_Test <- nrow(Threshold_max[AUCmaxTest >= 0.7])
Nsignif_Train <- nrow(Threshold_max[AUCmaxTrain >= 0.7])

Threshold_max[, length(which(AUCmaxTest >= 0.7)) / length(AUCmaxTest), by = c("Category","CompartmentName")]

Threshold_max[, Nsignif := length(which(AUCmaxTest >= 0.7)), by = c("Param", "Category")]
Comp3 <- unique(Threshold_max[Nsignif == 3, .(Category, Description_detailed, Name, Scale)])
Comp2 <- unique(Threshold_max[Nsignif == 2, .(Category, Name, Scale)])
Comp1 <- unique(Threshold_max[Nsignif == 1, .(Category, Name, Scale)])

Multiscale <- unique(Titles[!Category %in% c("Temperature", "Hydrology")][which(duplicated(Name)),Name])
Multiscale <- unique(Titles[Name %in% Multiscale, Name_Origine])
#####


## Difference between literature and observed thresholds
#####
Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))
Threshold_max <- as.data.table(merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], Titles))[
  Reference_Threshold != "",][, Reference_Threshold := as.numeric(Reference_Threshold)]
Threshold_max <- merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], 
                       Thresholds[,.(Name_Origine, minTHR, maxTHR)])

DiffThresh <- Threshold_max[, DiffThresh := THRmaxTest - Reference_Threshold][
  , Above := length(which(DiffThresh > 0)), by = Category][
    , Equal := length(which(DiffThresh == 0)), by = Category][
      , Below := length(which(DiffThresh < 0)), by = Category]

DiffThresh <- unique(melt(setDT(DiffThresh[,.(Category,Above, Equal, Below)]), id.vars = "Category"))


ggplot(DiffThresh, aes(x = Category, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = data.frame(Above = "#619CFF", Equal = "grey", Below = "firebrick")) + 
  labs(y = "Threshold comparison")

DiffThresh <- Threshold_max[, DiffThresh := THRmaxTest - Reference_Threshold][
  , Above := length(which(DiffThresh > 0)), by = CompartmentName][
    , Equal := length(which(DiffThresh == 0)), by = CompartmentName][
      , Below := length(which(DiffThresh < 0)), by = CompartmentName]

DiffThresh <- unique(melt(setDT(DiffThresh[,.(CompartmentName,Above, Equal, Below)]), id.vars = "CompartmentName"))

ggplot(DiffThresh, aes(x = CompartmentName, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = data.frame(Above = "#619CFF", Equal = "grey", Below = "firebrick")) + 
  labs(y = "Threshold comparison")


#####

grid.arrange(PdiffThreshFull, PdiffThreshLim, nrow = 1)


grid.arrange(arrangeGrob(Pradar, top = "(i) Proportion of significant threshold", padding = unit(0, "line")), 
             arrangeGrob(t, top = textGrob("(ii) Summary", vjust = -0.5), padding = unit(1, "line")), 
             arrangeGrob(arrangeGrob(PcatMean, PcatThresh, heights = c(2,4.5)), 
                         top = textGrob("(iii) Mean AUC and scaled threshold by category and community", vjust = 1), 
                         padding = unit(2, "line"), legd_Group, nrow = 1, widths = c(4,1)), 
             arrangeGrob(Pscales, top = textGrob("(iv) Scales significance comparison", vjust = -0.5)), 
             #arrangeGrob(PdiffThresh, top =  textGrob("(v) Difference with literature threshold by category", vjust = -0.5)),
             #arrangeGrob(PpieTmin, PpieTmax, ncol = 2, top = textGrob("(iv) Community occurence as:"), padding = unit(0, "line")), 
             nrow =4, heights = c(3.5,2,4,4))


#####
# Plot AUC range
CatColors <- c("springgreen4", "yellowgreen","#00BFC4", "#DF536B", "#619CFF","#F5C710", "tomato4"); names(CatColors)<- levels(Threshold_max$Category)

Plots <- lapply(levels(Threshold_max$CompartmentName), function(comm){
  
ret <- ggplot(Threshold_max[CompartmentName == comm], aes(x = reorder(Param, -as.numeric(AUCmaxTest)), y = AUCmaxTest, fill = Category)) +
  geom_col(position = position_dodge(0.6), width = 0.9) + scale_fill_manual(values = CatColors) +
  geom_hline(yintercept = c(0.6, 0.7, 0.8), color = c("darkgrey", "ivory4", "black"), lty = 1) + 
  coord_cartesian(ylim=c(0.4, 1)) +  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1))  +
  theme_classic() + labs(x = "", y = "F1-score")  + ggtitle(comm) + 
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))

if(comm == "Fish"){
  ret <- ret + 
    geom_text(aes(x = uniqueN(Threshold_max$Param) - 10, y = 0.6, label = "Poor", vjust = 1, hjust = 0), color = "darkgrey") +
    geom_text(aes(x = uniqueN(Threshold_max$Param) - 10, y = 0.7, label = "Correct", vjust = 1, hjust = 0), color = "ivory4") +
    geom_text(aes(x = uniqueN(Threshold_max$Param) - 10, y = 0.8, label = "Good", vjust = 1, hjust = 0), color = "black") 
}

return(ret)

})

toplot <- Threshold_max[AUCmaxTest > 0.8, Name_Origine]

PiePlot <- Threshold_max[, .(Ncorrect = length(which(AUCmaxTest >= 0.7)),
                             count = .N), by = c("Category", "CompartmentName")]
#PiePlot[, Ncorrect := Ncorrect/count]
ForPlot <- unique(PiePlot[,.(Ncorrect, Category, CompartmentName)])
ForPlot[, Total := sum(Ncorrect), by = "CompartmentName"][
  , FacetName := do.call(paste, c(.SD, sep = "\nN = ")),.SDcols = c("CompartmentName", "Total")][
  , FacetName := factor(ForPlot$FacetName, 
                 unique(ForPlot[order(CompartmentName),.(FacetName, CompartmentName)])[["FacetName"]])
  ]

PiePlot <- ggplot(ForPlot, aes(x="", y=Ncorrect, fill=Category)) + 
  geom_bar(stat="identity", width=1, position = "fill") + coord_polar("y") + theme_void() + 
  scale_fill_manual(values = CatColors ) + 
  facet_grid(. ~ FacetName) + 
  ggtitle("Distribution of correct models (F1score > 0.7) by category") + 
  theme(legend.position = "none", plot.title = element_text(size = 10))

Plots[[4]] <- PiePlot

ggarrange(plotlist = Plots, nrow = 4, common.legend = T, legend = "right")

Threshold_max[, AUCmax_cat := AUCmaxTest[which.max(AUCmaxTest)], by = Param][,
  Param_comm := paste0(Param, CompartmentName)]
Gplot_cat <- ggplot(Threshold_max, aes(x = reorder(Param_comm, -as.numeric(AUCmaxTest)), y = as.numeric(AUCmaxTest), fill = Category)) +
  geom_col(position = position_dodge(0.6), width = 0.9) + scale_fill_manual(values = CatColors) +
  geom_hline(yintercept = c(0.6, 0.7, 0.8), color = c("darkgrey", "ivory4", "black"), lty = 1) + 
  coord_cartesian(ylim=c(0.4, 1)) +  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1))  +
  theme_classic() + labs(x = "", y = "F1-score")  + ggtitle("Range by category") + 
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))

Gplot_comp <- ggplot(Threshold_max, aes(x = reorder(Param_comm, -as.numeric(AUCmaxTest)), y = as.numeric(AUCmaxTest), fill = CompartmentName)) +
  geom_col(position = position_dodge(0.6), width = 0.9) + scale_fill_manual(values = CompartmentCols) +
  geom_hline(yintercept = c(0.6, 0.7, 0.8), color = c("darkgrey", "ivory4", "black"), lty = 1) + 
  coord_cartesian(ylim=c(0.4, 1)) +  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1))  +
  theme_classic() + labs(x = "", y = "F1-score")  + ggtitle("Range by community") + 
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) #+ 
  #geom_text(aes(x = uniqueN(Threshold_max$Param_comm) - 10, y = 0.8, label = "Poor", vjust = 1, hjust = 0), color = "darkgrey") +
  #geom_text(aes(x = uniqueN(Threshold_max$Param_comm) - 10, y = 0.7, label = "Correct", vjust = 1, hjust = 0), color = "ivory4") +
  #geom_text(aes(x = uniqueN(Threshold_max$Param_comm) - 10, y = 0.8, label = "Good", vjust = 1, hjust = 0), color = "black") 


ggarrange(ggarrange(Gplot_comp,Pradar,  nrow = 2, heights = c(1,0.8)),
          ggarrange(Gplot_cat, PiePlot, common.legend = T, nrow = 2, legend = "right"), nrow = 2)
          

#####

#Comparison AUC and crossing
#####
load("variable_thresholds")
load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
load(paste0("ThresholdPlots_list_OtherMetrics", gsub("AUC_threshold_", "", toremove)))


Threshold_plot <- Threshold_plot[which(!unlist(lapply(Threshold_plot, is.null)))]
Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))

Threshold_max <- as.data.table(merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], Titles))
Threshold_max[, Category := factor(Category, levels = c("Natural Environment", "Physico-Chemistry",  "Temperature","Hydrology","Hydromorphology", 
                                                        "Flow Barriers", "Anthropic Environment"))]
Threshold_max <- Threshold_max[order(Category)]
Thresholds[, minTHR := apply(Thresholds[, !"Name_Origine", with = F],1,min)][
  , maxTHR := apply(Thresholds[, !"Name_Origine", with = F],1,max)]
Titles <- fread("SumUp_Seuil.csv")

Threshold_max <- merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], 
                       Thresholds[,.(Name_Origine, minTHR, maxTHR)])

Threshold_max[, Abs_THRmax := max(THRmaxTest), by = Param][, Abs_THRmin := min(THRmaxTest), by = Param][
  , scaleTHRmaxTest := (THRmaxTest-minTHR)/(maxTHR-minTHR), by = Param][
    #, scaleTHRmaxTrain := (THRmaxTrain-minTHR)/(maxTHR-minTHR), by = Param][
    is.na(scaleTHRmaxTest), c("scaleTHRmaxTest","scaleTHRmaxTrain") := list(0,0)][
      , Signif := 1][AUCmaxTest < 0.7, Signif := 0][, Signif := factor(Signif, levels = c(0,1))][
        , Nmax := length(which(THRmaxTest == Abs_THRmax)), by = Param][
          , Nmin := length(which(THRmaxTest == Abs_THRmin)), by = Param]

Threshold_max[, isTmax := 0][THRmaxTest == Abs_THRmax & Nmax == 1 & Signif == 1, isTmax := 1][, OccMax := sum(isTmax), by = CompartmentName][
  , isTmin := 0][THRmaxTest == Abs_THRmin & Nmin == 1 & Signif == 1, isTmin := 1][, OccMin := sum(isTmin), by = CompartmentName]



OtherMetrics_max <- do.call(rbind,lapply(1:length(OtherMetrics_plot), function(x){
  if(is.data.table(OtherMetrics_plot[x][[1]][[2]])){
    return(OtherMetrics_plot[x][[1]][[2]][, Param := names(OtherMetrics_plot)[x]])}}))


OtherMetrics_max <- as.data.table(merge(OtherMetrics_max[, Name_Origine := gsub(toremove, "", Param)], Titles))
OtherMetrics_max[, Category := factor(Category, levels = c("Natural Environment", "Physico-Chemistry",  "Temperature","Hydrology","Hydromorphology", 
                                                        "Flow Barriers", "Anthropic Environment"))]
OtherMetrics_max <- OtherMetrics_max[order(Category)]
Thresholds[, minTHR := apply(Thresholds[, !"Name_Origine", with = F],1,min)][
  , maxTHR := apply(Thresholds[, !"Name_Origine", with = F],1,max)]
Titles <- fread("SumUp_Seuil.csv")

OtherMetrics_max <- merge(OtherMetrics_max[, Name_Origine := gsub(toremove, "", Param)], 
                       Thresholds[,.(Name_Origine, minTHR, maxTHR)])

OtherMetrics_max[, Abs_THRmax_inter := max(THRintersect), by = Param][, Abs_THRmin_inter := min(THRintersect), by = Param][
  , scaleTHRintersect := (THRintersect-minTHR)/(maxTHR-minTHR), by = Param][
    , scaleTHRintersect := (THRintersect-minTHR)/(maxTHR-minTHR), by = Param][
    is.na(scaleTHRintersect), c("scaleTHRintersect","scaleTHRmaxTrain") := list(0,0)][
      , Signif := 1][AccuIntersect < 0.7, Signif := 0][, Signif := factor(Signif, levels = c(0,1))][
        , Nmax := length(which(THRintersect == Abs_THRmax_inter)), by = Param][
          , Nmin := length(which(THRintersect == Abs_THRmin_inter)), by = Param]


Thresh_comparison <- merge(Threshold_max[,.(CompartmentName, Param, Description, Category, THRmaxTest, scaleTHRmaxTest)],
                           OtherMetrics_max, by = c("CompartmentName", "Param"))[
                             , scaleDiff_Thr := scaleTHRmaxTest - scaleTHRintersect][
                               , Diff_Thr := THRmaxTest - THRintersect]
mean(Thresh_comparison$Diff_Thr)
mean(Thresh_comparison$scaleDiff_Thr)
