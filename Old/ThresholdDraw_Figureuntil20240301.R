invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid", "ggradar"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("../MetricsCatalogue.csv")[which(TokeepThreshold)]
Titles <- fread("SumUp_Seuil.csv")
Titles[, Reference_Threshold := as.numeric(gsub(",",".", Reference_Threshold))]

files <- list.files("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils/Genomig_Output")

load("../HYMOBIO_FULLDATA_202405.RData")

toremove <- "AUC_threshold_2405_Quant95_5step_"

files <- grep(toremove, files, value = T)
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

CompartmentCols <- c("Fish" = "cornflowerblue", "Macroinvertebrate" = "tan2", "Diatom" = "chartreuse4")


# Plot difference in threshold
load("variable_thresholds")
Thresholds[, minTHR := apply(Thresholds[, !"Name_Origine", with = F],1,min)][
  , maxTHR := apply(Thresholds[, !"Name_Origine", with = F],1,max)]

Table <- do.call(rbind,lapply(files, function(param){
  
  print(param)
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_draw <- Thresh_draw[unlist(lapply(Thresh_draw, is.list))]
  Thresh_names <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]][,"Compartment"])})
  Thresh_thresholds <- lapply(Thresh_draw, function(X){return(unique(X[["Thresh_AUCval"]]["Threshold"]))})
  
  Thresh_AUC <- as.data.table(do.call(rbind, lapply(1:length(Thresh_draw), function(X){
    
    auc <- Thresh_draw[[X]][["Thresh_AUCval"]]
    other <- Thresh_draw[[X]][["Thresh_OtherMetrics"]]
    Names <- Thresh_names[[X]]
    Thresh <- Thresh_thresholds[[X]]
    
    if(length(other) > 0){
      other <- do.call(rbind,lapply(1:length(other), function(XComp) {
        
        Comp <- other[[XComp]]
        
        if(length(Comp) !=0){
          
          Comp_Test <- as.data.table(do.call(rbind, lapply(Comp, function(comp){
            
            comp <- comp[["ConfusionTest"]]
            
            return(c(comp$overall["Accuracy"],
                     comp$byClass[c("Sensitivity", "Specificity", "Precision", "Recall", "F1")]))})))
          
          Comp_Test <-  Comp_Test[, lapply(.SD, quantile, prob =  .5, na.rm = T)][,
                                    Threshold := Thresh][,Compartment := Names[XComp]]
          
          return(Comp_Test)
        }
      })
      )
      
      other <- merge(other, data.table(Compartment = c("B_FISH", "B_INV", "B_DIA"),
                                       CompartmentName = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
      
      Table <- merge(auc[,c("AUC_valTest", "Threshold", "Compartment")], other, by = c("Compartment", "Threshold"))
      
      return(Table)}
  })))
  
  if(length(Thresh_AUC)!=0){
    return(Thresh_AUC[
      , Name_Origine := gsub(toremove, "",param)])
  }
})
)
cols <- setdiff(names(Table),c("Compartment","CompartmentName","Name_Origine"))
Table <- Table[, (cols) := lapply(.SD, function(X) round(as.numeric(X), 2)),
               .SDcols = cols]
setnames(Table, "CompartmentName", "Community")

Table <- merge(Table,Thresholds[,.(Name_Origine, minTHR, maxTHR)])

Table <- merge(Table, Titles[,.(Description, Name_Origine, Unit, Category, Scale)])

#####
ListBest09 <- unique(Table[AUC_valTest >= 0.9 & F1 >= 0.7, .(Name_Origine, Community)])
save(ListBest09, file = "Best_AUC09")
ListBest08 <- unique(Table[AUC_valTest >= 0.8 & AUC_valTest < 0.9 & F1 >= 0.7, .(Name_Origine, Community)])
save(ListBest08, file = "Best_AUC08")
ListBest07 <- unique(Table[AUC_valTest >= 0.7  & F1 >= 0.7,])
save(ListBest07, file = "Best_AUC07")

unique(Table[AUC_valTest >= 0.7 & F1 >= 0.7, .(Name_Origine, Community)])[,.N, by = Community]

## Pies of occurence as first/last
#####
Threshold_max <- Table[, AUCmaxTest := max(AUC_valTest), by = c("Name_Origine","Compartment")][
  ,THRmaxTest := Threshold[which.max(AUC_valTest)], by = Name_Origine]
Threshold_max <- unique(Threshold_max[AUC_valTest >= 0.7 & F1 >= 0.7,
                          .(Name_Origine, Community, minTHR, maxTHR, Description, Category, Threshold, AUC_valTest,
                                  AUCmaxTest)])
Threshold_max[, mean(AUC_valTest), by = Community]
Threshold_max[, MaxComm := Community[which.max(AUC_valTest)], by = Name_Origine]
unique(Threshold_max[,.(MaxComm, Name_Origine)])[,.N, by = MaxComm]



setwd("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils/Tabebuia_Output")
PI_bio_merge <- do.call(rbind, lapply(list.files(pattern = "ImpactProb_"), function(Param){
  load(Param); return(PI_bio_merge)}))
PI_bio_merge <- melt(PI_bio_merge, measure.vars = c("B_FISH", "B_INV", "B_DIA"), 
                     variable.name = "Compartment", value.name = "Pimpact")

setnames(Threshold_max, "Name_Origine", "Param")
PI_bio_merge <- merge(PI_bio_merge, unique(Threshold_max[,.(Param, Compartment, AUCmaxTest)]),
                      by = c("Param", "Compartment"))
PI_bio_merge <- dcast(PI_bio_merge, Param + ID + Year + AUCmaxTest ~ Compartment, value.var = "Pimpact")

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
CatColors <- c("springgreen4", "yellowgreen","#DF536B", "#00BFC4", "#619CFF","#F5C710", "tomato4"); names(CatColors)<- levels(Threshold_max$Category)

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

### Specific variable drawing
CompartmentCols <- c("Fish" = "cornflowerblue", "Macroinvertebrate" = "tan2", "Diatom" = "chartreuse4")
MetricShape <- c("Accuracy" = 16, "Sensitivity" = 1, "Recall" = 10, "Specificity" = 6, "Precision" = 8, "F1" = 19)

param <- grep("HIGHFARM13", files, value = T)

load(paste0(getwd(),"/Genomig_Output/",param))

Thresh_draw <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]])})

  Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
  
  category <- Cat[Names == param, Category]
  
    Thresh_plot <- as.data.table(do.call(rbind,lapply(Thresh_plot, function(x) {
      
      x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
      #x <- melt(x, id.vars = c("Threshold", "Compartment"), variable.name = Data)
      
      x <- as.data.table(x)[, grep(c("AUC_val|Low|Up"), colnames(x), value = T) := 
                              lapply(.SD, function(col){return(round(as.numeric(col), 2))}), 
                            .SDcols = grep(c("AUC_val|Low|Up"), colnames(x), value = T)]
      
      x <- merge(x, data.table(Compartment = c("B_FISH", "B_INV", "B_DIA"),
                               CompartmentName = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
      return(x)})))[
        , CompartmentName := factor(CompartmentName, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
          , N := length(which(!is.na(AUC_valTest))), by = CompartmentName]
    
    Thresh_plot[, AlphaTest := "S"][AUC_valTest < 0.7, AlphaTest := "NS"][
      , AlphaTrain := "S"][AUC_valTrain < 0.7, AlphaTrain := "NS"]
    
    Thresh_max <- unique(Thresh_plot[, c("AUCmaxTest","THRmaxTest","AUCmaxTrain","THRmaxTrain") 
                                     := list(max(AUC_valTest, na.rm = T),Threshold[which.max(AUC_valTest)],
                                             max(AUC_valTrain, na.rm = T),Threshold[which.max(AUC_valTrain)]), by = CompartmentName][
                                               ,.(CompartmentName, AUCmaxTest, THRmaxTest, AUCmaxTrain, THRmaxTrain)])[,Param := param]
    Thresh_max[Thresh_max==-Inf] <- NA
    
    if(!is.na(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold]) && 
       as.numeric(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold]) > 
       quantile(RFD_all[, gsub(toremove, "", param), with = F], probs = 0.975, na.rm = T) && 
       !any(Thresh_max$THRmaxTest == as.numeric(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold]))){
      Thresh_plot <- Thresh_plot[Threshold != as.numeric(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold])]
    }
    
    ifelse(any(nchar(Thresh_plot$Threshold) > 4) & diff(range(as.numeric(Thresh_plot$Threshold), na.rm = T))/10 < 0.01, 
           makeBreaks <- T, makeBreaks <- F)
    Thresh_plot[, Threshold := as.numeric(Threshold)]
    tonum <- c("THRmaxTrain", "THRmaxTest")
    Thresh_max[, (tonum) := lapply(.SD, as.numeric), .SDcols = tonum]
    
    plot <- ggplot(data = Thresh_plot, aes(x = Threshold, y = AUC_valTest, color = CompartmentName, fill = CompartmentName)) + 
      geom_hline(yintercept = c(0.7, 0.8), color = "darkgrey", lty = 1) + 
      geom_segment(data = Thresh_max, aes(x = THRmaxTest, xend = THRmaxTest, y = 0, yend = AUCmaxTest), lty = 1, linewidth = 0.7) +
      #geom_segment(data = Thresh_max, aes(x = THRmaxTrain, xend = THRmaxTrain, y = 0, yend = AUCmaxTrain), lty = 3, linewidth = 0.7) +
      ggtitle(Titles[Name_Origine == gsub(toremove, "", param), Description]) + 
      labs(subtitle = Titles[Name_Origine == gsub(toremove, "", param), Scale]) + 
      theme_classic() + facet_grid(. ~ CompartmentName) + labs(x = "",y = "AUC") + 
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
            axis.text.x = element_text(angle = 45, vjust = 0.6), axis.text=element_text(size=8),
            axis.title=element_text(size=10), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 9)) +
      scale_color_manual(values = CompartmentCols)
    
    if(makeBreaks){ plot <- plot + 
      scale_x_continuous(breaks = seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T), length.out = 10),
                         labels = format(seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T),
                                             length.out = 10), digits = 2, scientific = T))
    }
    
    if(any(Thresh_plot$N > 1)){
      plot <- plot + geom_ribbon(aes(ymax = UpTest, ymin = LowTest), colour = NA, alpha = 0.1) + 
        geom_point(data = Thresh_plot[N > 1], aes(x = Threshold, y = AUC_valTest, alpha = AlphaTest)) + 
        geom_ribbon(aes(ymax = UpTrain, ymin = LowTrain), colour = NA, alpha = 0.1) + 
        geom_point(data = Thresh_plot[N > 1], aes(x = Threshold, y = AUC_valTrain, alpha = AlphaTrain), shape = 24, size = 0.6) +
        scale_alpha_discrete(range = c(0.35,1)) +
        scale_fill_manual(values = CompartmentCols)
    } 
    
    if(any(Thresh_plot$N == 1)) {
      plot <- plot + geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = AUC_valTest, alpha = AlphaTest)) + 
        geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = AUC_valTrain, alpha = AlphaTrain), shape = 24) + 
        geom_errorbar(data = Thresh_plot[N == 1],aes(ymax = UpTest, ymin = LowTest), width = 0.2) + 
        geom_errorbar(data = Thresh_plot[N == 1],aes(ymax = UpTrain, ymin = LowTrain), width = 0.2) +
        scale_alpha_discrete(range = c(0.35,1)) +
        scale_fill_manual(values = CompartmentCols)
    }
    
    if(!is.na(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold])){
      if(as.numeric(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold]) > 
         max(Thresh_plot$Threshold, na.rm = T) + (max(Thresh_plot$Threshold, na.rm = T)-min(Thresh_plot$Threshold, na.rm = T))/2){
        plot <- plot + labs(tag = "Litterature threshold\nout of range (higher)") +
          theme(plot.tag = element_text(hjust = 0.5, color = "firebrick", size = 9), 
                plot.tag.position = c(0.85, 0.99))
      } else {plot <- plot + 
        geom_vline(aes(xintercept = as.numeric(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold])), 
                   color = "firebrick", lty = 2, linewidth = 0.9, alpha = .6)}
      
    }
    
  plotAUC <- plot
  
  Dens_plot <- RFD_all[, gsub(toremove, "", param), with = F]
  Dens_plot <- Dens_plot[complete.cases(Dens_plot)][
    get(gsub(toremove, "", param))>= quantile(Dens_plot, probs = 0.025, na.rm = T) &
      get(gsub(toremove, "", param))<= quantile(Dens_plot, probs = 0.975, na.rm = T)]
  
  Pdens <- ggplot(Dens_plot, aes(x = !!sym(gsub(toremove, "", param)))) + geom_density()+
    theme_classic() +  
    labs(x = paste0("(",Catalog[NameOrigin == gsub(toremove, "", param), Unit],") "), y = "density") + 
    geom_segment(data = Thresh_max, 
                 position = position_dodge2(width = (max(density(Dens_plot[[1]])$x)-min(density(Dens_plot[[1]])$x))/80),
                 aes(x = THRmaxTest, xend = THRmaxTest, color = CompartmentName, 
                     y = 0, yend = max(density(Dens_plot[[1]])$y)), lty = 1, linewidth = 1) +
    theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
    scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1)) +
    scale_color_manual(values = CompartmentCols)
  
  
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_names <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]][,"Compartment"])})
  Thresh_thresholds <- lapply(Thresh_draw, function(X){return(unique(X[["Thresh_AUCval"]]["Threshold"]))})
  Thresh_draw <- lapply(Thresh_draw, function(X){return(X[["Thresh_OtherMetrics"]])})
  
    Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.list))]
    Thresh_names <- Thresh_names[unlist(lapply(Thresh_draw, is.list))]
    Thresh_thresholds <- Thresh_thresholds[unlist(lapply(Thresh_draw, is.list))]
    
    category <- Cat[Names == param, Category]
    
    IsShp <- IsCol <- "none"

      Thresh_plot <- as.data.table(do.call(rbind,lapply(1:length(Thresh_plot), function(Xcomp) {
        
        Comp <- Thresh_plot[[Xcomp]]
        Names <- Thresh_names[[Xcomp]]
        Thresh <- Thresh_thresholds[[Xcomp]]
        
        if(length(Comp) !=0){
          Comp_Test <- as.data.table(do.call(rbind, lapply(1:length(Comp), function(xcomp){
            
            x <- Comp[[xcomp]]
            name <- Names[[xcomp]]
            
            ret <- data.table(do.call(rbind, lapply(x, function(petitX){return(
              c(petitX$ConfusionTest$overall["Accuracy"],
                petitX$ConfusionTest$byClass[c("Sensitivity", "Specificity", "Precision", "Recall", "F1")])
            )})))
            ret <-  ret[, lapply(.SD, quantile, prob = c(0.05, .5, 0.95), na.rm = T), 
                        .SDcols = colnames(ret)][, p := c("Low", "Med", "High")]
            
            ret <- dcast(melt(ret, id.vars = "p", variable.name = "Metric", value.name = "Value"),
                         Metric ~ p , value.var = "Value")[, Compartment := name][, Threshold := Thresh]
            
            return(ret)
          })
          ))
          Comp_Test <- merge(Comp_Test, data.table(Compartment = c("B_FISH", "B_INV", "B_DIA"),
                                                   CompartmentName = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
          return(Comp_Test)
        } 
      })))
  
      
      Thresh_plot[, CompartmentName := factor(CompartmentName, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
        , N := length(which(!is.na(Med))), by = c("CompartmentName", "Metric")]
      
      
      Thresh_max <- unique(dcast(Thresh_plot[!is.na(Threshold)], 
                                 CompartmentName + Threshold ~ Metric, value.var = "Med")[
                                   , c("THRintersect", "AccuIntersect") := list(Threshold[which.min(abs(Sensitivity - Specificity))],
                                                                                Accuracy[which.min(abs(Sensitivity - Specificity))]),
                                   by = CompartmentName][,.(CompartmentName, THRintersect, AccuIntersect)])
      
      Thresh_max[THRintersect ==-Inf] <- NA
      
      ifelse(any(nchar(Thresh_plot$Threshold) > 4) & diff(range(as.numeric(Thresh_plot$Threshold), na.rm = T))/10 < 0.01, 
             makeBreaks <- T, makeBreaks <- F)
      Thresh_plot[, Threshold := as.numeric(Threshold)]
      tonum <- c("THRintersect", "AccuIntersect")
      Thresh_max[, (tonum) := lapply(.SD, as.numeric), .SDcols = tonum]
      
      Thresh_plot <- Thresh_plot[Metric %in% c("Specificity", "Recall", "F1")]
      
      plot <- ggplot(data = Thresh_plot, aes(x = Threshold, y = Med, color = CompartmentName, fill = CompartmentName)) + 
        geom_segment(data = Thresh_max, aes(x = THRintersect, xend = THRintersect, y = -0.05, yend = AccuIntersect), lty = 2, linewidth = 0.7) +
         theme_classic() + facet_grid(. ~ CompartmentName) + labs(x = "",y = "Other Metrics", color = "Community") + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), axis.text.x = element_text(angle = 45, vjust = 0.6),
              axis.text=element_text(size=8), axis.title=element_text(size=10), plot.title = element_text(size = 10), 
              plot.subtitle = element_text(size = 9), strip.background = element_blank(),
              strip.text.x = element_blank(), legend.position = "none") +
        guides(colour = IsCol, shape = IsShp, fill = "none") + 
        scale_color_manual(values = CompartmentCols) + scale_fill_manual(values = CompartmentCols) 
      
      if(makeBreaks){ plot <- plot + 
        scale_x_continuous(breaks = seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T), length.out = 10),
                           labels = format(seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T), length.out = 10),
                                           digits = 2, scientific = T))
      }
      
      
        plot <- plot + geom_ribbon(aes(ymax = High, ymin = Low, group = Metric), colour = NA, alpha = 0.1) + 
          geom_point(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med, group = Metric, shape = Metric)) + 
          geom_line(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med, group = Metric)) + 
          scale_shape_manual(values = MetricShape, guide = IsShp) + labs(shape = "Model Metrics")
        
        
  tg <- text_grob(category, size = 10)
  
  sg <- text_grob("", size = 1)
  
  lt <- list(tg, sg)
  heights <- do.call(unit.c, lapply(lt, function(.g) 1.5*grobHeight(.g)))
  titles <- gtable::gtable_matrix('title', 
                                  grobs = matrix(lt, ncol=1), 
                                  widths = unit(2,'npc'),
                                  heights = heights)
  Legd <- arrangeGrob(rasterGrob(image_read("Legend_allPlots.png")))
  Toplot <- arrangeGrob(plotAUC, plot, Pdens, Legd, nrow = 4, heights = c(2.5,1.5,1,0.5))
 
  grid.arrange(tg, sg, Toplot, 
               heights = unit.c(unit.c(grobHeight(tg) + 1.2*unit(0.5, "line"), 
                                       grobHeight(sg) + unit(0.5, "line"), 
                                       unit(1,"null"))))
 
  