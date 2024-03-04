invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid", "ggradar"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]
Titles <- fread("SumUp_Seuil.csv")

files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")

load("HYMOBIO_FULLDATA_202401.RData")

toremove <- "AUC_threshold_2402_Quant95_5step_"

files <- grep(toremove, files, value = T)
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

CompartmentCols <- c("Fish" = "cornflowerblue", "Macroinvertebrate" = "tan2", "Diatom" = "chartreuse4")

#####
Threshold_plot <- lapply(files, function(param){
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
  
  category <- Cat[Names == param, Category]
  ifelse(which(Cat[Category == category, Names] == param) %in% c(1,4),
         Ytitle <- c("AUC", "density"), Ytitle <- c("", ""))
  
  if(length(Thresh_plot) > 0){
    
  Thresh_plot <- as.data.table(do.call(rbind,lapply(Thresh_plot, function(x) {
    x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
    #x <- melt(x, id.vars = c("Threshold", "Compartment"), variable.name = Data)
    x <- as.data.table(x)[, grep(c("AUC_val|Low|Up|Threshold"), colnames(x), value = T) := 
                            lapply(.SD, function(col){return(round(as.numeric(col), 2))}), 
                          .SDcols = grep(c("AUC_val|Low|Up|Threshold"), colnames(x), value = T)]
    x$Compartment = c("Fish", "Macroinvertebrate", "Diatom"); return(x)})))[
    , Compartment := factor(Compartment, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
    , N := length(which(!is.na(AUC_valTest))), by = Compartment]
  
  Thresh_plot[, AlphaTest := "S"][AUC_valTest < 0.7, AlphaTest := "NS"][
    , AlphaTrain := "S"][AUC_valTrain < 0.7, AlphaTrain := "NS"]
  
  Thresh_max <- unique(Thresh_plot[, c("AUCmaxTest","THRmaxTest","AUCmaxTrain","THRmaxTrain") 
                                   := list(max(AUC_valTest, na.rm = T),Threshold[which.max(AUC_valTest)],
                                           max(AUC_valTrain, na.rm = T),Threshold[which.max(AUC_valTrain)]), by = Compartment][
    ,.(Compartment, AUCmaxTest, THRmaxTest, AUCmaxTrain, THRmaxTrain)])[,Param := param]
  Thresh_max[Thresh_max==-Inf] <- NA
  
  plot <- ggplot(data = Thresh_plot, aes(x = Threshold, y = AUC_valTest, color = Compartment, fill = Compartment)) + 
    geom_hline(yintercept = 0.7, color = "darkgrey", lty = 1) + 
    geom_segment(data = Thresh_max, aes(x = THRmaxTest, xend = THRmaxTest, y = 0, yend = AUCmaxTest), lty = 1, linewidth = 0.7) +
    #geom_segment(data = Thresh_max, aes(x = THRmaxTrain, xend = THRmaxTrain, y = 0, yend = AUCmaxTrain), lty = 3, linewidth = 0.7) +
    ggtitle(Titles[Name_Origine == gsub(toremove, "", param), Description]) + 
    labs(subtitle = Titles[Name_Origine == gsub(toremove, "", param), Scale]) + 
    theme_classic() + facet_grid(. ~ Compartment) + labs(x = "",y = Ytitle[1]) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
          axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
          axis.title=element_text(size=10), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 9)) +
   scale_color_manual(values = CompartmentCols)
    
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
    plot <- plot + 
      geom_vline(aes(xintercept = as.numeric(Catalog[NameOrigin == gsub(toremove, "", param), LittThreshold])), 
                 color = "firebrick", lty = 2, linewidth = 0.9, alpha = .6)
  }
  
  
  } else { 
    plot <- ggplot(data = data.frame(Thresh = 0, AUC = 0, 
                  Compartment = factor(c("Fish", "Macroinvertebrate", "Diatom"),
                  levels = c("Fish", "Macroinvertebrate", "Diatom"))), 
                  aes(x = Thresh, y = AUC)) + 
    ggtitle(gsub("^[^_]*_|_.*", "", gsub(toremove, "", param))) + theme_classic() + facet_grid(. ~ Compartment) + 
    labs(x = "", y = Ytitle) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
          axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
          axis.title=element_text(size=10), plot.title = element_text(size = 10))
    
    Thresh_max <- NA
  }
  
  Dens_plot <- RFD_all[, gsub(toremove, "", param), with = F]
  Dens_plot <- Dens_plot[complete.cases(Dens_plot)][
    get(gsub(toremove, "", param))>= quantile(Dens_plot, probs = 0.025, na.rm = T) &
    get(gsub(toremove, "", param))<= quantile(Dens_plot, probs = 0.975, na.rm = T)]
  
  Pdens <- ggplot(Dens_plot, aes(x = !!sym(gsub(toremove, "", param)))) + geom_density()+
    theme_classic() +  
    labs(x = paste0("(",Catalog[NameOrigin == gsub(toremove, "", param), Unit],") "), y = Ytitle[2]) + 
    geom_segment(data = Thresh_max, 
                 position = position_dodge2(width = (max(density(Dens_plot[[1]])$x)-min(density(Dens_plot[[1]])$x))/80),
                 aes(x = THRmaxTest, xend = THRmaxTest, color = Compartment, 
                 y = 0, yend = max(density(Dens_plot[[1]])$y)), lty = 1, linewidth = 1) +
    theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
    #annotate("rect", xmin = quantile(Dens_plot[[1]], 0.025), xmax = quantile(Dens_plot[[1]], 0.975), ymin = 0, ymax = max(density(Dens_plot[[1]])$y), alpha = .3,fill = "grey") +
    scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1)) +
    scale_color_manual(values = CompartmentCols)
  
  #plot <- plot + geom_vline(xintercept = quantile(Dens_plot[[1]], c(0.025, 0.975)), colour = "grey", linetype="dotted" ) +
    #annotate("rect", xmin = quantile(Dens_plot[[1]], 0.025), xmax = quantile(Dens_plot[[1]], 0.975), ymin = layer_scales(plot)$y$range$range[1], ymax = layer_scales(plot)$y$range$range[2], alpha = .1,fill = "grey") 
  
  Plot <- arrangeGrob(plot,Pdens, nrow = 2, heights = c(3.5,1))
  
  return(list(Plot, Thresh_max))
})

names(Threshold_plot) <- files

save(Threshold_plot, file = paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
#####

# Printing plots
#####
load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
Titles <- fread("SumUp_Seuil.csv")
files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
files <- grep(toremove, files, value = T)

Threshold_plot <- lapply(Threshold_plot, function(param) {return(param[[1]])})
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
             Titles[,.(Name_Origine, Description, Category)])

lapply(unique(Cat$Category), function(cat){
  
  Threshold_cat <- Threshold_plot[Cat[Category == cat, Names]]
  
  SplitPlotCategory <- split(1:nrow(Cat[Category == cat]), ceiling(1:nrow(Cat[Category == cat])/6))
  
lapply(1:length(SplitPlotCategory), function(j){
    
    #Image <- grid.arrange(Threshold_cat[[j]])
     
    tg <- text_grob(cat, size = 12)
    
    ifelse(grepl("Quant95", toremove), sub1 <- "95% of data", sub1 <- "")
    ifelse(grepl("5step", toremove), sub2 <- "5% step", sub2 <- "")
    ifelse(grepl("10step", toremove), sub3 <- "10% step", sub3 <- "")
    
    sg <- text_grob(paste("50 sites limit ", sub1, sub2, sub3, sep = "  "), size = 10)
    
    lt <- list(tg, sg)
    heights <- do.call(unit.c, lapply(lt, function(.g) 1.5*grobHeight(.g)))
    titles <- gtable::gtable_matrix('title', 
                                    grobs = matrix(lt, ncol=1), 
                                    widths = unit(2,'npc'),
                                    heights = heights)
    gridedPlots <-  grid.arrange(grobs = Threshold_cat[SplitPlotCategory[[j]]],
                          ncol = 3, nrow = 2, common.legend = T, legend = "right")

    png(paste0("Threshold_plot_",gsub("AUC_threshold_", "", toremove), cat, j,".png"), width = 900, height = 800,)
   
    grid.arrange(tg, sg, gridedPlots, 
                 heights = unit.c(unit.c(grobHeight(tg) + 1.2*unit(0.5, "line"), 
                                         grobHeight(sg) + unit(0.5, "line"), 
                                         unit(1,"null"))))
    
    dev.off()
    
  })
  
})

#####



# Plot difference in threshold
load("variable_thresholds")
Thresholds[, minTHR := apply(Thresholds[, !"Name_Origine", with = F],1,min)][
  ,maxTHR := apply(Thresholds[, !"Name_Origine", with = F],1,max)]
load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
Titles <- fread("SumUp_Seuil.csv")
Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))

Threshold_max <- merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], 
                       Thresholds[,.(Name_Origine, minTHR, maxTHR)])
Threshold_max[, Abs_THRmax := max(THRmaxTest), by = Param][, minTHRmax := min(THRmaxTest), by = Param][
  , scaleTHRmaxTest := (THRmaxTest-minTHR)/(maxTHR-minTHR), by = Param][
  , scaleTHRmaxTrain := (THRmaxTrain-minTHR)/(maxTHR-minTHR), by = Param][
  is.na(scaleTHRmaxTest), c("scaleTHRmaxTest","scaleTHRmaxTrain") := list(0,0)][
  , Signif := 1][AUCmaxTest < 0.7, Signif := 0][, Signif := factor(Signif, levels = c(0,1))][
  , Nmax := length(which(THRmaxTest == Abs_THRmax)), by = Param][
  , Nmin := length(which(THRmaxTest == minTHRmax)), by = Param]

Threshold_max[, isTmax := 0][THRmaxTest == Abs_THRmax & Nmax == 1 & Signif == 1, isTmax := 1][, OccMax := sum(isTmax), by = Compartment][
  , isTmin := 0][THRmaxTest == minTHRmax & Nmin == 1 & Signif == 1, isTmin := 1][, OccMin := sum(isTmin), by = Compartment]

PpieTmin_t <- unique(Threshold_max[,.(OccMin, Compartment)][order(OccMin)])
PpieTmin_t[,Compartment := factor(Compartment, levels = PpieTmin_t[order(OccMin, decreasing = T),Compartment])][, prop := OccMin/sum(OccMin)][
  ,Ypos := cumsum(prop)- 0.5*prop]
PpieTmin <- ggplot(PpieTmin_t, aes(x="", y=prop, fill=Compartment)) +
  geom_bar(stat="identity", width=1) + coord_polar("y") + 
  geom_text(aes(y = Ypos, label = OccMin), color = "white") + theme_void() + 
  scale_fill_manual(values = CompartmentCols) +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank(), axis.title = element_blank(), plot.title = element_text(size =9, hjust = 1)) +
  ggtitle("Lowest threshold")

PpieTmax_t <- unique(Threshold_max[,.(OccMax, Compartment)][order(OccMax)])
PpieTmax_t[,Compartment := factor(Compartment, levels = PpieTmax_t[order(OccMax, decreasing = T),Compartment])][, prop := OccMax/sum(OccMax)][
  ,Ypos := cumsum(prop)- 0.5*prop]
PpieTmax <- ggplot(PpieTmax_t, aes(x="", y=prop, fill=Compartment)) +
  geom_bar(stat="identity", width=1) + coord_polar("y") + 
  geom_text(aes(y = Ypos, label = OccMax), color = "white") + theme_void() + 
  scale_fill_manual(values = CompartmentCols) +
  theme(legend.position = "none",axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank(), axis.title = element_blank(),  plot.title = element_text(size =9, hjust = 1)) +
  ggtitle("Higest threshold") 

AUCmeans <- unique(Threshold_max[, c("MTest","MTrain", "NTest", "NTrain", "MthreshTest", "MthreshTrain") := 
            list(round(mean(AUCmaxTest, na.rm = T),2), length(which(AUCmaxTest>= 0.7)), 
                 round(mean(AUCmaxTrain, na.rm = T),2), length(which(AUCmaxTrain>= 0.7)),
                 round(mean(scaleTHRmaxTest), 2), round(mean(scaleTHRmaxTrain), 2)),
            by = "Compartment"][,.(Compartment, MTest, MTrain, NTest, NTrain, MthreshTest, MthreshTrain)])
setnames(AUCmeans, c("Compartment", "MTest","MTrain", "NTest", "NTrain", "MthreshTest", "MthreshTrain"), 
   c("Community", "Mean AUC", "Mean\ntrain AUC", "N signif.", 
     "N signif.\n(train)", "Average\nthreshold", "Average\ntrain. threshold"))
table2 <- merge(AUCmeans, as.data.table(CompartmentCols, keep.rownames = T), 
                by.x = "Community", by.y = "rn")[
                , Community := factor(Community,levels = levels(AUCmeans$Community))][order(Community)]

t <- tableGrob(AUCmeans[,c("Community", "Mean AUC", "N signif.", "Average\nthreshold"), with = F],
               rows = NULL, theme = ttheme_default(
  core=list(bg_params = list(fill=c(table2$CompartmentCols, rep("grey90", nrow(table2)*(ncol(table2)-1))),
                             alpha = .5)), base_size = 8))

AUCmeansCat <- merge(Threshold_max, Titles[,.(Name_Origine, Category)], by = "Name_Origine")
AUCmeansCat <- unique(AUCmeansCat[, c("MTest","MthreshTest") := list(round(mean(AUCmaxTest, na.rm = T),2), round(mean(scaleTHRmaxTest), 2)),
                                 by = c("Compartment", "Category")][,.(Compartment,Category, MTest, MthreshTest)])
PcatMean <- dcast(AUCmeansCat[,.(Compartment, Category, MTest)], Compartment ~ Category, value.var = "MTest")
PcatThresh <- dcast(AUCmeansCat[,.(Compartment, Category, MthreshTest)], Compartment ~ Category, value.var = "MthreshTest")

PcatMean <- ggplot(AUCmeansCat, aes(x = Category, y = MTest, fill = Compartment)) +
  geom_col(position = position_dodge(0.6), width = 0.5) + scale_fill_manual(values = CompartmentCols) +
  theme_classic() + labs(x = "", y = "AUC")  + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"), axis.text.x=element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = 0.7, color = "darkgrey", lty = 1) 
PcatThresh <- ggplot(AUCmeansCat, aes(x = Category, y = MthreshTest, fill = Compartment)) +
  geom_col(position = position_dodge(0.6), width = 0.5) +  theme_classic() + 
  theme(plot.margin = unit(c(0, 0, 1, 0), "lines"),legend.position = "none") +
  scale_fill_manual(values = CompartmentCols) +
  labs(x = "", y = "Threshold")  


Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))
Threshold_max <- as.data.table(merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], Titles))
Threshold_max[, Category := factor(Category, levels = c("Natural Environment", "Hydromorphology", "Hydrology", "Temperature",
                                                        "Anthropic Environment", "Flow Barriers"))]
Threshold_max <- Threshold_max[order(Category)]

Threshold_max[, Ncomp := length(which(AUCmaxTest>=0.7)), by = c("Param", "Compartment")][
  , NcompCat := sum(Ncomp)/uniqueN(Param), by = c("Category", "Compartment")]

Pradar <- dcast(unique(Threshold_max[,.(Category, Compartment, NcompCat)]), Compartment ~ Category, value.var = "NcompCat")
Pradar <- ggradar(background.circle.colour = "white", Pradar, legend.position = "none", plot.extent.x.sf = 2,
        axis.label.size = 3.5, grid.label.size = 3.5, group.line.width = 1, group.point.size = 3, fill = T,
        fill.alpha = 0.1, group.colours = CompartmentCols, label.gridline.min = F, label.gridline.mid = F,
        label.gridline.max = F)

grid.arrange(arrangeGrob(Pradar, top = "(i) Proportion of significant threshold", padding = unit(0, "line")), 
             arrangeGrob(t, top = "(ii) Summary", padding = unit(1.5, "line")), 
             arrangeGrob(PcatMean, PcatThresh, heights = c(2,4.5), 
                         top = "(iii) Mean AUC and scaled threshold by category and community", padding = unit(2, "line")), 
             arrangeGrob(Pscales, top = "(iv) Scales significance comparison", padding = unit(2, "line")), 
             
             #arrangeGrob(PpieTmin, PpieTmax, ncol = 2, top = textGrob("(iv) Community occurence as:"), padding = unit(0, "line")), 
             nrow =4, heights = c(3.5,2,4,4))


## Summary table
toremove <- "AUC_threshold_2402_Quant95_5step_"
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

Threshold_max[, length(which(AUCmaxTest >= 0.7)) / length(AUCmaxTest), by = c("Category","Compartment")]

Threshold_max[, Nsignif := length(which(AUCmaxTest >= 0.7)), by = c("Name", "Category")]
Comp3 <- unique(Threshold_max[Nsignif == 3, .(Category, Name, Scale)])
Comp2 <- unique(Threshold_max[Nsignif == 2, .(Category, Name, Scale)])
Comp1 <- unique(Threshold_max[Nsignif == 1, .(Category, Name, Scale)])

Multiscale <- unique(Titles[!Category %in% c("Temperature", "Hydrology")][which(duplicated(Name)),Name])
Multiscale <- unique(Titles[Name %in% Multiscale, Name_Origine])

Threshold_max <- merge(Threshold_max[, Name_Origine := gsub(toremove, "", Param)], 
                       Thresholds[,.(Name_Origine, minTHR, maxTHR)])
Threshold_max[, Abs_THRmax := max(THRmaxTest), by = Param][, minTHRmax := min(THRmaxTest), by = Param][
  , scaleTHRmaxTest := (THRmaxTest-minTHR)/(maxTHR-minTHR), by = Param][
    , scaleTHRmaxTrain := (THRmaxTrain-minTHR)/(maxTHR-minTHR), by = Param][
      is.na(scaleTHRmaxTest), c("scaleTHRmaxTest","scaleTHRmaxTrain") := list(0,0)]
ScaleComp <- Threshold_max[Name_Origine %in% Multiscale,
              .(Name_Origine, Name, Description, Category, Scale, AUCmaxTest, scaleTHRmaxTest)]
ScaleComp[, c("meanAUC", "meanTHR") := list(mean(AUCmaxTest), mean(scaleTHRmaxTest)), by = Name_Origine]
ScaleComp <- unique(ScaleComp[,.(Description, Category, Name, Scale, meanAUC, meanTHR)])[
  grep("Hydro area|Proximal", Scale), Scale := "Local"][
  grep("Upstream subsystem|Downstream stream", Scale), Scale := "Watershed"][
    , Signif := ifelse(meanAUC >= 0.7, "Signif", "Nsignif")]

Pscales <- ggplot(ScaleComp, aes(x = Description, y = meanAUC, fill = Scale, alpha = Signif)) +
  geom_col(position = position_dodge(0.6), width = 0.5) + scale_alpha_discrete(range = c(0.5,1)) +
  scale_fill_manual(values = c("Local" = "burlywood3",  "Watershed" = "brown")) +
  theme_classic() + labs(x = "", y = "AUC")  +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x", 
             labeller = labeller(Category = setNames(str_wrap(unique(ScaleComp$Category), width = 10),unique(ScaleComp$Category)))) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"), axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_hline(yintercept = 0.7, color = "darkgrey", lty = 1) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) 



### Plots for first ppt

CatColors <- c("#00BA38", "#00BFC4", "#619CFF","#DF536B", "#F5C710", "tomato4"); names(CatColors)<- levels(Threshold_max$Category)

Stats_cat <- unique(Threshold_max[, NcompParam := length(which(AUCmax>=0.7)), by = c("Compartment", "Category")][
  ,.(NcompParam, Compartment, Category)])

Breaks  <- unique(as.data.table(merge(Threshold_max[, Name_Origine := gsub("AUC_threshold_", "", Param)], Titles))[
  , Brk := Param[1], by = Category][,c("Brk", "Category")])[, Category := factor(Category, levels = levels(Threshold_max$Category))][
    order(Category)]

PNcomp <- ggplot(Threshold_max, aes(x=Param, y = Ncomp, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = CatColors) +
  scale_x_discrete(limits = Threshold_max$Param, breaks = Breaks[["Brk"]], labels = Breaks[["Category"]]) +
  theme_classic() + theme( axis.text.x = element_text(angle = 75, vjust = 0.5)) + 
  labs(x = "Categories", y = "N compartments") 


cols <-ggplot_build(PNcomp)
labs <- rep(" ", uniqueN(Stats_cat$Category)); names(labs) <- unique(Stats_cat$Category)

PlotStats <- ggplot(Stats_cat, aes(x=Compartment, y = NcompParam, fill = Compartment)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Category, labeller = labeller(Category = labs))  +
  theme_classic() + theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 7)) + 
  labs(x = "Categories", y = "N impacting variables") 

g <- ggplot_gtable(ggplot_build(PlotStats))
stripr <- which(grepl('strip-t', g$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- CatColors[k]
  k <- k+1
}
grid::grid.draw(g)


Threshold_max[, Nsuccess := .N/3, by = Ncomp][, propNsuccess := Nsuccess/uniqueN(Param)]


tdia <- tableGrob(Threshold_max[Compartment == "DIA" & AUCmax >=0.7, .(AUCmax, Description, Category)], rows = NULL, 
                  theme = ttheme_default(
  core=list(padding=unit(c(6, 2), "mm"), bg_params = list(alpha = .5), fg_params = list(hjust = 0, x = 0.1)), base_size = 8))

PpieA_t <- unique(Threshold_max[,.(NdiffA, Compartment)][order(NdiffA)])
PpieA_t[,Compartment := factor(Compartment, levels = PpieA_t[order(NdiffA, decreasing = T),Compartment])][, prop := NdiffA/sum(NdiffA)][
  ,Ypos := cumsum(prop)- 0.5*prop]

Ppie_success <- unique(Threshold_max[, .(Ncomp, Nsuccess, propNsuccess)])[order(Ncomp)]
Ppie_success[, Ncomp := factor(Ncomp, levels = Ppie_success[order(Ncomp, decreasing = T), Ncomp])][
  ,Ypos := cumsum(propNsuccess)- 0.5*propNsuccess]
ggplot(Ppie_success, aes(x="", y=propNsuccess, fill=Ncomp)) +
  geom_bar(stat="identity", width=1) + coord_polar("y")+ 
  geom_text(aes(y = Ypos, label = Nsuccess), color = "white") +
  theme_void() +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.title = element_text(size = 8),
        panel.grid  = element_blank(), axis.title = element_blank(), plot.title = element_text(size =12)) +
  ggtitle("Variables responsiveness") + guides(fill = guide_legend(title = "Responding\ncompartments"))

ggarrange(PNcomp, tdia, ncol = 2, widths = c(2.5,1))


# Old
######

Prel <- ggplot(data = Threshold_max, aes(x = Param, y = scaleTHRmax, group = Compartment)) + 
  geom_point(aes(colour = Compartment, shape = Signif), size = 2) +
  scale_shape_manual(values=c(1,19) , guide = "none") +
  scale_x_discrete(breaks = Breaks[["Brk"]], labels = Breaks[["Category"]]) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 75, vjust = 0.5)) + 
  labs(x = "Categories", y = "Scaled threshold") + ggtitle("Relative threshold value")


Threshold_max[,Abs_AUCmax := max(AUCmax), by = Param][, isAmax := 0][AUCmax == Abs_AUCmax, isAmax := 1][
  , NdiffA := sum(isAmax), by = Compartment]


Pauc <- ggplot(data = Threshold_max, aes(x = Param, y = AUCmax, group = Compartment)) + 
  geom_point(aes(colour = Compartment, shape = Signif), size = 2) +
  scale_shape_manual(values=c(1,19) , guide = "none") +
  scale_x_discrete(breaks = Breaks[["Brk"]], labels = Breaks[["Category"]]) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 75, vjust = 0.5)) + 
  labs(x = "Categories", y = "") + ggtitle("Maximum AUC on gradient")

