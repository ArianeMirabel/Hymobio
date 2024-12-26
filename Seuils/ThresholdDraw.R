#devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid", "devtools","scales", "magick"),function(pk){
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

CommunityCols <- c("Fish" = "cornflowerblue", "Macroinvertebrate" = "tan2", "Diatom" = "chartreuse4")
MetricShape <- c("Accuracy" = 16, "Sensitivity" = 1, "Recall" = 10, "Specificity" = 6, "Precision" = 8, "F1" = 19)


# Result Table

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
                                   Threshold := Thresh][, Compartment := Names[XComp]]
          
          return(Comp_Test)
        }
      })
      )
      
      other <- merge(other, data.table(Compartment = c("B_FISH", "B_INV", "B_DIA"),
                                       Community = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
      
      Table <- merge(auc[,c("AUC_valTest", "Threshold", "Compartment")], other, by = c("Compartment", "Threshold"))
      
      return(Table)}
  })))
  
  
  if(length(Thresh_AUC)!=0){
    
  Thresh_AUC[, DegradationDirection := Catalog[Name_Origine == gsub(toremove, "", param), DegradationDirection]]
    
    if(nrow(Thresh_AUC[AUC_valTest >= 0.7 & F1 >= 0.7]) > 0){
      ifelse(length(unique(Thresh_AUC[AUC_valTest >= 0.7 & F1 >= 0.7, Community])) > 1,
             Resp_community <- substr(unique(
               Thresh_AUC[AUC_valTest >= 0.7 & F1 >= 0.7 & !Community %in% Thresh_AUC[AUC_valTest == max(AUC_valTest),Community], Community]
             ), start=1, stop=4),
             Resp_community <- c("")
      )
      Resp_community <- paste(Resp_community, collapse = " ")
      
      #Thresh_AUC <- Thresh_AUC[AUC_valTest >= 0.7 & F1 >= 0.7]
      
      Thresh_AUC[, THRintersect := Threshold[which.min(abs(Sensitivity - Specificity))], by = Compartment]
      Thresh_AUC[, Threshold_detect := Threshold[which.max(AUC_valTest)[1]], by = Compartment]
      
      Thresh_AUC[, Threshold_detect := ifelse(DegradationDirection == 1,min(Threshold),
                                       ifelse(DegradationDirection == -1,max(Threshold), mean(Threshold)))]
      
      Thresh_AUC[, THRintersect_detect := ifelse(DegradationDirection == 1,min(THRintersect),
                                          ifelse(DegradationDirection == -1,max(THRintersect), mean(THRintersect)))]
      
      Main_comm <- unique(substr(Thresh_AUC[AUC_valTest == max(AUC_valTest), Community], start = 1, stop = 4))
      Main_comm <- paste(Main_comm, collapse = " ")
      
      Thresh_AUC <- Thresh_AUC[, Resp_comm := Resp_community][, Name_Origine := gsub(toremove, "",param)][
        , Main_comm := Main_comm]#[AUC_valTest == max(AUC_valTest)]
      
      #cols <- c("Accuracy","Sensitivity", "Specificity", "Precision", "Recall","F1")
      #Thresh_AUC[, (cols) := as.list(apply(Thresh_AUC[AUC_valTest == max(AUC_valTest), cols, with = F],2,mean))]
         
    } else { Thresh_AUC[, c("THRintersect","Threshold_detect", "THRintersect_detect", "Resp_comm", 
                    "Main_comm") := NA][, Name_Origine := gsub(toremove, "",param)]
    }
      return(Thresh_AUC)
  } 

}))

cols <- names(which(unlist(lapply(Table, is.numeric))))
Table <- Table[, (cols) := lapply(.SD, function(X) round(as.numeric(X), 2)),
               .SDcols = cols]
#setnames(Table, c("Community","Resp_comm"), c("Main_comm", "Responsive_comm"))
Table <- merge(Table, Titles[,.(Description, Name_Origine, Unit, Category, Scale)])

write.table(Table, file = "FinalTable.csv", sep = ";", row.names = F)
save(Table, file = "FinalTable")
###

#####
Table <- fread("FinalTable.csv")
files <- unique(paste0(toremove, Table[AUC_valTest >= 0.7 & F1 >= 0.7, Name_Origine]))

Threshold_plot <- lapply(files, function(param){
  
  print(param)
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_AUC <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]])})
  
  if(any(unlist(lapply(Thresh_AUC, is.data.frame)))){
    
    Thresh_AUC <- Thresh_AUC[unlist(lapply(Thresh_AUC, is.data.frame))]
    
    Thresh_names <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]][,"Compartment"])})
    Thresh_thresholds <- lapply(Thresh_draw, function(X){return(unique(X[["Thresh_AUCval"]]["Threshold"]))})
    Thresh_other <- lapply(Thresh_draw, function(X){return(X[["Thresh_OtherMetrics"]])})
    
    category <- Cat[Names == param, Category]
    ifelse(which(Cat[Category == category, Names] == param) %in% seq(1, length(Cat[Category == category, Names]), by = 2),
           Ytitle <- c("AUC", "density","Other Metrics") , Ytitle <- c("", "", ""))
    
    IsShp <- IsCol <- "none"
    if(which(Cat[Category == category, Names] == param) %in% 
       seq(min(3, length(Cat[Category == category, Names])), length(Cat[Category == category, Names]), by = min(4, length(Cat[Category == category, Names])))){
      IsCol <- "none"#IsCol <- guide_legend(title.position="top", title.hjust = 0.5)
      #sapply(sapply(seq(3, 14, by = 4), function(x) return(seq(x,length.out = 2))),c),
    }
    
    if(which(Cat[Category == category, Names] == param) %in% 
       seq(min(4, length(Cat[Category == category, Names])), length(Cat[Category == category, Names]), by = min(4, length(Cat[Category == category, Names])))){
      IsShp <- "none"#,IsShp <- guide_legend(title.position="top", title.hjust = 0.5))
    }
    
    if(length(Thresh_AUC) > 0){
      
    # For AUC plot
      Thresh_AUC <- as.data.table(do.call(rbind,lapply(Thresh_AUC, function(x) {
        
        x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
        
        x <- as.data.table(x)[, grep(c("AUC_val|Low|Up"), colnames(x), value = T) := 
                                lapply(.SD, function(col){return(round(as.numeric(col), 2))}), 
                              .SDcols = grep(c("AUC_val|Low|Up"), colnames(x), value = T)]
        
        x <- merge(x, data.table(Compartment = c("B_FISH", "B_INV", "B_DIA"),
                                 Community = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
        return(x)})))[
          , Community := factor(Community, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
            , N := length(which(!is.na(AUC_valTest))), by = Community]
      
      Thresh_AUC[, AlphaTest := "S"][AUC_valTest < 0.7, AlphaTest := "NS"][
        , AlphaTrain := "S"][AUC_valTrain < 0.7, AlphaTrain := "NS"]
      
      Thresh_AUC_max <- unique(Thresh_AUC[, c("AUCmaxTest","THRmaxTest","AUCmaxTrain","THRmaxTrain") 
                                       := list(max(AUC_valTest, na.rm = T),Threshold[which.max(AUC_valTest)],
                                               max(AUC_valTrain, na.rm = T),Threshold[which.max(AUC_valTrain)]), by = Community][
                                                 ,.(Community, AUCmaxTest, THRmaxTest, AUCmaxTrain, THRmaxTrain)])[,Param := param]
      Thresh_AUC_max[Thresh_AUC_max==-Inf] <- NA
      
      if(!is.na(Catalog[Name_Origine == gsub(toremove, "", param), LittThreshold]) && 
         as.numeric(Catalog[Name_Origine == gsub(toremove, "", param), LittThreshold]) > 
         quantile(RFD_all[, gsub(toremove, "", param), with = F], probs = 0.975, na.rm = T) && 
         !any(Thresh_AUC_max$THRmaxTest == as.numeric(Catalog[Name_Origine == gsub(toremove, "", param), LittThreshold]))){
        Thresh_AUC <- Thresh_AUC[Threshold != as.numeric(Catalog[Name_Origine == gsub(toremove, "", param), LittThreshold])]
      }
      
      ifelse(any(nchar(Thresh_AUC$Threshold) > 4) & diff(range(as.numeric(Thresh_AUC$Threshold), na.rm = T))/10 < 0.01, 
             makeBreaks <- T, makeBreaks <- F)
      Thresh_AUC[, Threshold := as.numeric(Threshold)]
      tonum <- c("THRmaxTrain", "THRmaxTest")
      Thresh_AUC_max[, (tonum) := lapply(.SD, as.numeric), .SDcols = tonum]
      
      # For other metrics
      Thresh_other <- as.data.table(do.call(rbind,lapply(1:length(Thresh_other), function(Xcomp) {
        
        Comp <- Thresh_other[[Xcomp]]
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
                                                   Community = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
          return(Comp_Test)
        } 
      })))
      
      Thresh_other[, Community := factor(Community, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
        , N := length(which(!is.na(Med))), by = c("Community", "Metric")]
      
      ReturnTable_other <- unique(dcast(Thresh_other[!is.na(Threshold)], 
                                  Community + Threshold ~ Metric, value.var = "Med"))
      ReturnTable_other$Param <- gsub(toremove,"", param)
      
      Thresh_other_max <- unique(dcast(Thresh_other[!is.na(Threshold)], 
                                 Community + Threshold ~ Metric, value.var = "Med")[
                                   , c("THRintersect", "AccuIntersect") := list(Threshold[which.min(abs(Sensitivity - Specificity))],
                                                                                Accuracy[which.min(abs(Sensitivity - Specificity))]),
                                   by = Community][,.(Community, THRintersect, AccuIntersect)])
      
      Thresh_other_max[THRintersect ==-Inf] <- NA
      
      Thresh_other[, Threshold := as.numeric(Threshold)]
      tonum <- c("THRintersect", "AccuIntersect")
      Thresh_other_max[, (tonum) := lapply(.SD, as.numeric), .SDcols = tonum]
      
      Thresh_other <- Thresh_other[Metric %in% c("Specificity", "Recall", "F1")]
      
      # Plot AUC
      plot_AUC <- ggplot(data = Thresh_AUC, aes(x = Threshold, y = AUC_valTest, color = Community, fill = Community)) + 
        geom_hline(yintercept = c(0.7, 0.8), color = "darkgrey", lty = 1) + 
        geom_segment(data = Thresh_AUC_max, aes(x = THRmaxTest, xend = THRmaxTest, y = 0, yend = AUCmaxTest), lty = 1, linewidth = 0.7) +
        ggtitle(Titles[Name_Origine == gsub(toremove, "", param), Description]) + 
        labs(subtitle = Titles[Name_Origine == gsub(toremove, "", param), Scale]) + 
        theme_classic() + facet_grid(. ~ Community) + labs(x = "",y = Ytitle[1]) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
              axis.text.x = element_text(angle = 45, vjust = 0.6), axis.text=element_text(size=8),
              axis.title=element_text(size=10), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 9)) +
        scale_color_manual(values = CommunityCols)
      
      if(makeBreaks){plot_AUC <- plot_AUC + 
        scale_x_continuous(breaks = seq(min(Thresh_AUC$Threshold, na.rm = T), max(Thresh_AUC$Threshold, na.rm = T), length.out = 10),
                           labels = format(seq(min(Thresh_AUC$Threshold, na.rm = T), max(Thresh_AUC$Threshold, na.rm = T),
                                               length.out = 10), digits = 2, scientific = T))
      }
      
      if(any(Thresh_AUC$N > 1)){
        plot_AUC <- plot_AUC + geom_ribbon(aes(ymax = UpTest, ymin = LowTest), colour = NA, alpha = 0.1) + 
          geom_point(data = Thresh_AUC[N > 1], aes(x = Threshold, y = AUC_valTest, alpha = AlphaTest)) + 
          geom_ribbon(aes(ymax = UpTrain, ymin = LowTrain), colour = NA, alpha = 0.1) + 
          geom_point(data = Thresh_AUC[N > 1], aes(x = Threshold, y = AUC_valTrain, alpha = AlphaTrain), shape = 24, size = 0.6) +
          scale_alpha_discrete(range = c(0.35,1)) +
          scale_fill_manual(values = CommunityCols)
      } 
      
      if(any(Thresh_AUC$N == 1)) {
        plot_AUC <- plot_AUC + geom_point(data = Thresh_AUC[N == 1], aes(x = Threshold, y = AUC_valTest, alpha = AlphaTest)) + 
          geom_point(data = Thresh_AUC[N == 1], aes(x = Threshold, y = AUC_valTrain, alpha = AlphaTrain), shape = 24) + 
          geom_errorbar(data = Thresh_AUC[N == 1],aes(ymax = UpTest, ymin = LowTest), width = 0.2) + 
          geom_errorbar(data = Thresh_AUC[N == 1],aes(ymax = UpTrain, ymin = LowTrain), width = 0.2) +
          scale_alpha_discrete(range = c(0.35,1)) +
          scale_fill_manual(values = CommunityCols)
      }
      
      if(!is.na(Catalog[Name_Origine == gsub(toremove, "", param), LittThreshold])){
        if(as.numeric(Catalog[Name_Origine == gsub(toremove, "", param), LittThreshold]) > 
           max(Thresh_AUC$Threshold, na.rm = T) + (max(Thresh_AUC$Threshold, na.rm = T)-min(Thresh_AUC$Threshold, na.rm = T))/2){
          plot_AUC <- plot_AUC + labs(tag = "Litterature threshold\nout of range (higher)") +
            theme(plot.tag = element_text(hjust = 0.5, color = "firebrick", size = 9), 
                  plot.tag.position = c(0.85, 0.99))
        } else {plot_AUC <- plot_AUC + 
          geom_vline(aes(xintercept = as.numeric(Catalog[Name_Origine == gsub(toremove, "", param), LittThreshold])), 
                     color = "firebrick", lty = 2, linewidth = 0.9, alpha = .6)}
        
      }
      
      # Plot other
      plot_other <- ggplot(data = Thresh_other, aes(x = Threshold, y = Med, color = Community, fill = Community)) + 
        geom_segment(data = Thresh_other_max, aes(x = THRintersect, xend = THRintersect, y = -0.05, yend = AccuIntersect), lty = 2, linewidth = 0.7) +
        theme_classic() + facet_grid(. ~ Community) + labs(x = "",y = Ytitle[3], color = "Community") + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), axis.text.x = element_text(angle = 45, vjust = 0.6),
              axis.text=element_text(size=8), axis.title=element_text(size=10), plot.title = element_text(size = 10), 
              plot.subtitle = element_text(size = 9), strip.background = element_blank(),
              strip.text.x = element_blank(), legend.position = "none") +
        guides(colour = IsCol, shape = IsShp, fill = "none") + 
        scale_color_manual(values = CommunityCols) + scale_fill_manual(values = CommunityCols) 
      
      if(makeBreaks){ plot_other <- plot_other + 
        scale_x_continuous(breaks = seq(min(Thresh_AUC$Threshold, na.rm = T), max(Thresh_AUC$Threshold, na.rm = T), length.out = 10),
                           labels = format(seq(min(Thresh_AUC$Threshold, na.rm = T), max(Thresh_AUC$Threshold, na.rm = T), length.out = 10),
                                           digits = 2, scientific = T))
      }
      
      if(any(Thresh_other$N > 1)){
        plot_other <- plot_other + geom_ribbon(aes(ymax = High, ymin = Low, group = Metric), colour = NA, alpha = 0.1) + 
          geom_point(data = Thresh_other[N > 1], aes(x = Threshold, y = Med, group = Metric, shape = Metric)) + 
          geom_line(data = Thresh_other[N > 1], aes(x = Threshold, y = Med, group = Metric)) + 
          scale_shape_manual(values = MetricShape, guide = IsShp) + labs(shape = "Model Metrics")
      } 
      
      if(any(Thresh_other$N == 1)) {
        plot_other <- plot_other + geom_point(data = Thresh_other[N == 1], aes(x = Threshold, y = Med, group = Metric)) + 
          geom_point(data = Thresh_other[N == 1], aes(x = Threshold, y = Med, group = Metric, shape = Metric)) + 
          geom_errorbar(data = Thresh_other[N == 1],aes(ymax = High, ymin = Low), width = 0.2) + 
          scale_fill_manual(values = CommunityCols) + scale_shape_manual(values = MetricShape)
      }
      
      ReturnTable <- merge(Thresh_AUC_max, Thresh_other_max)[, Thresh_diff := THRintersect - THRmaxTest][,
                     c("Rect_min", "Rect_max") := list(min(THRintersect, THRmaxTest), max(THRintersect, THRmaxTest)), by = Community]
      
      if(any(ReturnTable$Thresh_diff != 0)){
        plot_other <- plot_other + 
          geom_rect(data = ReturnTable, inherit.aes = FALSE, aes(xmin = Rect_min, xmax = Rect_max, ymin = -0.05, ymax = AccuIntersect),
                    alpha = 0.1) 
        
        plot_AUC <- plot_AUC + 
          geom_rect(data = ReturnTable, inherit.aes = FALSE, alpha = 0.1,  
                    aes(xmin = Rect_min, xmax = Rect_max,ymin = 0, ymax = AUCmaxTest))
      }
      
      
    } else { 
      plot_AUC <- ggplot(data = data.frame(Thresh = 0, AUC = 0, 
                                       Community = factor(c("Fish", "Macroinvertebrate", "Diatom"),
                                                                levels = c("Fish", "Macroinvertebrate", "Diatom"))), 
                     aes(x = Thresh, y = AUC)) + 
        ggtitle(gsub("^[^_]*_|_.*", "", gsub(toremove, "", param))) + theme_classic() + facet_grid(. ~ Community) + 
        labs(x = "", y = Ytitle) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
              axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
              axis.title=element_text(size=10), plot.title = element_text(size = 10))
      
      plot_other <- ggplot(data = data.frame(Thresh = 0, AUC = 0, 
                                       Community = factor(c("Fish", "Macroinvertebrate", "Diatom"),
                                                                levels = c("Fish", "Macroinvertebrate", "Diatom"))), 
                     aes(x = Thresh, y = AUC)) + 
       theme_classic() + facet_grid(. ~ Community) + 
        labs(x = "", y = Ytitle[3]) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
              axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
              axis.title=element_text(size=10), plot.title = element_text(size = 10))
      
      
      ReturnTable <- NA
    }
    
    Dens_plot <- RFD_all[, gsub(toremove, "", param), with = F]
    Dens_plot <- Dens_plot[complete.cases(Dens_plot)][
      get(gsub(toremove, "", param))>= quantile(Dens_plot, probs = 0.025, na.rm = T) &
        get(gsub(toremove, "", param))<= quantile(Dens_plot, probs = 0.975, na.rm = T)]
    
    Pdens <- ggplot(Dens_plot, aes(x = !!sym(gsub(toremove, "", param)))) + geom_density()+
      theme_classic() +  
      labs(x = paste0("(",Catalog[Name_Origine == gsub(toremove, "", param), Unit],") "), y = Ytitle[2]) + 
      geom_segment(data = Thresh_AUC_max, 
                   position = position_dodge2(width = (max(density(Dens_plot[[1]])$x)-min(density(Dens_plot[[1]])$x))/80),
                   aes(x = THRmaxTest, xend = THRmaxTest, color = Community, 
                       y = 0, yend = max(density(Dens_plot[[1]])$y)), lty = 1, linewidth = 1) +
      theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
      scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1)) +
      scale_color_manual(values = CommunityCols)
    
    Plot <- arrangeGrob(plot_AUC, plot_other, Pdens,nrow = 3, heights = c(3.5,2.5,1))
    
    return(list(Plot, ReturnTable))
    
  }
  })

names(Threshold_plot) <- files

#save(Threshold_plot, file = "ThresholdPlots_Final_2411")
#####


#####
Auclim <- "7"
Titles <- fread("SumUp_Seuil.csv")

Threshold_plot <- Threshold_plot[which(!unlist(lapply(Threshold_plot, is.null)))]
Threshold_plot <- lapply(Threshold_plot, function(param) {return(param[[1]])})

Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

Legd <- arrangeGrob(rasterGrob(image_read("Legend_allPlots.png")))

SplitPlotCategory <- split(names(Threshold_plot), ceiling(1:length(Threshold_plot)/4))
  
  
lapply(1:length(SplitPlotCategory), function(j){
    
    print(j)
    
    tg <- text_grob(paste0("AUC > 0.",Auclim), size = 12)
    
    #sg <- text_grob("50 sites limit  95% of data  5% step", size = 10)
    sg <- text_grob("", size = 10)
    
    lt <- list(tg, sg)
    heights <- do.call(unit.c, lapply(lt, function(.g) 1.5*grobHeight(.g)))
    titles <- gtable::gtable_matrix('title', 
                                    grobs = matrix(lt, ncol=1), 
                                    widths = unit(2,'npc'),
                                    heights = heights)
    Toplot <- Threshold_plot[SplitPlotCategory[[j]]]
    #for(i in 1:length(Toplot)){
      #Toplot[[i]][["grobs"]][[1]][[1]][[1]] <- grid.arrange(top = textGrob(Cat[Names == names(Toplot)[i], Category], gp=grid::gpar(fontsize=11)),
                                                           # Toplot[[i]][["grobs"]][[1]][[1]][[1]])}
    
  Toplot[[5]] <- Legd
    gridedPlots <-  grid.arrange(grobs = Toplot,layout_matrix = rbind(c(1, 2), c(3, 4), c(5,5)),heights = c(1,1,0.2)) #, ncol = 2, nrow = 2)
    #gridePlots <- grid.arrange(gridedPlots, Legd, nrow = 2, heights = c(1,0.3))
    
    png(paste0("Best0", Auclim,"_",gsub("AUC_threshold_2405_Q95_5stp_", "", toremove), j,".png"), width = 900, height = 1000)
    
    grid.arrange(tg, sg, gridedPlots, 
                 heights = unit.c(unit.c(grobHeight(tg) + 1.2*unit(0.5, "line"), 
                                         grobHeight(sg) + unit(0.5, "line"), 
                                         unit(1,"null"))))
    
    dev.off()
    
  })
  
## Drawing single graph
Var <- Catalog[grep("PCB153", Name_Origine),]
grid.arrange(Threshold_plot[[grep("PCB153", names(Threshold_plot), value = T)]])


#####

F1score_plot <- lapply(files, function(param){
  
  print(param)
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_names <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]][,"Compartment"])})
  Thresh_thresholds <- lapply(Thresh_draw, function(X){return(unique(X[["Thresh_AUCval"]]["Threshold"]))})
  Thresh_draw <- lapply(Thresh_draw, function(X){return(X[["Thresh_OtherMetrics"]])})
  
  if(any(unlist(lapply(Thresh_draw, is.list)))){
    
    Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.list))]
    Thresh_names <- Thresh_names[unlist(lapply(Thresh_draw, is.list))]
    Thresh_thresholds <- Thresh_thresholds[unlist(lapply(Thresh_draw, is.list))]
    
    category <- Cat[Names == param, Category]
    ifelse(which(Cat[Category == category, Names] == param) %in% c(1,4),
           Ytitle <- c("", "density"), Ytitle <- c("", ""))
    
    ifelse(which(Cat[Category == category, Names] == param) %in% c(3,6),
           Islegd <- "none", Islegd <- "none")
    
    if(length(Thresh_plot) > 0){
      
      Thresh_plot <- as.data.table(do.call(rbind,lapply(1:length(Thresh_plot), function(Xcomp) {
        
        Comp <- Thresh_plot[[Xcomp]]
        Names <- Thresh_names[[Xcomp]]
        Thresh <- Thresh_thresholds[[Xcomp]]
        
        if(length(Comp) !=0){
          Comp_Test <- as.data.table(do.call(rbind, lapply(1:length(Comp), function(xcomp){
            
            x <- Comp[[xcomp]]
            name <- Names[[xcomp]]
            
            #ret[,F1 := 2*((Accuracy * Sensitivity)/(Precision + Sensitivity))]
            #[, Precision := TP / (TP + FP)][, Recall := TP / (TP + FN)][ , F1 := 2*(Precision * Recall)/(Precision + Recall)]#
            
            #ret <-  as.data.table(t(ret[, quantile(F1, prob = c(0.05, .5, 0.95), na.rm = T)]))
            
            #ret <- as.data.table(t(ret[, lapply(.SD, quantile, prob = c(0.05, .5, 0.95), na.rm = T), .SDcols = c("Precision", "Recall", "F1")]), keep.rownames = T)
            #names(ret) <- c("Metric", "Low","Med","High")
            #names(ret) <- c("Low","Med","High")
            
            ret[, Compartment := name][, Threshold := Thresh]
            
            return(unique(ret[,.(Compartment, Threshold, Low, Med, High)]))
          })
          ))
          Comp_Test <- merge(Comp_Test, data.table(Compartment = c("B_FISH", "B_INV", "B_DIA"),
                                                   Community = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
          return(Comp_Test)
        } 
      })))
    }
    
    if(length(Thresh_plot) != 0){
      
      Thresh_plot[, Community := factor(Community, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
        , N := length(which(!is.na(Med))), by = Community]
      
      #Thresh_max <- dcast(Thresh_plot[!is.na(Threshold)], Community + Threshold ~ Metric, value.var = "Med")
      #Thresh_max[is.na(Thresh_max)] <- 0
      #Thresh_max <- unique(Thresh_max[, c("THRintersect", "AccuIntersect") := 
      #             list(Threshold[which.min(abs(Precision - Recall))], F1[which.min(abs(Precision - Recall))]),
      #               by = Community][,.(Community, THRintersect, AccuIntersect)])
      
      Thresh_max <- unique(Thresh_plot[!is.na(Threshold)][, c("AccuIntersect","THRintersect") 
                                                          := list(max(Med, na.rm = T),Threshold[which.max(Med)]), by = Community][
                                                            ,.(Community, AccuIntersect, THRintersect)])[,Param := param]
      
      Thresh_max[THRintersect ==-Inf] <- NA
      
      ifelse(any(nchar(Thresh_plot$Threshold) > 4) & diff(range(as.numeric(Thresh_plot$Threshold), na.rm = T))/10 < 0.01, 
             makeBreaks <- T, makeBreaks <- F)
      Thresh_plot[, Threshold := as.numeric(Threshold)]
      tonum <- c("THRintersect", "AccuIntersect")
      Thresh_max[, (tonum) := lapply(.SD, as.numeric), .SDcols = tonum]
      
      plot <- ggplot(data = Thresh_plot, aes(x = Threshold, y = Med, color = Community, fill = Community)) + 
        geom_segment(data = Thresh_max, aes(x = THRintersect, xend = THRintersect, y = -0.05, yend = AccuIntersect), lty = 1, linewidth = 0.7) +
        # ggtitle(Titles[Name_Origine == gsub(toremove, "", param), Description]) + 
        #labs(subtitle = Titles[Name_Origine == gsub(toremove, "", param), Scale]) + 
        theme_classic() + facet_grid(. ~ Community) + labs(x = "F1-score",y = Ytitle[1]) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = "none", 
              axis.text.x = element_text(angle = 45, vjust = 0.6), axis.text=element_text(size=8),
              axis.title=element_text(size=10), plot.title = element_text(size = 10), 
              plot.subtitle = element_text(size = 9), strip.background = element_blank(),
              strip.text.x = element_blank()) +
        scale_color_manual(values = CommunityCols)
      
      if(makeBreaks){ plot <- plot + 
        scale_x_continuous(breaks = seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T), length.out = 10),
                           labels = format(seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T), length.out = 10),
                                           digits = 2, scientific = T))
      }
      
      if(any(Thresh_plot$N > 1)){
        plot <- plot + geom_ribbon(aes(ymax = High, ymin = Low), colour = NA, alpha = 0.1) + 
          geom_point(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med)) + 
          geom_line(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med)) + 
          scale_fill_manual(values = CommunityCols) +
          guides(color = "none", fill = "none")
      } 
      
      if(any(Thresh_plot$N == 1)) {
        plot <- plot + geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = Med)) + 
          geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = Med)) + 
          geom_errorbar(data = Thresh_plot[N == 1],aes(ymax = High, ymin = Low), width = 0.2) + 
          scale_fill_manual(values = CommunityCols) 
      }
      
      
    } else { 
      plot <- ggplot(data = data.frame(Thresh = 0, AUC = 0, 
                                       Community = factor(c("Fish", "Macroinvertebrate", "Diatom"),
                                                                levels = c("Fish", "Macroinvertebrate", "Diatom"))), 
                     aes(x = Thresh, y = AUC)) + 
        #ggtitle(gsub("^[^_]*_|_.*", "", gsub(toremove, "", param))) + 
        theme_classic() + facet_grid(. ~ Community) + 
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
      labs(x = paste0("(",Catalog[Name_Origine == gsub(toremove, "", param), Unit],") "), y = Ytitle[2]) +
      theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
      scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1)) +
      scale_color_manual(values = CommunityCols)
    
    Plot <- arrangeGrob(plot,Pdens, nrow = 2, heights = c(3.5,1))
    
    #grid.arrange(Plot)
    return(Plot)
    
  }
})

names(F1score_plot) <- files

save(F1score_plot, file = paste0("ThresholdPlots_list_F1score", gsub("AUC_threshold_", "", toremove)))

