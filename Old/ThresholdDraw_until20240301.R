invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid", "ggradar","scales", "magick"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("../MetricsCatalogue.csv")[which(TokeepThreshold)]
Titles <- fread("SumUp_Seuil.csv")
Titles[, Reference_Threshold := as.numeric(gsub(",",".", Reference_Threshold))]

files <- list.files("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils/Genomig_Output")

load("../HYMOBIO_FULLDATA_202405.RData")

#tomatch <-
#tomatch <- gsub(tomatch, "", grep(tomatch, files, value = T))
toremove <- "AUC_threshold_2405_Quant95_5step_"#"AUC_threshold_2405_Q95_5stp_SMOTE"

files <- grep(toremove, files, value = T)
#files <- grep(paste0(tomatch, collapse = "|"), files, value = T)
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

CompartmentCols <- c("Fish" = "cornflowerblue", "Macroinvertebrate" = "tan2", "Diatom" = "chartreuse4")

#####
Table <- fread("FinalTable2.csv")
files <- paste0(toremove, Table$Name_Origine)
Threshold_plot <- lapply(files, function(param){
  
  print(param)
  
  version <- "new"
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  if(version == "new"){Thresh_draw <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]])})}
  
  if(any(unlist(lapply(Thresh_draw, is.data.frame)))){
    
    Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
    
    category <- Cat[Names == param, Category]
    
    ifelse(which(Cat[Category == category, Names] == param) %in% seq(1, length(Cat[Category == category, Names]), by = 2),
           Ytitle <- c("AUC", "density"), Ytitle <- c("", ""))
    
    if(length(Thresh_plot) > 0){
      
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
        theme_classic() + facet_grid(. ~ CompartmentName) + labs(x = "",y = Ytitle[1]) + 
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
      
      
    } else { 
      plot <- ggplot(data = data.frame(Thresh = 0, AUC = 0, 
                                       CompartmentName = factor(c("Fish", "Macroinvertebrate", "Diatom"),
                                                                levels = c("Fish", "Macroinvertebrate", "Diatom"))), 
                     aes(x = Thresh, y = AUC)) + 
        ggtitle(gsub("^[^_]*_|_.*", "", gsub(toremove, "", param))) + theme_classic() + facet_grid(. ~ CompartmentName) + 
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
                   aes(x = THRmaxTest, xend = THRmaxTest, color = CompartmentName, 
                       y = 0, yend = max(density(Dens_plot[[1]])$y)), lty = 1, linewidth = 1) +
      theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
      #annotate("rect", xmin = quantile(Dens_plot[[1]], 0.025), xmax = quantile(Dens_plot[[1]], 0.975), ymin = 0, ymax = max(density(Dens_plot[[1]])$y), alpha = .3,fill = "grey") +
      scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1)) +
      scale_color_manual(values = CompartmentCols)
    
    #plot <- plot + geom_vline(xintercept = quantile(Dens_plot[[1]], c(0.025, 0.975)), colour = "grey", linetype="dotted" ) +
    #annotate("rect", xmin = quantile(Dens_plot[[1]], 0.025), xmax = quantile(Dens_plot[[1]], 0.975), ymin = layer_scales(plot)$y$range$range[1], ymax = layer_scales(plot)$y$range$range[2], alpha = .1,fill = "grey") 
    
    Plot <- arrangeGrob(plot,Pdens, nrow = 2, heights = c(3.5,1))
    
    return(list(Plot, Thresh_max))
    
  }
  })

names(Threshold_plot) <- files

save(Threshold_plot, file = paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
#####


#####
MetricShape <- c("Accuracy" = 16, "Sensitivity" = 1, "Recall" = 10, "Specificity" = 6, "Precision" = 8, "F1" = 19)
OtherMetrics_plot <- lapply(files, function(param){
  
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
    ifelse(which(Cat[Category == category, Names] == param) %in% seq(1, length(Cat[Category == category, Names]), by = 2),
           Ytitle <- c("Other Metrics", "density"), Ytitle <- c("", ""))
    
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
            
    
    if(length(Thresh_plot) > 0){
      
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
    }
      
      if(length(Thresh_plot) != 0){
      
      Thresh_plot[, CompartmentName := factor(CompartmentName, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
            , N := length(which(!is.na(Med))), by = c("CompartmentName", "Metric")]
      
        ReturnTable <- unique(dcast(Thresh_plot[!is.na(Threshold)], 
                       CompartmentName + Threshold ~ Metric, value.var = "Med"))
        ReturnTable$Param <- gsub(toremove,"", param)
        
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
        # ggtitle(Titles[Name_Origine == gsub(toremove, "", param), Description]) + 
        #labs(subtitle = Titles[Name_Origine == gsub(toremove, "", param), Scale]) + 
        theme_classic() + facet_grid(. ~ CompartmentName) + labs(x = "",y = Ytitle[1], color = "Community") + 
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
      
      if(any(Thresh_plot$N > 1)){
        plot <- plot + geom_ribbon(aes(ymax = High, ymin = Low, group = Metric), colour = NA, alpha = 0.1) + 
          geom_point(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med, group = Metric, shape = Metric)) + 
          geom_line(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med, group = Metric)) + 
          scale_shape_manual(values = MetricShape, guide = IsShp) + labs(shape = "Model Metrics")
      } 
      
      if(any(Thresh_plot$N == 1)) {
        plot <- plot + geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = Med, group = Metric)) + 
          geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = Med, group = Metric, shape = Metric)) + 
          geom_errorbar(data = Thresh_plot[N == 1],aes(ymax = High, ymin = Low), width = 0.2) + 
          scale_fill_manual(values = CompartmentCols) + scale_shape_manual(values = MetricShape)
      }
      
      
    } else { 
      plot <- ggplot(data = data.frame(Thresh = 0, AUC = 0, 
                                       CompartmentName = factor(c("Fish", "Macroinvertebrate", "Diatom"),
                                                                levels = c("Fish", "Macroinvertebrate", "Diatom"))), 
                     aes(x = Thresh, y = AUC)) + 
        #ggtitle(gsub("^[^_]*_|_.*", "", gsub(toremove, "", param))) + 
        theme_classic() + facet_grid(. ~ CompartmentName) + 
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
      theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
       scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1)) +
      scale_color_manual(values = CompartmentCols)
    
     Plot <- arrangeGrob(plot,Pdens, nrow = 2, heights = c(3.5,1))
    
     #grid.arrange(Plot)#list(Plot, ReturnTable)
   return(Plot)
    
  }
})

names(OtherMetrics_plot) <- files

save(OtherMetrics_plot, file = paste0("ThresholdPlots_list_OtherMetrics", gsub("AUC_threshold_", "", toremove)))
#####

#####
load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
load(paste0("ThresholdPlots_list", gsub("Other_Metrics_", "", toremove)))

Titles <- fread("SumUp_Seuil.csv")
files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
files <- grep(toremove, files, value = T)

Threshold_plot <- Threshold_plot[which(!unlist(lapply(Threshold_plot, is.null)))]
Threshold_plot <- lapply(Threshold_plot, function(param) {return(param[[1]])})
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

Final <- lapply(intersect(names(OtherMetrics_plot), names(Threshold_plot)), function(param){
  return(arrangeGrob(Threshold_plot[[param]],
                     #F1score_plot[[param]][["grobs"]][[1]],
                     OtherMetrics_plot[[param]][["grobs"]][[1]], nrow = 2, heights = c(2,1)))
})
names(Final) <- intersect(names(OtherMetrics_plot), names(Threshold_plot))

Legd <- arrangeGrob(rasterGrob(image_read("Legend_allPlots.png")))

lapply(unique(Cat$Category), function(cat){
  print(cat)
  
  SplitPlotCategory <- split(Cat[Category == cat,Names], ceiling(1:nrow(Cat[Category == cat])/4))
  
  lapply(1:length(SplitPlotCategory), function(j){
    
    print(j)
    
    #Image <- grid.arrange(Threshold_cat[[j]])
    
    tg <- text_grob(cat, size = 12)
    
    sg <- text_grob("50 sites limit  95% of data  5% step", size = 10)
    
    lt <- list(tg, sg)
    heights <- do.call(unit.c, lapply(lt, function(.g) 1.5*grobHeight(.g)))
    titles <- gtable::gtable_matrix('title', 
                                    grobs = matrix(lt, ncol=1), 
                                    widths = unit(2,'npc'),
                                    heights = heights)
    Toplot <- Final[SplitPlotCategory[[j]]]
    Toplot[[5]] <- Legd
    gridedPlots <-  grid.arrange(grobs = Toplot,
                                 layout_matrix = rbind(c(1, 2), c(3, 4), c(5,5)),heights = c(1,1,0.2)) #, ncol = 2, nrow = 2)
    #gridePlots <- grid.arrange(gridedPlots, Legd, nrow = 2, heights = c(1,0.3))
    
    png(paste0("AllMetrics_",gsub("AUC_threshold_2405_Q95_5stp_", "", toremove), cat, j,".png"), width = 900, height = 1000)
    
    grid.arrange(tg, sg, gridedPlots, 
                 heights = unit.c(unit.c(grobHeight(tg) + 1.2*unit(0.5, "line"), 
                                         grobHeight(sg) + unit(0.5, "line"), 
                                         unit(1,"null"))))
    
    dev.off()
    
  })
  
})

#####

load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
load(paste0("ThresholdPlots_list_OtherMetrics", gsub("AUC_threshold_", "", toremove)))

Titles <- fread("SumUp_Seuil.csv")

Threshold_plot <- Threshold_plot[which(!unlist(lapply(Threshold_plot, is.null)))]
Threshold_plot <- lapply(Threshold_plot, function(param) {return(param[[1]])})

Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

Final <- lapply(intersect(names(OtherMetrics_plot), names(Threshold_plot)), function(param){
  return(arrangeGrob(Threshold_plot[[param]],
                     #F1score_plot[[param]][["grobs"]][[1]],
                     OtherMetrics_plot[[param]][["grobs"]][[1]],
                     nrow = 2, heights = c(2,1)))
})
names(Final) <- intersect(names(OtherMetrics_plot), names(Threshold_plot))

Legd <- arrangeGrob(rasterGrob(image_read("Legend_allPlots.png")))


SplitPlotCategory <- split(names(Final), ceiling(1:length(Final)/6))
  
  lapply(1:length(SplitPlotCategory), function(j){
    
    print(j)
    
    #Image <- grid.arrange(Threshold_cat[[j]])
    
    tg <- text_grob("AUC > 0.7 & F1 > 0.7", size = 12)
    
    sg <- text_grob("50 sites limit  95% of data  5% step", size = 10)
    
    lt <- list(tg, sg)
    heights <- do.call(unit.c, lapply(lt, function(.g) 1.5*grobHeight(.g)))
    titles <- gtable::gtable_matrix('title', 
                                    grobs = matrix(lt, ncol=1), 
                                    widths = unit(2,'npc'),
                                    heights = heights)
    Toplot <- Final[SplitPlotCategory[[j]]]
    for(i in 1:length(Toplot)){
      Toplot[[i]][["grobs"]][[1]][[1]][[1]] <- grid.arrange(top = textGrob(Cat[Names == names(Toplot)[i], Category], gp=grid::gpar(fontsize=11), vjust = 1),
                                                            Toplot[[i]][["grobs"]][[1]][[1]][[1]])}
    
  Toplot[[7]] <- Legd
    gridedPlots <-  grid.arrange(grobs = Toplot,
                                 layout_matrix = rbind(c(1, 2), c(3, 4), c(5,6), c(7,7)),heights = c(1,1,1,0.2)) #, ncol = 2, nrow = 2)
    #gridePlots <- grid.arrange(gridedPlots, Legd, nrow = 2, heights = c(1,0.3))
    
    png(paste0("Best07_",gsub("AUC_threshold_2405_Q95_5stp_", "", toremove), j,".png"), width = 900, height = 1500)
    
    grid.arrange(tg, sg, gridedPlots, 
                 heights = unit.c(unit.c(grobHeight(tg) + 1.2*unit(0.5, "line"), 
                                         grobHeight(sg) + unit(0.5, "line"), 
                                         unit(1,"null"))))
    
    dev.off()
    
  })
  


lapply(unique(Cat$Category), function(cat){
  print(cat)
  
  SplitPlotCategory <- split(Cat[Category == cat,Names], ceiling(1:nrow(Cat[Category == cat])/4))
  
  lapply(1:length(SplitPlotCategory), function(j){
    
    print(j)
    
    #Image <- grid.arrange(Threshold_cat[[j]])
    
    tg <- text_grob(cat, size = 12)
    
    sg <- text_grob("", size = 10)
    
    lt <- list(tg, sg)
    heights <- do.call(unit.c, lapply(lt, function(.g) 1.5*grobHeight(.g)))
    titles <- gtable::gtable_matrix('title', 
                                    grobs = matrix(lt, ncol=1), 
                                    widths = unit(2,'npc'),
                                    heights = heights)
    Toplot <- Final[SplitPlotCategory[[j]]]
    Toplot[[5]] <- Legd
    gridedPlots <-  grid.arrange(grobs = Toplot,
                                 layout_matrix = rbind(c(1, 2), c(3, 4), c(5,5)),heights = c(1,1,0.2)) #, ncol = 2, nrow = 2)
    #gridePlots <- grid.arrange(gridedPlots, Legd, nrow = 2, heights = c(1,0.3))
    
    png(paste0("AllMetrics_",gsub("AUC_threshold_2405_Q95_5stp_", "", toremove), cat, j,".png"), width = 900, height = 1000)
    
    grid.arrange(tg, sg, gridedPlots, 
                 heights = unit.c(unit.c(grobHeight(tg) + 1.2*unit(0.5, "line"), 
                                         grobHeight(sg) + unit(0.5, "line"), 
                                         unit(1,"null"))))
    
    dev.off()
    
  })
  
})


#Table
#####
toremove <- "AUC_threshold_2405_Quant95_5step_"
files <- list.files("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
files <- grep(toremove, files, value = T)

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
  return(Thresh_AUC[which.max(AUC_valTest)[1]][
    , Name_Origine := gsub(toremove, "",param)][, Compartment := NULL])
    }
  })
)
    
cols <- setdiff(names(Table),c("CompartmentName","Name_Origine"))
Table <- Table[, (cols) := lapply(.SD, function(X) round(as.numeric(X), 2)),
      .SDcols = cols]
setnames(Table, "CompartmentName", "Community")
Table <- merge(Table, Titles[,.(Description, Name_Origine, Unit, Category, Scale)])

Comm_formater <- formatter("span", style = x ~ style(color = ifelse(x == "Fish", "cornflowerblue", 
                                                                    ifelse(x == "Diatom", "green", "gold"))))
Cat_formater <- formatter("span", style = x ~ style(color = ifelse(x == "Natural Environment", "green", 
                          ifelse(x == "Physico-Chemistry", "yellowgreen", ifelse(x == "Temperature", "#DF536B",
                          ifelse(x == "Hydrology", "#00BFC4", ifelse(x == "Hydromorphology", "#619CFF",
                          ifelse(x == "Flow Barriers", "#F5C710", "brown"))))))))


TableForm <- formattable(Table[order(-AUC_valTest,-F1),.(Community, Category, Description, Scale,
                 Unit, Threshold, AUC_valTest, Accuracy, Sensitivity, Specificity,
                 Precision, Recall, F1)], 
                 list(Community = Comm_formater, Category = Cat_formater))

write.table(Table, file = "FinalTable.csv", sep = ";", row.names = F)

saveWidget(as.htmlwidget(TableForm, width = '40%'), file = file.path(getwd(), "FinalTable.html"))




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
                                                   CompartmentName = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
          return(Comp_Test)
        } 
      })))
    }
    
    if(length(Thresh_plot) != 0){
      
      Thresh_plot[, CompartmentName := factor(CompartmentName, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
        , N := length(which(!is.na(Med))), by = CompartmentName]
      
      #Thresh_max <- dcast(Thresh_plot[!is.na(Threshold)], CompartmentName + Threshold ~ Metric, value.var = "Med")
      #Thresh_max[is.na(Thresh_max)] <- 0
      #Thresh_max <- unique(Thresh_max[, c("THRintersect", "AccuIntersect") := 
      #             list(Threshold[which.min(abs(Precision - Recall))], F1[which.min(abs(Precision - Recall))]),
      #               by = CompartmentName][,.(CompartmentName, THRintersect, AccuIntersect)])
      
      Thresh_max <- unique(Thresh_plot[!is.na(Threshold)][, c("AccuIntersect","THRintersect") 
                                                          := list(max(Med, na.rm = T),Threshold[which.max(Med)]), by = CompartmentName][
                                                            ,.(CompartmentName, AccuIntersect, THRintersect)])[,Param := param]
      
      Thresh_max[THRintersect ==-Inf] <- NA
      
      ifelse(any(nchar(Thresh_plot$Threshold) > 4) & diff(range(as.numeric(Thresh_plot$Threshold), na.rm = T))/10 < 0.01, 
             makeBreaks <- T, makeBreaks <- F)
      Thresh_plot[, Threshold := as.numeric(Threshold)]
      tonum <- c("THRintersect", "AccuIntersect")
      Thresh_max[, (tonum) := lapply(.SD, as.numeric), .SDcols = tonum]
      
      plot <- ggplot(data = Thresh_plot, aes(x = Threshold, y = Med, color = CompartmentName, fill = CompartmentName)) + 
        geom_segment(data = Thresh_max, aes(x = THRintersect, xend = THRintersect, y = -0.05, yend = AccuIntersect), lty = 1, linewidth = 0.7) +
        # ggtitle(Titles[Name_Origine == gsub(toremove, "", param), Description]) + 
        #labs(subtitle = Titles[Name_Origine == gsub(toremove, "", param), Scale]) + 
        theme_classic() + facet_grid(. ~ CompartmentName) + labs(x = "F1-score",y = Ytitle[1]) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = "none", 
              axis.text.x = element_text(angle = 45, vjust = 0.6), axis.text=element_text(size=8),
              axis.title=element_text(size=10), plot.title = element_text(size = 10), 
              plot.subtitle = element_text(size = 9), strip.background = element_blank(),
              strip.text.x = element_blank()) +
        scale_color_manual(values = CompartmentCols)
      
      if(makeBreaks){ plot <- plot + 
        scale_x_continuous(breaks = seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T), length.out = 10),
                           labels = format(seq(min(Thresh_plot$Threshold, na.rm = T), max(Thresh_plot$Threshold, na.rm = T), length.out = 10),
                                           digits = 2, scientific = T))
      }
      
      if(any(Thresh_plot$N > 1)){
        plot <- plot + geom_ribbon(aes(ymax = High, ymin = Low), colour = NA, alpha = 0.1) + 
          geom_point(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med)) + 
          geom_line(data = Thresh_plot[N > 1], aes(x = Threshold, y = Med)) + 
          scale_fill_manual(values = CompartmentCols) +
          guides(color = "none", fill = "none")
      } 
      
      if(any(Thresh_plot$N == 1)) {
        plot <- plot + geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = Med)) + 
          geom_point(data = Thresh_plot[N == 1], aes(x = Threshold, y = Med)) + 
          geom_errorbar(data = Thresh_plot[N == 1],aes(ymax = High, ymin = Low), width = 0.2) + 
          scale_fill_manual(values = CompartmentCols) 
      }
      
      
    } else { 
      plot <- ggplot(data = data.frame(Thresh = 0, AUC = 0, 
                                       CompartmentName = factor(c("Fish", "Macroinvertebrate", "Diatom"),
                                                                levels = c("Fish", "Macroinvertebrate", "Diatom"))), 
                     aes(x = Thresh, y = AUC)) + 
        #ggtitle(gsub("^[^_]*_|_.*", "", gsub(toremove, "", param))) + 
        theme_classic() + facet_grid(. ~ CompartmentName) + 
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
      theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
      scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1)) +
      scale_color_manual(values = CompartmentCols)
    
    Plot <- arrangeGrob(plot,Pdens, nrow = 2, heights = c(3.5,1))
    
    #grid.arrange(Plot)
    return(Plot)
    
  }
})

names(F1score_plot) <- files

save(F1score_plot, file = paste0("ThresholdPlots_list_F1score", gsub("AUC_threshold_", "", toremove)))



# Printing plots
#####
load(paste0("ThresholdPlots_list", gsub("AUC_threshold_", "", toremove)))
Titles <- fread("SumUp_Seuil.csv")
files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
files <- grep(toremove, files, value = T)

Threshold_plot <- OtherMetrics_plot

Threshold_plot <- Threshold_plot[which(!unlist(lapply(Threshold_plot, is.null)))]
Threshold_plot <- lapply(Threshold_plot, function(param) {return(param[[1]])})
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

lapply(unique(Cat$Category), function(cat){
  print(cat)
  
  Threshold_cat <- Threshold_plot[Cat[Category == cat, Names]]
  
  SplitPlotCategory <- split(1:nrow(Cat[Category == cat]), ceiling(1:nrow(Cat[Category == cat])/6))
  
  lapply(1:length(SplitPlotCategory), function(j){
    
    print(j)
    
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
    
    png(paste0("Test", cat, j,".png"), width = 900, height = 800,)
    
    grid.arrange(tg, sg, gridedPlots, 
                 heights = unit.c(unit.c(grobHeight(tg) + 1.2*unit(0.5, "line"), 
                                         grobHeight(sg) + unit(0.5, "line"), 
                                         unit(1,"null"))))
    
    dev.off()
    
  })
  
})

#####
