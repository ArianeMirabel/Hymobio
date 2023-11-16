invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]

files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
files <- files[which(gsub("AUC_threshold_", "", files) %in% Catalog$NameOrigin)]

legend <- NULL

Threshold_plot <- lapply(files, function(param){
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
 
  ifelse(which(files == param) %in% seq(1, length(files), by = 3), Ytitle <- "AUC", Ytitle <- "")
  
  if(length(Thresh_plot) > 0){
    
  Thresh_plot <- as.data.table(do.call(rbind,lapply(Thresh_plot, function(x) {
    x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
    x$Compartment = c("FISH", "INV", "DIA"); return(x)})))[
    , Compartment := factor(Compartment, levels = c("FISH", "INV", "DIA"))][, N := .N, by = Compartment]
  
  Thresh_max <- unique(Thresh_plot[, AUCmax := Threshold[which.max(AUC_val)], by = Compartment][,.(Compartment, AUCmax)])
  
  plot <- ggplot(data = Thresh_plot, aes(x = Threshold, y = AUC_val, color = Compartment, fill = Compartment)) + 
    geom_hline(yintercept = 0.7, color = "darkgrey", lty = 2) + 
    geom_vline(data = Thresh_max, aes(xintercept = AUCmax, color = Compartment), lty = 2, show.legend = F) +
    ggtitle(gsub("^[^_]*_|_.*", "", gsub("AUC_threshold_", "", param))) + 
    theme_classic() + facet_grid(. ~ Compartment) + 
    labs(x = paste0("Threshold (",Catalog[NameOrigin == gsub("AUC_threshold_", "", param), Unit],") "),
        y = Ytitle) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
          axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
          axis.title=element_text(size=10))
    
  if(any(Thresh_plot$N > 1)){
  plot <- plot + geom_line(data = Thresh_plot[N > 1]) + 
    geom_ribbon(aes(ymax = Up, ymin = Low), color = "white", alpha = 0.1) 
    
  } 
  
  if(any(Thresh_plot$N == 1)) {
  plot <- plot + geom_point(data = Thresh_plot[N == 1]) + geom_errorbar(aes(ymax = Up, ymin = Low), width = 0.2) + 
    scale_x_continuous(breaks = unique(Thresh_plot$Threshold),
      labels = round(unique(Thresh_plot$Threshold), 1), limits = round(c(unique(Thresh_plot$Threshold)-1, unique(Thresh_plot$Threshold)+1)))
  }
      
  } else {
    
    plot <- ggplot(data = data.frame(Thresh = 0, AUC = 0, Compartment = c("FISH", "INV", "DIA")), 
                   aes(x = Thresh, y = AUC)) + 
      ggtitle(gsub("^[^_]*_|_.*", "", gsub("AUC_threshold_", "", param))) + theme_classic() + facet_grid(. ~ Compartment) + 
      labs(x = paste0("Threshold (",Catalog[NameOrigin == gsub("AUC_threshold_", "", param), Unit],") "),
           y = Ytitle) + 
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
            axis.text.x = element_blank(), axis.text.y = element_blank(), 
            axis.title=element_text(size=10))
  }
  
    return(plot)
})

names(Threshold_plot) <- files
Cat <- data.table(Names = files, Cat = gsub("_.*", "", gsub("AUC_threshold_", "", files)),
                  Param = gsub("^[^_]*_|_.*", "", gsub("AUC_threshold_", "", files)))
Cat[grep("BATHY_MAX", Names), Param := "BATHY MAX"][grep("BATHY_MIN", Names), Param := "BATHY MIN"]

lapply(unique(Cat$Cat), function(cat){
  
  Threshold_cat <- Threshold_plot[Cat[Cat == cat, Names]]
  
  SplitPlotCategory <- split(1:nrow(Cat[Cat == cat]), ceiling(1:nrow(Cat[Cat == cat])/6))

lapply(1:length(SplitPlotCategory), function(j){
  
  Image <- ggarrange(plotlist = Threshold_cat[SplitPlotCategory[[j]]], ncol = 3, nrow = 2, 
                     common.legend = T, legend = "right")
  
  png(paste0("Threshold_plot_", cat, j,".png"))
  
  print(annotate_figure(Image, top=text_grob(cat, size = 10)))
  
  dev.off()
  
})
})

  


