invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]
Titles <- fread("MetricsTitles.csv")

files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
files <- files[which(gsub("AUC_threshold_", "", files) %in% Catalog$NameOrigin)]

load("HYMOBIO_FULLDATA_20231129.RData")

Threshold_plot_withlitt <- lapply(files, function(param){
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
  
  ifelse(which(files == param) %in% seq(1, length(files), by = 3), Ytitle <- "AUC", Ytitle <- "")
  
  if(length(Thresh_plot) > 0){
    
  Thresh_plot <- as.data.table(do.call(rbind,lapply(Thresh_plot, function(x) {
    x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
    x <- as.data.table(x)[, c("AUC_val","Low","Up","Threshold") := lapply(.SD, function(col){return(round(as.numeric(col), 2))}), .SDcols = c("AUC_val","Low","Up","Threshold")]
    x$Compartment = c("FISH", "INV", "DIA"); return(x)})))[
    , Compartment := factor(Compartment, levels = c("FISH", "INV", "DIA"))][, N := length(which(!is.na(AUC_val))), by = Compartment]
  
  Thresh_max <- unique(Thresh_plot[, c("AUCmax","THRmax") := list(max(AUC_val, na.rm = T),Threshold[which.max(AUC_val)]), by = Compartment][
    ,.(Compartment, AUCmax, THRmax)])[,Param := param]
  Thresh_max[Thresh_max==-Inf] <- NA
  
  plot <- ggplot(data = Thresh_plot, aes(x = Threshold, y = AUC_val, color = Compartment, fill = Compartment)) + 
    geom_hline(yintercept = 0.7, color = "darkgrey", lty = 1) + 
    geom_segment(data = Thresh_max, aes(x = THRmax, xend = THRmax, y = 0, yend = AUCmax), lty = 1, linewidth = 0.7) +
    ggtitle(Titles[Name_Origine == gsub("AUC_threshold_", "", param), Description]) + 
    theme_classic() + facet_grid(. ~ Compartment) + labs(x = "",y = Ytitle) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
          axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
          axis.title=element_text(size=10), plot.title = element_text(size = 10))
    
  if(any(Thresh_plot$N > 1)){
  plot <- plot + geom_line(data = Thresh_plot[N > 1], aes(x = Threshold, y = AUC_val)) + 
    geom_ribbon(aes(ymax = Up, ymin = Low), color = "white", alpha = 0.1) 
  } 
  
  if(any(Thresh_plot$N == 1)) {
  plot <- plot + geom_point(data = Thresh_plot[N == 1]) + geom_errorbar(data = Thresh_plot[N == 1],aes(ymax = Up, ymin = Low), width = 0.2) 
  }
  
  if(!is.na(Catalog[NameOrigin == gsub("AUC_threshold_", "", param), LittThreshold])){
    plot <- plot + 
      geom_vline(aes(xintercept = as.numeric(Catalog[NameOrigin == gsub("AUC_threshold_", "", param), LittThreshold])), 
                 color = "royalblue", lty = 2, linewidth = 1.1, alpha = .6)
  }
  
  Dens_plot <- RFD_all[, gsub("AUC_threshold_", "", param), with = F]
  Dens_plot <- Dens_plot[complete.cases(Dens_plot)]
  
  Pdens <- ggplot(Dens_plot, aes(x = !!sym(gsub("AUC_threshold_", "", param)))) + geom_density()+
    theme_classic() + labs(x = paste0("(",Catalog[NameOrigin == gsub("AUC_threshold_", "", param), Unit],") "), y = "") + 
    theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.position = 'none', 
          axis.text.x = element_text(angle = 90), axis.text=element_text(size=8)) +
    annotate("rect", xmin = quantile(Dens_plot[[1]], 0.025), xmax = quantile(Dens_plot[[1]], 0.975),
             ymin = 0, ymax = max(density(Dens_plot[[1]])$y), alpha = .3,fill = "grey") +
  scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1))
  
  grid.arrange(grobs = list(plot,Pdens,  grid.rect(gp=gpar(col="white")), grid.rect(gp=gpar(col="white"))), nrow = 2, 
               layout_matrix = rbind(c(1,1,1),c(2,3,4)), heights = c(3,1))
  
  gb2 <- ggplotGrob(Pdens)
  gb1 <- ggplotGrob(plot)
  gb2$grobs[[13]]$widths <-  gb1$grobs[[13]]$widths
  
  ggdraw(gridExtra::gtable_rbind(g1,g2,g3))
  
  
  } else {
    return(NA)
  }
  
    return(list(plot, Thresh_max))
})

names(Threshold_plot_withlitt) <- files
Threshold_plot_withlitt <- Threshold_plot_withlitt[which(!unlist(lapply(Threshold_plot_withlitt, is.na)))]
Threshold_plot <- Threshold_plot_withlitt 

Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))

Threshold_plot <- lapply(Threshold_plot, function(param) {return(param[[1]])})
Cat <- left_join(data.table(Names = files, Cat = gsub("_.*", "", gsub("AUC_threshold_", "", files)),
                        Name_Origine = gsub("AUC_threshold_", "", files)),
             Titles)

lapply(unique(Cat$Cat), function(cat){
  
  Threshold_cat <- Threshold_plot[Cat[Cat == cat, Names]]
  
  SplitPlotCategory <- split(1:nrow(Cat[Cat == cat]), ceiling(1:nrow(Cat[Cat == cat])/6))

lapply(1:length(SplitPlotCategory), function(j){
  
  Image <- ggarrange(plotlist = Threshold_cat[SplitPlotCategory[[j]]], ncol = 3, nrow = 2, 
                     common.legend = T, legend = "right")
  
  png(paste0("Threshold_plot_withlitt", cat, j,".png"), width = 800, height = 400,)
  
  print(annotate_figure(Image, top=text_grob(unique(Cat[Cat == cat, Category]), size = 10)))
  
  dev.off()
  
})
})

# Plot difference in threshold
Threshold_max[, Abs_THRmax := THRmax/max(THRmax), by = Param][
  , Signif := 1][AUCmax < 0.7, Signif := 0][, Signif := factor(Signif, levels = c(0,1))]

Breaks  <- unique(as.data.table(merge(Threshold_max[, Name_Origine := gsub("AUC_threshold_", "", Param)], Titles))[
  , Brk := Param[1], by = Category][,c("Brk", "Category")])

ggplot(data = Threshold_max, aes(x = Param, y = Abs_THRmax, group = Compartment, color = Compartment)) + 
  geom_line() +
  scale_x_discrete(breaks = Breaks[["Brk"]], labels = Breaks[["Category"]]) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5))

Threshold_max[,Parammax := max(THRmax), by = Param][, ismax := 1][Parammax == THRmax, ismax := 0][
  , Ndiff := sum(ismax), by = Compartment]

ggplot(Threshold_max, aes(x="", y=Ndiff, fill=Compartment)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)

# Plot difference in AUC
AUCmeans <- Threshold_max[, round(mean(AUCmax, na.rm = T),2), by = "Compartment"]
setnames(AUCmeans, c("V1", "Compartment"), c("Av. AUC", "Comp."))
AUCmeans <- as.data.table(gsub("DIA", "Diatom", gsub("INV", "Macroinv.", gsub("FISH", "Fish", as.matrix(AUCmeans)))))

P <- ggplot(data = Threshold_max, aes(x = Param, y = AUCmax, group = Compartment, color = Compartment)) + 
  geom_line() +
  scale_x_discrete(breaks = Breaks[["Brk"]], labels = Breaks[["Category"]]) +
  theme_classic() + theme(axis.text.x = element_text(angle = 75, vjust = 0.5), legend.position = "none")


Pb <- ggplot_build(P)
scale_cols <- unique(Pb$data[[1]][,c("colour", "group")])
table2 <- merge(AUCmeans[, group := as.integer(rownames(AUCmeans))], scale_cols, by = "group", sort = FALSE)

t <- tableGrob(AUCmeans[,c("Comp.", "Av. AUC")], rows = NULL, theme = ttheme_default(
  core=list(bg_params = list(fill=c(table2$colour, rep("grey90", nrow(table2))),
                             alpha = .5)), base_size = 9))
grid.arrange(P, t, nrow = 1, ncol = 2, widths=c(2,.5))


AUCmeans <- merge(Threshold_max[, Name_Origine:= gsub("AUC_threshold_", "", Param)], Titles)
AUCmeans <- AUCmeans[, AUCmax := round(AUCmax,2)][complete.cases(AUCmeans)]
AUCmeans[, Compartment := gsub("DIA", "Diatom", gsub("INV", "Macroinv.", gsub("FISH", "Fish", Compartment)))]
AUCmeans[, Compartment := factor(Compartment, levels = c("Fish", "Macroinv.", "Diatom"))]

ggplot(AUCmeans, aes(x = Category, y = AUCmax, fill = Compartment)) + geom_hline(yintercept = 0.7, color = "darkgrey", lty = 1) +
  geom_boxplot() + labs(x = "Category", y = "Mean AUC", fill = "Compartment") +
  theme_classic()


