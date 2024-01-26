invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]
Titles <- fread("SumUp_Seuil.csv")

files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")

load("HYMOBIO_FULLDATA_20231129.RData")

#files <- grep(paste0(sub("AUC_threshold_50sites", "",grep("50sites", files, value = T)), collapse = "|"), files, value = T)
files <- grep(paste0(sub(".*50sites *(.*?) *_WD.*", "\\1", grep("_WD", files, value = T)), collapse = "|"),
              files, value = T)
files <- files[order(sub("50sites", "", files))]

files <- files[which(gsub("50sites", "",gsub("AUC_threshold_", "", files)) %in% Catalog$NameOrigin)]

#####
Threshold_plot <- lapply(files, function(param){
  
  N <- "50"
  if(!grepl("50", param)) {N <- "100"}
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  param <- gsub("50sites", "", param)
  
  Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
  
  Ytitle <- c("AUC", "density")
  #if(!which(files == param) %in% seq(1, length(files), by = 3){Ytitle <- c("","")}
  
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
    labs(subtitle = paste0(Titles[Name_Origine == gsub("AUC_threshold_", "", param), Scale], "\n", N, " sites limit & 10% step"))+ 
    theme_classic() + facet_grid(. ~ Compartment) + labs(x = "",y = Ytitle[1]) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
          axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
          axis.title=element_text(size=10), plot.title = element_text(size = 10), plot.subtitle = element_text(size = 9)) 
    
  if(any(Thresh_plot$N > 1)){
  plot <- plot + geom_ribbon(aes(ymax = Up, ymin = Low), color = "white", alpha = 0.1) + 
    geom_point(data = Thresh_plot[N > 1], aes(x = Threshold, y = AUC_val)) 
  } 
  
  if(any(Thresh_plot$N == 1)) {
  plot <- plot + geom_point(data = Thresh_plot[N == 1]) + geom_errorbar(data = Thresh_plot[N == 1],aes(ymax = Up, ymin = Low), width = 0.2) 
  }
  
  if(!is.na(Catalog[NameOrigin == gsub("AUC_threshold_", "", param), LittThreshold])){
    plot <- plot + 
      geom_vline(aes(xintercept = as.numeric(Catalog[NameOrigin == gsub("AUC_threshold_", "", param), LittThreshold])), 
                 color = "royalblue", lty = 2, linewidth = 1.1, alpha = .6)
  }
  
  
  } else { 
    plot <- ggplot(data = data.frame(Thresh = 0, AUC = 0, Compartment = c("FISH", "INV", "DIA")), 
                          aes(x = Thresh, y = AUC)) + 
    ggtitle(gsub("^[^_]*_|_.*", "", gsub("AUC_threshold_", "", param))) + theme_classic() + facet_grid(. ~ Compartment) + 
    labs(x = "", y = Ytitle) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
          axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
          axis.title=element_text(size=10), plot.title = element_text(size = 10))
    
    Thresh_max <- NA
  }
  
  Dens_plot <- RFD_all[, gsub("AUC_threshold_", "", param), with = F]
  Dens_plot <- Dens_plot[complete.cases(Dens_plot)]
  
  Pdens <- ggplot(Dens_plot, aes(x = !!sym(gsub("AUC_threshold_", "", param)))) + geom_density()+
    theme_classic() + 
    labs(x = paste0("(",Catalog[NameOrigin == gsub("AUC_threshold_", "", param), Unit],") "), y = Ytitle[2]) + 
    theme(plot.margin = unit(c(0,0.5,1,0.5), "lines"), legend.position = 'none', axis.text=element_text(size=8)) +
    annotate("rect", xmin = quantile(Dens_plot[[1]], 0.025), xmax = quantile(Dens_plot[[1]], 0.975),
             ymin = 0, ymax = max(density(Dens_plot[[1]])$y), alpha = .3,fill = "grey") +
    scale_y_continuous(breaks = round(c(min(density(Dens_plot[[1]])$y), max(density(Dens_plot[[1]])$y)), 1))
  
  plot <- plot + geom_vline(xintercept = quantile(Dens_plot[[1]], c(0.025, 0.975)),
                            colour = "grey", linetype="dotted" ) +
    annotate("rect", xmin = quantile(Dens_plot[[1]], 0.025), xmax = quantile(Dens_plot[[1]], 0.975),
             ymin = layer_scales(plot)$y$range$range[1], ymax = layer_scales(plot)$y$range$range[2],
             alpha = .1,fill = "grey") 
  
  Plot <- arrangeGrob(plot,Pdens, nrow = 2, heights = c(3.5,1))
  
  return(list(Plot, Thresh_max))
})

names(Threshold_plot) <- files

save(Threshold_plot, file = "ThresholdPlots_list")
#####

# Printing plots
#####
load("ThresholdPlots_list")
Threshold_plot <- lapply(Threshold_plot, function(param) {return(param[[1]])})
Cat <- left_join(data.table(Names = files, Cat = gsub("50sites" ,"", gsub("_.*", "", gsub("AUC_threshold_", "", files))),
                        Name_Origine = gsub("50sites" ,"", gsub("AUC_threshold_", "", files))),
             Titles[,.(Name_Origine, Description, Category)])

lapply(unique(Cat$Category), function(cat){
  
  Threshold_cat <- Threshold_plot[Cat[Category == cat, Names]]
  
  SplitPlotCategory <- split(1:nrow(Cat[Category == cat]), ceiling(1:nrow(Cat[Category == cat])/6))

lapply(1:length(SplitPlotCategory), function(j){
  
  Image <- grid.arrange(grobs = Threshold_cat[SplitPlotCategory[[j]]], ncol = 3, nrow = 2, 
                     common.legend = T, legend = "right")
  
  png(paste0("Threshold_plot", cat, j,".png"), width = 800, height = 800,)
  
  print(annotate_figure(Image, top=text_grob(cat, size = 10)))
  
  dev.off()
  
})
})


lapply(unique(Cat$Category), function(cat){
  
  Threshold_cat <- Threshold_plot[Cat[Category == cat, Names]]
  
  SplitPlotCategory <- split(1:nrow(Cat[Category == cat]), ceiling(1:nrow(Cat[Category == cat])/6))
  
lapply(1:length(SplitPlotCategory), function(j){
    
    #Image <- grid.arrange(Threshold_cat[[j]])
    
    Image <- grid.arrange(grobs = Threshold_cat[SplitPlotCategory[[j]]], ncol = 3, nrow = 2, common.legend = T, legend = "right")
    
    #png(paste0("Threshold_plot", cat, j,".png"), width = 800, height = 800,)
    
    print(annotate_figure(Image, top=text_grob(cat, size = 10)))
    
    #dev.off()
    
  })
  
})

grid.arrange(grobs = Threshold_plot, ncol = 2, nrow = 2, 
                      common.legend = T, legend = "right")


#####



# Plot difference in threshold
load("ThresholdPlots_list")
Titles <- fread("MetricsTitles.csv")
Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))

Threshold_max[, Abs_THRmax := max(THRmax), by = Param][, minTHRmax := min(THRmax), by = Param][
  , scaleTHRmax := scale(THRmax), by = Param][is.na(scaleTHRmax), scaleTHRmax := 0][
  , Signif := 1][AUCmax < 0.7, Signif := 0][, Signif := factor(Signif, levels = c(0,1))][
  , Nmax := length(which(THRmax==Abs_THRmax)), by = Param][, Nmin := length(which(THRmax==minTHRmax)), by = Param]

Threshold_max[, isTmax := 0][THRmax == Abs_THRmax & Nmax == 1, isTmax := 1][, OccMax := sum(isTmax), by = Compartment][
  , isTmin := 0][THRmax == minTHRmax & Nmin == 1, isTmin := 1][, OccMin := sum(isTmin), by = Compartment]

PpieTmax_t <- unique(Threshold_max[,.(OccMax, Compartment)][order(OccMax)])
PpieTmax_t[,Compartment := factor(Compartment, levels = PpieTmax_t[order(OccMax, decreasing = T),Compartment])][, prop := OccMax/sum(OccMax)][
  ,Ypos := cumsum(prop)- 0.5*prop]
PpieTmax <- ggplot(PpieTmax_t, aes(x="", y=prop, fill=Compartment)) +
  geom_bar(stat="identity", width=1) + coord_polar("y") + 
  coord_polar("y") + geom_text(aes(y = Ypos, label = OccMax), color = "white") +
  theme_void() + 
  scale_fill_manual(values = c("FISH" = "#F8766D", "INV" = "#00BA38", "DIA" = "#619CFF")) +
  theme(legend.position = "none",axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank(), axis.title = element_blank(), plot.title = element_text(size =9)) +
  ggtitle("(i) Occurence higest threshold") +
  geom_text(aes(y = Ypos, label = OccMax), color = "white")


PpieTmin_t <- unique(Threshold_max[,.(OccMin, Compartment)][order(OccMin)])
PpieTmin_t[,Compartment := factor(Compartment, levels = PpieTmin_t[order(OccMin, decreasing = T),Compartment])][, prop := OccMin/sum(OccMin)][
  ,Ypos := cumsum(prop)- 0.5*prop]
PpieTmin <- ggplot(PpieTmin_t, aes(x="", y=prop, fill=Compartment)) +
  geom_bar(stat="identity", width=1) + 
  coord_polar("y") + geom_text(aes(y = Ypos, label = OccMin), color = "white") +
  theme_void() + 
  scale_fill_manual(values = c("FISH" = "#F8766D", "INV" = "#00BA38", "DIA" = "#619CFF")) +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank(), axis.title = element_blank(), plot.title = element_text(size =9)) +
  ggtitle("(ii) Occurence lowest threshold")

grid.arrange(grobs = list(Prel, PpieTmax, PpieTmin), common.legend = T, legend = "bottom", 
             layout_matrix = rbind(c(1, 1, 1, 2),
                                   c(1, 1, 1, 3),c(1, 1, 1, NA)))

Threshold_max[,Abs_AUCmax := max(AUCmax), by = Param][, isAmax := 0][AUCmax == Abs_AUCmax, isAmax := 1][
  , NdiffA := sum(isAmax), by = Compartment]


Pauc <- ggplot(data = Threshold_max, aes(x = Param, y = AUCmax, group = Compartment)) + 
  geom_point(aes(colour = Compartment, shape = Signif), size = 2) +
  scale_shape_manual(values=c(1,19) , guide = "none") +
  scale_x_discrete(breaks = Breaks[["Brk"]], labels = Breaks[["Category"]]) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 75, vjust = 0.5)) + 
  labs(x = "Categories", y = "") + ggtitle("Maximum AUC on gradient.")


PpieA_t <- unique(Threshold_max[,.(NdiffA, Compartment)][order(NdiffA)])
PpieA_t[,Compartment := factor(Compartment, levels = PpieA_t[order(NdiffA, decreasing = T),Compartment])][, prop := NdiffA/sum(NdiffA)][
  ,Ypos := cumsum(prop)- 0.5*prop]
PpieA <- ggplot(PpieA_t, aes(x="", y=prop, fill=Compartment))  +
  geom_bar(stat="identity", width=1) + coord_polar("y")+ 
  geom_text(aes(y = Ypos, label = NdiffA), color = "white") +
  theme_void() + 
  scale_fill_manual(values = c("FISH" = "#F8766D", "INV" = "#00BA38", "DIA" = "#619CFF")) +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid  = element_blank(), axis.title = element_blank(), plot.title = element_text(size =9)) +
  ggtitle("(iii) Occurence highest AUC")

AUCmeans <- Threshold_max[, round(mean(AUCmax, na.rm = T),2), by = "Compartment"]
setnames(AUCmeans, c("V1", "Compartment"), c("Av. AUC", "Comp."))
AUCmeans <- as.data.table(gsub("DIA", "Diatom", gsub("INV", "Macroinv.", gsub("FISH", "Fish", as.matrix(AUCmeans)))))
Pb <- ggplot_build(PpieTmax)
scale_cols <- unique(Pb$data[[1]][,c("fill", "group")])
table2 <- merge(AUCmeans[, group := as.integer(rownames(AUCmeans))], scale_cols, by = "group", sort = FALSE)

t <- tableGrob(AUCmeans[,c("Comp.", "Av. AUC")], rows = NULL, theme = ttheme_default(
  core=list(bg_params = list(fill=c(table2$fill, rep("grey90", nrow(table2))),
                             alpha = .5)), base_size = 9))

grid.arrange(grobs = list(PpieTmax, PpieTmin, PpieA, t), common.legend = T, legend = "bottom", ncol = 2, nrow = 2)




# Plot by compartment
load("ThresholdPlots_list")
Titles <- fread("MetricsTitles.csv")
Threshold_max <- do.call(rbind,lapply(Threshold_plot, function(param) {if(length(param[[2]])!=1){return(param[[2]])}}))
Threshold_max <- as.data.table(merge(Threshold_max[, Name_Origine := gsub("AUC_threshold_", "", Param)], Titles))
Threshold_max[, Category := factor(Category, levels = c("Natural Env.", "Hydromorphology", "Hydrology", "Temperature",
                                                        "Anthropic Env.", "Obstacles"))]
Threshold_max <- Threshold_max[order(Category)]
CatColors <- c("#00BA38", "#00BFC4", "#619CFF","#DF536B", "#F5C710", "tomato4"); names(CatColors)<- levels(Threshold_max$Category)

Threshold_max[, Ncomp := length(which(AUCmax>=0.7)), by = Param][, NcompCat := sum(Ncomp), by = Category]

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
  geom_bar(stat = "identity", position = "dodge")+
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

