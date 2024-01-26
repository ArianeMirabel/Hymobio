invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel", "dplyr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]

load("AMOBIO_WP3_2_FISH_INV_DIA_ENTROPIE_20231123.Rdata")
RFD_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN)[, H_Ncrue_S1_ALL := NULL]

load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/TaxonomicBiodiv_Janv24_Diatom")
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/TaxonomicBiodiv_Diatom_SuppTaxo")
Index_Taxo_Diatom <- merge(Index_Taxo_Diatom, SuppTaxo_Diatom, by = c("Group", "CdStation", "ID_AMOBIO_START", "Year"))
save(Index_Taxo_Diatom, file = "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/TaxonomicBiodiv_Janv24_Diatom")


grep("B_FISH", colnames(RFD_all), value = T)
directory <- "C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/Indices/Genomig_Output/"
lapply(paste0(directory,grep("Janv24",list.files(directory), value =T)), load, .GlobalEnv)


load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Hydro/HydroIndex_36125all_AM_20232911")
RFD_all <- left_join(RFD_all, 
                         Hydrolaps[, c("ID_AMOBIO_START", "Samp_date", grep("Ncrues|Netiage", colnames(Hydrolaps), value = T)), with = F],
                         by = c("ID" = "ID_AMOBIO_START", "DATE_OPERATION" = "Samp_date"))
colnames(RFD_all) <- sub(" ", "_", colnames(RFD_all))
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/FinalMAJ_Naiades.Alric.Traits/SuppTaxo_Diatom")
RFD_all <- left_join(RFD_all, 
                     SuppTaxo_Diatom[, c("ID_AMOBIO_START", "Date_PrelBio", grep("B_DIA", colnames(SuppTaxo_Diatom), value = T)), with = F],
                     by = c("ID" = "ID_AMOBIO_START", "DATE_OPERATION" = "Date_PrelBio"))
RFD_all <- RFD_all[, grep(paste0("B_FISH|B_INV|B_DIA|",paste(Catalog$NameOrigin, collapse = "|")), colnames(RFD_all), value = T), with = F]
RFD_all <- RFD_all[, lapply(.SD, function(X) {X[X == "-Inf"] <- NA; return (X)})]
save(RFD_all, file = "HYMOBIO_FULLDATA_20240112.RData")
rm(list= c("MATRIX_AMOBIO_WP3_CLEAN", "Hydrolaps"));


Nslice <- 3

load("HYMOBIO_FULLDATA_20240112.RData")
#Params <- paste(Catalog$Category, Catalog$Name, sep = "_")
Params <- grep(paste(Catalog$NameOrigin, collapse = "|"), colnames(RFD_all), value = "T")
Params <- split(Params, ceiling(seq_along(Params)/(length(Params)%/%9)))[[Nslice]]

#Params <- setdiff(Params, gsub("AUC_threshold_", "",list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")))

for(param in Params){
  
  print(param)
  
    RFD_all <- RFD_all[, grep(paste0("B_DIA",param), colnames(RFD_all), value = T), with = F]
    Thresholds <- seq(from = range(RFD_all[,..param], na.rm = T)[1], to = range(RFD_all[,..param], na.rm = T)[2], length.out = 10)

if(!is.na(Catalog[NameOrigin == param,LittThreshold])){
  Thresholds[which.min(abs(Thresholds-as.numeric(Catalog[NameOrigin == param,LittThreshold])))] <- Catalog[NameOrigin == param,LittThreshold]
}

Thresh_draw <- lapply(Thresholds, function(Thresh){
         
         pb <- txtProgressBar(min = 0, max = 10, style = 3)
         setTxtProgressBar(pb, which(Thresholds == Thresh))
         
Thresh_AUC <- lapply(c("B_FISH", "B_INV", "B_DIA"), function(comp){
  
  Rfd_comp <- RFD_all[, grep(paste0(comp,"|State|",param), colnames(RFD_all), value = T), with = F]
  Rfd_comp <- Rfd_comp[complete.cases(Rfd_comp)]
  Rfd_comp <- Rfd_comp[get(param) > Thresh, State := "bad"][get(param) <= Thresh, State := "good"][
  , State := as.factor(State)]
  
 if(!any(c(nrow(Rfd_comp[State == "good"]), nrow(Rfd_comp[State == "bad"])) < 50) & uniqueN(Rfd_comp$State) > 1){
      
      set.seed(1)
      Rfd_comp[, Cross_grp := sample(x=1:10, nrow(Rfd_comp),replace=TRUE)]
     
      AUC_10Cross <- lapply(1:10, function(grp){
       
        Cross_train <- Rfd_comp[Cross_grp != grp]
        Cross_test <- Rfd_comp[Cross_grp == grp]
        
        if(uniqueN(Cross_train$State) > 1){
           
        rf_cross <- tryCatch({randomForest(State ~ .,
                                           data = Cross_train[, grep(paste0(comp,"|State"), colnames(Rfd_comp), value = T), with = F],
                                           scale.permutation.importance = TRUE,
                                           respect.unordered.factors = TRUE,
                                           mtry = 3)},
                             warning = function(w) {print(paste(w, "Cross val on threshold", Thresh, "for", param))},
                             error = function(e) {print(paste(e, "Cross val on threshold", Thresh, "for", param))})
        
        AUC_cross <- performance(prediction(
          as.numeric(predict(rf_cross, 
          Cross_test)), Cross_test$State),"auc")
        
        return(AUC_cross@y.values[[1]])
        
        } else {return(NA)}
        
      }) 
      
    return(data.frame(AUC_val = mean(do.call(c,AUC_10Cross), na.rm = T),
                      Low = quantile(do.call(c,AUC_10Cross), probs = 0.05, na.rm = T),
                      Up = quantile(do.call(c,AUC_10Cross), probs = 0.95, na.rm = T),
                      Threshold = Thresh,
                      Compartment = comp))
      
    } else {print(paste("Not enough data for", which(Thresholds == Thresh)));return(NA)}
    
})

return(do.call(rbind, Thresh_AUC)) 

})


save(Thresh_draw, file = paste0("AUC_threshold_50sites", param))
}


load("AUC_threshold_NC_ALT2_L4_test")

Thresh_plot <- do.call(rbind, lapply(Thresh_draw, function(x) return(x[[1]])))

Thresh_envelope <- as.data.frame(do.call(rbind, lapply(
  Thresh_draw[which(unlist(lapply(Thresh_draw, length)) == 2)],
  function(x) return(x[[2]]))))

ggplot(data = Thresh_plot[complete.cases(Thresh_plot),], aes(x = Threshold, y = AUC_val, color = AUC_type)) + 
  geom_line() + theme_classic() + ggtitle("Riparian vegetation score") + ylab("AUC") + 
  geom_ribbon(data = Thresh_envelope[complete.cases(Thresh_envelope),], aes(x = Threshold, y= AUC_val, ymax = Up, ymin = Low), 
              color = "white", fill = "olivedrab", alpha = 0.3)
  
  

## Plot density distributions

Titles <- fread("MetricsTitles.csv")

Density_plots <- lapply(Params, function(param){
  
Dens_plot <- RFD_all[, ..param]
Dens_plot <-Dens_plot[complete.cases(Dens_plot)]

Pdens <- ggplot(Dens_plot, aes(x = !!sym(param))) + geom_density()+
  ggtitle(Titles[Name_Origine == param, Description]) + 
        theme_classic() +
        labs(x = paste0("(",Catalog[NameOrigin == param, Unit],") ")) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
              axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
              axis.title=element_text(size=10), plot.title = element_text(size = 10))

return(Pdens)
})

names(Density_plots) <- Params

lapply(unique(Titles$Category), function(cat){
  
  Density_cat <- Density_plots[Titles[Category == cat, Name_Origine]]
  
  SplitPlotCategory <- split(1:length(Density_cat), ceiling(1:length(Density_cat)/6))
  
  lapply(1:length(SplitPlotCategory), function(j){
    
    Image <- ggarrange(plotlist = Density_cat[SplitPlotCategory[[j]]], ncol = 3, nrow = 2, 
                       common.legend = T, legend = "right")
    
    png(paste0("Density_plot", cat, j,".png"), width = 800, height = 400,)
    
    print(annotate_figure(Image, top = text_grob(cat, size = 10)))
    
    dev.off()
    
  })
})


