invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

load("AMOBIO_WP3_2_FISH_INV_DIA_ENTROPIE_20231123.Rdata")
Carhyce_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN)[, H_Ncrue_S1_ALL := NULL]
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Hydro/HydroIndex_36125all_AM_20232911")
Carhyce_all <- left_join(Carhyce_all, 
                         Hydrolaps[, c("ID_AMOBIO_START", "Samp_date", grep("Ncrues|Netiage", colnames(Hydrolaps), value = T)), with = F],
                         by = c("ID" = "ID_AMOBIO_START", "DATE_OPERATION" = "Samp_date"))
save(Carhyce_all, file = "HYMOBIO_FULLDATA_20231129")
rm(list= c("MATRIX_AMOBIO_WP3_CLEAN", "Hydrolaps"))


Nslice <- 1

load("HYMOBIO_FULLDATA_20231129")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]
Params <- paste(Catalog$Category, Catalog$Name, sep = "_")
Params <- grep(paste(Params, collapse = "|"), colnames(Carhyce_all), value = "T")
Params <- split(Params, ceiling(seq_along(Params)/(length(Params)%/%9)))[[Nslice]]

for(param in Params){
  
  print(param)
  
Carhyce_all <- Carhyce_all[, grep(paste0("B_FISH|B_INV|B_DIA|",param), colnames(Carhyce_all), value = T), with = F]


Thresholds <- seq(from = range(Carhyce_all[,..param], na.rm = T)[1], to = range(Carhyce_all[,..param], na.rm = T)[2], length.out = 10)

if(!is.na(Catalog[NameOrigin == param,LittThreshold])){
  Thresholds[which.min(abs(Thresholds-Catalog[NameOrigin == param,LittThreshold]))] <- Catalog[NameOrigin == param,LittThreshold]
}

Thresh_draw <- lapply(Thresholds, function(Thresh){
         
         pb <- txtProgressBar(min = 0, max = 10, style = 3)
         setTxtProgressBar(pb, which(Thresholds == Thresh))
         
Rfd <- Carhyce_all[get(param) > Thresh, State := "bad"][get(param) <= Thresh, State := "good"][
  , State := as.factor(State)]

Thresh_AUC <- lapply(c("B_FISH", "B_INV", "B_DIA"), function(comp){
  
  Rfd_comp <- Rfd[, grep(paste0(comp,"|State|",param), colnames(Rfd), value = T), with = F]
  Rfd_comp <- Rfd_comp[complete.cases(Rfd_comp)]
  
  if(uniqueN(Rfd_comp$State) > 1){
    
    if(!any(c(nrow(Rfd_comp[State == "good"]), nrow(Rfd_comp[State == "bad"])) < 100)){
      
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
      
    } else {return(NA)}
  } else {return(NA)}
    
})

return(do.call(rbind, Thresh_AUC)) 

})


save(Thresh_draw, file = paste0("AUC_threshold_", param))
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
  
  
  


