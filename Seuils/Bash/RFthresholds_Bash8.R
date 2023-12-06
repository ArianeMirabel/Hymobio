#!/usr/bin/Rscript 
invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

Nslice <- 8

load("HYMOBIO_FULLDATA_20231129.RData")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]
Params <- grep(paste(Catalog$NameOrigin, collapse = "|"), colnames(RFD_all), value = "T")
Params <- split(Params, ceiling(seq_along(Params)/(length(Params)%/%9)))[[Nslice]]

for(param in Params){
  
  RFD <- RFD_all[, grep(paste0("B_FISH|B_INV|B_DIA|",param), colnames(RFD_all), value = T), with = F]
  Thresholds <- seq(from = range(RFD[,..param], na.rm = T)[1], to = range(RFD[,..param], na.rm = T)[2], length.out = 10)
  
  if(!is.na(Catalog[NameOrigin == param,LittThreshold])){
    Thresholds[which.min(abs(Thresholds-as.numeric(Catalog[NameOrigin == param,LittThreshold])))] <- Catalog[NameOrigin == param,LittThreshold]
  }
  
  Thresh_draw <- lapply(Thresholds, function(Thresh){
    
    Rfd <- RFD[get(param) > Thresh, State := "bad"][get(param) <= Thresh, State := "good"][
      , State := as.factor(State)]
    
    Thresh_AUC <- lapply(c("B_FISH", "B_INV", "B_DIA"), function(comp){
      
      Rfd_comp <- Rfd[, grep(paste0(comp,"|State|",param), colnames(Rfd), value = T), with = F]
      Rfd_comp <- Rfd_comp[complete.cases(Rfd_comp)]
      
      if(!any(c(nrow(Rfd_comp[State == "good"]), nrow(Rfd_comp[State == "bad"])) < 100) & uniqueN(Rfd_comp$State) > 1){
        
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
      
    })
    
    return(do.call(rbind, Thresh_AUC)) 
    
  })
  
  save(Thresh_draw, file = paste0("AUC_threshold_", param))
}

