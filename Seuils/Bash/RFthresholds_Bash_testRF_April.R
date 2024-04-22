#!/usr/bin/Rscript 
invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel", "dplyr", "caret","randomForestSRC"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

load("HYMOBIO_FULLDATA_202403.RData")

Catalog <- fread("MetricsCatalogue.csv")[which(TokeepThreshold)]

param <- "NC_CatchStruct5_L5"

RFD <- RFD_all[, c("ID",grep(paste0("B_FISH|B_INV|B_DIA|",param), colnames(RFD_all), value = T)), with = F]
  Thresholds <- seq(from = quantile(RFD[,..param], probs = 0.025, na.rm = T), 
                    to = quantile(RFD[,..param], probs = 0.975, na.rm = T), 
                    length.out = 20)

  if(!is.na(Catalog[NameOrigin == param,LittThreshold])){
    Thresholds[which.min(abs(Thresholds-as.numeric(Catalog[NameOrigin == param,LittThreshold])))] <- Catalog[NameOrigin == param,LittThreshold]
    }

  RFD[, ThreshGroup := floor((get(param)-min(get(param), na.rm = T)) / 
                               (max(get(param), na.rm = T)-min(get(param), na.rm = T))*10)]
  
  Thresh_draw <- lapply(Thresholds, function(Thresh){
    
    #Thresh <- Thresholds[1]
    
     Rfd <- RFD[get(param) > Thresh, State := "bad"][get(param) <= Thresh, State := "good"][
       , State := as.factor(State)]
  
     Thresh_AUC <- lapply(c("B_FISH", "B_INV","B_DIA"), function(comp){
  
       Rfd_comp <- Rfd[, c("ID", grep(paste0(comp,"|State|",param), colnames(Rfd), value = T)), with = F]
       Rfd_comp <- Rfd_comp[complete.cases(Rfd_comp)]
  
       if(!any(c(nrow(Rfd_comp[State == "good"]), nrow(Rfd_comp[State == "bad"])) < 50) & uniqueN(Rfd_comp$State) > 1){
      
      set.seed(1)
         
      AUC_10Cross <- lapply(1:10, function(cross){
           
           Train_grp <- sample(unique(Rfd_comp$ID), round(0.25*uniqueN(Rfd_comp$ID)))
           
           Cross_train <- Rfd_comp[ID %in% Train_grp]
           Cross_test <- Rfd_comp[!ID %in% Train_grp]
           
        if(uniqueN(Cross_train$State) > 1){
           
        rf_cross <- tryCatch({randomForest(State ~ .,
                                           data = Cross_train[, grep(paste0(comp,"|State"), colnames(Rfd_comp), value = T), with = F],
                                           scale.permutation.importance = TRUE, ntree=150,
                                           mtry = 100, nodesize = nrow(Cross_train)*0.005)},
                             warning = function(w) {print(paste(w, "Cross val on threshold", Thresh, "for", param))},
                             error = function(e) {print(paste(e, "Cross val on threshold", Thresh, "for", param))})
        
        # Seems that reducing nodesize only helps reducing overfitting. See Segel et al. 2004, higher nodesize = sometimes better results
        
        AUC_crossTest <- performance(prediction(
          as.numeric(predict(rf_cross, 
                             Cross_test)), Cross_test$State),"auc")
        
        AUC_crossTrain <- performance(prediction(
          as.numeric(predict(rf_cross, 
                             Cross_train)), Cross_train$State),"auc")
        
        return(data.frame(Test = AUC_crossTest@y.values[[1]], Train = AUC_crossTrain@y.values[[1]]))
        
        } else {return(NA)}
        
      }) 
      
      AUC_10CrossTest <- lapply(AUC_10Cross, function(cross){return(cross["Test"])})
      AUC_10CrossTrain <- lapply(AUC_10Cross, function(cross){return(cross["Train"])})
      
      return(data.frame(AUC_valTest = mean(unlist(AUC_10CrossTest), na.rm = T),
                        LowTest = quantile(unlist(AUC_10CrossTest), probs = 0.05, na.rm = T),
                        UpTest = quantile(unlist(AUC_10CrossTest), probs = 0.95, na.rm = T),
                        AUC_valTrain = mean(unlist(AUC_10CrossTrain), na.rm = T),
                        LowTrain = quantile(unlist(AUC_10CrossTrain), probs = 0.05, na.rm = T),
                        UpTrain = quantile(unlist(AUC_10CrossTrain), probs = 0.95, na.rm = T),
                        Threshold = Thresh,
                        Compartment = comp))
      
    } else {return(NA)}
    
})
     return(do.call(rbind, Thresh_AUC)) 

})


save(Thresh_draw, file = paste0("AUC_threshold_RFtuneTest_Node0005", param))

