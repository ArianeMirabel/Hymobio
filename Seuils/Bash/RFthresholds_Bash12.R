#!/usr/bin/Rscript 
invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel", "dplyr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

Nslice <- 12

load("HYMOBIO_FULLDATA_202403.RData")

Catalog <- fread("MetricsCatalogue.csv")[which(TokeepThreshold)]
Params <- grep(paste(Catalog$NameOrigin, collapse = "|"), colnames(RFD_all), value = "T")
Params <- split(Params, ceiling(seq_along(Params)/(length(Params)%/%15)))[[Nslice]]

for(param in Params){
  
  RFD <- RFD_all[, c("ID",grep(paste0("B_FISH|B_INV|B_DIA|",param), colnames(RFD_all), value = T)), with = F]
  Thresholds <- seq(from = quantile(RFD[,..param], probs = 0.025, na.rm = T), 
                    to = quantile(RFD[,..param], probs = 0.975, na.rm = T), 
                    length.out = 20)
  
  if(!is.na(Catalog[NameOrigin == param,LittThreshold])){
    Thresholds[which.min(abs(Thresholds-as.numeric(Catalog[NameOrigin == param,LittThreshold])))] <- Catalog[NameOrigin == param,LittThreshold]
  }
  
  Thresh_draw <- lapply(Thresholds, function(Thresh){
    
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
                                               scale.permutation.importance = TRUE,
                                               respect.unordered.factors = TRUE, ntree = 150,
                                               mtry = 100, nodesize = nrow(Cross_train)*0.005                                               )},
                                 warning = function(w) {print(paste(w, "Cross val on threshold", Thresh, "for", param))},
                                 error = function(e) {print(paste(e, "Cross val on threshold", Thresh, "for", param))})
            AUC_crossTest <- performance(prediction(
              as.numeric(predict(rf_cross, 
                                 Cross_test)), Cross_test$State),"auc")
            
            AUC_crossTrain <- performance(prediction(
              as.numeric(predict(rf_cross, 
                                 Cross_train)), Cross_train$State),"auc")
            
            Metric_Test <- confusionMatrix(predict(rf_cross, Cross_test), Cross_test$State)
            
            Metric_Train <- confusionMatrix(predict(rf_cross, Cross_train), Cross_train$State)
            
            return(list(AUCvalues = data.frame(Test = AUC_crossTest@y.values[[1]], Train = AUC_crossTrain@y.values[[1]]),
                        OtherMetrics = list(ConfusionTest = Metric_Test, ConfusionTrain = Metric_Train),
                        Model = rf_cross))
            
          } else {return(NA)}
          
        }) 
        
        AUC_10CrossTest <- unlist(lapply(AUC_10Cross, function(cross){return(cross[["AUCvalues"]]["Test"])}))
        AUC_10CrossTrain <- unlist(lapply(AUC_10Cross, function(cross){return(cross[["AUCvalues"]]["Train"])}))
        
        return(list(AUCval = data.frame(AUC_valTest = mean(AUC_10CrossTest, na.rm = T),
                                        LowTest = quantile(AUC_10CrossTest, probs = 0.05, na.rm = T),
                                        UpTest = quantile(AUC_10CrossTest, probs = 0.95, na.rm = T),
                                        AUC_valTrain = mean(AUC_10CrossTrain, na.rm = T),
                                        LowTrain = quantile(AUC_10CrossTrain, probs = 0.05, na.rm = T),
                                        UpTrain = quantile(AUC_10CrossTrain, probs = 0.95, na.rm = T),
                                        Threshold = Thresh,
                                        Compartment = comp),
                    OtherMetrics = lapply(AUC_10Cross, function(cross){return(cross[["OtherMetrics"]])}),
                    Models = lapply(AUC_10Cross, function(cross){return(cross[["Model"]])})
        ))
        
      } else {return(NA)}
      
    })
    
    Thresh_AUCval <- do.call(rbind,lapply(Thresh_AUC, function(T) return(T$AUCval)))
    Thresh_OtherMetrics <- lapply(Thresh_AUC, function(T) return(T$OtherMetrics))
    Thresh_Models <- lapply(Thresh_AUC, function(T) return(T$Models))
    
    return(list(Thresh_AUCval = Thresh_AUCval,
                Thresh_OtherMetrics = Thresh_OtherMetrics,
                Thresh_Models = Thresh_Models))
    
  })
  
  save(Thresh_draw, file = paste0("AUC_threshold_2404_Quant95_5step_", param))
}
