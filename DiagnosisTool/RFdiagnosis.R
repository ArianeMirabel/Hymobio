#!/usr/bin/Rscript 
invisible(lapply(c("data.table", "randomForest", "ROCR"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/DiagnosisTool")

load("../HYMOBIO_FULLDATA_202403.RData")

Nslice <- 1

Catalog <- fread("../MetricsCatalogue.csv")[which(TokeepDiagnosis)]
Catalog[, FinalThreshold := as.numeric(gsub(",",".", FinalThreshold))]
Params <- grep(paste(Catalog$NameOrigin, collapse = "|"), colnames(RFD_all), value = "T")
Params <- split(Params, ceiling(seq_along(Params)/(length(Params)%/%9)))[[Nslice]]

param <- "HM_RIPARIAN2_S1_3Y";comp<-"B_FISH"

for(param in Params){
  
  RFD <- RFD_all[, grep(paste0("ID|B_FISH|B_INV|B_DIA|",param), colnames(RFD_all), value = T), with = F]
  
  RFD[get(param) <= Catalog[NameOrigin == param, FinalThreshold], State := "good"][is.na(State), State := "bad"][
    , State := as.factor(State)]
  
  Param_diag <- lapply(c("B_FISH", "B_INV","B_DIA"), function(comp){
  
       Rfd_comp <- RFD[, c("ID",grep(paste0(comp,"|State|",param), colnames(RFD), value = T)), with = F]
       Rfd_comp <- Rfd_comp[complete.cases(Rfd_comp)]
  
       if(!any(c(nrow(Rfd_comp[State == "good"]), nrow(Rfd_comp[State == "bad"])) < 50) & uniqueN(Rfd_comp$State) > 1){
      
      set.seed(1)
      AUC_10Cross <- lapply(1:10, function(cross){
        
        Train_grp <- sample(unique(Rfd_comp$ID), round(0.2*uniqueN(Rfd_comp$ID)))
       
        Cross_train <- Rfd_comp[ID %in% Train_grp]
        Cross_test <- Rfd_comp[!ID %in% Train_grp]
        
        if(uniqueN(Cross_train$State) > 1){
           
        rf_cross <- tryCatch({randomForest(State ~ .,
                                           data = Cross_train[, grep(paste0(comp,"|State"), colnames(Rfd_comp), value = T), with = F],
                                           scale.permutation.importance = TRUE,
                                           respect.unordered.factors = TRUE,
                                           mtry = 3)},
                             warning = function(w) {print(paste(w, "Cross val for", param))},
                             error = function(e) {print(paste(e, "Cross val for", param))})
        
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
      
      return(data.frame(AUC_valTest = mean(do.call(c,AUC_10CrossTest), na.rm = T),
                        LowTest = quantile(do.call(c,AUC_10CrossTest), probs = 0.05, na.rm = T),
                        UpTest = quantile(do.call(c,AUC_10CrossTest), probs = 0.95, na.rm = T),
                        AUC_valTrain = mean(do.call(c,AUC_10CrossTrain), na.rm = T),
                        LowTrain = quantile(do.call(c,AUC_10CrossTrain), probs = 0.05, na.rm = T),
                        UpTrain = quantile(do.call(c,AUC_10CrossTrain), probs = 0.95, na.rm = T),
                        Threshold = Thresh,
                        Compartment = comp))
      
    } else {print(paste("Not enough data for", which(Thresholds == Thresh)));return(NA)}
    
})


save(Param_diag, file = paste0("Diagnosis_", param))
}

