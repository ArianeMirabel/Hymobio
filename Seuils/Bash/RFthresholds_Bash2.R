#!/usr/bin/Rscript 
invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

load("AMOBIO_WP3_2_FISH_INV_DIA_20230817.Rdata")
RFD_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN)
rm("MATRIX_AMOBIO_WP3_CLEAN")

Nslice <- 2

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]
Params <- paste(Catalog$Category, Catalog$Name, sep = "_")
Params <- grep(paste(Params, collapse = "|"), colnames(RFD_all), value = "T")
Params <- split(Params, ceiling(seq_along(Params)/(length(Params)%/%9)))[[Nslice]]


for(param in Params){
  
  RFD <- RFD_all[, grep(paste0("B_FISH|B_INV|B_DIA|",param), colnames(RFD_all), value = T), with = F]
  Thresholds <- seq(from = range(RFD[,..param], na.rm = T)[1], to = range(RFD[,..param], na.rm = T)[2], length.out = 10)
  
  Thresh_draw <- lapply(Thresholds, function(Thresh){
    
    Rfd <- RFD[get(param) > Thresh, State := "bad"][get(param) <= Thresh, State := "good"][
      , State := as.factor(State)]
    
    Thresh_AUC <- lapply(c("FISH", "INV", "DIA"), function(comp){
      
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

