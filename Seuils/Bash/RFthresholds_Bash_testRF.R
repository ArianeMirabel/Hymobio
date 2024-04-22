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
        
        layout(matrix(c(1,2),nrow=1),
               width=c(4,1)) 
        par(mar=c(5,4,4,0)) #No margin on the right side
        plot(rf_cross, log="y")
        par(mar=c(5,0,4,2)) #No margin on the left side
        plot(c(0,1),type="n", axes=F, xlab="", ylab="")
        legend("top", colnames(rf_cross$err.rate),col=1:4,cex=0.8,fill=1:4)
        
        
        Metric_Test <- confusionMatrix(predict(rf_cross, Cross_test), Cross_test$State)
        
        Metric_Train <- confusionMatrix(predict(rf_cross, Cross_train), Cross_train$State)
        
        AUC_crossTest <- performance(prediction(
          as.numeric(predict(rf_cross, 
                             Cross_test)), Cross_test$State),"auc")
        
        AUC_crossTrain <- performance(prediction(
          as.numeric(predict(rf_cross, 
                             Cross_train)), Cross_train$State),"auc")
        
        return(list(Test_metric = Metric_Test, Train_metric = Metric_Train, 
                    Test = AUC_crossTest@y.values[[1]], Train = AUC_crossTrain@y.values[[1]]))
        
        } else {return(NA)}
        
      }) 
      
      names(AUC_10Cross) <- paste0("Cross_", 1:10)
      
      Metric_CrossTest <- c(quantile(unlist(lapply(AUC_10Cross, function(cross){
        ret <- cross[["Test_metric"]][["overall"]]["Accuracy"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Test_metric"]][["byClass"]]["Sensitivity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Test_metric"]][["byClass"]]["Specificity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T)
       )
      names(Metric_CrossTest) <- as.vector(outer(c("Low", "Med", "Upp"),
                                 paste0(c("Accuracy", "Sensitivity", "Specificity"), "_Test"), paste, sep ="_"))
      
      Metric_CrossTrain <- c(quantile(unlist(lapply(AUC_10Cross, function(cross){
        ret <- cross[["Train_metric"]][["overall"]]["Accuracy"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Train_metric"]][["byClass"]]["Sensitivity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Train_metric"]][["byClass"]]["Specificity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T)
      )
      names(Metric_CrossTrain) <- as.vector(outer(c("Low", "Med", "Upp"),
                                  paste0(c("Accuracy", "Sensitivity", "Specificity"), "_Train"), paste, sep ="_"))
      
      AUC_10CrossTest <- lapply(AUC_10Cross, function(cross){return(cross["Test"])})
      AUC_10CrossTrain <- lapply(AUC_10Cross, function(cross){return(cross["Train"])})
      
      AUC_CrossAll <- data.frame(AUC_valTest = mean(unlist(AUC_10CrossTest), na.rm = T),
                        LowTest = quantile(unlist(AUC_10CrossTest), probs = 0.05, na.rm = T),
                        UpTest = quantile(unlist(AUC_10CrossTest), probs = 0.95, na.rm = T),
                        AUC_valTrain = mean(unlist(AUC_10CrossTrain), na.rm = T),
                        LowTrain = quantile(unlist(AUC_10CrossTrain), probs = 0.05, na.rm = T),
                        UpTrain = quantile(unlist(AUC_10CrossTrain), probs = 0.95, na.rm = T),
                        Threshold = Thresh,
                        Compartment = comp)
      
      return(list(c(Metric_CrossTest, Metric_CrossTrain, AUC_CrossAll,
                        Threshold = Thresh,
                        Compartment = comp), AUC_10Cross))
      
    } else {print(paste("Not enough data for", which(Thresholds == Thresh)));return(NA)}
    
})
     names(Thresh_AUC) <- c("FISH", "INV","DIA")
     
return(list(Metrics = do.call(rbind, lapply(Thresh_AUC, function(Comp){return(Comp[[1]])})),
            Models =  lapply(Thresh_AUC, function(Comp){return(Comp[[2]])})))

})


save(Thresh_draw, file = paste0("AUC_threshold_RFtuneTest_Node0005", param))

