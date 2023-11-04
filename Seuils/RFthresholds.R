invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))


setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")


# Plot variables density distribution
#####
load("AMOBIO_WP3_2_FISH_INV_DIA_20230817.Rdata")
Carhyce_all <- as.data.table(MATRIX_AMOBIO_WP3_CLEAN)
Carhyce_all <- Carhyce_all[, grep("B_FISH|B_MAC|B_DIATOM|RIPARIAN3", colnames(Carhyce_all), value = T), with = F]
#Carhyce_all$NAs <- apply(Carhyce_all[,grep("B_FISH|B_MAC|B_DIATOM", colnames(Carhyce_all), value = T), with = F], 1, function(r){any(!is.na(r))}) 
#Carhyce_all <- Carhyce_all[which(NAs)][, NAs := NULL][!is.na(HM_RIPARIAN3_S1_3Y)]
Carhyce_all <- Carhyce_all[complete.cases(Carhyce_all)]

rm("MATRIX_AMOBIO_WP3_CLEAN")
Carhyce_test <- Carhyce_all[1:1000]

Thresholds <- seq(from=range(Carhyce_test$HM_RIPARIAN3_S1_3Y)[1], to = range(Carhyce_test$HM_RIPARIAN3_S1_3Y)[2], length.out = 10)

Thresh_draw <- lapply(Thresholds, function(Thresh){
         
         pb <- txtProgressBar(min = 0, max = 10, style = 3)
         
         setTxtProgressBar(pb, which(Thresholds == Thresh))
         
Rfd <- Carhyce_test[HM_RIPARIAN3_S1_3Y > Thresh, Riparian_state := "bad"][
  HM_RIPARIAN3_S1_3Y <= Thresh, Riparian_state := "good"][
  ,Riparian_state := as.factor(Riparian_state)]

if(uniqueN(Rfd$Riparian_state) == 1){return(NA)
  } else {
    
    if(!any(c(nrow(Rfd[Riparian_state == "good"]),
             nrow(Rfd[Riparian_state == "bad"])) < 100)){
      
      set.seed(1)
      train_index <- sample(c(TRUE, FALSE), nrow(Rfd), replace=TRUE, prob=c(0.7,0.3))
      Rfd_train <- Rfd[train_index,][,!"HM_RIPARIAN3_S1_3Y", with = F]
      Rfd_test <- Rfd[!train_index,][,!"HM_RIPARIAN3_S1_3Y", with = F]
      
      rf_train <- tryCatch({randomForest(Riparian_state ~ .,
                                   data = Rfd_train,
                                   scale.permutation.importance = TRUE,
                                   respect.unordered.factors = TRUE,
                                   mtry = 3)},
                     warning = function(w) {print(paste(w, "On threshold", Thresh))},
                     error = function(e) {print(paste(e, "On threshold", Thresh))})
      
      
      rf_all <- tryCatch({randomForest(Riparian_state ~ .,
                                       data = Rfd[,!"HM_RIPARIAN3_S1_3Y", with = F],
                                       scale.permutation.importance = TRUE,
                                       respect.unordered.factors = TRUE,
                                       mtry = 3)},
                         warning = function(w) {print(paste(w, "On threshold", Thresh))},
                         error = function(e) {print(paste(e, "On threshold", Thresh))})
 
      AUC_all <- performance(prediction(
        as.data.table(rf_all$votes[,2]), Rfd$Riparian_state), "auc")
      
      AUC_test_trainRf <- performance(prediction(
        as.numeric(predict(rf_train, 
        Rfd_test[,!"Riparian_state", with = F])), Rfd_test$Riparian_state),"auc")
      
      AUC_test_allRf <- performance(prediction(
        as.numeric(predict(rf_all, 
        Rfd_test[,!"Riparian_state", with = F])), Rfd_test$Riparian_state),"auc")
      
      #####
      # 10-fold cross-validation
      set.seed(1)
      Rfd[, Cross_grp := sample(x=1:10, nrow(Rfd),replace=TRUE)]
     
      AUC_10Cross <- lapply(1:10, function(grp){
       
        Cross_train <- Rfd[Cross_grp != grp]
        Cross_test <- Rfd[Cross_grp == grp]
        
        if(uniqueN(Cross_train$Riparian_state) > 1){
           
        rf_cross <- tryCatch({randomForest(Riparian_state ~ .,
                                           data = Cross_train,
                                           scale.permutation.importance = TRUE,
                                           respect.unordered.factors = TRUE,
                                           mtry = 3)},
                             warning = function(w) {print(paste(w, "Cross val on threshold", Thresh))},
                             error = function(e) {print(paste(e, "Cross val on threshold", Thresh))})
        
        AUC_cross <- performance(prediction(
          as.numeric(predict(rf_cross, 
          Cross_test[,!"Riparian_state", with = F])), Cross_test$Riparian_state),"auc")
        
        return(AUC_cross@y.values[[1]])
        
        } else {return(NA)}
        
      }) 
      
      AUC_10Cross <- mean(do.call(c,AUC_10Cross), na.rm = T)
      
    #####
      
      
      return(data.frame(
        AUC_val = c(AUC_all@y.values[[1]], AUC_test_trainRf@y.values[[1]], 
                    AUC_test_allRf@y.values[[1]], AUC_10Cross),
        AUC_type = c("All_Rf", "Test_trainRf", "Test_allRf", "Mean_10crossValidation"),
        Threshold = Thresh))
      
    } else {return(NA)}
    
} 
       })

Thresh_plot <- do.call(rbind, Thresh_draw)

ggplot(data = Thresh_plot[complete.cases(Thresh_plot),], aes(x = Threshold, y = AUC_val, color = AUC_type)) + 
  geom_line() + theme_classic() + ggtitle("Riparian vegetation score") + ylab("AUC")
  


