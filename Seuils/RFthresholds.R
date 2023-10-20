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

Thresholds <- seq(from=range(Carhyce_all$HM_RIPARIAN3_S1_3Y)[1], to = range(Carhyce_all$HM_RIPARIAN3_S1_3Y)[2], length.out = 10)
Thresh_draw <- lapply(Thresholds, function(Thresh){
         
         pb <- txtProgressBar(min = 0, max = 10, style = 3)
         
Rfd <- Carhyce_test[HM_RIPARIAN3_S1_3Y > Thresh, Riparian_state := "bad"][HM_RIPARIAN3_S1_3Y <= Thresh, Riparian_state := "good"][
  ,Riparian_state := as.factor(Riparian_state)]

if(uniqueN(Rfd$Riparian_state) == 1){return(NA)
  } else {

rf <- tryCatch({randomForest(Riparian_state ~ .,
             data = Rfd[,!"HM_RIPARIAN3_S1_3Y", with = F],
             scale.permutation.importance = TRUE,
             respect.unordered.factors = TRUE,
             mtry = 3)},
      warning = function(w) {print(paste(w, "On threshold", Thresh))},
      error = function(e) {print(paste(e, "On threshold", Thresh))})

AUC <- performance(prediction(as.data.table(rf$votes[,2]), Rfd$Riparian_state), "auc")

setTxtProgressBar(pb, which(Thresholds == Thresh))

return(AUC@y.values[[1]])
  
} 
       })

Thresh_plot <- data.table(AUC = do.call(rbind, Thresh_draw),
                             Threshold = Thresholds)

ggplot(data = Thresh_plot, aes(x= Threshold, y = AUC)) + 
  geom_line() + theme_classic() + ggtitle("Riparian vegetation score")
  


