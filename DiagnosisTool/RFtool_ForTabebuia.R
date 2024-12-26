invisible(lapply(c("data.table", "randomForest", "microbenchmark", "parallel"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))


setwd("D:/Mes Donnees/Hymobio/DataBase_treatment")

toremove <- "AUC_threshold_2405_Quant95_5step_"#"AUC_threshold_2405_Q95_5stp_SMOTE"

Catalog <- fread("MetricsCatalogue.csv")[which(TokeepThreshold)]

invisible(lapply(c("data.table", "randomForest"),function(pk){
    if(!pk %in% row.names(installed.packages())){install.packages(pk)}
    library(pk,character.only=T)}))
  
load("HYMOBIO_FULLDATA_202405.RData")
  
toremove <- "AUC_threshold_2405_Quant95_5step_"#"AUC_threshold_2405_Q95_5stp_SMOTE"
  
load("For_ImpactProba")

mclapply(Threshold_detect$Param, function(param){
  
  
  RFD <- RFD_all[, c("ID", "Year", grep(paste0("|B_FISH|B_INV|B_DIA|",param), colnames(RFD_all), value = T)), with = F]
 
  Thresh <- Threshold_detect[Param == param, Threshold]
  
  DegradationDirection <- as.character(Catalog[NameOrigin == param, DegradationDirection])

  Rfd <- RFD[, State := ifelse(DegradationDirection == "1", "good", "bad")][get(param) > Thresh, 
             State := ifelse(DegradationDirection == "1", "bad", "good")][, State := as.factor(State)]
  
  PI_bio <- lapply(c("B_FISH", "B_INV","B_DIA"), function(comp){
    
    Rfd_comp <- Rfd[, c("ID", "Year", grep(paste0(comp,"|State|",param), colnames(Rfd), value = T)), with = F]
    Rfd_comp <- Rfd_comp[complete.cases(Rfd_comp)]
    
    rf <- tryCatch({randomForest(State ~ .,
                                 data = Rfd_comp[, grep(paste0(comp,"|State"), colnames(Rfd_comp), value = T), with = F],
                                 scale.permutation.importance = TRUE,
                                 respect.unordered.factors = TRUE, ntree=500,
                                 mtry = 100, nodesize = nrow(Rfd_comp)*0.005)
            }, warning = function(w) {print(paste(w, "Cross val on threshold", Thresh, "for", param))},
            error = function(e) {print(paste(e, "Cross val on threshold", Thresh, "for", param))})
            
    PI <- cbind(data.table(predict(rf, Rfd_comp, type = "prob")), Rfd_comp[, .(ID, Year)])[
      , c("Comp","Param") := list(comp, param)]
    #setnames(PI, c("bad", "good"), paste0(comp, c("_bad", "_good")))
    
    return(PI)
  })
  
  PI_bio_merge <- do.call(rbind, PI_bio)
  PI_bio_merge <- dcast(PI_bio_merge, formula = ID + Year + Param ~ Comp, value.var = "bad", fun.aggregate = mean)
  
  save(PI_bio_merge, file = paste0("ImpactProb_", param))
}, mc.cores=4) 
        

stopCluster(cl)       
          
          
      



#####
Threshold_detect <- lapply(files, function(param){
  
  print(param)
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  DegradationDirection <- as.character(Catalog[NameOrigin == gsub(toremove, "", param), DegradationDirection])
  
  Thresh_draw <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]])})
  
  if(any(unlist(lapply(Thresh_draw, is.data.frame)))){
    
    Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
    
    if(length(Thresh_plot) > 0){
      
      Thresh_plot <- as.data.table(do.call(rbind,lapply(Thresh_plot, function(x) {
        
        x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
        
        x <- as.data.table(x)[, grep(c("AUC_val|Low|Up"), colnames(x), value = T) := 
                                lapply(.SD, function(col){return(round(as.numeric(col), 2))}), 
                              .SDcols = grep(c("AUC_val|Low|Up"), colnames(x), value = T)]
        
        return(x)})))
      
      
      ThreshAUC_max <- unique(Thresh_plot[, c("AUCmaxTest","THRmaxTest") 
                                          := list(max(AUC_valTest, na.rm = T),Threshold[which.max(AUC_valTest)]),
                                          by = Compartment][,
                                                            .(Compartment, AUCmaxTest, THRmaxTest)])[,Param := gsub(toremove, "", param)]
      
    }
  }
  
  
  load(paste0(getwd(),"/Genomig_Output/",param))
  
  Thresh_names <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]][,"Compartment"])})
  Thresh_thresholds <- lapply(Thresh_draw, function(X){return(unique(X[["Thresh_AUCval"]]["Threshold"]))})
  Thresh_draw <- lapply(Thresh_draw, function(X){return(X[["Thresh_OtherMetrics"]])})
  
  if(any(unlist(lapply(Thresh_draw, is.list)))){
    
    Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.list))]
    Thresh_names <- Thresh_names[unlist(lapply(Thresh_draw, is.list))]
    Thresh_thresholds <- Thresh_thresholds[unlist(lapply(Thresh_draw, is.list))]
    
    if(length(Thresh_plot) > 0){
      
      Thresh_plot <- as.data.table(do.call(rbind,lapply(1:length(Thresh_plot), function(Xcomp) {
        
        Comp <- Thresh_plot[[Xcomp]]
        Names <- Thresh_names[[Xcomp]]
        Thresh <- Thresh_thresholds[[Xcomp]]
        
        if(length(Comp) !=0){
          Comp_Test <- as.data.table(do.call(rbind, lapply(1:length(Comp), function(xcomp){
            
            x <- Comp[[xcomp]]
            name <- Names[[xcomp]]
            
            ret <- data.table(do.call(rbind, lapply(x, function(petitX){return(
              c(petitX$ConfusionTest$overall["Accuracy"],
                petitX$ConfusionTest$byClass[c("Sensitivity", "Specificity", "Precision", "Recall", "F1")])
            )})))
            ret <-  ret[, lapply(.SD, quantile, prob = c(0.05, .5, 0.95), na.rm = T), 
                        .SDcols = colnames(ret)][, p := c("Low", "Med", "High")]
            
            ret <- dcast(melt(ret, id.vars = "p", variable.name = "Metric", value.name = "Value"),
                         Metric ~ p , value.var = "Value")[, Compartment := name][, Threshold := Thresh]
            
            return(ret)
          })
          ))
          return(Comp_Test)
        } 
      })))
    }
    
    if(length(Thresh_plot) != 0){
      
      ThreshInt_max <- dcast(Thresh_plot[!is.na(Threshold)], 
                             Compartment + Threshold ~ Metric, value.var = "Med")
      
      ThreshInt_max <- unique(ThreshInt_max[, c("THRintersect", "THRF1", "F1max") := 
                                              list(Threshold[which.min(abs(Sensitivity - Specificity))],
                                                   Threshold[which.max(F1)], max(F1, na.rm = T)),
                                            by = Compartment][,.(Compartment, THRintersect, THRF1, F1max)])
    }
  }
  
  ifelse(any(grep("Quant95_5step_HM_", param)), AUClim <- 0.6, AUClim <- 0.7)
  
  if(exists("ThreshInt_max") & exists("ThreshAUC_max")){
    
    if(any(ThreshAUC_max$AUCmaxTest >= AUClim)){
      
      Thresh_max <- merge(ThreshInt_max, ThreshAUC_max)[AUCmaxTest >= AUClim]
      
      cols <- c("THRintersect", "THRmaxTest")
      Thresh_max[, (cols) := lapply(.SD, switch(DegradationDirection, "1" = min,"-1" = max, mean)), .SDcols = cols]
      
      Thresh_max <- unique(Thresh_max[, Threshold :=  apply(.SD, 1, switch(DegradationDirection, "1" = min,
                                                                           "-1" = max, mean)), .SDcols = cols][
                                                                             , .(Threshold, THRintersect, THRF1, THRmaxTest, Param)])
      
      return(Thresh_max)       
      
    }}
  
  
})

Threshold_detect <- do.call(rbind,Threshold_detect)

save(Threshold_detect, file = "For_ImpactProba")