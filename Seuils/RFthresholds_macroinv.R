#!/usr/bin/Rscript 
invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "randomForest", "ROCR", "doParallel", "dplyr", "caret"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

Nslice <- 1

load("../HYMOBIO_FULLDATA_202403.RData")

Catalog <- fread("../MetricsCatalogue.csv")[which(TokeepThreshold)]
Params <- grep(paste(Catalog$NameOrigin, collapse = "|"), colnames(RFD_all), value = "T")
Remaining <- setdiff(Params, gsub("AUC_threshold_2402_Quant95_5step_","",
                     grep("2402_Quant95_5step",list.files("Genomig_Output"), value = T)))
save(Remaining, file = "RemainingVariables")

load("RemainingVariables")
Params <- Remaining
Params <- split(Params, ceiling(seq_along(Params)/(length(Params)%/%9)))[[Nslice]]

param <- "NC_CatchStruct5_L5"; comp<-"B_INV"; cross <-1

for(param in Params){
  
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
           #Cross_train <- Cross_train[sample(1:nrow(Cross_train), 1000, replace=F)]
           Cross_test <- Rfd_comp[!ID %in% Train_grp]
           
        if(uniqueN(Cross_train$State) > 1){
           
        rf_cross <- tryCatch({randomForest(State ~ .,
                                           data = Cross_train[, grep(paste0(comp,"|State"), colnames(Rfd_comp), value = T), with = F],
                                           scale.permutation.importance = TRUE,
                                           respect.unordered.factors = TRUE,
                                           mtry = 3, nodesize = 50)},
                             warning = function(w) {print(paste(w, "Cross val on threshold", Thresh, "for", param))},
                             error = function(e) {print(paste(e, "Cross val on threshold", Thresh, "for", param))})
        
        Metric_Test <- confusionMatrix(predict(rf_cross, Cross_test), Cross_test$State)
        
        Metric_Train <- confusionMatrix(predict(rf_cross, Cross_train), Cross_train$State)
        
        AUC_crossTest <- performance(prediction(
          as.numeric(predict(rf_cross, Cross_test)), Cross_test$State),"auc")
        
        AUC_crossTrain <- performance(prediction(
          as.numeric(predict(rf_cross, Cross_train)), Cross_train$State),"auc")
        
        return(list(Test = Metric_Test, Train = Metric_Train, 
                    AUC_Test = AUC_crossTest@y.values[[1]], AUC_Train = AUC_crossTrain@y.values[[1]]))
        
        } else {return(NA)}
        
      }) 
      
      names(AUC_10Cross) <- paste0("Cross_", 1:10)
      
      Metric_CrossTest <- c(quantile(unlist(lapply(AUC_10Cross, function(cross){
        ret <- cross[["Test"]][["overall"]]["Accuracy"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Test"]][["byClass"]]["Sensitivity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Test"]][["byClass"]]["Specificity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T)
       )
      names(Metric_CrossTest) <- as.vector(outer(c("Low", "Med", "Upp"),
                                 paste0(c("Accuracy", "Sensitivity", "Specificity"), "_Test"), paste, sep ="_"))
      
      Metric_CrossTrain <- c(quantile(unlist(lapply(AUC_10Cross, function(cross){
        ret <- cross[["Train"]][["overall"]]["Accuracy"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Train"]][["byClass"]]["Sensitivity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T),
        quantile(unlist(lapply(AUC_10Cross, function(cross){
          ret <- cross[["Train"]][["byClass"]]["Specificity"]})), probs = c(0.05, 0.5, 0.95), na.rm = T)
      )
      names(Metric_CrossTrain) <- as.vector(outer(c("Low", "Med", "Upp"),
                                  paste0(c("Accuracy", "Sensitivity", "Specificity"), "_Train"), paste, sep ="_"))
      
      
      return(list(c(Metric_CrossTest, Metric_CrossTrain,
                        Threshold = Thresh,
                        Compartment = comp), AUC_10Cross))
      
    } else {print(paste("Not enough data for", which(Thresholds == Thresh)));return(NA)}
    
})
     names(Thresh_AUC) <- c("FISH", "INV","DIA")
     
return(list(Metrics = do.call(rbind, lapply(Thresh_AUC, function(Comp){return(Comp[[1]])})),
            Models =  lapply(Thresh_AUC, function(Comp){return(Comp[[2]])})))

})


save(Thresh_draw, file = paste0("AUC_threshold_RFtuneTest", param))
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
  
#####
# Get threshold ranges
load("HYMOBIO_FULLDATA_202401.RData")

Catalog <- fread("MetricsCatalogue.csv")[which(Tokeep)]
Params <- grep(paste(Catalog$NameOrigin, collapse = "|"), colnames(RFD_all), value = "T")

Thresholds <- as.data.table(do.call(rbind,lapply(Params, function(param){
  
  RFD <- RFD_all[, ..param]
  ret <- tryCatch({seq(from = quantile(RFD[,..param], probs = 0.025, na.rm = T), 
                       to = quantile(RFD[,..param], probs = 0.975, na.rm = T), 
                       length.out = 20)},
           error = function(e) {print(e); return(NA)})
  
})))[, Name_Origine := Params]

save(Thresholds, file = "variable_thresholds")
  

## Plot density distributions

Titles <- fread("MetricsTitles.csv")

Density_plots <- lapply(Params, function(param){
  
Dens_plot <- RFD_all[, ..param]
Dens_plot <-Dens_plot[complete.cases(Dens_plot)]

Pdens <- ggplot(Dens_plot, aes(x = !!sym(param))) + geom_density()+
  ggtitle(Titles[Name_Origine == param, Description]) + 
        theme_classic() +
        labs(x = paste0("(",Catalog[NameOrigin == param, Unit],") ")) + 
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), legend.position = 'none', 
              axis.text.x = element_text(angle = 90), axis.text=element_text(size=8),
              axis.title=element_text(size=10), plot.title = element_text(size = 10))

return(Pdens)
})

names(Density_plots) <- Params

lapply(unique(Titles$Category), function(cat){
  
  Density_cat <- Density_plots[Titles[Category == cat, Name_Origine]]
  
  SplitPlotCategory <- split(1:length(Density_cat), ceiling(1:length(Density_cat)/6))
  
  lapply(1:length(SplitPlotCategory), function(j){
    
    Image <- ggarrange(plotlist = Density_cat[SplitPlotCategory[[j]]], ncol = 3, nrow = 2, 
                       common.legend = T, legend = "right")
    
    png(paste0("Density_plot", cat, j,".png"), width = 800, height = 400,)
    
    print(annotate_figure(Image, top = text_grob(cat, size = 10)))
    
    dev.off()
    
  })
})


