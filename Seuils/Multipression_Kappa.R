invisible(lapply(c("data.table", "psych", "sf", "ggplot2","ggpubr", "scales"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils")

load("../HYMOBIO_FULLDATA_202405.RData")

files <- list.files("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
Catalog <- fread("../MetricsCatalogue.csv")#[DegradationDirection != 0]
toremove <- "AUC_threshold_2405_Quant95_5step_"
Titles <- fread("SumUp_Seuil.csv")
Titles[, Reference_Threshold := as.numeric(gsub(",",".", Reference_Threshold))]


files <- grep(paste(Catalog$NameOrigin, collapse = "|"), grep(toremove, files, value = T), value = T)
files <- grep(toremove, files, value = T)
files <- gsub(toremove, "", files)

RFD <- RFD_all[,c("ID", "Year",grep(paste(files,collapse = "|"),colnames(RFD_all),value = T)),with = F]
#####
load("FinalTable")

files <- Table[AUC_valTest >= 0.7 & F1 >= 0.7 & DegradationDirection != 0, Name_Origine]

Qualificate_Samples <- function(RFD_sample, Files){
  
  RFD_return <- lapply(Files, function(param){
    print(param)
  
  load(paste0(getwd(),"/Genomig_Output/", toremove, param))
    
  Thresh_draw <- lapply(Thresh_draw, function(X){return(X[["Thresh_AUCval"]])})
    
  if(any(unlist(lapply(Thresh_draw, is.data.frame)))){
    
    #Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
    
  #Thresh_plot <- as.data.table(do.call(rbind,lapply(Thresh_plot, function(x) {
    
    #x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
    
    #x <- as.data.table(x)[, grep(c("AUC_val|Low|Up|Threshold"), colnames(x), value = T) := 
    #                        lapply(.SD, function(col){return(round(as.numeric(col), 2))}), 
    #                      .SDcols = grep(c("AUC_val|Low|Up|Threshold"), colnames(x), value = T)]
    
    #x <- merge(x, data.table(Compartment = c("B_FISH", "B_INV", "B_DIA"),
    #                         CompartmentName = c("Fish", "Macroinvertebrate", "Diatom")), all.y = T)
   # return(x)})))[
    #  , CompartmentName := factor(CompartmentName, levels = c("Fish", "Macroinvertebrate", "Diatom"))][
     #   , N := length(which(!is.na(AUC_valTest))), by = CompartmentName]
  
 # Thresh_max <- unique(Thresh_plot[, c("AUCmaxTest","THRmaxTest") := list(max(AUC_valTest, na.rm = T),
  #                                                                        Threshold[which.max(AUC_valTest)]), by = CompartmentName][
   #                                                                         ,.(CompartmentName, AUCmaxTest, THRmaxTest)])
  #Thresh_max[Thresh_max==-Inf] <- NA
  
    Thresh_max <- Table[Name_Origine == param] 
    
  #if(any(Thresh_max$AUCmaxTest >= 0.5)){
    if(Thresh_max[,AUC_valTest] >= 0.7){
    Dir <- Catalog[NameOrigin == param, DegradationDirection]
    
    Threshold <- Thresh_max[, ifelse(Dir > 0, min(Threshold, THRintersect), max(Threshold, THRintersect))]
    
    ret <- RFD_sample[, param, with = F]
    ret[!is.na(get(param)), paste0("State_", param) := "good"]
    #ret[!is.na(get(param)), paste0("State_", param) := max(Thresh_max$AUC_valTest)]
    
    if(Dir > 0){ 
      ret[get(param) >= Threshold, paste0("State_", param) := "bad"]
    } else {ret[get(param) <= Threshold, paste0("State_", param) := "bad"]}
    
    #if(Dir > 0){ 
     # ret[get(param) >= Threshold, paste0("State_", param) := -max(Thresh_max$AUCmaxTest)]
    #} else {ret[get(param) <= Threshold, paste0("State_", param) := -max(Thresh_max$AUCmaxTest)]}
    
  return(ret)
  }
    #else {ret <- RFD_sample[, param, with = F][!is.na(get(param)), paste0("State_", param) := NA][
    #, paste0("AUC_", param) := max(Thresh_max$AUCmaxTest)]}
  
  
 # return(list(DirMat = ret,
  #            AUC = data.table(AUC = max(Thresh_max$AUCmaxTest), Var = param)))
  }
  
  })
  
  return(cbind(RFD_sample[,.(ID, Year)],do.call(cbind, RFD_return)))
  
 # return(list(DirMat = cbind(RFD_sample[,.(ID, Year)],
  #                           do.call(cbind, lapply(RFD_return, function(param){return(param[["DirMat"]])}))),
   #           AUC = do.call(rbind, lapply(RFD_return, function(param){return(param[["AUC"]])}))))
}


BioState <- Qualificate_Samples(RFD_sample = RFD[, .SD[which.max(Year)], ,by = "ID"], Files = files)
save(BioState, file = "BioState_save")

BioState_Kappa <- as.data.table(t(apply(combn(grep("State", colnames(BioState), value = T), 2), 2, function(x) {
  res <- BioState[,x, with = F][complete.cases(BioState[,x, with = F])]
  if(nrow(res)!= 0){
    return(c(gsub("State_", "", x), round(
               cohen.kappa(BioState[,x, with = F][complete.cases(BioState[,x, with = F])])[["confid"]]["unweighted kappa",], 2)))
  } else {return(c(gsub("State_", "", x), "lower" = NA, "estimate" = NA, "upper" = NA))}
})))


BioState_Kappa <- BioState_Kappa[complete.cases(BioState_Kappa)]

Pkappa <- BioState_Kappa[, lapply(.SD, function(x){mean(as.numeric(x))}), .SDcols = c("lower", "estimate", "upper"), by = V1]
setnames(Pkappa, "V1", "Name_Origine")
Pkappa <- merge(Pkappa, Titles[,.(Name_Origine, Description, Category, Scale)])
Pkappa[, Category := gsub("Environment", "Env.", Category)]
Pkappa[, Category := factor(Category, levels = c("Natural Env.", "Hydromorphology", "Hydrology", "Temperature",
                            "Anthropic Env.", "Flow Barriers"))]
Pkappa[grep("Hydro area|Proximal|USRA", Scale), Scale := "Local"][
    Scale == "Downstream stream to the mouth", Scale := "Downstream"][
      grep("(edge)", Scale), Scale := "Upstream"][
        grep("Downstream stream|Upstream subsystem", Scale), Scale := "Large"]

CatColors <- c("Natural Env." = "#00BA38", "Hydromorphology" = "#00BFC4", "Hydrology" = "#619CFF",
               "Temperature" = "#DF536B", "Anthropic Env." = "#F5C710", "Flow Barriers"="tomato4") 

ggplot(Pkappa, aes(x= Description, y = estimate, color = Category)) +
  geom_crossbar(aes(ymin = lower, ymax = upper)) +
  geom_text(data = Pkappa[Description %in% Pkappa[duplicated(Description), Description]],
             aes(x = Description, y = estimate, label = Scale), vjust = 1.4, size = 3) +
  theme_classic() + labs(x = "", y = "Mean Cohen Kappa") + 
  facet_grid(. ~ Category, scales = "free_x", space = "free_x", 
  labeller = labeller(Category = setNames(str_wrap(unique(Pkappa$Category), width = 10),unique(Pkappa$Category)))) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) +
  scale_color_manual(values = CatColors[which(names(CatColors) %in% Pkappa$Category)]) +
  geom_hline(yintercept = 0, color = "darkgrey", lty = 1) 


#####
load("D:/Mes Donnees/Hymobio/DataBase_treatment/BD_Bio/AllStation_ArianeMarch2023")
load("BioState_save")
Titles <- fread("SumUp_Seuil.csv")
BioStateMap <- merge(BioState, unique(AllStations[,.(CdStation, CoordX_L93, CoordY_L93)]), 
                  all.x = T, all.y = F, by.x = "ID", by.y = "CdStation")
BioStateMap <- BioStateMap[, c(c("ID", "CoordX_L93", "CoordY_L93"), 
                  grep(paste(paste0("State_|",Catalog[DegradationDirection != 0, NameOrigin]), collapse = "|"), 
                  colnames(BioStateMap), value = T)), with = F]

  
BioStateMap_bads <- BioStateMap[, "PC_Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), 
                                         .SDcols = patterns('State_PC_'), by = 'ID'][
                               , "NC_Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), 
                                         .SDcols = patterns('State_NC_'), by = 'ID'][
                               , "POE_Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), 
                                         .SDcols = patterns('State_POE_'), by = 'ID'][
                               , "T_Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), 
                                         .SDcols = patterns('State_T_'), by = 'ID'][
                               #, "HM_Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), 
                                #         .SDcols = patterns('State_HM_'), by = 'ID'][
                               , "P_Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), 
                                         .SDcols = patterns('State_P_'), by = 'ID'][
                               #, "H_Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), 
                                #         .SDcols = patterns('State_H_'), by = 'ID'][
                               , "Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))),
                                         .SDcols = patterns('State'), by = 'ID']

#BioStateMap_bads <- BioStateMap[, "PC_Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State_PC_'), by = 'ID'][
#  , "NC_Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State_NC_'), by = 'ID'][
#  , "POE_Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State_POE_'), by = 'ID'][
#  , "T_Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State_T_'), by = 'ID'][
#  , "HM_Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State_HM_'), by = 'ID'][
#  , "P_Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State_P_'), by = 'ID'][
#  , "H_Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State_H_'), by = 'ID'][
#  , "Bads" := -sum(.SD, na.rm = T),.SDcols = patterns('State'), by = 'ID']
                                

BioStateMap_toplot <- BioStateMap[, "Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))),
 .SDcols = patterns('State'), by = 'ID']

#BioStateMap_toplot <- BioStateMap[, "Bads" := -sum(.SD, na.rm = T), .SDcols = patterns('State'), by = 'ID']

BioStateMap_tofacet <- BioStateMap_bads[ , c("ID", "CoordX_L93", "CoordY_L93", 
                                           grep("Bads", colnames(BioStateMap_bads), value = T)), with = F]

BioStateMap_tofacet <- melt(setDT(BioStateMap_tofacet), measure = grep("Bads", colnames(BioStateMap_tofacet), value = T), 
     value.name = "Bads",
     variable.name = "Category")
BioStateMap_tofacet <- unique(BioStateMap_tofacet[, Bads := mean(Bads), 
                        by = c("ID", "CoordX_L93", "CoordY_L93", "Category")])
BioStateMap_tofacet[,Category := gsub("Bads","All variables",gsub("_Bads","",Category))]

#BioStateMap_tofacet <- unique(BioStateMap[, "Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), .SDcols = patterns('State'), by = 'ID'][
#  ,.(ID, CoordX_L93, CoordY_L93, Bads)])

#BioStateMap <- cbind(BioStateMap, data.frame(GlobalState = apply(BioStateMap,1,
#                        function(R){c("Bads" = length(which(R=="bad")), "Goods" = length(which(R=="good")))})))

FrenchOutline <- st_read("BassinHydrographique_TOPAGE_UNION_20240301.shp")
st_as_sf(AllStations, coords = c("CoordX_L93", "CoordY_L93"))

BioStateMap_tofacet[Category == "PC", Category := "Physico-chemistry"][Category == "NC", Category := "Natural Env."][
  Category == "POE", Category := "Flow barriers"][Category == "T", Category := "Temperature"][
  Category == "HM", Category := "Hydromorphology"][Category == "P", Category := "Anthropic Env."][
  Category == "H", Category := "Hydrology"]

Plot_cat <- 
  
  ggplot(data = FrenchOutline) +  
  geom_sf(fill="white", color="gray") +
  geom_point(data = BioStateMap_tofacet[Category != 'All variables'], aes(x = CoordX_L93, y = CoordY_L93, color = Bads), size = 1) +
  scale_colour_distiller(#low = "cornflowerblue", mid = "snow2",high = "firebrick",  midpoint = 0.5, oob=squish, 
    na.value = NA, palette = "YlOrRd", direction = 1,
    name = "Biological state\n", breaks = c(0,0.5,1), limits = c(0,1), 
    labels = c("Least impacted", "", "Most impacted")) + 
  facet_wrap(Category ~ ., ncol = 4) + theme_void() +
  ggtitle("(B)") + theme(plot.title = element_text(hjust = 0))


Plot_All <-
  
  ggplot(data = FrenchOutline) +  
  geom_sf( fill="white", color="gray") +
  geom_point(data = BioStateMap_toplot, aes(x = CoordX_L93, y = CoordY_L93, colour = Bads), size = 1) +
  scale_colour_distiller(#low = "cornflowerblue", mid = "snow2",high = "firebrick",  midpoint = 0.5, oob=squish, 
                        na.value = NA, palette = "YlOrRd", direction = 1,
                        name = "Biological state\n", breaks = c(0,0.5,1), limits = c(0,1), 
                        labels = c("Least impacted", "", "Most impacted")) + 
  theme_void() + theme(plot.title = element_text(hjust = 0)) +
  ggtitle("(A) Main biological state regarding detected thresholds", subtitle = "All variables")

ggarrange(Plot_All,Plot_cat,  nrow = 2, common.legend = T, legend = "right")
