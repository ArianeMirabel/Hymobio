invisible(lapply(c("data.table", "psych", "sf", "ggplot2"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils")

load("../HYMOBIO_FULLDATA_202403.RData")

files <- list.files("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/Seuils/Genomig_Output")
Catalog <- fread("../MetricsCatalogue.csv")#[DegradationDirection != 0]
toremove <- "AUC_threshold_2402_Quant95_5step_"
Titles <- fread("SumUp_Seuil.csv")
Titles[, Reference_Threshold := as.numeric(gsub(",",".", Reference_Threshold))]


files <- grep(paste(Catalog$NameOrigin, collapse = "|"), grep(toremove, files, value = T), value = T)
files <- grep(toremove, files, value = T)
files <- gsub(toremove, "", files)

RFD <- RFD_all[,c("ID", "Year",grep(paste(files,collapse = "|"),colnames(RFD_all),value = T)),with = F]
#####

Qualificate_Samples <- function(RFD_sample, Files){
  
  RFD_return <- lapply(Files, function(param){
  
  load(paste0(getwd(),"/Genomig_Output/", toremove, param))
  
  Thresh_plot <- Thresh_draw[unlist(lapply(Thresh_draw, is.data.frame))]
  
  Thresh_plot <- as.data.table(do.call(rbind,lapply(Thresh_plot, function(x) {
    x$Threshold <- unique(x$Threshold)[!is.na(unique(x$Threshold))]
    x <- as.data.table(x)[, grep(c("AUC_val|Low|Up|Threshold"), colnames(x), value = T) := 
                            lapply(.SD, function(col){return(round(as.numeric(col), 2))}), 
                          .SDcols = grep(c("AUC_val|Low|Up|Threshold"), colnames(x), value = T)]
    x$Compartment = c("Fish", "Macroinvertebrate", "Diatom"); return(x)})))
  
  Thresh_max <- unique(Thresh_plot[, c("AUCmaxTest","THRmaxTest") := list(max(AUC_valTest, na.rm = T),
                                                                          Threshold[which.max(AUC_valTest)]), by = Compartment][
                                                                            ,.(Compartment, AUCmaxTest, THRmaxTest)])
  Thresh_max[Thresh_max==-Inf] <- NA
  
  if(any(Thresh_max$AUCmaxTest >= 0.7)){
    Dir <- Catalog[NameOrigin == param, DegradationDirection]
    
    Threshold <- Thresh_max[, ifelse(Dir > 0, min(THRmaxTest), max(THRmaxTest))]
    
    ret <- RFD_sample[, param, with = F]
    ret[!is.na(get(param)), paste0("State_", param) := "good"]
    if(Dir > 0){ 
      ret[get(param) >= Threshold, paste0("State_", param) := "bad"]
    } else {ret[get(param) <= Threshold, paste0("State_", param) := "bad"]}
    
  } else {ret <- RFD_sample[, param, with = F][!is.na(get(param)), paste0("State_", param) := NA]}
  
  return(ret)
  
  })
  
  return(cbind(RFD_sample[,.(ID, Year)],do.call(cbind, RFD_return))) 
}

BioState <- Qualificate_Samples(RFD, files)


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
load("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio/AllStation_ArianeMarch2023")
BioStateMap <- merge(BioState, unique(AllStations[,.(CdStation, CoordX_L93, CoordY_L93)]), 
                  all.x = T, all.y = F, by.x = "ID", by.y = "CdStation")
BioStateMap <- BioStateMap[, c(c("ID", "CoordX_L93", "CoordY_L93"), 
                  grep(paste(paste0("State_",Catalog[DegradationDirection != 0, NameOrigin]), collapse = "|"), 
                  colnames(BioStateMap), value = T)), with = F]

BioStateMap <- unique(BioStateMap[, "Bads" := length(which(.SD == "bad"))/length(which(!is.na(.SD))), .SDcols = patterns('State'), by = 'ID'][
  ,.(ID, CoordX_L93, CoordY_L93, Bads)])


#BioStateMap <- cbind(BioStateMap, data.frame(GlobalState = apply(BioStateMap,1,
#                        function(R){c("Bads" = length(which(R=="bad")), "Goods" = length(which(R=="good")))})))

FrenchOutline <- st_read("BassinHydrographique_TOPAGE_UNION_20240301.shp")
st_as_sf(AllStations, coords = c("CoordX_L93", "CoordY_L93"))


ggplot() +
  geom_sf(data = FrenchOutline, fill="white", color="black") +
  theme_void() +
  geom_point(data = BioStateMap, aes(x = CoordX_L93, y = CoordY_L93, color = Bads), size = 1) +
  scale_color_gradient(low = "cornflowerblue", high ="firebrick", breaks=c(0.3,0.7),
                     labels = c("Least impacted", "Impacted"), name = "Biological state") +
  ggtitle("Main biological state regarding detected thresholds", subtitle = "(All years and variables)")





