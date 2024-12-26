invisible(lapply(c("data.table", "ggplot2", "stringr", "ggpubr", "dplyr", "gridExtra","grid", "ggradar"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils")

Catalog <- fread("../MetricsCatalogue.csv")[which(TokeepThreshold)]
Titles <- fread("SumUp_Seuil.csv")
Titles[, Reference_Threshold := as.numeric(gsub(",",".", Reference_Threshold))]

files <- list.files("D:/Mes Donnees/Hymobio/DataBase_treatment/Seuils/Genomig_Output")

load("../HYMOBIO_FULLDATA_202405.RData")

toremove <- "AUC_threshold_2405_Quant95_5step_"

files <- grep(toremove, files, value = T)
Cat <- left_join(data.table(Names = files, Cat = gsub(toremove, "", files), Name_Origine = gsub(toremove, "", files)),
                 Titles[,.(Name_Origine, Description, Category)])

CommunityCols <- c("Fish" = "cornflowerblue", "Macroinvertebrate" = "tan2", "Diatom" = "chartreuse4")
CategoryCols <- c("Natural Environment" = "yellowgreen","Physico-Chemistry" = "lightgoldenrod","Temperature" = "#DF536B", 
            "Hydrology" = "#619CFF", "Hydromorphology" = "lightcyan2","Flow Barriers" = "antiquewhite4", "Anthropic Environment"="tomato4")

#CategoryCols_ForArticle <- c("Natural Environment" = "springgreen4", "Physico-chemistry" = "yellowgreen", "Temperature" = "#DF536B", 
                 # "Hydrology" = "#00BFC4", "Hydromorphology" = "#619CFF","Flow Barriers" = "#F5C710", "Anthropic Environment" = "tomato4"
SelectTHR <- function(x, Dir) {
  return(ifelse(length(x)==1, x, switch(unique(Dir),
                                        "1" = min(x),
                                        "-1" = max(x),
                                        "0" = mean(x))))
}

#####
load("FinalTable")
load("variable_thresholds")
Thresholds[, minTHR := apply(Thresholds[, !"Name_Origine", with = F],1,min)][
  , maxTHR := apply(Thresholds[, !"Name_Origine", with = F],1,max)]

Table <- merge(Table,Thresholds[,.(Name_Origine, minTHR, maxTHR)])[, DegradationDirection := as.character(DegradationDirection)][
  , Threshold := as.numeric(Threshold)]

#Overall number of significant models
nrow(unique(Table[,.(Name_Origine, Community)]))
unique(setdiff(
  do.call(paste0, unique(Table[AUC_valTest >= 0.7 & AUC_valTest < 0.8 & F1 >= 0.7,.(Name_Origine, Community)])),
  do.call(paste0, unique(Table[AUC_valTest >= 0.8 & AUC_valTest < 0.9 & F1 >= 0.7,.(Name_Origine, Community)]))))
#Percentage of significantly impacting variables
uniqueN(Table[AUC_valTest >= 0.7  & F1 >= 0.7,Name_Origine])/uniqueN(Table$Name_Origine)

#Number of significant variables by community
#unique(Table[AUC_valTest >= 0.7 & F1 >= 0.7, .(Name_Origine, Community)])[,.N, by = Community]

# Get max AUC and corresponding F1 and Threshold by community
Threshold_max <- Table[AUC_valTest >= 0.7 & F1 >= 0.7, c("AUCmaxTest","F1maxAUC") := 
                         list(max(AUC_valTest[which(F1>=0.7)]), 
                              max(F1[AUC_valTest == max(AUC_valTest[which(F1>=0.7)])])),
                       by = c("Name_Origine","Community")]
#Threshold_max <- Threshold_max[complete.cases(Threshold_max)]
Threshold_max <- Threshold_max[!is.na(AUCmaxTest),THRmaxAUC := SelectTHR(Threshold[AUC_valTest == AUCmaxTest & F1 == F1maxAUC], DegradationDirection), by = c("Name_Origine","Community")]



#Mean AUC by community, all models
#Threshold_max[, mean(AUC_valTest), by = Community]

Threshold_max <- unique(Threshold_max[#AUCmaxTest >= 0.7 & F1maxAUC >= 0.7
  ,.(Name_Origine, Community, minTHR, maxTHR, Description, Category, AUCmaxTest, F1maxAUC, THRmaxAUC)])

Threshold_max[, Category := factor(Category, levels = names(CategoryCols))]
Threshold_max[, Community := factor(Community, levels = names(CommunityCols))]

#Mean AUC by community, only significant
#Threshold_max[, mean(AUCmaxTest), by = Community]

Threshold_max[, MaxComm := Community[which.max(AUCmaxTest)], by = Name_Origine]
#Most responsive communities
#Nmaincom <- unique(Threshold_max[!is.na(AUCmaxTest),.(MaxComm, Name_Origine)])[,.N, by = MaxComm][,Nperc := N/sum(N)]
#table(Threshold_max[!is.na(AUCmaxTest),.(MaxComm, Category)])


# Description of significant models 
#unique(Table[AUC_valTest >= 0.9 & F1 >= 0.7, .(Name_Origine, Community)])
#unique(Table[AUC_valTest >= 0.8 & AUC_valTest < 0.9 & F1 >= 0.7, .(Name_Origine, Community)])
#unique(Table[AUC_valTest >= 0.7 & F1 >= 0.7, .(Name_Origine, Community)])

#N_respCom <- table(unique(Table[AUC_valTest >= 0.7 & F1 >= 0.7, .(Name_Origine, Community)])[,Name_Origine])
#N_respCom

#unique(setdiff(
  #do.call(paste, unique(Table[AUC_valTest >= 0.7  & F1 >= 0.7,.(Category,Description)])),
  #do.call(paste, unique(Table[AUC_valTest >= 0.8 & F1 >= 0.7,.(Category,Description)]))))

CommunityStats <- unique(Table[AUC_valTest >= 0.7  & F1 >= 0.7, .(Name_Origine, Community)])
uniqueN(CommunityStats[Community == "Diatom", Name_Origine])/uniqueN(CommunityStats$Name_Origine)

## Radar graphic
#####
Ncat <- Table[, uniqueN(Name_Origine), by = Category]
setnames(Ncat, "V1", "Ncat")
PradarComp <- merge(Threshold_max[, Nsignif_byCat := length(which(AUCmaxTest >= 0.7 & F1maxAUC >= 0.7)),
                                  by = c("Category", "Community")][, .(Nsignif_byCat, Category, Community)],
                    Ncat)[, Nsignif_byCat := Nsignif_byCat/Ncat]

PradarComp <- dcast(unique(PradarComp), Community ~ Category, value.var = "Nsignif_byCat")
PradarComp[is.na(PradarComp)] <- 0

PlotradarComp <- ggradar(background.circle.colour = "white", PradarComp, legend.position = "none", plot.extent.x.sf = 2,
        axis.label.size = 3.5, grid.label.size = 3.5, group.line.width = 1, group.point.size = 3, fill = T,
        fill.alpha = 0.1, group.colours = CommunityCols, label.gridline.min = F, label.gridline.mid = F,
        label.gridline.max = F) +
  theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5, face = "bold")) +
  ggtitle("Proportion of significant models by category\n") 
#####


## Summary table
#####
CommunityStats <- unique(Threshold_max[AUCmaxTest >= 0.7 & F1maxAUC >= 0.7, .(Name_Origine, Category, Description, Community, AUCmaxTest)])
CommunityStats[,Nsignif_community := .N, by = "Name_Origine"]
uniqueN(CommunityStats[Nsignif_community == 3, Name_Origine])
table(CommunityStats[Nsignif_community == 2, Community])
table(CommunityStats[Nsignif_community == 3, Category])

uniqueN(CommunityStats[Nsignif_community == 2, Name_Origine])

table(unique(CommunityStats[Nsignif_community == 3, .(Name_Origine, Category)])[,Category])

write.table(Threshold_max[!is.na(F1maxAUC)][order(-AUCmaxTest)], file = "tablemax.csv", sep = ";", row.names = F)
#####


# Plot AUC range
#####
Plots <- lapply(levels(Threshold_max$Community), function(comm){
  
return(ggplot(Threshold_max[AUCmaxTest >= 0.7 & F1maxAUC >= 0.7 & Community == comm], 
              aes(x = reorder(Name_Origine, -as.numeric(AUCmaxTest)), y = AUCmaxTest, fill = Category)) +
         geom_col(position = position_dodge(0.6), width = 0.9) + scale_fill_manual(values = CategoryCols) +
         geom_hline(yintercept = c(0.6, 0.7, 0.8), color = c("darkgrey", "ivory4", "black"), lty = 1) + 
         coord_cartesian(ylim=c(0.4, 1)) +  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1))  +
         theme_classic() + labs(x = "", y = "AUCmax")  + ggtitle(comm) + 
         theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)))

})

#PiePlot[, Ncorrect := Ncorrect/count]
ForPlot <- Threshold_max[!is.na(AUCmaxTest), .N, by = c("Category", "Community")]
ForPlot[, Total := sum(N), by = "Community"][
  , FacetName := do.call(paste, c(.SD, sep = "\nN = ")),.SDcols = c("Community", "Total")][
  , FacetName := factor(ForPlot$FacetName, 
                 unique(ForPlot[order(Community),.(FacetName, Community)])[["FacetName"]])
  ]

PiePlot <- ggplot(ForPlot, aes(x="", y=N, fill=Category)) + 
  geom_bar(stat="identity", width=1, position = "fill") + coord_polar("y") + theme_void() + 
  scale_fill_manual(values = CategoryCols) + 
  facet_grid(. ~ FacetName) + 
  ggtitle("Distribution of significant models by community\n") + 
  theme(legend.position = "right", plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))

Plots[[4]] <- PiePlot

ggarrange(plotlist = Plots, nrow = 4, common.legend = T, legend = "right")

Threshold_max[, AUCmax_cat := AUCmaxTest[which.max(AUCmaxTest)], by = Name_Origine][,
  Name_comm := paste0(Name_Origine, Community)]
Gplot_cat <- ggplot(Threshold_max[AUCmaxTest >= 0.7 & F1maxAUC >= 0.7], 
                    aes(x = reorder(Name_comm, -as.numeric(AUCmaxTest)), y = as.numeric(AUCmaxTest), fill = Category)) +
  geom_col(position = position_dodge(0.6), width = 0.9) + scale_fill_manual(values = CategoryCols) +
  geom_hline(yintercept = c(0.7, 0.8), color = c( "ivory4", "black"), lty = 1) + 
  coord_cartesian(ylim=c(0.4, 1)) +  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1))  +
  theme_classic() + labs(x = "", y = "AUC")  + ggtitle("Range by category") + 
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))

Gplot_comp <- ggplot(Threshold_max[AUCmaxTest >= 0.7 & F1maxAUC >= 0.7], 
                     aes(x = reorder(Name_comm, -as.numeric(AUCmaxTest)), y = as.numeric(AUCmaxTest), fill = Community)) +
  geom_col(position = position_dodge(0.6), width = 0.9) + scale_fill_manual(values = CommunityCols) +
  geom_hline(yintercept = c(0.7, 0.8), color = c("ivory4", "black"), lty = 1) + 
  coord_cartesian(ylim=c(0.4, 1)) +  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1))  +
  theme_classic() + labs(x = "", y = "AUC")  + ggtitle("Range by community") + 
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) #+ 
  #geom_text(aes(x = uniqueN(Threshold_max$Param_comm) - 10, y = 0.8, label = "Poor", vjust = 1, hjust = 0), color = "darkgrey") +
  #geom_text(aes(x = uniqueN(Threshold_max$Param_comm) - 10, y = 0.7, label = "Correct", vjust = 1, hjust = 0), color = "ivory4") +
  #geom_text(aes(x = uniqueN(Threshold_max$Param_comm) - 10, y = 0.8, label = "Good", vjust = 1, hjust = 0), color = "black") 

ggarrange(Gplot_cat, Gplot_comp, nrow = 2, labels = c("a)", "b)"))

ggarrange(ggarrange(Gplot_comp,PlotradarComp,  nrow = 2, heights = c(1,0.8)),
          ggarrange(Gplot_cat, PiePlot, common.legend = T, nrow = 2, legend = "right"), nrow = 2)
          
png("../../Redaction/Seuils/Fig4.png", width = 600, height = 1000)

ggarrange(ggarrange(Gplot_comp,PlotradarComp, common.legend = T, nrow = 2, legend = "right", heights = c(1,1.3)),
          NULL,
          ggarrange(Gplot_cat, PiePlot, common.legend = T, nrow = 2, legend = "right"), nrow = 3, heights = c(1,0.1,1))

dev.off()
#####


## Difference between literature and observed thresholds
#####
DiffThresh <- merge(Threshold_max, Titles[,.(Name_Origine, Reference_Threshold)])[
  Reference_Threshold != "",][, Reference_Threshold := as.numeric(Reference_Threshold)][
    , THRmaxAUC := as.numeric(THRmaxAUC)]

DiffThresh <- DiffThresh[!is.na(AUCmaxTest)][, DiffThresh := THRmaxAUC - Reference_Threshold][,
               PosThresh := ifelse(DiffThresh < 0, "Below", ifelse(DiffThresh > 0, "Above", "Equal"))][
               , PosThresh := factor(PosThresh, levels = c("Above", "Equal", "Below"))]

DiffThresh_tab <-table(DiffThresh[, .(Category, PosThresh)])
round(colSums(DiffThresh_tab)/sum(colSums(DiffThresh_tab)), 2)

DiffThresh <- unique(melt(setDT(DiffThresh[,.(Category,Above, Equal, Below)]), id.vars = "Category"))


Pdiff <- ggplot(DiffThresh[, .N, by = c("Category", "PosThresh")], aes(x = Category, y = N, fill = PosThresh)) + 
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = data.frame(Above = "#619CFF", Equal = "grey", Below = "firebrick")) + 
  labs(y = "Threshold comparison" , fill = "Position to threshold") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0)) 

+
  ggtitle("(a)") 
#####


#Multiple scales comparions
#####
Multiscale <- unique(Titles[!Category %in% c("Temperature", "Hydrology")][which(duplicated(Name)),Name])
Multiscale <- unique(Titles[Name %in% Multiscale, .(Name_Origine, Scale)])
Threshold_max <- left_join(Threshold_max, Multiscale)[!is.na(Scale)]

ScaleComp <- Threshold_max[grep("Upstream", Scale), Scale := "Upstream"][
  grep("Downstream", Scale), Scale := "Downstream"][grep("Downstream", Scale), Scale := "Downstream"][
    grep("USRA buffer 100", Scale), Scale := "USRA 100"][grep("USRA buffer 30", Scale), Scale := "USRA 30"][
        grep("USRA floodplain", Scale), Scale := "USRA floodplain"][grep("Proximal", Scale), Scale := "Proximal linear"][
          grep("Hydro", Scale), Scale := "Hydro area"][grep("Station", Scale), Scale := "Station"]#[
            #AUCmaxTest >= 0.7 & F1maxAUC >= 0.7]
ScaleComp[, Scale := factor(Scale, 
        levels = c("Station", "Proximal linear", "USRA floodplain", "USRA 30", "USRA 100", "Hydro area", "Connected river segments", "Downstream", "Upstream"))]


ScaleComp <- left_join(Table, Multiscale)[!is.na(Scale)]

ScaleComp <- ScaleComp[grep("Upstream", Scale), Scale := "Upstream"][
  grep("Downstream", Scale), Scale := "Downstream"][grep("Downstream", Scale), Scale := "Downstream"][
    grep("USRA buffer 100", Scale), Scale := "USRA 100"][grep("USRA buffer 30", Scale), Scale := "USRA 30"][
      grep("USRA floodplain", Scale), Scale := "USRA floodplain"][grep("Proximal", Scale), Scale := "Proximal linear"][
        grep("Hydro", Scale), Scale := "Hydro area"][grep("Station", Scale), Scale := "Station"]#[
#AUCmaxTest >= 0.7 & F1maxAUC >= 0.7]
ScaleComp[, Scale := factor(Scale, 
                            levels = c("Station", "Proximal linear", "USRA floodplain", "USRA 30", "USRA 100", "Hydro area", "Connected river segments", "Downstream", "Upstream"))]



ScaleComp <- ScaleComp[!Scale %in% c("All serie")][, meanAUC := mean(AUC_valTest), by = Name_Origine]
ScaleComp <- unique(ScaleComp[, .(Scale, Name_Origine, meanAUC)])[!is.na(Scale)]

pairwise.wilcox.test(ScaleComp$meanAUC, ScaleComp$Scale, p.adjust.method="none")

Pscale <- ggplot(ScaleComp, aes(x = Scale, y = meanAUC, fill = Scale)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge()) +
  stat_summary(fun.data = mean_se, geom = "errorbar", linewidth = 0.5, color = "darkgray") +
  scale_fill_brewer(palette = "YlOrRd") + theme_classic() + labs(y = "Mean AUC test", x = "Scale") +
  theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        plot.title = element_text(hjust = 0)) + ggtitle("(b)")

grid.arrange(Pdiff, Pscale, nrow = 1)
#####


# Figures hydromorphology
#####
Threshold_maxHM <- Threshold_max[Category == "Hydromorphology" & !is.na(AUCmaxTest)][
  , Name_comm := paste0(Description, Community)]
Threshold_maxHM <- merge(Threshold_maxHM, Titles[,.(Name_Origine, Unit)])

ggplot(Threshold_maxHM, aes(x = reorder(Name_comm, -as.numeric(AUCmaxTest)), y = as.numeric(AUCmaxTest), fill = Community)) +
  geom_col(position = position_dodge(0.6), width = 0.9) + scale_fill_manual(values = CommunityCols) +
  geom_hline(yintercept = c(0.7, 0.8), color = c("ivory4", "black"), lty = 1) + 
  coord_cartesian(ylim=c(0.4, 1)) +  scale_y_continuous(breaks = seq(0.4, 1, by = 0.1))  +
  theme_classic() + labs(x = "", y = "AUC")  + ggtitle("Hydromorphology variables") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(labels= Threshold_maxHM[order(-as.numeric(AUCmaxTest)), Description])


HM_var <- unique(Threshold_maxHM[, .(Description, Community, AUCmaxTest, F1maxAUC, THRmaxAUC, Unit)])
HM_var[, THRmaxAUC := round(THRmaxAUC, 2)]
write.table(HM_var, file = "HMvariables.csv", sep=";", row.names = F)





