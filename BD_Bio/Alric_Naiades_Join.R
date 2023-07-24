invisible(lapply(c("data.table", "rgeos", "sp", "ggplot2", "viridis", "gridExtra"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("../../DataBase_treatment/BD_Bio")

load("Naiades_Fish"); load("Naiades_Macroinvertebrate"); load("Naiades_Diatom")
Nai_Station <- as.data.table(read.csv("../../../../DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/stations.csv", 
               sep=";"))[,.(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, 
                            CodeRegion, LbRegion)][!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")),]

source("Alric_Fish.R"); source("Alric_Macroinvertebrate.R"); source("Alric_Diatom.R")


#
Abs_NaiFish <- data.table(NomLatinAppelTaxon = setdiff(Nai_Abd_Fish$NomLatinAppelTaxon, Alr_traits_Fish$sciname))[
  , tax.GENUS := sub(" .*", "", NomLatinAppelTaxon)]

Corr_Fish <- Correspondance_Fish(Alric_traits = Alr_traits_Fish, Alric_Abd = Alr_Abd_Fish, 
                               Naiades_Abd = Nai_Abd_Fish, Naiades_Station = Nai_Station, Abs_Nai = Abs_NaiFish)
save(Corr_Fish, file = "Correspondance_Fish")
load("Correspondance_Fish")

names(Corr_Fish) <- unique(Abs_NaiFish$tax.GENUS)
No_match_Fish <- merge(data.table(tax.GENUS = names(Corr_Fish)[which(unlist(lapply(Corr_Fish, length))==0)]),
                  Abs_NaiFish)
Matched_Fish <- Corr_Fish[which(unlist(lapply(Corr_Fish, length))!=0)]


#
load("Naiades_Macroinvertebrate")

Abs_Nai_Macroinv <- data.table(NomLatinAppelTaxon = setdiff(Nai_Abd_Macroinv$NomLatinAppelTaxon, Alr_Names_Macroinv$TAXON))[
  , tax.GENUS := sub(" .*", "", NomLatinAppelTaxon)]

Corr_Macroinv <- Correspondance_Macroinv(Alric_traits = Alr_Names_Macroinv, Alric_Abd = Alr_Abd_Macroinv, 
                              Naiades_Abd = Nai_Abd_Macroinv, Naiades_Station = Nai_Station, Abs_Nai = Abs_Nai_Macroinv)
save(Corr_Macroinv, file = "Correspondance_Macroinvertebrate2")
load("Correspondance_Macroinvertebrate")

names(Corr_Macroinv) <- unique(Abs_Nai_Macroinv$tax.GENUS)
No_match_Macroinv <- merge(data.table(tax.GENUS = names(Corr_Macroinv)[which(unlist(lapply(Corr_Macroinv, length))==0)]),
                       Abs_Nai_Macroinv)
Matched_Macroinv <- Corr_Macroinv[which(unlist(lapply(Corr_Macroinv, length))!=0)]

#
Abs_Nai_Diatom <- data.table(NomLatinAppelTaxon = setdiff(Nai_Abd_Diatom$NomLatinAppelTaxon, Alr_traits_Diatom$SCIENTIFIC_NAME))[
  , tax.GENUS := sub(" .*", "", NomLatinAppelTaxon)]

Corr_Diatom <- Correspondance_Diatom(Alric_traits = Alr_traits_Diatom, Alric_Abd = Alr_Abd_Diatom, 
                                Naiades_Abd = Nai_Abd_Diatom, Naiades_Station = Nai_Station, Abs_Nai = Abs_Nai_Diatom)
save(Corr_Diatom, file = "Correspondance_Diatom")

load("Correspondance_Diatom")
names(Corr_Diatom) <- unique(Abs_Nai_Diatom$tax.GENUS)
No_match_Diatom <- merge(data.table(tax.GENUS = names(Corr_Diatom)[which(unlist(lapply(Corr_Diatom, length))==0)]),
                       Abs_Nai_Diatom)
Matched_Diatom <- Corr_Diatom[which(unlist(lapply(Corr_Diatom, length))!=0)]


##

Nai_Abd <- Nai_Abd_Diatom; Matched <- Matched_Diatom; Abs_Nai <- Abs_Nai_Diatom; group <- "diatoms"

NomLatin_Join <- function(Nai_Abd, Matched, Abs_Nai, group){
  
  invisible(ifelse(group == "fish", NomRef <- "sciname", ifelse(group == "macroinv", NomRef <- "TAXON",
                                                                NomRef <- "SCIENTIFIC_NAME")))
  
  Nai_Abd_matched <- Nai_Abd[, NomLatin_Join := NomLatinAppelTaxon]
  
for(i in 1:length(Matched)){
  Gen <- Matched[[i]]
  if(is.character(Gen)){
    Nai_Abd_matched[NomLatinAppelTaxon %in% Abs_Nai[tax.GENUS == Gen["Nai.GENUS"], NomLatinAppelTaxon], 
                    NomLatin_Join := Gen["Alr.SPECIES"]]
  }
  else{ if(!is.data.frame(Gen)){print(paste("Problem with", Gen))}
    Nai_Abd_matched[Gen, on = .(CdStationMesureEauxSurface, NomLatinAppelTaxon), NomLatin_Join := get(paste0("i.",NomRef))]
  }
}
  return(Nai_Abd_matched)
}

Fish_Joined <- NomLatin_Join(Nai_Abd = Nai_Abd_Fish, Matched = Matched_Fish, Abs_Nai = Abs_NaiFish, group = "fish")
Fish_Joined[NomLatinAppelTaxon %in% No_match_Fish$NomLatinAppelTaxon, NomLatin_Join := "Unmatched"]
write.table(Fish_Joined, file = "Correspondance_Fish.csv", sep=";", row.names = F)
save(Fish_Joined, file = "Correspondance_Fish")

lsNaiades_F <- rbind(do.call(rbind, lapply(Matched_Fish, function(sp){
  if(is.character(sp)){
    ret <- t(data.frame(sp)); colnames(ret) <- c("NomLatinAppelTaxon_Naiades", "NomLatin_Join")
    return(ret)} 
  else { ret <- unique(sp[,.(NomLatinAppelTaxon, sciname)]); colnames(ret) <- c("NomLatinAppelTaxon_Naiades", "NomLatin_Join")
  return(ret)}
})), data.table(NomLatinAppelTaxon_Naiades = No_match_Fish$NomLatinAppelTaxon)[, NomLatin_Join := "Unmatched"])[
  , NomLatinAppelTaxon_Naiades := as.character(NomLatinAppelTaxon_Naiades)][base::order(NomLatinAppelTaxon_Naiades)]
write.table(lsNaiades_F, file = "Resume_Fish.csv", sep=";", row.names = F)



Macroinv_Joined <- NomLatin_Join(Nai_Abd = Nai_Abd_Macroinv, Matched = Matched_Macroinv, Abs_Nai = Abs_Nai_Macroinv, group = "macroinv")
Macroinv_Joined[NomLatinAppelTaxon %in% No_match_Macroinv$NomLatinAppelTaxon, NomLatin_Join := "Unmatched"]
write.table(unique(Macroinv_Joined[,.(NomLatinAppelTaxon, NomLatin_Join)]), file = "Correspondance_Macroinvertebrate.csv", sep=";", row.names = F)
save(Macroinv_Joined, file = "Correspondance_Macroinvertebrate")

lsNaiades_Macroinv <- rbind(do.call(rbind, lapply(Matched_Macroinv, function(sp){
  if(is.character(sp)){
    #ret <- t(data.frame(sp)); colnames(ret) <- c("NomLatinAppelTaxon_Naiades", "NomLatin_Join")
    return(Alr_Names_Macroinv[grep(sp, TAXON), TAXON])} 
  else { ret <- unique(sp[,.(NomLatinAppelTaxon, TAXON)]); colnames(ret) <- c("NomLatinAppelTaxon_Naiades", "NomLatin_Join");return(ret)}
})), data.table(NomLatinAppelTaxon_Naiades = No_match_Macroinv$NomLatinAppelTaxon)[, NomLatin_Join := "Unmatched"])[
  , NomLatinAppelTaxon_Naiades := as.character(NomLatinAppelTaxon_Naiades)][base::order(NomLatinAppelTaxon_Naiades)]
write.table(lsNaiades_Macroinv, file = "Resume_macroinvertebrate.csv", sep=";", row.names = F)





Diatom_Joined <- NomLatin_Join(Nai_Abd = Nai_Abd_Diatom, Matched = Matched_Diatom, Abs_Nai = Abs_Nai_Diatom, group = "diatomees")
Diatom_Joined[NomLatinAppelTaxon %in% No_match_Diatom$NomLatinAppelTaxon, NomLatin_Join := "Unmatched"]

lsNaiades_Diatom <- do.call(rbind, lapply(Matched_Diatom, function(sp){
  if(is.character(sp)){
    ret <- t(data.frame(sp)); colnames(ret) <- c("NomLatinAppelTaxon_Naiades", "NomLatin_Join")
    return(ret)} 
  else { 
    ret <- unique(sp[,.(NomLatinAppelTaxon, SCIENTIFIC_NAME)]); colnames(ret) <- c("NomLatinAppelTaxon_Naiades", "NomLatin_Join")
    return(ret)}
}))[, NomLatinAppelTaxon_Naiades := as.character(NomLatinAppelTaxon_Naiades)][base::order(NomLatinAppelTaxon_Naiades)]
lsNaiades_Diatom <- lsNaiades_Diatom[!grep("var.|abnormal|ssp.|subsp.|f.|morphotype", NomLatinAppelTaxon_Naiades)][,
                               tax.GENUS := sub(" .*", "", NomLatinAppelTaxon_Naiades)]
lsNaiades_Diatom <- lsNaiades_Diatom[NomLatinAppelTaxon_Naiades != tax.GENUS]
write.table(lsNaiades_Diatom, file = "Resume_diatoms.csv", sep=";", row.names = F)


Plot_props <- function(Group_joined, No_match, Abs_Nai, group){
  
  invisible(ifelse(group == "fish", Title <- "Fish", ifelse(group == "macroinv", Title <- "Macroinvertabrate",
                                                            Title <- "Diatoms")))
  
  Group_joined <- Group_joined[MnTypTaxRep == "NbrTax"]
  
  Ntax <- sum(Group_joined$RsTaxRep)
  NonMatching <- Group_joined[,Matching := "Matching"][NomLatinAppelTaxon %in% Abs_Nai[,NomLatinAppelTaxon], Matching := "Replaced"][
    NomLatinAppelTaxon %in% No_match$NomLatinAppelTaxon, Matching := "Unmatched"][ , Prop := sum(RsTaxRep)/Ntax, by = Matching][
      ,.(Matching, NomLatinAppelTaxon, Prop, RsTaxRep)]
  
  p1 <- ggplot(unique(NonMatching[,.(Matching, Prop)]), aes(x = "", y = Prop, fill = Matching)) +
    geom_col(color = "black") + coord_polar(theta = "y") +
    scale_fill_manual(values = c("Matching" = "olivedrab", 
                                 "Replaced" = "darkorange",
                                 "Unmatched" = "darkred")) +
    theme_void()+ theme(legend.position = "bottom", plot.margin = unit(c(0,0,3,0), "lines"),legend.text=element_text(size=8)) +
    guides(fill = guide_legend(title="")) + ggtitle("(a) Proportions of taxa sampled for abundance")
  
  Nnomatch <- sum(NonMatching[Matching != "Matching",RsTaxRep])
  NonMatching.rel <- NonMatching[Matching != "Matching"][,Prop := NULL][Matching == "Replaced", 
                     Matching := NomLatinAppelTaxon][, NomLatinAppelTaxon := NULL][
                      ,Prop.rel := sum(RsTaxRep)/Nnomatch, by = Matching] 
  if(group != "diatoms"){
    p2 <- ggplot(unique(NonMatching.rel[,.(Matching, Prop.rel)]),aes(x = "", y = Prop.rel, fill = Matching)) +
    geom_col(color = "black") + coord_polar(theta = "y") + 
    scale_fill_viridis(discrete = TRUE, option = "inferno") + guides(fill = guide_legend(title="")) +
    theme_void()+ theme(legend.position = "bottom", plot.margin = unit(c(0,0,0,0), "lines"),
                        legend.text=element_text(size=8),legend.key.size = unit(0.6, 'lines')) +
    ggtitle("(b) Distribution of replaced and unmatched samples")}

  else {
    NonMatching.rel <- NonMatching[Matching != "Matching"][,Prop := NULL][, Genus := sub(" .*", "", NomLatinAppelTaxon)][
      Matching == "Replaced", Matching := Genus][, NomLatinAppelTaxon := NULL][, Genus := NULL][ ,
      Prop.rel := sum(RsTaxRep)/Nnomatch, by = Matching]
    
    
      p2 <- ggplot(unique(NonMatching.rel[,.(Matching, Prop.rel)]),aes(x = "", y = Prop.rel, fill = Matching)) +
      geom_col(color = "black") + coord_polar(theta = "y") + 
      scale_fill_viridis(discrete = TRUE, option = "inferno") + guides(fill = guide_legend(title="")) +
      theme_void()+ theme(legend.position = "bottom", plot.margin = unit(c(0,0,0,0), "lines"),
                          legend.text=element_text(size=8),legend.key.size = unit(0.8, 'lines')) +
      ggtitle("Distribution of replaced and unmatched samples\nby genus")
  }
  
  grid.arrange(p1, p2, ncol = 2, widths = c(1/3, 2/3), top = Title)
}

Plot_props(Group_joined = Fish_Joined, No_match = No_match_Fish, Abs_Nai = Abs_Nai_Fish, group = "fish")

Plot_props(Group_joined = Macroinv_Joined, No_match = No_match_Macroinv, Abs_Nai = Abs_Nai_Macroinv, group = "macroinv")

Plot_props(Group_joined = Diatom_Joined, No_match = No_match_Diatom, Abs_Nai = Abs_Nai_Diatom, group = "diatoms")




