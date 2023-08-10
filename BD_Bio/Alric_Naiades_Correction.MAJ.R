invisible(lapply(c("data.table", "rgeos", "sp", "dplyr", "gridExtra", "lubridate", "stringr", "ggplot2", "ggpubr"),function(pk){
  if(!pk %in% row.names(installed.packages())){install.packages(pk)}
  library(pk,character.only=T)}))

setwd("C:/Users/armirabel/Documents/INRAE/Hymobio/DataBase_treatment/BD_Bio")

#####
Resume_graphique <- function(Alric, Naiades, Group){
  
  par(mfrow=c(2,2), mar=c(1,1,4,1))
  Naiplot <- Nai_Station[CdStationMesureEauxSurface %in% Naiades$CdStationMesureEauxSurface]
  coordinates(Naiplot) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
  plot(Naiplot, cex = 0.1, pch = 20, main = "All Naiades", cex.main = 0.8)
  mtext(paste0(min(Naiades$year), " - ", max(Naiades$year), "\n", uniqueN(Naiplot$CdStationMesureEauxSurface), " Stations"), side = 3, cex = 0.6, outer = F)
  
  Alrplot <- Alric
  coordinates(Alrplot) <- ~ x + y
  plot(Alrplot, main = "All Alric", cex = 0.1, pch = 20, cex.main = 0.8)
  mtext(paste0(min(Alrplot$year), " - ", max(Alrplot$year), "\n", uniqueN(Alric$cd_site), " Stations"), side = 3, cex = 0.6, outer = F)
  
  diff_Alr_Nai <- Naiades[CdStationMesureEauxSurface %in% setdiff(Naiades$CdStationMesureEauxSurface, Alric$cd_site)][year %in% Alric$year]
  diff_Alr_Nai <- Nai_Station[CdStationMesureEauxSurface %in% diff_Alr_Nai$CdStationMesureEauxSurface]
  coordinates(diff_Alr_Nai) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
  plot(diff_Alr_Nai, main = "In Naiades not in Alric", cex = 0.1, pch = 20, cex.main = 0.8)
  mtext(paste0(min(Alric$year), " - ", max(Alric$year), " (same period)\n", uniqueN(diff_Alr_Nai$CdStationMesureEauxSurface), " Stations"), side = 3, cex = 0.6, outer = F)
  
  diff_Nai_Alr <- Alric[cd_site %in% setdiff(Alric$cd_site, Naiades$CdStationMesureEauxSurface)]
  coordinates(diff_Nai_Alr) <- ~ x + y
  plot(diff_Nai_Alr, main = "In Alric not in Naiades (any time)", cex = 0.1, pch = 20, cex.main = 0.8)
  mtext(paste0(min(Alric$year), " - ", max(Alric$year), "\n", uniqueN(diff_Nai_Alr$cd_site), " Stations"), side = 3, cex = 0.6, outer = F)
  
  mtext(Group, side = 3, line = -1, outer = T)
  mtext(paste(length(unique(intersect(Alric$cd_site, Naiades$CdStationMesureEauxSurface)))," common stations"), 
        side = 3, line = -2, col = "red", outer = T, cex = 0.7)
  
  Naiades$Common <- "Naiades only"
  Stack <- unique(unique(Naiades[,.(CdStationMesureEauxSurface, year, Common)])[
    CdStationMesureEauxSurface %in% intersect(Naiades$CdStationMesureEauxSurface, Alric$cd_site), Common := "Common stations"][
      ,.(year, Common)][, stack := .N, by = c("year","Common")])
  
  
  Pstack <- ggplot(Stack, aes(x=year, y=stack, fill=Common)) + geom_area() + theme_minimal() + 
    labs(x = "Year", y = "N stations") + ggtitle(Group) +
    scale_fill_manual(values = c("Common stations" = "olivedrab",  "Naiades only" = "darkred")) + guides(fill = guide_legend(title=""))
  
  
  PieContinue <- unique(unique(unique(Naiades[CdStationMesureEauxSurface %in% intersect(Naiades$CdStationMesureEauxSurface, 
                               Alric$cd_site), .(CdStationMesureEauxSurface, year)])[, Continue := "Maintained"][,MaxYear := max(year),                                                                                                                                                                  by = "CdStationMesureEauxSurface"][MaxYear <= 2017, Continue := "Stopped"])[, Eff := .N, by = "Continue"][,.(Continue, Eff)])
  
  Ppie <- ggplot(PieContinue,aes(x = "", y=Eff, fill = Continue)) +
    geom_col(color = "black") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Maintained" = "olivedrab",  "Stopped" = "darkred")) +
    theme_void() + guides(fill = guide_legend(title="")) + ggtitle(Group) 
  
  return(list(stack = Pstack, pie = Ppie))
}


AlricF <- Alr_Abd_Fish
setnames(AlricF , c("X_L93", "Y_L93","CODE_STATION"), c("x", "y", "cd_site"))
AlricF[, year := format(as.Date(DATE, "%d/%m/%Y"), "%Y")]

Resume_graphique(Alric = Alr_Abd_Diatom, Naiades = Nai_Abd_Diatom, Group = "Diatoms")
#####


# Macroinvertebrates

#####
load("Naiades_Macroinvertebrate"); source("C:/Users/armirabel/Documents/INRAE/Hymobio/Hymobio_GitHub/BD_Bio/Alric_Macroinvertebrate.R")

Nai_Station <- as.data.table(read.csv("C:/Users/armirabel/Documents/DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/stations.csv", 
                                      sep=";"))[,.(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, CodeRegion, LbRegion)]
Stations_metropole <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")), CdStationMesureEauxSurface]
Nai_Station <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")),]

Corr_philippe <- data.table(read.table("Correspondance_Macroinvertebrate.csv", sep = ";", quote = "\"", header = T))[,
                 .(NomLatinAppelTaxon_Naiades, NomLatin_Join, Propositions)]
setnames(Corr_philippe, "NomLatinAppelTaxon_Naiades", "NomLatinAppelTaxon")

# Format the initial Naiades base and make replacements from Philippe's suggestions
Nai_Abd_Macroinv_maj <- Nai_Abd_Macroinv[, "NomLatin_Join" := NomLatinAppelTaxon]

for(i in which(sapply(Nai_Abd_Macroinv_maj, class ) == "factor")) Nai_Abd_Macroinv_maj[[i]] = as.character(Nai_Abd_Macroinv_maj[[i]])

Nai_Abd_Macroinv_maj <- unique(Nai_Abd_Macroinv_maj[, RsTaxRep := sum(RsTaxRep), by = c("CdStationMesureEauxSurface", "year", "NomLatinAppelTaxon")])
Nai_Abd_Macroinv_maj[Corr_philippe[Propositions == "OK",], on = .(NomLatinAppelTaxon), NomLatin_Join := i.NomLatin_Join]
Nai_Abd_Macroinv_maj[Corr_philippe[!grep("Ne pas tenir compte|OK|Trop large",Propositions)], on = .(NomLatinAppelTaxon), NomLatin_Join := i.Propositions]
Nai_Abd_Macroinv_maj <- Nai_Abd_Macroinv_maj[!NomLatin_Join %in% Corr_philippe[grep("Ne pas tenir compte|Trop large",Propositions), NomLatinAppelTaxon]]

# integrate additional traits profiles, attribute a species code
Traits_complement <- fread("C:/Users/armirabel/Documents/DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/BIOLOGICAL_TRAITS_TREATMENTS/1_MACROINVERTEBRATES/2_BASE_AMOBIO/ComplementSpeciesTraits_Macroinvertebrate.csv")[
  ,-c("Identifiant", "Nom_variable", "Libelle", "Attributs", "Description", "Fuzzy")][, Nom_final := gsub("INV_", "INV.", Nom_final)]
CdToAdd <- dcast(melt(unique(Nai_Abd_Macroinv_maj[NomLatin_Join %in% names(Traits_complement)[-1], 
                 .(NomLatin_Join, CdAppelTaxon)]), id.vars = "NomLatin_Join"), variable ~ NomLatin_Join, fun.aggregate = function(x) {min(as.numeric(x))})
setnames(CdToAdd, "variable", "Nom_final")
CdToAdd <- CdToAdd[, names(Traits_complement), with = F]
Traits_complement <- rbind(Traits_complement, CdToAdd)

Alr_Names_Macroinv <- Alr_Names_Macroinv[, Ncode := .N, by = "TAXON"][Ncode > 1, CODE_TAXON := intersect(CODE_TAXON, Alr_Traits_Macroinv$code_taxon_Alric2021), by = "TAXON"]
Alr_Names_Macroinv[TAXON == "Helodes", TAXON := "Elodes"]
Alr_Names_Macroinv <- unique(Alr_Names_Macroinv)

setnames(Alr_Names_Macroinv, "CODE_TAXON", "CdAppelTaxon")
Alr_Names_Macroinv <- rbind(Alr_Names_Macroinv[,.(CdAppelTaxon, TAXON, tax.GENUS, FAMILLE)], dcast(melt(
  Traits_complement[Nom_final %in% c("CdAppelTaxon", "tax.GENUS", "FAMILLE")], id.vars = "Nom_final", variable.name = "TAXON"), 
  TAXON ~ Nom_final))
setnames(Alr_Names_Macroinv, c("CdAppelTaxon","TAXON"), c("CdAppelTaxon_Join","NomLatin_Join"))

for(i in which(sapply(Alr_Names_Macroinv, class ) == "factor")) Alr_Names_Macroinv[[i]] = as.character(Alr_Names_Macroinv[[i]])
                   
# Get the code taxon for Naiades entries
Nai_Abd_Macroinv_maj <- left_join(Nai_Abd_Macroinv_maj, Alr_Names_Macroinv[CdAppelTaxon_Join %in% Alr_Traits_Macroinv$code_taxon_Alric2021], 
                                  by = "NomLatin_Join")[,
  .(CdStationMesureEauxSurface, DateDebutOperationPrelBio, year, tax.GENUS, FAMILLE, CdAppelTaxon, NomLatinAppelTaxon, NomLatin_Join, CdAppelTaxon_Join, RsTaxRep)]
Nai_Abd_Macroinv_maj[Alr_Names_Macroinv[!CdAppelTaxon_Join %in% Alr_Traits_Macroinv$code_taxon_Alric2021], on = "NomLatin_Join", 
                     c("tax.GENUS", "FAMILLE") := list(i.tax.GENUS, i.FAMILLE)]

# Attribute a species from same genus to undetermined samples, based on dendritic & geographical distance + relative abundance
load("dendritic_distance_table_biological_data.Rdata")
DistanceMatrix <- as.data.table(distance_table_biological_data)[COMPARTIMENT_START == "MACROINVERTEBRATE" & COMPARTIMENT_TARGET == "MACROINVERTEBRATE"]

Torep <- unique(Nai_Abd_Macroinv_maj[!CdAppelTaxon_Join %in% c(unlist(CdToAdd[,-1]), Alr_Traits_Macroinv$code_taxon_Alric2021),
                                     .(tax.GENUS, NomLatin_Join, CdAppelTaxon_Join, CdAppelTaxon)])

Nomatch <- data.table()
for(i.species in 1:nrow(Torep)){
  species <- Torep[i.species,]
  
  Cd_TAXON <- intersect(Nai_Abd_Macroinv_maj[tax.GENUS == species$tax.GENUS, CdAppelTaxon_Join], 
                        Alr_Traits_Macroinv$code_taxon_Alric2021)
  
  if(length(Cd_TAXON) != 0) {
  Nai_torep <- Nai_Abd_Macroinv_maj[NomLatin_Join %in% species$NomLatin_Join]
  Ref_Nai <- Nai_Abd_Macroinv_maj[CdAppelTaxon_Join %in% Cd_TAXON, ][, Abce_ref := sum(RsTaxRep), 
               by = c("CdStationMesureEauxSurface", "CdAppelTaxon")]
  Name_ref <- unique(Ref_Nai[,.(NomLatin_Join, CdAppelTaxon_Join)])
    
  idref <- unique(DistanceMatrix[ID_STATION_START %in% Nai_torep$CdStationMesureEauxSurface & ID_STATION_TARGET %in% Ref_Nai$CdStationMesureEauxSurface][, 
                   closest := as.numeric(ID_STATION_TARGET[which.min(distance_m)]), by = "ID_STATION_START"][,.(ID_STATION_START, closest)])
    
  Ref_closest_dendri <- right_join(Nai_torep, idref, by = c("CdStationMesureEauxSurface" = "ID_STATION_START"))
    
    if(any(!Nai_torep$CdStationMesureEauxSurface %in% Ref_closest_dendri$CdStationMesureEauxSurface)){
      Nai_torep_XY <- Nai_Station[CdStationMesureEauxSurface %in% 
                                    setdiff(Nai_torep$CdStationMesureEauxSurface,Ref_closest_dendri$CdStationMesureEauxSurface)]
      coordinates(Nai_torep_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
      
      Ref_Nai_XY <- Nai_Station[CdStationMesureEauxSurface %in% 
                                  setdiff(Ref_Nai$CdStationMesureEauxSurface, Ref_closest_dendri$CdStationMesureEauxSurface)] 
      coordinates(Ref_Nai_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
      
      Ref_closest <- data.frame(closest = Ref_Nai_XY[apply(gDistance(Nai_torep_XY, Ref_Nai_XY, byid=TRUE), 2, which.min),
                                                          ]$CdStationMesureEauxSurface)
      Ref_closest <- cbind(data.table(CdStationMesureEauxSurface = Nai_torep_XY$CdStationMesureEauxSurface), Ref_closest)
      Ref_closest <- right_join(Nai_torep, Ref_closest, by = "CdStationMesureEauxSurface")
      
      Nai_torep <- rbindlist(list(Ref_closest, Ref_closest_dendri), use.names = T)
    } else { Nai_torep <- Ref_closest_dendri; Nai_torep$closest <- as.character(Nai_torep$closest)}
    
    Ref_Nai <- dcast(Ref_Nai, CdStationMesureEauxSurface ~ CdAppelTaxon_Join, value.var = "Abce_ref", fun.aggregate =  sum )
    setnames(Ref_Nai, "CdStationMesureEauxSurface","closest")
    
    Nai_torep <- left_join(Nai_torep, Ref_Nai, by = "closest")
    
    Nai_torep[, c("CdAppelTaxon_Join", "Ninv") := list(names(.SD)[max.col(.SD, ties.method = "first")], 
              length(which(.SD!=0))), .SDcols = Cd_TAXON, by = c("year", "CdStationMesureEauxSurface")]
    
    Nai_torep[Ninv > 1, CdAppelTaxon_Join := sample(names(.SD),1, prob = colSums(.SD)), .SDcols = Cd_TAXON, by = c("year", "CdStationMesureEauxSurface")]
    Nai_torep[Name_ref, on = .(CdAppelTaxon_Join), NomLatin_Join := i.NomLatin_Join]
    
    Nai_Abd_Macroinv_maj[Nai_torep[,.(CdStationMesureEauxSurface, NomLatinAppelTaxon, NomLatin_Join, CdAppelTaxon_Join)], 
                       on = .(CdStationMesureEauxSurface, NomLatinAppelTaxon), 
                       c("CdAppelTaxon_Join", "NomLatin_Join") := list(i.CdAppelTaxon_Join, i.NomLatin_Join)]
      } else { Nomatch <- c(Nomatch, species$NomLatin_Join)}
}

Nai_Abd_Macroinv_maj[NomLatin_Join %in% unlist(Nomatch), NomLatin_Join := NA][NomLatin_Join == "Unmatched", CdAppelTaxon_Join := NA]

# Complete with additional information
setnames(Nai_Abd_Macroinv_maj, c("CdStationMesureEauxSurface", "DateDebutOperationPrelBio", "year", "tax.GENUS", "FAMILLE", "RsTaxRep"),
         c("CdStation", "Date_PrelBio", "Year", "Genus", "Family", "Abundance"))
Nai_Abd_Macroinv_maj[ , BDsource := "Naiades"][, Group := "Macroinvertebrate"][, IDoperation := paste(Group, CdStation, Year, sep = "_")][,MAJ := "Ariane_March2023"]
Nai_Abd_Macroinv_maj[, CdStation := as.character(CdStation)][, Date_PrelBio := as.character(Date_PrelBio)]

Nai_Abd_Macroinv_maj[CdStation %in% Alr_Abd_Macroinv$cd_site, BDsource := "Naiades+Alric"]
# Get additional samplings from Alric
Supp_MacroinvertebratesAlr <- Alr_Abd_Macroinv[cd_site %in% setdiff(Alr_Abd_Macroinv$cd_site, Nai_Abd_Macroinv_maj$CdStation)]
Supp_MacroinvertebratesAlr <- melt(Supp_MacroinvertebratesAlr, id.vars = c("cd_site", "cd_opecont", "date_opecont", "year", "x", "y", "her1",
     "typo_nationale", "rang_corr", "agence"), measure.vars = grep("X",colnames(Supp_MacroinvertebratesAlr), value = T),
     variable.name = "CdAppelTaxon_Join", value.name = "RsTaxRep")[, CdAppelTaxon_Join := gsub("X","", CdAppelTaxon_Join)][RsTaxRep != 0]
Supp_MacroinvertebratesAlr <- left_join(Supp_MacroinvertebratesAlr, Alr_Names_Macroinv, by = "CdAppelTaxon_Join")
Supp_MacroinvertebratesAlr[!CdAppelTaxon_Join %in% Alr_Traits_Macroinv$code_taxon_Alric2021, CdAppelTaxon_Join := NA]

setnames(Supp_MacroinvertebratesAlr, c("cd_site", "date_opecont", "year", "tax.GENUS", "FAMILLE",  "RsTaxRep"), 
         c("CdStation", "Date_PrelBio", "Year", "Genus", "Family", "Abundance"))
Supp_MacroinvertebratesAlr$NomLatinAppelTaxon <- Supp_MacroinvertebratesAlr$NomLatin_Join
Supp_MacroinvertebratesAlr <- Supp_MacroinvertebratesAlr[ ,.(CdStation, Date_PrelBio, Year, Genus, Family, NomLatinAppelTaxon, NomLatin_Join, CdAppelTaxon_Join, Abundance)]
Supp_MacroinvertebratesAlr[, CdAppelTaxon := CdAppelTaxon_Join][ , BDsource := "Alric"][, Group := "Macroinvertebrate"][, IDoperation := paste(Group, CdStation, Year, sep = "_")][,MAJ := "Ariane_March2023"]
Supp_MacroinvertebratesAlr[NomLatin_Join %in% Nai_Abd_Macroinv_maj[is.na(NomLatin_Join), NomLatinAppelTaxon], 
                           c("NomLatin_Join", "CdAppelTaxon_Join") := NA]

AllInv_Macroinvertebrate <- rbind(Nai_Abd_Macroinv_maj,Supp_MacroinvertebratesAlr)[ ,.(Group, MAJ, BDsource, IDoperation, CdStation, Date_PrelBio, Year, Genus, Family,
               NomLatinAppelTaxon, NomLatin_Join, CdAppelTaxon, CdAppelTaxon_Join, Abundance)][order(IDoperation)]

# Remove stations with single misidentified shellfish
SingleSp <- unique(AllInv_Macroinvertebrate[, Nsp := uniqueN(NomLatinAppelTaxon), by = c("CdStation","Year")][Nsp == 1, IDoperation])
AllInv_Macroinvertebrate <- AllInv_Macroinvertebrate[!IDoperation %in% SingleSp & !grepl("Procambarus|Orconectes|Pacifastacus", NomLatinAppelTaxon)]

AllInv_Macroinvertebrate[grep("/", Date_PrelBio), Samp_date := as.Date(Date_PrelBio, format = "%d/%m/%Y")]
AllInv_Macroinvertebrate[grep("-", Date_PrelBio), Samp_date := as.Date(Date_PrelBio, format = "%Y-%m-%d")]
AllInv_Macroinvertebrate[, Date_PrelBio := Samp_date][, Samp_date := NULL]

save(AllInv_Macroinvertebrate, file = "FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Inventories")
write.csv2(AllInv_Macroinvertebrate, row.names = F, file = "FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Inventories.csv")

Traits_Macroinvertebrate <- Alr_Traits_Macroinv
setnames(Traits_Macroinvertebrate, "code_taxon_Alric2021","CdAppelTaxon")

AllInv_Macroinvertebrate[!is.na(CdAppelTaxon), CdAppelTaxon := as.numeric(CdAppelTaxon)]

Traits_complement[is.na(Traits_complement)] <- 0
Traits_complement <- dcast(melt(Traits_complement, id.vars = "Nom_final"), variable ~ Nom_final)
setnames(Traits_complement, "variable","NomLatinAppelTaxon")
Traits_complement <- Traits_complement[,intersect(names(Traits_Macroinvertebrate), names(Traits_complement)), with = F][, 
                   setdiff(names(Traits_Macroinvertebrate), names(Traits_complement)) := NA]
Traits_Macroinvertebrate <- rbind(Traits_Macroinvertebrate, Traits_complement)

save(Traits_Macroinvertebrate, file = "FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Traits")
write.csv2(Traits_Macroinvertebrate, row.names = F, file = "FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Traits.csv")

#####


# Fish

#####
load("Naiades_Fish"); source("C:/Users/armirabel/Documents/INRAE/Hymobio/Hymobio_GitHub/BD_Bio/Alric_Fish.R")
Nai_Station <- as.data.table(read.csv("C:/Users/armirabel/Documents/DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/stations.csv", 
                             sep=";"))[,.(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, CodeRegion, LbRegion)]
Stations_metropole <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")), CdStationMesureEauxSurface]
Nai_Station <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")),]

Corr_Jerome <- data.table(read.table("Correspondance_Fish.csv", sep = ";", quote = "\"", header = T))[
  , .(NomLatinAppelTaxon_Naiades, NomLatin_Join, Propositions)]
setnames(Corr_Jerome, "NomLatinAppelTaxon_Naiades", "NomLatinAppelTaxon")

Nai_Abd_Fish_maj <- Nai_Abd_Fish[, NomLatin_Join := NomLatinAppelTaxon]
Nai_Abd_Fish_maj <- Nai_Abd_Fish_maj[!NomLatinAppelTaxon %in% Corr_Jerome[Propositions == "Discard", NomLatinAppelTaxon]]
Nai_Abd_Fish_maj[Corr_Jerome[Propositions == "OK", .(NomLatinAppelTaxon, NomLatin_Join)], 
                 on = "NomLatinAppelTaxon", NomLatin_Join := i.NomLatin_Join]
Nai_Abd_Fish_maj[grep("yprini", NomLatinAppelTaxon), NomLatinAppelTaxon := "Cyprinidae"]

Nai_Abd_Fish_maj[NomLatinAppelTaxon %in% Corr_Jerome[grep("Unmatched|Too large", Propositions),NomLatinAppelTaxon], 
                 NomLatin_Join := "Unmatched"]

load("dendritic_distance_table_biological_data.Rdata")
DistanceMatrix <- as.data.table(distance_table_biological_data)[COMPARTIMENT_START == "FISH" & COMPARTIMENT_TARGET == "FISH"]

for(species in c("Carassius", "Abramis")){
  
Cd_TAXON <- Alr_traits_Fish[tax.GENUS == species,sciname]

if(species == "Abramis"){Nai_torep <- Nai_Abd_Fish_maj[grep(species, NomLatin_Join)][NomLatin_Join != Cd_TAXON]
} else {Nai_torep <- Nai_Abd_Fish_maj[NomLatin_Join == species]}
       
Ref_Nai <- Nai_Abd_Fish_maj[NomLatinAppelTaxon %in% Cd_TAXON, ][,Abce_ref := sum(RsTaxRep), 
            by = c("CdStationMesureEauxSurface", "NomLatinAppelTaxon")]

idref <- unique(DistanceMatrix[ID_STATION_START %in% Nai_torep$CdStationMesureEauxSurface & ID_STATION_TARGET %in% Ref_Nai$CdStationMesureEauxSurface][, 
         closest := as.character(ID_STATION_TARGET[which.min(distance_m)]), by = "ID_STATION_START"][,.(ID_STATION_START, closest)])

Ref_closest_dendri <- right_join(Nai_torep, idref, by = c("CdStationMesureEauxSurface" = "ID_STATION_START"))

if(any(!Nai_torep$CdStationMesureEauxSurface %in% Ref_closest_dendri$CdStationMesureEauxSurface)){
  Nai_torep_XY <- Nai_Station[CdStationMesureEauxSurface %in% 
                                setdiff(Nai_torep$CdStationMesureEauxSurface,Ref_closest_dendri$CdStationMesureEauxSurface)]
  coordinates(Nai_torep_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
  
  Ref_Nai_XY <- Nai_Station[CdStationMesureEauxSurface %in% 
                              setdiff(Ref_Nai$CdStationMesureEauxSurface, Ref_closest_dendri$CdStationMesureEauxSurface)] 
  coordinates(Ref_Nai_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
  
  Ref_closest <- data.frame(closest = Ref_Nai_XY[apply(gDistance(Nai_torep_XY, Ref_Nai_XY, byid=TRUE), 2, which.min),
                                                 ]$CdStationMesureEauxSurface)
  Ref_closest <- cbind(data.table(CdStationMesureEauxSurface = Nai_torep_XY$CdStationMesureEauxSurface), Ref_closest)
  Ref_closest <- right_join(Nai_torep, Ref_closest, by = "CdStationMesureEauxSurface")
  
  Nai_torep <- rbindlist(list(Ref_closest, Ref_closest_dendri), use.names = T)
} else { Nai_torep <- Ref_closest_dendri; Nai_torep$closest <- as.character(Nai_torep$closest)}

Ref_Nai <- dcast(Ref_Nai, CdStationMesureEauxSurface ~ NomLatinAppelTaxon, value.var = "Abce_ref", fun.aggregate =  sum )
setnames(Ref_Nai, "CdStationMesureEauxSurface","closest")

Nai_torep <- left_join(Nai_torep, Ref_Nai, by = "closest")

Nai_torep[, c("NomLatin_Join", "Ninv") := list(names(.SD)[max.col(.SD, ties.method = "first")], 
               length(which(.SD!=0))), .SDcols = Cd_TAXON, by = c("year", "CdStationMesureEauxSurface")]

Nai_torep[Ninv > 1, NomLatin_Join := sample(names(.SD),1, prob = colSums(.SD)), .SDcols = Cd_TAXON, by = c("year", "CdStationMesureEauxSurface")]

Nai_Abd_Fish_maj[Nai_torep[,.(CdStationMesureEauxSurface, NomLatinAppelTaxon, NomLatin_Join)], 
                     on = .(CdStationMesureEauxSurface, NomLatinAppelTaxon), "NomLatin_Join" := i.NomLatin_Join]
}

Nai_Abd_Fish_maj[NomLatinAppelTaxon %in% setdiff(Nai_Abd_Fish_maj$NomLatin_Join, Alr_traits_Fish$sciname), NomLatin_Join := "Unmatched"]

Nai_Abd_Fish_maj[ , BDsource := "Naiades"][, Group := "Fish"][, IDoperation := paste(Group, CdStationMesureEauxSurface, year, sep = "_")][
  , MAJ := "Ariane_April2023"][, DateDebutOperationPrelBio := as.character(DateDebutOperationPrelBio)]

Nai_Abd_Fish_maj <- left_join(Nai_Abd_Fish_maj, Alr_traits_Fish, by = c("NomLatin_Join" = "sciname"))

setnames(Nai_Abd_Fish_maj, c("CdStationMesureEauxSurface", "DateDebutOperationPrelBio", "year", "tax.GENUS", "tax.FAMILY", "Code_TAXON", "RsTaxRep"),
         c("CdStation", "Date_PrelBio", "Year", "Genus", "Family", "CdAppelTaxon_Join","Abundance"))

Nai_Abd_Fish_maj[CdStation %in% Alr_Abd_Fish$CODE_STATION, BDsource := "Naiades+Alric"]

Supp_FishAlr <- Alr_Abd_Fish[CODE_STATION %in% setdiff(Alr_Abd_Fish$CODE_STATION, Nai_Abd_Fish_maj$CdStation)]
Supp_FishAlr <- melt(Supp_FishAlr, id.vars = c("CODE_STATION", "DATE", "year", "X_L93", "Y_L93"), 
        measure.vars = setdiff(colnames(Supp_FishAlr)[24:ncol(Supp_FishAlr)],"year"), variable.name = "CdAppelTaxon", 
        value.name = "RsTaxRep")[RsTaxRep != 0]
Supp_FishAlr[ , BDsource := "Alric"][, Group := "Fish"][, IDoperation := paste(Group, CODE_STATION, year, sep = "_")][, MAJ := "Ariane_March2023"]
Supp_FishAlr <- left_join(Supp_FishAlr, Alr_traits_Fish, by = c("CdAppelTaxon" = "Code_TAXON"))[, NomLatin_Join := sciname][, CdAppelTaxon_Join := CdAppelTaxon]
setnames(Supp_FishAlr, c("CODE_STATION", "DATE", "year", "tax.GENUS", "tax.FAMILY", "CdAppelTaxon", "sciname", "RsTaxRep"), 
    c("CdStation", "Date_PrelBio", "Year", "Genus", "Family", "CdAppelTaxon", "NomLatinAppelTaxon", "Abundance"))
Supp_FishAlr[!CdAppelTaxon %in% Alr_traits_Fish$Code_TAXON, "CdAppelTaxon_Join" := NA]

AllInv_Fish <- rbind(Nai_Abd_Fish_maj[ ,.(Group, MAJ, BDsource, IDoperation, CdStation, Date_PrelBio, Year, Genus, Family, 
                                          NomLatinAppelTaxon, CdAppelTaxon, NomLatin_Join, CdAppelTaxon_Join, Abundance)],
                 Supp_FishAlr[ ,.(Group, MAJ, BDsource, IDoperation, CdStation, Date_PrelBio, Year, Genus, Family, NomLatinAppelTaxon,
                                  CdAppelTaxon, NomLatin_Join, CdAppelTaxon_Join, Abundance)])[order(IDoperation)]

AllInv_Fish[grep("/", Date_PrelBio), Samp_date := as.Date(Date_PrelBio, format = "%d/%m/%Y")]
AllInv_Fish[grep("-", Date_PrelBio), Samp_date := as.Date(Date_PrelBio, format = "%Y-%m-%d")]
AllInv_Fish[, Date_PrelBio := Samp_date][, Samp_date := NULL]

save(AllInv_Fish, file = "FinalMAJ_Naiades.Alric.Traits/Fish_Inventories")
write.csv2(AllInv_Fish, row.names = F, file = "FinalMAJ_Naiades.Alric.Traits/Fish_Inventories.csv")

Traits_Fish <- Alr_traits_Fish
setnames(Traits_Fish, "Code_TAXON","CdAppelTaxon")
save(Traits_Fish, file = "FinalMAJ_Naiades.Alric.Traits/Fish_Traits")
write.csv2(Traits_Fish, row.names = F, file = "FinalMAJ_Naiades.Alric.Traits/Fish_Traits.csv")

#####

# Diatoms

#####
Nai_Station <- as.data.table(read.csv("C:/Users/armirabel/Documents/DB/1_MISE_A_JOUR_DONNEES_BIOLOGIQUE/HYDROBIOLOGIE NAIADE (ne pas modifier)/stations.csv", 
                                      sep=";"))[,.(CdStationMesureEauxSurface, CoordXStationMesureEauxSurface, CoordYStationMesureEauxSurface, CodeRegion, LbRegion)]
Stations_metropole <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")), CdStationMesureEauxSurface]
Nai_Station <- Nai_Station[!(LbRegion %in% c("", "La Réunion", "Martinique", "Guyane")),]

load("Naiades_Diatom"); source("C:/Users/armirabel/Documents/INRAE/Hymobio/Hymobio_GitHub/BD_Bio/Alric_Diatom.R")

Nai_Abd_Diatom_maj <- Nai_Abd_Diatom[, NomLatin_Simple := sub("^(\\S*\\s+\\S+).*", "\\1", NomLatinAppelTaxon)][, NomLatin_Join := NomLatin_Simple]
Nai_Abd_Diatom_maj[NomLatin_Join == "Pseudostaurosira (as", NomLatin_Join := "Fragilaria parasitica"]
Nai_Abd_Diatom_maj <- Nai_Abd_Diatom_maj[NomLatin_Join != "Groupe Diatomées"][!grep("Terpsi",NomLatin_Join)]

Nai_Abd_Diatom_maj[, Replacement_type := "None"]

Correspondances_Diatom <- fread("Correspondance_Diatom.csv")

Nai_Abd_Diatom_maj[Correspondances_Diatom, on = .(NomLatin_Simple), 
                   c("NomLatin_Join", "Replacement_type") := list(i.NomLatin_Join, i.Replacement_type)]
Nai_Abd_Diatom_maj[, Genus := sub(" .*", "", NomLatin_Simple)]
Nai_Abd_Diatom_maj[grep("Most common", NomLatin_Join), Torep := sub("^\\S+\\s+\\S+.", "", NomLatin_Join)][!is.na(Torep) & Torep != "", Genus := Torep][,Torep := NULL]

load("dendritic_distance_table_biological_data.Rdata")
DistanceMatrix <- as.data.table(distance_table_biological_data)[COMPARTIMENT_START == "DIATOM" & COMPARTIMENT_TARGET == "DIATOM"]

ToRep <- unique(Nai_Abd_Diatom_maj[grep("Most common", NomLatin_Join), .(NomLatin_Simple,Genus)])

Nomatch <- data.frame()
for(i.species in 1:nrow(ToRep)){
 
  species <- ToRep[i.species,]
  
  Cd_TAXON <- intersect(Nai_Abd_Diatom_maj[Genus == species$Genus & sub(" .*", "", NomLatin_Join) == species$Genus, NomLatin_Join],
                        Alr_Names_Diatom$SCIENTIFIC_NAME)
  
  if(length(Cd_TAXON) != 0) {
    
    Nai_torep <- Nai_Abd_Diatom_maj[NomLatin_Simple == species$NomLatin_Simple]
    Ref_Nai <- Nai_Abd_Diatom_maj[NomLatin_Simple %in% Cd_TAXON, ][, Abce_ref := sum(RsTaxRep), 
                                  by = c("CdStationMesureEauxSurface", "NomLatin_Simple")]
    
    idref <- unique(DistanceMatrix[ID_STATION_START %in% Nai_torep$CdStationMesureEauxSurface & ID_STATION_TARGET %in% Ref_Nai$CdStationMesureEauxSurface][, 
             closest := ID_STATION_TARGET[which.min(distance_m)], by = "ID_STATION_START"][,.(ID_STATION_START, closest)])
    setnames(idref, "ID_STATION_START", "CdStationMesureEauxSurface")
    
    Nai_torep[, closest := NA_character_][idref, on = "CdStationMesureEauxSurface", closest := i.closest]
   
  if(anyNA(Nai_torep$closest)){
    Nai_torep_XY <- Nai_Station[CdStationMesureEauxSurface %in% Nai_torep[is.na(closest), CdStationMesureEauxSurface]]
    coordinates(Nai_torep_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
    
    Ref_Nai_XY <- Nai_Station[CdStationMesureEauxSurface %in% Ref_Nai$CdStationMesureEauxSurface] 
    coordinates(Ref_Nai_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
    
    Ref_closest <- data.frame(closest = Ref_Nai_XY[apply(gDistance(Nai_torep_XY, Ref_Nai_XY, byid=TRUE), 2, which.min),
    ]$CdStationMesureEauxSurface)
    Ref_closest <- cbind(data.table(CdStationMesureEauxSurface = Nai_torep_XY$CdStationMesureEauxSurface), Ref_closest)
    
    Nai_torep[Ref_closest, on = "CdStationMesureEauxSurface", closest := i.closest]
    
  } 
    Ref_Nai <- dcast(Ref_Nai, CdStationMesureEauxSurface ~ NomLatin_Simple, value.var = "Abce_ref", fun.aggregate =  sum )
    setnames(Ref_Nai, "CdStationMesureEauxSurface","closest")
    
    Nai_torep <- left_join(Nai_torep, Ref_Nai, by = "closest")
    
    Nai_torep[, c("NomLatin_Join", "Ninv") := list(names(.SD)[max.col(.SD, ties.method = "first")], 
                length(which(.SD!=0))), .SDcols = Cd_TAXON, by = c("year", "CdStationMesureEauxSurface")]
    
    Nai_torep[Ninv > 1, NomLatin_Join := sample(names(.SD),1, prob = colSums(.SD)), .SDcols = Cd_TAXON, by = c("year", "CdStationMesureEauxSurface")]
     
    Nai_Abd_Diatom_maj[Nai_torep[,.(CdStationMesureEauxSurface, NomLatinAppelTaxon, NomLatin_Join)], 
                         on = .(CdStationMesureEauxSurface, NomLatinAppelTaxon), NomLatin_Join := i.NomLatin_Join]
  } else {Nomatch <- c(Nomatch, species$NomLatin_Simple)}

}

Nai_Abd_Diatom_maj[NomLatin_Simple %in% unlist(Nomatch), NomLatin_Join := "Unmatched"]
Nai_Abd_Diatom_maj[NomLatin_Join == "Unmatched", NomLatin_Simple := gsub(" diatom", "", NomLatin_Simple)]
Nai_Abd_Diatom_maj[, NomLatin_Taxo := NomLatin_Simple]

ToRep <- unique(Nai_Abd_Diatom_maj[NomLatin_Join == "Unmatched" & NomLatin_Simple == Genus, NomLatin_Taxo])

for(species in ToRep){
  Ref_Nai <- Nai_Abd_Diatom_maj[Genus == species & NomLatin_Simple != species][, Abce_ref := sum(RsTaxRep), 
            by = c("CdStationMesureEauxSurface", "NomLatin_Simple")][,.(CdStationMesureEauxSurface, NomLatin_Simple, Abce_ref)]
  
  if(nrow(Ref_Nai) != 0){
  
    Nai_torep <- Nai_Abd_Diatom_maj[NomLatin_Simple == species][
      CdStationMesureEauxSurface %in% Ref_Nai$CdStationMesureEauxSurface, closest := CdStationMesureEauxSurface]
    
    if(anyNA(Nai_torep$closest)){
      
    idref <- unique(DistanceMatrix[ID_STATION_START %in% Nai_torep[is.na(closest), CdStationMesureEauxSurface]][, 
                    closest := ID_STATION_TARGET[which.min(distance_m)], by = "ID_STATION_START"][,.(ID_STATION_START, closest)])
    setnames(idref, "ID_STATION_START", "CdStationMesureEauxSurface")
    
    Nai_torep[idref, on = "CdStationMesureEauxSurface", closest := i.closest]
    
    }
    
    if(anyNA(Nai_torep$closest)){
      Nai_torep_XY <- Nai_Station[!CdStationMesureEauxSurface %in% Nai_torep[is.na(closest), CdStationMesureEauxSurface]]
      coordinates(Nai_torep_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
      
      Ref_Nai_XY <- Nai_Station[CdStationMesureEauxSurface %in% Nai_torep[is.na(closest), CdStationMesureEauxSurface]]
      coordinates(Ref_Nai_XY) <- ~ CoordXStationMesureEauxSurface + CoordYStationMesureEauxSurface
      
      Ref_closest <- as.data.table(Nai_torep_XY[apply(gDistance(Ref_Nai_XY, Nai_torep_XY, byid=TRUE), 2, which.min),])
      setnames(Ref_closest, "CdStationMesureEauxSurface", "closest")
      Ref_closest$CdStationMesureEauxSurface <- Ref_Nai_XY$CdStationMesureEauxSurface
      
      Nai_torep[Ref_closest, on = "CdStationMesureEauxSurface", closest := i.closest]
      
    }
    
    Cd_TAXON <- unique(Ref_Nai$NomLatin_Simple)
    Ref_Nai <- dcast(Ref_Nai, CdStationMesureEauxSurface ~ NomLatin_Simple, value.var = "Abce_ref", fun.aggregate =  sum )
    setnames(Ref_Nai, "CdStationMesureEauxSurface", "closest")
    
    Nai_torep <- left_join(Nai_torep, Ref_Nai, by = "closest")
    
    Nai_torep[, c("NomLatin_Taxo", "Ninv") := list(names(.SD)[max.col(.SD, ties.method = "first")], 
                length(which(.SD!=0))), .SDcols = Cd_TAXON, by = RefOperationPrelBio]
    
    Nai_torep[Ninv > 1, NomLatin_Taxo := sample(names(.SD),1, prob = colSums(.SD)), .SDcols = Cd_TAXON, by = RefOperationPrelBio]
    
    Nai_Abd_Diatom_maj[Nai_torep[,.(CdStationMesureEauxSurface, NomLatinAppelTaxon, NomLatin_Taxo)], 
                       on = .(CdStationMesureEauxSurface, NomLatinAppelTaxon), NomLatin_Taxo := i.NomLatin_Taxo]
    
  }
  
}
       
Nai_Abd_Diatom_maj <- Nai_Abd_Diatom_maj[,.(CdStationMesureEauxSurface, DateDebutOperationPrelBio, year, Genus, CdAppelTaxon,
                                            NomLatinAppelTaxon, NomLatin_Simple, NomLatin_Join, NomLatin_Taxo, RsTaxRep)]

setnames(Nai_Abd_Diatom_maj, c("CdStationMesureEauxSurface", "DateDebutOperationPrelBio", "year", "Genus","RsTaxRep"),
         c("CdStation", "Date_PrelBio", "Year", "Genus", "Abundance"))
Nai_Abd_Diatom_maj[ , BDsource := "Naiades"][, Group := "Diatom"][, IDoperation := paste(Group, CdStation, Year, sep = "_")][,MAJ := "Ariane_April2023"]
Nai_Abd_Diatom_maj[, CdStation := as.character(CdStation)][, Date_PrelBio := as.character(Date_PrelBio)]

Nai_Abd_Diatom_maj[CdStation %in% Alr_Abd_Diatom$cd_site, BDsource := "Naiades+Alric"]

Supp_DiatomAlr <- Alr_Abd_Diatom[cd_site %in% setdiff(Alr_Abd_Diatom$cd_site, Nai_Abd_Diatom_maj$CdStation)]

Supp_DiatomAlr <- melt(Supp_DiatomAlr, id.vars = c("cd_site", "cd_opecont", "date_opecont", "year", "x", "y", "her1", "typo_nationale", "rang_corr", "agence"),
                       measure.vars = intersect(colnames(Supp_DiatomAlr),unique(Alr_Names_Diatom$CODE_TAXON)),
                  variable.name = "CdAppelTaxon", value.name = "RsTaxRep")[RsTaxRep != 0]
Supp_DiatomAlr <- left_join(Supp_DiatomAlr, Alr_Names_Diatom, by = c("CdAppelTaxon" = "CODE_TAXON"))
setnames(Supp_DiatomAlr, c("cd_site", "date_opecont", "year", "GENUS", "SCIENTIFIC_NAME", "RsTaxRep"), 
         c("CdStation", "Date_PrelBio", "Year", "Genus", "NomLatin_Join",  "Abundance"))
Supp_DiatomAlr[, c("NomLatinAppelTaxon", "NomLatin_Simple", "NomLatin_Taxo") := NomLatin_Join]
Supp_DiatomAlr <- Supp_DiatomAlr[ ,.(CdStation, Date_PrelBio, Year, Genus, CdAppelTaxon, NomLatinAppelTaxon, NomLatin_Simple, NomLatin_Join, NomLatin_Taxo, Abundance)]

Supp_DiatomAlr[ , BDsource := "Alric"][, Group := "Diatom"][, IDoperation := paste(Group, CdStation, Year, sep = "_")][,MAJ := "Ariane_July2023"]

AllInv_Diatom <- rbind(Nai_Abd_Diatom_maj,Supp_DiatomAlr)

AllInv_Diatom[grep("/", Date_PrelBio), SampDate := as.Date(Date_PrelBio, format = "%d/%m/%Y")]
AllInv_Diatom[grep("-", Date_PrelBio), SampDate := as.Date(Date_PrelBio, format = "%Y-%m-%d")]
AllInv_Diatom[, Date_PrelBio := SampDate][, SampDate := NULL]

AllInv_Diatom <- merge(AllInv_Diatom, unique(merge(Alr_Traits_Diatom[, "code"], 
                      Alr_Names_Diatom[, .(CODE_TAXON, SCIENTIFIC_NAME)], by.x = "code", by.y = "CODE_TAXON")[
                      , GenCode := code[1], by = SCIENTIFIC_NAME][, code := NULL]), 
                      by.x = "NomLatin_Join", by.y = "SCIENTIFIC_NAME", all.x = T)
setnames(AllInv_Diatom, "GenCode", "CdAppelTaxon_Join")
AllInv_Diatom <- AllInv_Diatom[ ,.(Group, MAJ, BDsource, IDoperation, CdStation, Date_PrelBio, Year, Genus, 
     NomLatinAppelTaxon, NomLatin_Simple, NomLatin_Join, NomLatin_Taxo, CdAppelTaxon, CdAppelTaxon_Join, Abundance)][order(IDoperation)]

save(AllInv_Diatom, file = "FinalMAJ_Naiades.Alric.Traits/Diatom_Inventories")
write.csv2(AllInv_Diatom, row.names = F, file = "FinalMAJ_Naiades.Alric.Traits/Diatom_Inventories.csv")

Traits_Diatom <- fread("FinalMAJ_Naiades.Alric.Traits/Diatom_Traits.csv")
traits <- setdiff(colnames(Traits_Diatom), c("code", "type", "taxo", "refuni", "fam", "genre", "terato", "tax_ibd"))
Traits_Diatom[, (traits) := lapply(.SD, function(x){as.numeric(gsub(",",".",x))}), .SDcols = traits]
setnames(Traits_Diatom, "code", "CdAppelTaxon")
  
save(Traits_Diatom, file = "FinalMAJ_Naiades.Alric.Traits/Diatom_Traits")
write.csv2(Traits_Diatom, row.names = F, file = "FinalMAJ_Naiades.Alric.Traits/Diatom_Traits.csv")
#####


# Overview 
load("FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Traits")
load("FinalMAJ_Naiades.Alric.Traits/Macroinvertebrate_Inventories")

load("FinalMAJ_Naiades.Alric.Traits/Diatom_Traits")
load("FinalMAJ_Naiades.Alric.Traits/Diatom_Inventories")
Correspondances_Diatom <- fread("Correspondance_Diatom.csv")

Correspondances <-  AllInv_Diatom[,.(BDsource, NomLatin_Simple, NomLatin_Join, Abundance)][, Replacement_type := "None"]
Correspondances[Correspondances_Diatom[,.(NomLatin_Simple, Replacement_type)], 
                on = "NomLatin_Simple", Replacement_type := i.Replacement_type]

load("FinalMAJ_Naiades.Alric.Traits/Fish_Traits")
load("FinalMAJ_Naiades.Alric.Traits/Fish_Inventories")
Corr_Jerome <- data.table(read.table("Correspondance_Fish.csv", sep = ";", quote = "\"", header = T))[
  , .(NomLatinAppelTaxon_Naiades, NomLatin_Join, Propositions)]
setnames(Corr_Jerome, "NomLatinAppelTaxon_Naiades", "NomLatinAppelTaxon")


#  Macroinv
#####
SumUp <- as.data.frame(t(data.frame("Discarded" = sum(AllInv_Macroinvertebrate[NomLatin_Join == "Unmatched", Abundance]),
        "No_Traits" = sum(AllInv_Macroinvertebrate[NomLatin_Join != "Unmatched" & CdAppelTaxon_Join == "Unmatched", Abundance]),
        "Substitute_Traits" = sum(AllInv_Macroinvertebrate[NomLatinAppelTaxon != NomLatin_Join & NomLatin_Join != "Unmatched", Abundance]),
        "Matching" = sum(AllInv_Macroinvertebrate[NomLatinAppelTaxon == NomLatin_Join, Abundance]),
        
        "Naiades" = sum(AllInv_Macroinvertebrate[BDsource == "Naiades", Abundance]),
        "Alric" = sum(AllInv_Macroinvertebrate[BDsource == "Alric", Abundance]),
        "Both" = sum(AllInv_Macroinvertebrate[BDsource == "Naiades+Alric", Abundance]))))
colnames(SumUp) <- "Eff"
SumUp$groups <- factor(rownames(SumUp), levels = c("Discarded", "Matching", "No_Traits","Substitute_Traits", "Both", "Naiades", "Alric"))

p1m <-ggplot(SumUp[1:4,],aes(x = "", y = Eff, fill = groups)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_viridis_d() + theme_void() + 
  guides(fill = guide_legend(title="")) + ggtitle("Macroinvertebrates")

p2m <- ggplot(SumUp[5:7,],aes(x = "", y = Eff, fill = groups)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_viridis_d() + theme_void() + 
  guides(fill = guide_legend(title="")) 
#####
# Diatom
#####

SumUp <- as.data.frame(t(data.frame("Substitute_Traits" = sum(AllInv_Diatom[NomLatin_Simple %in% Correspondances_Diatom[grep("Holotype|common|Basionym", Replacement_type)
                                                                          , NomLatin_Simple] & NomLatin_Join != "Unmatched", Abundance]),
                                    "No_Traits" = sum(AllInv_Diatom[!NomLatin_Simple %in% Correspondances_Diatom[Replacement_type == "Too large", NomLatin_Simple] 
                                                      & NomLatin_Join == "Unmatched", Abundance]),
                                    "Matching" = sum(Correspondances[Replacement_type == "None" | NomLatin_Simple %in% Correspondances_Diatom[grep("Synonym|Correction", Replacement_type), NomLatin_Simple],
                                                                     Abundance]),
                                    "Discarded" = sum(AllInv_Diatom[NomLatin_Simple %in% Correspondances_Diatom[Replacement_type == "Too large", NomLatin_Simple], Abundance]),
                                     
                                    "Naiades" = sum(AllInv_Diatom[BDsource == "Naiades", Abundance]),
                                    "Alric" = sum(AllInv_Diatom[BDsource == "Alric", Abundance]),
                                    "Both" = sum(AllInv_Diatom[BDsource == "Naiades+Alric", Abundance]))))
colnames(SumUp) <- "Eff"
SumUp$groups <-  factor(rownames(SumUp), levels = c("Discarded", "Matching", "Substitute_Traits", "No_Traits","Both", "Naiades", "Alric"))

p1d <-ggplot(SumUp[1:4,],aes(x = "", y = Eff, fill = groups)) + 
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_viridis_d() + theme_void() + 
  guides(fill = guide_legend(title="")) + ggtitle("Diatoms")

p2d <- ggplot(SumUp[5:7,],aes(x = "", y = Eff, fill = groups)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_viridis_d() + theme_void() + 
  guides(fill = guide_legend(title="")) 
#####
# Fish
######

SumUp <- as.data.frame(t(data.frame("Substitute_Traits" = sum(AllInv_Fish[NomLatinAppelTaxon %in%
                                                    Corr_Jerome[Propositions == "OK", NomLatinAppelTaxon], Abundance]),
                                    "No_Traits" = sum(AllInv_Fish[NomLatinAppelTaxon %in%
                                                      Corr_Jerome[Propositions == "Unmatched", NomLatinAppelTaxon], Abundance]),
                                    "Discarded" = sum(AllInv_Fish[NomLatinAppelTaxon %in%
                                                      Corr_Jerome[Propositions == "Too large", NomLatinAppelTaxon], Abundance]),
                                    "Matching" = sum(AllInv_Fish[!NomLatinAppelTaxon %in% Corr_Jerome$NomLatinAppelTaxon, Abundance]),
                                    
                                    "Naiades" = sum(AllInv_Fish[BDsource == "Naiades", Abundance]),
                                    "Alric" = sum(AllInv_Fish[BDsource == "Alric", Abundance]),
                                    "Both" = sum(AllInv_Fish[BDsource == "Naiades+Alric", Abundance]))))
colnames(SumUp) <- "Eff"
SumUp$groups <-  factor(rownames(SumUp), levels = c("Discarded", "Matching", "Substitute_Traits", "No_Traits","Both", "Naiades", "Alric"))

p1f <-ggplot(SumUp[1:4,],aes(x = "", y = Eff, fill = groups)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_viridis_d() + theme_void() + 
  guides(fill = guide_legend(title="")) + ggtitle("Fishes")

p2f <- ggplot(SumUp[5:7,],aes(x = "", y = Eff, fill = groups)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_viridis_d() + theme_void() + 
  guides(fill = guide_legend(title="")) 
#####

ggarrange(p1d, p2d, p1m, p2m, p1f, p2f, ncol= 2, nrow = 3, widths = c(4,3))


# Single species inventories
SingleCheck <- rbind(unique(AllInv_Diatom[, Nsp := .N, by = c("CdStation", "Year", "Group")][
  , .(Group, CdStation, Year, Nsp)]),
  unique(AllInv_Macroinvertebrate[, Nsp := .N, by = c("CdStation", "Year", "Group")][, .(Group, CdStation, Year, Nsp)]),
  unique(AllInv_Fish[, Nsp := .N, by = c("CdStation", "Year", "Group")][, .(Group, CdStation, Year, Nsp)])
  )
Effectifs <- unique(SingleCheck[,Ngrp := .N, by = c("Nsp", "Group")][,.(Ngrp,Nsp, Group)])[,Pgrp := round(Ngrp/sum(Ngrp),3), by = Group]
Effectifs[Pgrp > 0.01]

AllInv_Fish[CdStation %in% SingleCheck[Nsp == 1 & Group == "Fish", CdStation]][Nsp == 1,.(Group, BDsource, CdStation, Year, NomLatinAppelTaxon, Abundance)]



SingleCheck <- rbind(unique(AllInv_Diatom[, Nsp := uniqueN(NomLatinAppelTaxon), by = c("CdStation", "Group")][
  , .(Group, CdStation, Year, Nsp)]),
  unique(AllInv_Macroinvertebrate[, Nsp := uniqueN(NomLatinAppelTaxon), by = c("CdStation", "Group")][, .(Group, CdStation, Year, Nsp)]),
  unique(AllInv_Fish[, Nsp := uniqueN(NomLatinAppelTaxon), by = c("CdStation", "Group")][, .(Group, CdStation, Year, Nsp)])
)
Effectifs <- unique(SingleCheck[, Ngrp := .N, by = c("Nsp", "Group")][,.(Ngrp, Nsp, Group)])[,Pgrp := round(Ngrp/sum(Ngrp),3), by = Group]
Effectifs[Pgrp > 0.01]
Effectifs[Nsp ==1]

AllInv_Fish[CdStation %in% SingleCheck[Nsp == 1 & Group == "Fish", CdStation]][,.(Group, BDsource, CdStation, Year, NomLatinAppelTaxon, Abundance)]

SingleSpeciesStations <- SingleCheck[Nsp == 1 & Group == "Fish", CdStation]
save(SingleSpeciesStations, file = "../BD_Hydro/CheckSingleSpeciesStations")




