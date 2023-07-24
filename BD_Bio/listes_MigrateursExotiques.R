
load("Fishes_MAJ")
listMig <- fread("listes_MA_EE/poissons_amphihalins.csv")
Inv_Migrateurs <- AllFish[NomLatin %in% listMig$scientificName]
listMig[listMig$scientificName %in% AllFish$NomLatin, match := "Exact"]

Inv_Migrateurs_unmatched <- listMig[scientificName %in% setdiff(listMig[,scientificName],AllFish[,NomLatin])][, Genus := sub(" .*", "", scientificName)]
Inv_Migrateurs <- rbind(Inv_Migrateurs, AllFish[grep(paste0(Inv_Migrateurs_unmatched$Genus, collapse = "|"), NomLatin)])
listMig[scientificName %in% grep(paste0(sub(" .*", "", unique(grep(paste0(Inv_Migrateurs_unmatched$Genus, collapse = "|"),
                                                                   AllFish$NomLatin, value = T))), collapse = "|"), listMig$scientificName, value = TRUE) & is.na(match), match := "Genus"]
Inv_Migrateurs <- rbind(Inv_Migrateurs, AllFish[Family %in% str_to_title(listMig[is.na(match),familyName]),])
listMig[familyName %in% intersect(AllFish$Family, str_to_title(listMig[is.na(match),familyName])), match := "Family"]
listMig[is.na(match), match := "Unmatched"]
write.table(listMig, file = "listes_MA_EE/matching_migrateurs.csv", sep = ";", row.names = F)


listExotiq <- fread("listes_MA_EE/poissons_exotiques_envahissants.csv")
Inv_Exotiq <- AllFish[NomLatin %in% listExotiq$scientificName]
listExotiq[listExotiq$scientificName %in% AllFish$NomLatin, match := "Exact"]
Inv_Exotiq_unmatched <- listExotiq[scientificName %in% setdiff(listExotiq[,scientificName],AllFish[,NomLatin])][, Genus := sub(" .*", "", scientificName)]
Inv_Exotiq <- rbind(Inv_Exotiq, AllFish[grep(paste0(Inv_Exotiq_unmatched$Genus, collapse = "|"), NomLatin)])
listExotiq[scientificName %in% grep(paste0(sub(" .*", "", unique(grep(paste0(Inv_Exotiq_unmatched$Genus, collapse = "|"),
                                                                      AllFish$NomLatin, value = T))), collapse = "|"), listExotiq$scientificName, value = TRUE) & is.na(match), match := "Genus"]
listExotiq[is.na(match), match := "Unmatched"]
write.table(listExotiq, file = "listes_MA_EE/matching_exotics.csv", sep = ";", row.names = F)
