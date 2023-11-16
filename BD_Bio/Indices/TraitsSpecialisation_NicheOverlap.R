# Food specialisation
tsi <-  filter(traitm1, 
               TRAIT %in% c("food")) %>%
  left_join(group_by(traitm1, TRAIT) %>%
              summarise(k = n_distinct(CATEGORIE)),
            by = "TRAIT")                               %>%
  group_by(CODE_TAXON, TRAIT)                           %>%
  summarise(TSI = (sum(VALEUR^2) - 1 / unique(k)) /
              (1 - 1 / unique(k)))

specialisation <- data_entree1%>%
  select(CODE_OPERATION,CODE_TAXON,RESULTAT)                      %>%
  left_join(tsi,by = c("CODE_TAXON"))     %>%
  filter(!is.na(TRAIT))                                          %>%
  left_join(abondance, by = "CODE_OPERATION")                    %>%
  group_by(CODE_OPERATION, TRAIT)                                %>%
  summarise(CSI = sum(TSI * log(RESULTAT + 1)) / 
              unique(logAbTot))                                %>%
  mutate(TRAIT = paste("CSI", TRAIT, sep = "_"))                 %>%
  spread(key = "TRAIT", value = "CSI")


# Niches overlap
calc_piankaDist <- function(m) {
  designdist(m, method = "J/sqrt(A*B)", terms = "quadratic")
} 

distances <- filter(traitm1,
                    TRAIT %in% c("food"),
                    CODE_TAXON %in% unique(data_entree1$CODE_TAXON))%>%
  split(.$TRAIT)                                             %>%
  lapply(select, -TRAIT)                                     %>%
  lapply(spread, key = "CATEGORIE", 
         value = "VALEUR", fill = 0)                         %>%
  lapply(function(df) {
    taxa <- df$CODE_TAXON
    
    df <- as.matrix(df[, -1]) %>%
      apply(MARGIN = 2, as.numeric)
    rownames(df) <- taxa
    
    return(df)
  })                                                         %>%
  lapply(calc_piankaDist)                                    %>%
  lapply(as.matrix)

calc_recouvrement <- function(taxa, dist) {
  tx <- unique(taxa) %>%
    as.character() %>%
    (function(x) {x[x %in% colnames(dist)]})
  
  mean(as.dist(dist[tx, tx]))
}

recouvrement <- lapply(names(distances),
                       function(i) {
                         d <- distances[[i]]
                         
                         data_entree1%>%
                           select(CODE_OPERATION,CODE_TAXON,RESULTAT)   %>%
                           group_by(CODE_OPERATION)                %>%
                           summarise(pianka = 
                                       calc_recouvrement(taxa = CODE_TAXON,
                                                         dist = d)) %>%
                           mutate(TRAIT = i)
                       })                              %>%
  bind_rows()                                        %>%
  mutate(TRAIT = paste("SIMILARITE", TRAIT, sep = "_")) %>%
  spread(key = "TRAIT", value = "pianka")
