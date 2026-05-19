
#====================================================================================
# ANALYSE DES INCOHERENCES GPS — EACI MALI 2017
# Script complet et definitif
#====================================================================================
# Objectif 1 : Validation géographique systématique
# Objectif 2 : Taxonomie des erreurs GPS
# Objectif 3 : Détection des doublons spatiaux et distinction des origines
# Objectif 4 : Quantifier surface déclarée vs emprise GPS
# Objectif 5 : Mesurer l'impact des corrections sur les estimations finales
#====================================================================================
# 0. INITIALISATION
#====================================================================================
rm(list = ls())

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(dbscan)
library(units)

#====================================================================================
# 1. IMPORTATION DES BASES DE DONNEES
#====================================================================================
chemin <- "C:/Users/bmd/Desktop/Expose_R_ZONGO"

# Fichier principal : coordonnées GPS + variables géographiques
donnees <- read.csv(file.path(chemin, "eaci_geovariables_2017.csv"),
                    stringsAsFactors = FALSE)
View(donnees)
# Fichier superficie déclarée (s7cq07 = superficie parcelle en hectares)
parcelles <- read.csv(file.path(chemin, "eaci17_s01p1.csv"),
                      stringsAsFactors = FALSE)

# Fichier production réelle (s7fq13d = production campagne 2017/18 en kg)
production <- read.csv(file.path(chemin, "eaci17_s7fp2.csv"),
                       stringsAsFactors = FALSE)

cat("=== CHARGEMENT DES DONNEES ===\n")
cat("Observations GPS         :", nrow(donnees), "\n")
cat("Observations parcelles   :", nrow(parcelles), "\n")
cat("Observations production  :", nrow(production), "\n")

#====================================================================================
# 2. NETTOYAGE ET PREPARATION DES COORDONNEES GPS
#====================================================================================
donnees$lat_dd_mod <- as.numeric(donnees$lat_dd_mod)
donnees$lon_dd_mod <- as.numeric(donnees$lon_dd_mod)

donnees <- donnees %>%
  filter(!is.na(lat_dd_mod), !is.na(lon_dd_mod),
         is.finite(lat_dd_mod), is.finite(lon_dd_mod))

cat("Observations après nettoyage initial :", nrow(donnees), "\n")

#====================================================================================
# 3. CONVERSION SPATIALE — TROIS PROJECTIONS PREPAREES
# EPSG:4326 = WGS84 (degrés décimaux)
# EPSG:3857 = Web Mercator (mètres) — utilisé pour distances
# EPSG:32629 = UTM Zone 29N (mètres) — projection locale Mali
#====================================================================================

# Projection WGS84 de base
points_wgs84 <- st_as_sf(donnees,
                         coords = c("lon_dd_mod", "lat_dd_mod"),
                         crs = 4326)

# Projection Web Mercator (mètres) — utilisée pour DBSCAN et distances
points_sf <- points_wgs84 %>% st_transform(3857)

# Projection UTM Zone 29N — projection locale adaptée au Mali
# Utilisée pour le piège 1 (vérification projection)
points_utm <- points_wgs84 %>% st_transform(32629)

# Contour du Mali dans les 3 projections
mali_wgs84 <- ne_states(country = "Mali", returnclass = "sf") %>%
  st_transform(4326)
mali       <- mali_wgs84 %>% st_transform(3857)
mali_utm   <- mali_wgs84 %>% st_transform(32629)

#====================================================================================
# ETAPE 1 — VALIDATION GEOGRAPHIQUE SYSTEMATIQUE
#====================================================================================

cat("\n=== ETAPE 1 — VALIDATION GEOGRAPHIQUE ===\n")

# --- 1A. Points hors frontières en EPSG:3857 ---
dans_mali_3857 <- lengths(st_intersects(points_sf, mali)) > 0
donnees$hors_mali <- ifelse(dans_mali_3857, 0, 1)
cat("Points hors Mali (EPSG:3857) :", sum(donnees$hors_mali), "\n")

# ---------------------------------------------------------------
# PIEGE 1 : Un point hors frontière peut être une projection
# mal définie, pas une vraie erreur GPS
# Solution : tester en UTM Zone 29N (projection locale Mali)
# Si le point rentre dans le Mali en UTM → erreur de projection
# ---------------------------------------------------------------
dans_mali_utm <- lengths(st_intersects(points_utm, mali_utm)) > 0
donnees$hors_mali_utm <- ifelse(dans_mali_utm, 0, 1)

# Erreur de projection : hors Mali en 3857 MAIS dans Mali en UTM
donnees$erreur_projection <- ifelse(
  donnees$hors_mali == 1 & donnees$hors_mali_utm == 0, 1, 0
)

# Hors Mali réel : hors Mali dans les DEUX projections
donnees$hors_mali_reel <- ifelse(
  donnees$hors_mali == 1 & donnees$hors_mali_utm == 1, 1, 0
)

cat("Erreurs de projection (rentrent en UTM) :",
    sum(donnees$erreur_projection), "\n")
cat("Hors Mali réel (hors dans les 2 proj.)  :",
    sum(donnees$hors_mali_reel), "\n")

# --- 1B. Distance à la frontière (proxy zones aquatiques) ---
frontiere <- st_union(st_boundary(mali))
distance_min <- st_distance(points_sf, frontiere)
donnees$proche_frontiere <- as.numeric(distance_min)
donnees$dans_eau_suspect <- ifelse(donnees$proche_frontiere > 50000, 1, 0)
cat("Zones aquatiques (distance > 50km)      :",
    sum(donnees$dans_eau_suspect), "\n")

# --- 1C. Zones humides via TWI ---
seuil_twi <- quantile(donnees$twi, 0.90, na.rm = TRUE)
donnees$zone_humide_twi <- ifelse(donnees$twi > seuil_twi, 1, 0)
cat("Zones humides (TWI > Q90)               :",
    sum(donnees$zone_humide_twi, na.rm = TRUE), "\n")

# --- 1D. Zones non cultivées via NDVI ---
seuil_ndvi_bas <- quantile(donnees$ndvi_avg, 0.05, na.rm = TRUE)
donnees$zone_non_cultivee <- ifelse(donnees$ndvi_avg < seuil_ndvi_bas, 1, 0)
cat("Zones non cultivées (NDVI < Q5)         :",
    sum(donnees$zone_non_cultivee, na.rm = TRUE), "\n")

# --- 1E. Zones inhabitées via densité de population ---
seuil_pop <- quantile(donnees$popdensity, 0.05, na.rm = TRUE)
donnees$zone_inhabitee <- ifelse(donnees$popdensity < seuil_pop, 1, 0)
cat("Zones inhabitées (popdensity < Q5)      :",
    sum(donnees$zone_inhabitee, na.rm = TRUE), "\n")

#====================================================================================
# ETAPE 2 — TAXONOMIE DES ERREURS GPS
# Identification des causes probables
#====================================================================================

cat("\n=== ETAPE 2 — TAXONOMIE DES ERREURS GPS ===\n")

# --- 2A. Erreur de saisie manuelle ---
# Point hors frontière réel + NDVI normal → probable inversion lat/lon
donnees$erreur_saisie <- ifelse(
  donnees$hors_mali_reel == 1 &
    donnees$ndvi_avg > seuil_ndvi_bas,
  1, 0
)

# --- 2B. Erreur de projection ---
# Déjà calculée à l'étape 1 : donnees$erreur_projection

# --- 2C. Dérive GPS ---
# Dans le Mali + zone non cultivée + loin des routes
seuil_route <- quantile(donnees$dist_road, 0.90, na.rm = TRUE)
donnees$derive_gps <- ifelse(
  donnees$hors_mali_reel == 0 &
    donnees$zone_non_cultivee == 1 &
    donnees$dist_road > seuil_route,
  1, 0
)

# --- 2D. Résumé taxonomie ---
cat("Erreurs de projection                   :",
    sum(donnees$erreur_projection, na.rm = TRUE), "\n")
cat("Erreurs de saisie (inversion lat/lon)   :",
    sum(donnees$erreur_saisie, na.rm = TRUE), "\n")
cat("Dérives GPS                             :",
    sum(donnees$derive_gps, na.rm = TRUE), "\n")
cat("(Fraudes enquêteur calculées étape 3)\n")

#====================================================================================
# ETAPE 3 — DETECTION DES DOUBLONS SPATIAUX
#====================================================================================

cat("\n=== ETAPE 3 — DOUBLONS SPATIAUX ===\n")

# Charger le fichier s01
s01 <- read.csv(file.path(chemin, "eaci17_s01p1.csv"), 
                stringsAsFactors = FALSE)
# Étape 1 : Réduire s01 à une ligne par grappe
s01_reduit <- s01 %>%
  select(grappe, exploitation, passage) %>%
  group_by(grappe) %>%
  slice(1) %>%
  ungroup()

cat("Nb grappes uniques s01 :", nrow(s01_reduit), "\n")

# Étape 2 : Jointure avec donnees
donnees <- donnees %>%
  left_join(s01_reduit, by = "grappe")

cat("Nb lignes après jointure :", nrow(donnees), "\n")
# Doit afficher 953

# Étape 3 : Vérifier que passage existe bien
cat("passage" %in% names(donnees))
# Doit afficher TRUE

# Étape 2 : Relancer DBSCAN
donnees$cluster              <- NULL
donnees$doublon              <- NULL
donnees$origine_doublon      <- NULL
donnees$nb_grappes_cluster   <- NULL
donnees$nb_passages_cluster  <- NULL
donnees$nb_exploitations_cluster <- NULL

# Étape 2 : DBSCAN
coords <- st_coordinates(points_sf)
db <- dbscan(coords, eps = 500, minPts = 2)
donnees$cluster <- db$cluster
donnees$doublon <- ifelse(donnees$cluster > 0, 1, 0)

# Étape 3 : Distinction origines doublons
donnees <- donnees %>%
  group_by(cluster) %>%
  mutate(
    nb_grappes_cluster       = n_distinct(grappe),
    nb_passages_cluster      = n_distinct(passage),
    nb_exploitations_cluster = n_distinct(exploitation),
    origine_doublon = case_when(
      cluster == 0
      ~ "Non_doublon",
      nb_passages_cluster > 1
      ~ "Meme_parcelle_deux_passages_legitime",
      nb_grappes_cluster == 1 &
        nb_passages_cluster == 1 &
        nb_exploitations_cluster == 1
      ~ "Meme_saisie_fraude_probable",
      nb_grappes_cluster > 1
      ~ "Meme_zone_jittering_legitime",
      TRUE ~ "Indetermine"
    )
  ) %>%
  ungroup()

# Étape 4 : Résultat
cat("Résumé origines doublons :\n")
print(table(donnees$origine_doublon))
#====================================================================================
# ETAPE 4 — OUTLIERS SPATIAUX (DISTANCE AU CENTROIDE)
#====================================================================================

centre <- st_centroid(st_union(points_sf))
distances <- st_distance(points_sf, centre)
donnees$distance_centre <- as.numeric(distances)
seuil_dist <- quantile(donnees$distance_centre, 0.95, na.rm = TRUE)
donnees$suspect_distance <- ifelse(donnees$distance_centre > seuil_dist, 1, 0)

cat("Outliers spatiaux (Q95%) :", sum(donnees$suspect_distance), "\n")
#====================================================================================
# ETAPE 5 — SCORE GLOBAL D'ANOMALIE
#====================================================================================

donnees$score <- donnees$hors_mali_reel +
  donnees$dans_eau_suspect +
  donnees$zone_humide_twi +
  donnees$doublon +
  donnees$suspect_distance

donnees <- donnees %>%
  mutate(
    niveau = case_when(
      score == 0   ~ "Bon",
      score <= 2   ~ "A_verifier",
      TRUE         ~ "Suspect"
    )
  )

cat("\n=== SCORE GLOBAL ===\n")
print(table(donnees$niveau))
#====================================================================================
# ETAPE 5 — REGLES DE DECISION JUSTIFIEES
#====================================================================================

cat("\n=== REGLES DE DECISION ===\n")

donnees <- donnees %>%
  mutate(
    decision = case_when(
      
      # EXCLURE : hors territoire réel (vérifié dans 2 projections)
      # Justification : ne peut pas représenter une exploitation malienne
      hors_mali_reel == 1 & erreur_projection == 0
      ~ "Exclure_hors_territoire",
      
      # CORRIGER PROJECTION : hors 3857 mais dans UTM = erreur de projection
      # Justification : erreur technique, pas une vraie anomalie GPS
      erreur_projection == 1
      ~ "Corriger_projection",
      
      # EXCLURE : fraude probable AVEC score élevé
      # Justification : même saisie + anomalies multiples = non fiable
      origine_doublon == "Meme_saisie_fraude_probable" & score >= 2
      ~ "Exclure_fraude",
      
      # GARDER : même parcelle visitée 2 passages = légitime
      # Justification : PIEGE 2 — c'est une vraie double visite
      origine_doublon == "Meme_parcelle_deux_passages_legitime"
      ~ "Garder_double_visite_legitime",
      
      # FUSIONNER : doublons de zone jittering = garder un seul
      # Justification : coordonnées perturbées LSMS, pas une erreur
      origine_doublon == "Meme_zone_jittering_legitime"
      ~ "Fusionner_jittering",
      
      # VERIFIER MANUELLEMENT : zones aquatiques ou humides
      # Justification : besoin de confirmation terrain avant exclusion
      dans_eau_suspect == 1 | zone_humide_twi == 1
      ~ "Verifier_manuel",
      
      # IMPUTER : outlier isolé dans le Mali avec score faible
      # Justification : anomalie légère, on préfère imputer qu'exclure
      suspect_distance == 1 & hors_mali_reel == 0 & score == 1
      ~ "Imputer_coordonnees",
      
      # GARDER : aucune anomalie détectée
      TRUE ~ "Garder"
    )
  )

print(table(donnees$decision))

#====================================================================================
# OBJECTIF 4 — SURFACE CULTIVEE ESTIMEE VIA hybrid_V8
#====================================================================================

points_buffer_geom <- lapply(
  st_geometry(points_sf),
  function(geom) st_buffer(geom, dist = 100)
)
points_buffer_geom <- st_sfc(points_buffer_geom, crs = 3857)
donnees$surface_gps_m2 <- as.numeric(st_area(points_buffer_geom))
donnees$surface_gps_ha  <- donnees$surface_gps_m2 / 10000

# Corriger hybrid_V8 en proportion (diviser par 100)
donnees$hybrid_V8_prop <- donnees$hybrid_V8 / 100

# Recalculer la surface cultivée
donnees$surface_cultivee_ha <- donnees$surface_gps_ha * donnees$hybrid_V8_prop
cat("\n=== SURFACE CULTIVEE ESTIMEE (proxy hybrid_V8) ===\n")
cat("Surface GPS moyenne (ha)   :",
    round(mean(donnees$surface_gps_ha, na.rm=TRUE), 3), "\n")
cat("Surface cultivée moy. (ha) :",
    round(mean(donnees$surface_cultivee_ha, na.rm=TRUE), 3), "\n")
cat("Écart moyen (ha)           :",
    round(mean(donnees$surface_gps_ha - donnees$surface_cultivee_ha,
               na.rm=TRUE), 3), "\n")
cat("hybrid_V8 moyen            :",
    round(mean(donnees$hybrid_V8, na.rm=TRUE), 3), "\n")
library(ggplot2)
ggplot(donnees, aes(x = surface_cultivee_ha)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(donnees$surface_cultivee_ha, na.rm=TRUE),
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution des surfaces cultivees estimees",
       subtitle = "Estimee via hybrid_V8 x emprise GPS (buffer 100m)",
       x = "Surface cultivee (ha)", 
       y = "Frequence") + 
  theme_minimal()
#====================================================================================
# OBJECTIF 5 — IMPACT DES CORRECTIONS SUR LES ESTIMATIONS FINALES
#====================================================================================

cat("\n=== OBJECTIF 5 — IMPACT DES CORRECTIONS ===\n")

#====================================================================================
# OBJECTIF 5 — IMPACT DES CORRECTIONS SUR LES ESTIMATIONS FINALES
#====================================================================================

# Base avant correction
base_avant <- donnees

# Base après correction
base_apres <- donnees %>%
  filter(!decision %in% c("Exclure_hors_territoire", "Exclure_fraude"))

# Calcul des estimations
estim_avant <- data.frame(
  Etape              = "Avant correction",
  N                  = nrow(base_avant),
  Surface_moy_ha     = round(mean(base_avant$surface_cultivee_ha, na.rm=TRUE), 3),
  NDVI_moyen         = round(mean(base_avant$ndvi_avg, na.rm=TRUE), 4),
  hybrid_V8_moyen    = round(mean(base_avant$hybrid_V8, na.rm=TRUE), 2),
  Dist_moy_km        = round(mean(base_avant$distance_centre, na.rm=TRUE)/1000, 2)
)

estim_apres <- data.frame(
  Etape              = "Après correction",
  N                  = nrow(base_apres),
  Surface_moy_ha     = round(mean(base_apres$surface_cultivee_ha, na.rm=TRUE), 3),
  NDVI_moyen         = round(mean(base_apres$ndvi_avg, na.rm=TRUE), 4),
  hybrid_V8_moyen    = round(mean(base_apres$hybrid_V8, na.rm=TRUE), 2),
  Dist_moy_km        = round(mean(base_apres$distance_centre, na.rm=TRUE)/1000, 2)
)

comparaison <- rbind(estim_avant, estim_apres)
cat("\n=== IMPACT DES CORRECTIONS SUR LES ESTIMATIONS ===\n")
print(comparaison)

cat("\n--- Variations relatives ---\n")
cat("Surface cultivée :",
    round((estim_apres$Surface_moy_ha - estim_avant$Surface_moy_ha) /
            estim_avant$Surface_moy_ha * 100, 2), "%\n")
cat("NDVI moyen       :",
    round((estim_apres$NDVI_moyen - estim_avant$NDVI_moyen) /
            estim_avant$NDVI_moyen * 100, 2), "%\n")
cat("hybrid_V8 moyen  :",
    round((estim_apres$hybrid_V8_moyen - estim_avant$hybrid_V8_moyen) /
            estim_avant$hybrid_V8_moyen * 100, 2), "%\n")
#====================================================================================
# PIEGE 3 — VERIFICATION DU BIAIS DE SELECTION
# Comparer les exclus vs conservés sur plusieurs indicateurs
# Si les exclus ont des caractéristiques très différentes → biais documenté
#====================================================================================

cat("\n=== PIEGE 3 — VERIFICATION BIAIS DE SELECTION ===\n")

biais_check <- donnees %>%
  mutate(
    statut = ifelse(
      decision %in% c("Exclure_hors_territoire", "Exclure_fraude"),
      "Exclus",
      "Conserves"
    )
  ) %>%
  group_by(statut) %>%
  summarise(
    N                  = n(),
    NDVI_moyen         = round(mean(ndvi_avg, na.rm=TRUE), 4),
    TWI_moyen          = round(mean(twi, na.rm=TRUE), 3),
    Precipitation_moy  = round(mean(anntot_avg, na.rm=TRUE), 1),
    Dist_route_moy_km  = round(mean(dist_road, na.rm=TRUE), 2),
    Elevation_moy_m    = round(mean(srtm_1k, na.rm=TRUE), 1),
    .groups = "drop"
  )

print(biais_check)

# --- Interprétation automatique ---
ndvi_exclus    <- biais_check$NDVI_moyen[biais_check$statut == "Exclus"]
ndvi_conserves <- biais_check$NDVI_moyen[biais_check$statut == "Conserves"]

if (length(ndvi_exclus) > 0 && length(ndvi_conserves) > 0) {
  diff_ndvi <- abs(ndvi_exclus - ndvi_conserves) / ndvi_conserves * 100
  if (diff_ndvi < 10) {
    cat("\n✅ Biais faible : NDVI des exclus proche des conservés (",
        round(diff_ndvi, 1), "% de différence)\n")
    cat("   Les exclusions ne créent pas de biais géographique systématique.\n")
  } else {
    cat("\n⚠️  Biais détecté : NDVI des exclus diffère de",
        round(diff_ndvi, 1), "% par rapport aux conservés\n")
    cat("   Ce biais doit être documenté dans le rapport.\n")
  }
}

# --- Visualisation biais ---
biais_long <- biais_check %>%
  select(statut, NDVI_moyen, TWI_moyen, Dist_route_moy_km, Elevation_moy_m) %>%
  pivot_longer(cols = -statut,
               names_to = "Indicateur",
               values_to = "Valeur")

ggplot(biais_long, aes(x = statut, y = Valeur, fill = statut)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Indicateur, scales = "free_y") +
  scale_fill_manual(values = c("Conserves" = "steelblue",
                               "Exclus"    = "tomato")) +
  labs(title = "Vérification du biais de sélection",
       subtitle = "Comparaison exclus vs conservés sur indicateurs géographiques",
       x = "", y = "Valeur moyenne") +
  theme_minimal() +
  theme(legend.position = "none")

#====================================================================================
# ANALYSE DE SENSIBILITE AU SEUIL Q
#====================================================================================
seuils_test <- c(0.90, 0.95, 0.975)
resultats_sensibilite <- data.frame()

for (q in seuils_test) {
  seuil_q <- quantile(donnees$distance_centre, q, na.rm = TRUE)
  
  base_test <- donnees %>%
    filter(!decision %in% c("Exclure_hors_territoire", "Exclure_fraude")) %>%
    filter(distance_centre <= seuil_q)
  
  resultats_sensibilite <- rbind(resultats_sensibilite, data.frame(
    Quantile        = paste0("Q", q * 100, "%"),
    N_conserves     = nrow(base_test),
    N_exclus        = nrow(donnees) - nrow(base_test),
    Pct_exclus      = round((nrow(donnees) - nrow(base_test)) /
                              nrow(donnees) * 100, 1),
    Surface_moy_ha  = round(mean(base_test$surface_cultivee_ha,
                                 na.rm=TRUE), 3),
    NDVI_moyen      = round(mean(base_test$ndvi_avg, na.rm=TRUE), 4),
    hybrid_V8_moy   = round(mean(base_test$hybrid_V8, na.rm=TRUE), 2)
  ))
}

cat("\n=== ANALYSE DE SENSIBILITE ===\n")
print(resultats_sensibilite)

# Visualisation
ggplot(resultats_sensibilite,
       aes(x = Quantile, y = Surface_moy_ha, group = 1)) +
  geom_line(color = "darkblue", linewidth = 1.2) +
  geom_point(size = 4, color = "darkblue") +
  geom_text(aes(label = paste0(Pct_exclus, "% exclus")),
            vjust = -1.5, size = 3.5, color = "darkred", fontface = "bold") +
  labs(
    title    = "Robustesse des estimations au choix du seuil d'exclusion des outliers",
    subtitle = "Comparaison Q90%, Q95% et Q97.5% --- stabilite confirmee des estimations",
    x        = "Quantile utilise",
    y        = "Surface cultivee moy. (ha)",
    caption  = "Source : MLI_2017_EAC-I_v03, Banque mondiale | Traitement : ZONGO, ISEP2 S4"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption  = element_text(face = "italic", color = "gray50", hjust = 0.5)
  )
ggsave("C:/Users/bmd/Desktop/Expose_R_ZONGO/sensibilite.png", 
       width = 6,    # largeur
       height = 4,   # hauteur
       dpi = 150)
#==================================================================================
# VISUALISATIONS FINALES
#====================================================================================

# --- Carte des niveaux de qualité ---
points_sf$niveau   <- donnees$niveau
points_sf$decision <- donnees$decision

ggplot() +
  geom_sf(data = mali, fill = "grey95", color = "black") +
  geom_sf(data = points_sf, aes(color = niveau),
          size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("Bon"        = "green3",
                                "A_verifier" = "orange",
                                "Suspect"    = "red")) +
  labs(title = "Carte des anomalies GPS — EACI Mali 2017",
       color  = "Niveau") +
  theme_minimal()

# --- Histogramme distances ---
ggplot(donnees, aes(x = distance_centre / 1000)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "white") +
  geom_vline(xintercept = seuil_dist / 1000,
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = seuil_dist / 1000 + 30, y = 40,
           label = "Seuil Q95%", color = "red") +
  labs(title = "Distribution des distances au centroide",
       x = "Distance (km)", y = "Frequence") +
  theme_minimal()
# --- Diagramme décisions ---

ggplot(donnees, aes(x = decision, fill = decision)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Repartition des decisions de traitement",
       x = "Decision", y = "Nombre de points") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

#====================================================================================
# EXPORT DE LA BASE NETTOYEE
#====================================================================================
dossier_sortie <- "C:/Users/bmd/Desktop/Expose_R_ZONGO/Fichiers_a_Rendre/Base_netoyee"
chemin_complet <- file.path(dossier_sortie, "gps_mali_nettoye_definitif.csv")
write.csv(donnees, chemin_complet, row.names = FALSE)
cat("\nFichier exporté : gps_mali_nettoye_definitif.csv\n")
