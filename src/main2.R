#------------------------------------------------------------------------------
# Type: script
# Name: get_sentinel.R
# Author: Chris Reudenbach, creuden@gmail.com
# Description:  retrieves sentinel data
#               and exemplary defines AOI and calculates albedo
# Dependencies:
# Output: original sentinel tile
#         AOI window of this tile (research_area)
#         top of atmosphere albedo AOI
#         surface albedo AOI
# Copyright: GPL (>= 3)
#------------------------------------------------------------------------------

#--- laden der notwendigen Bibliotheken
# Achtung Pakete müssen evtl. manuell installiert werden
library(envimaR)
library(tmaptools)
library(motif)
library(spatialEco)

#--- Schalter für den Download der sentinel daten
get_sen = FALSE

#--- schalter ob digitalisiert werden muss falls er auf FALSE gesetzt ist werden die
# (zuvor erstellten und gesciherten Daten ) im else Teil der Verzeigung eingelesen
digitize = FALSE

## setzen des aktuellen Projektverzeichnisses (erstellt mit envimaR) als root_folder
#root_folder = find_rstudio_root_file()
root_folder = "~/edu/geoinfo/"

# einlesen des zuvor erstellten Setup-Skripts
source(file.path(root_folder, "src/functions/000_setup.R"))

# weitere Parameter

# größe der Aggregationsfesnter für die CD Auswertung
win_size = 26

#--- Download der Daten
# gui = TRUE ruft die GUI zur Kontrolle auf
if (get_sen){
   out_paths_3 <- sen2r(
    gui = T,
    param_list = "~/edu/geoinfo/data/harz.json",
    tmpdir = envrmt$path_tmp,
  )
}

#--- Einlesen der Daten aus den Verzeichnissen
# RGB stack der beiden Jahre
pred_stack_2019 = raster::stack(list.files(file.path(envrmt$path_data_lev1,"RGB843B"),pattern = "20190619",full.names = TRUE))
pred_stack_2020 = raster::stack(list.files(file.path(envrmt$path_data_lev1,"RGB843B"),pattern = "20200623",full.names = TRUE))

# Einlesen der Corinne Daten
# Für den Download ist ein Konto notwendig https://land.copernicus.eu/pan-european/corine-land-cover
# angenommen die Daten liegen im level0 Verzeichnis

if (!file.exists(file.path(envrmt$path_data_lev0,"u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))){
corine_eu = raster(file.path(envrmt$path_data_lev0,"u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))
tmp = projectRaster(pred_stack_2019[[1]],crs = crs(corine_eu))
corine_crop = raster::crop(corine_eu,tmp)
cronine_utm = projectRaster(corine_crop,crs = crs(pred_stack_2019))
corine_crop = raster::crop(cronine_utm,pred_stack_2019)
corine = resample(corine_crop,pred_stack_2019[[1]])
raster::writeRaster(corine,file.path(envrmt$path_data_lev0,"/corine.tif"),overwrite=TRUE)
} else{
corine = raster::raster(file.path(envrmt$path_data_lev0,"/corine.tif"))
}
mapview(corine) +pred_stack_2019[[1]]


# Stack-Loop über die Daten
for (pat in c("RGB432B","EVI","MSAVI2","NDVI","SAVI")){
  pred_stack_2019 = raster::stack(pred_stack_2019,raster::stack(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20190619",full.names = TRUE)))
  pred_stack_2020 = raster::stack(pred_stack_2020,raster::stack(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20200623",full.names = TRUE)))
}

# Zuweisen von leserlichen Namen auf die Datenebenen
names(pred_stack_2019) = c("nir","red","green","red","green","blue","EVI","MSAVI2","NDVI","SAVI")
names(pred_stack_2020) = c("nir","red","green","red","green","blue","EVI","MSAVI2","NDVI","SAVI")
saveRDS(pred_stack_2019,paste0(envrmt$path_data,"pred_stack_2019.rds"))
saveRDS(pred_stack_2020,paste0(envrmt$path_data,"pred_stack_2020.rds"))
pred_stack_2019 = readRDS(paste0(envrmt$path_data,"pred_stack_2019.rds"))
pred_stack_2020 = readRDS(paste0(envrmt$path_data,"pred_stack_2020.rds"))
##---- Kmeans Klassifikation der Daten----

## k-means über RStoolbox
# Verändern Sie den Wert von nclasses und vergleichen sie die Ergebnisse
prediction_kmeans_2019 = unsuperClass(pred_stack_2019, nClasses = 2,norm = TRUE, algorithm = "MacQueen" )
prediction_kmeans_2020 = unsuperClass(pred_stack_2020, nClasses = 2,norm = TRUE, algorithm = "MacQueen")

# Visualisierung Note mit dem + Zeichen können beliebige raster und vektor  Datenebnenen kombiniert werden
mapview(prediction_kmeans_2019$map,col.regions = mapviewPalette("mapviewTopoColors"), at = seq(0, 2, 1), legend = TRUE,alpha.regions = 0.5) +
mapview(prediction_kmeans_2020$map,col.regions = mapviewPalette("mapviewTopoColors"), at = seq(0, 2, 1), legend = FALSE,alpha.regions = 0.5)


#---- Digitalisierung der Trainingsdaten ----

if (digitize) {
# Für die überwachte Klassifikation benötigen wir Trainingsgebiete. Sie können Sie wie nachfolgend digitalisieren oder alternativ z.B. QGis verwenden

#--- 2019
# Kahlschlag
# Es wird das Falschfarbenkomosit in originaler Auflösung genutzt (maxpixels =  1693870)
# Bitte beachten Sie dass es (1) deutlich länger lädt und (2) Vegetation in Rot dargestellt wird.
# Die Kahlschäge sind jetzt grün
train_area <- mapview::viewRGB(pred_stack_2019, r = 1, g = 2, b = 3, maxpixels =  1693870) %>% mapedit::editMap()
# Hinzufügen der Attribute class (text) und id (integer)
clearcut <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "clearcut", id = 1)

# other: hier gilt es möglichst verteilt übers Bild möglichst alle nicht zu Kahlschlag  gehörenden Flächen zu erfassen.
train_area <- mapview::viewRGB(pred_stack_2019, r = 1, g = 2, b = 3) %>% mapedit::editMap()
other <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "other", id = 2)

# rbind  kopiert die beiden obigen Vektorobjekte in eine Datei
train_areas_2019 <- rbind(clearcut, other)

# Umprojizieren auf die Raster Datei
train_areas_2019 = sf::st_transform(train_areas_2019,crs = sf::st_crs(pred_stack_2019))


# Extraktion der Trainingsdaten für die digitalisierten Flächen
tDF_2019 = exactextractr::exact_extract(pred_stack_2019, train_areas_2019,  force_df = TRUE,
                                        include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "class")
#  auch hier wieder zusamenkopieren in eine Datei
tDF_2019 = dplyr::bind_rows(tDF_2019)

# Löschen von etwaigen Zeilen die NA (no data) Werte enthalten
tDF_2019 = tDF_2019[complete.cases(tDF_2019) ,]
tDF_2019 = tDF_2019[ ,rowSums(is.na(tDF_2019)) == 0]

# check der extrahierten Daten
summary(tDF_2019)
mapview(train_areas_2019)+pred_stack_2019[[1]]

# Abspeichern als R-internes Datenformat
# ist im Repo hinterlegt und kann desahlb (zeile drunter) eingeladen werden
saveRDS(tDF_2019, paste0(envrmt$path_data,"train_areas_2019.rds"))


# # ---- Das gleiche muss für 2020 wiederholt werden zum digitalisieren und extrahieren bitte ent-kommentieren ----

# # Kahlschlag
train_area <- mapview::viewRGB(pred_stack_2020, r = 4, g =5, b = 6,maxpixels =  1693870) %>% mapedit::editMap()
clearcut <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "clearcut", id = 1)
train_area <- mapview::viewRGB(pred_stack_2020, r = 4, g = 5, b = 6) %>% mapedit::editMap()
other <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "other", id = 2)
train_areas_2020 <- rbind(clearcut, other)
train_areas_2020 = sf::st_transform(train_areas_2020,crs = sf::st_crs(pred_stack_2020))
tDF_2020 = exactextractr::exact_extract(pred_stack_2020, train_areas_2020,  force_df = TRUE,
                                        include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "class")
tDF_2020 = dplyr::bind_rows(tDF_2020)
tDF_2020 = tDF_2020[  rowSums(is.na(tDF_2020)) == 0,]
saveRDS(tDF_2020, paste0(envrmt$path_data,"train_areas_2020.rds"))

} else {
  train_areas_2019 = readRDS(paste0(envrmt$path_data,"train_areas_2019.rds"))
  train_areas_2020 = readRDS(paste0(envrmt$path_data,"train_areas_2020.rds"))

}

mapview(train_areas_2020,fgb=F)

## ---- Überwachte  mit Klassifikation Random Forest ----

## Hier wird der Entscheidungsbaum Algorithmus Random Forest über das Utility Paket caret aufgerufen

# Setzen eines "seed" ermöglicht reproduzierbaren Zufall
set.seed(123)

# Zufälliges Ziehen von 25 der Daten
trainDat_2019 =  train_areas_2019[createDataPartition(train_areas_2019$class,list = FALSE,p = 0.25),]

#  Response-Variable (=Spalte "class") muss den Datentyp "factor" haben
trainDat_2019$class <- as.factor(trainDat_2019$class)

# Einstellungen Modelltraining: cross-validation, 10 Wiederholungen
ctrlh = trainControl(method = "cv",
                     number = 10,
                     savePredictions = TRUE)
#--- random forest model training

# Modelltraining
cv_model_2019 = train(trainDat_2019[,2:11], # in den Spalten 2 bis 20 stehen die Trainingsdaten (Prediktoren genannt)
                 trainDat_2019[,1],         # in der Spalte 1 steht die zu klassizierende Variable (Response genannt)
                 method = "rf",             # Methode hier rf für random forest
                 metric = "Kappa",          # Qualitäts/Performanzmaß KAppa
                 trControl = ctrlh,         # obig erzeugte Trainingssteuerung soll eingelsen werden
                 importance = TRUE)         # Die Bedeung der Variablen wird mit abgespeichert

# Klassifikation, Häufig auch Vorhersage genannt.
prediction_rf_2019  = raster::predict(pred_stack_2019 ,cv_model_2019, progress = "text")

#--- 2020
trainDat_2020 =  train_areas_2020[createDataPartition(train_areas_2020$class,list = FALSE,p = 0.25),]
trainDat_2020$class <- as.factor(trainDat_2020$class)
ctrlh = trainControl(method = "cv",
                     number = 10,
                     savePredictions = TRUE)
cv_model_2020 = train(trainDat_2020[,2:11],
                      trainDat_2020[,1],
                      method = "rf",
                      metric = "Kappa",
                      trControl = ctrlh,
                      importance = TRUE)
prediction_rf_2020  = raster::predict(pred_stack_2020 ,cv_model_2020, progress = "text")


# Visualisierung mit tmap qtm (interaktive = "view", statisch = "plot")
tmap_mode("view")
w1 = qtm(prediction_rf_2019)
w2 = qtm(prediction_kmeans_2019$map)
w3 = qtm(prediction_rf_2020)
w4 = qtm(prediction_kmeans_2020$map)
tmap_arrange(w1, w2, w3, w4)
tmap_mode("plot")
w1 = qtm(prediction_rf_2019)
w2 = qtm(prediction_kmeans_2019$map)
w3 = qtm(prediction_rf_2020)
w4 = qtm(prediction_kmeans_2020$map)
tmap_arrange(w1, w2, w3, w4, widths = c(.33, .66))


## Hier wird eine Maximum Likelihood Classification
#  mit superClass aus dem Paket RSToolbox durchgeführt

# umwandeln der Tabelle in das geforderte SpatialdataPoint Objekt
# Identifikation der Koordinaten
trainDat_2020 = sp::coordinates(trainDat_2020) = ~x+y
# Zuweisen der korrekten Projektion (hier aus dem zugehörigen Raster-Stack)
crs(trainDat_2020) = crs(pred_stack_2020)

# superClass method "mlc" berechnet die Statisik und klassifiziert in einem aufruf
prediction_mlc_2020       <- superClass(pred_stack_2020, trainData = trainDat_2020[,1:11], responseCol = "class",
                                  model = "mlc", tuneLength = 1, trainPartition = 0.5)
# Ergebnisse
prediction_mlc_2020

## Interaktive Karte
tmap_mode("view")
qtm(prediction_rf_2020 ) + qtm(prediction_mlc_2020$map)


# Differenzbild
comp1 = prediction_rf_2020 - prediction_rf_2019
qtm(comp1)

# Extraktion KLASSENNAMEN
categories = levels(prediction_rf_2019)[[1]]$value
categories

# Berechnung der Kontingenz Tabelle mit raster Basisfunktion
ct = raster::crosstab(prediction_rf_2020,prediction_rf_2019)
rownames(ct) = categories
colnames(ct) = categories
ct
# Das Visualisieren der Veränderungen von jeder Klasse zu jeder Klasse
# mit Bordmitteln ist etwas aufwendiger. Es werden dazu eininge Hilfsfunktionen
# und Schleifen (loops) benötigt, da jede KAtegorie mit jeder Vergleihen werden muss

# Hilfsfunktion vergleicht je zwei Raster
changefrom=function(r1,r2,i,j){r1==i & r2==j}

# Erzeugen aller vergleichsraster aus der Kontingenztabelle
# mit der erstellten Hilfsunktion Lapply ist eine funktion die über eine Liste looped
r = lapply(1:length(categories), function(i){lapply(1:length(categories), function(j){changefrom(prediction_rf_2019, prediction_rf_2020, i,j)})})

# Plotten der Raster  hierzu werden erneut alle Kategorien einzeln geplottet
# i und j sind hilfsvariablen um die korrekten Raster Layer ansprechen zu können.
# t ist eine Hilfvariable um eine Liste für die Ergebnisbilder hochzählen zu können

tmap_mode("view")
t=i=j=1 # setze zählvariablen auf 1
m=list() # erzeuge leere Liste für die Karten
for(cat1 in categories)  { # für jede Kategorie
  for(cat2 in categories)  { # mit jeder Kategorie
    # plotte die Karte
    m[[t]]  = tm_shape(st_as_stars(r[[i]][[j]])) +
    tm_raster(col = "layer",palette = "viridis",style = "cat",
              title = if(cat1==cat2) {
                paste("no changes " , unique(cat1,cat2))
              }
              else if (cat1!=cat2) {
                paste(cat1," -> ",cat2)
              })
    # zähle die innere schleife hoch
    j = j + 1
    # zähle die ergebnisliste hoch
    t = t + 1
  }
  # innere schleife abgearbeitet setze sie auf  den Anfang zurück
  j = 1
  # zähle die äußere Schleife hoch
  i = i + 1
}
# plotte die erzeugten Karten
tmap::tmap_arrange(m)

# #-- Das spatialEco Paket: Berechnung von Kappa für die win_size Umgebung
# r.kappa2 <- spatialEco::raster.change(prediction_rf_2019,
#                          prediction_rf_2020,
#                          stat="kappa",
#                          force.memory = TRUE)
# mapview(r.kappa2)

#-- Das package motif stellt ausgezeichnete Werkzeuge für die Change Detection Postanalyse zur Verfügung

# lsp_compare vergleicht zwei kategoriale Karten (change detction klassifikationen etc)
# miteinander und nutzt für die Ausgabe der Wahrscheinlichkeiten verschiedene
# Aggregierungsstufen und Wahrscheinlichkeitsschwellenwert eines zu in dist_fun bestimmten Entfernungsmasses
mrf_compare_2020_2019 = lsp_compare(st_as_stars(prediction_rf_2020), st_as_stars(prediction_rf_2019),
                        type = "cove", dist_fun = "jensen-shannon",
                        window = win_size, threshold = 0.9)


# Visualisierung der Gesamtwahrscheinlichkeiten
# logarithmische Skala
tmap_mode("view")
my_breaks = c(0, 0.001, 0.01, 0.1, 1.01)
tm_compare2 = tm_shape(mrf_compare_2020_2019) +
  tm_raster("dist",breaks = my_breaks, palette = "viridis")
tm_compare2

# Identifikation der Gebiete (Referenz ist win_size) die maximale Veränderungen aufweisen
# hier der logarithmisch skalierte Paramter "dist" soll größer 0.1 sein
lc_am_compare_sel = st_as_sf(mrf_compare_2020_2019) %>%
  subset(dist > 0.1)
# Sortierung des Ergebnis nach Größe
lc_am_compare_sel = lc_am_compare_sel[order(lc_am_compare_sel$dist,
                                            decreasing = TRUE), ]
# Extraktion der Top 3 Gebiete
lc_1 = lsp_extract(c(st_as_stars(prediction_rf_2019), st_as_stars(prediction_rf_2020)), window = win_size, id = 2446)
lc_2 = lsp_extract(c(st_as_stars(prediction_rf_2019), st_as_stars(prediction_rf_2020)), window = win_size, id = 2099)
lc_3 = lsp_extract(c(st_as_stars(prediction_rf_2019), st_as_stars(prediction_rf_2020)), window = win_size, id = 2040)

# Interaktive Visualisierung
tm_plot = function(x,bg){
  bg +
    tm_raster() +
    tm_shape(x) +
    tm_raster(style = "cat",
              title = "Land cover:") +
    tm_facets(ncol = 2) +
    tm_layout(legend.show = FALSE,
              panel.labels = c(2020, 2019))

}

tm1 = tm_plot(lc_1,tm_compare2)
tm2 = tm_plot(lc_2,tm_compare2)
tm3 = tm_plot(lc_3,tm_compare2)


#####

