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
win_size = 20

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

# Einlesen der Corine Daten
# Für den Download ist ein Konto notwendig https://land.copernicus.eu/pan-european/corine-land-cover
# angenommen die Daten liegen im level0 Verzeichnis

if (!file.exists(file.path(envrmt$path_data_lev0,"u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))){
corine_eu = raster(file.path(envrmt$path_data_lev0,"u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))
tmp = projectRaster(pred_stack_2019[[1]],crs = crs(corine_eu))
corine_crop = raster::crop(corine_eu,tmp)
corine_utm = projectRaster(corine_crop,crs = crs(pred_stack_2019))
corine = resample(corine_utm,pred_stack_2019[[1]])
raster::writeRaster(corine,file.path(envrmt$path_data_lev0,"/corine.tif"),overwrite=TRUE)
} else{
corine = raster::raster(file.path(envrmt$path_data_lev0,"/corine.tif"))
}
# Erstellen einer Wald-Maske
# Agro-forestry areas code=22, Broad-leaved forest code=23,
# Coniferous forest code=24, Mixed forest code=25
mask = reclassify(corine,c(-100,22,0,22,26,1,26,500,0))

# Stack-Loop über die Daten
for (pat in c("RGB432B","EVI","MSAVI2","NDVI","SAVI")){
  pred_stack_2019 = raster::stack(pred_stack_2019,raster::stack(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20190619",full.names = TRUE)))
  pred_stack_2020 = raster::stack(pred_stack_2020,raster::stack(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20200623",full.names = TRUE)))
}
# get rid of NA
pred_stack_2019 = reclassify(pred_stack_2019, cbind(NA, 0))
pred_stack_2020 = reclassify(pred_stack_2020, cbind(NA, 0))

# Zuweisen von leserlichen Namen auf die Datenebenen
names(pred_stack_2019) = c("nir","red_1","green_1","red_2","green_2","blue","EVI","MSAVI2","NDVI","SAVI")
names(pred_stack_2020) = c("nir","red_1","green_1","red_2","green_2","blue","EVI","MSAVI2","NDVI","SAVI")
saveRDS(pred_stack_2019,paste0(envrmt$path_data,"pred_stack_2019.rds"))
saveRDS(pred_stack_2020,paste0(envrmt$path_data,"pred_stack_2020.rds"))
pred_stack_2019 = readRDS(paste0(envrmt$path_data,"pred_stack_2019.rds"))
pred_stack_2020 = readRDS(paste0(envrmt$path_data,"pred_stack_2020.rds"))

##---- Kmeans Klassifikation der Daten----
## k-means über RStoolbox
# Verändern Sie den Wert von nclasses und vergleichen sie die Ergebnisse
nclasses = 2
prediction_kmeans_2019 = unsuperClass(pred_stack_2019,
                                      nClasses = nclasses,
                                      norm = TRUE,
                                      algorithm = "MacQueen" )
prediction_kmeans_2020 = unsuperClass(pred_stack_2020,
                                      nClasses = nclasses,
                                      norm = TRUE,
                                      algorithm = "MacQueen")

# Visualisierung mit Mapview: mit dem + Zeichen können beliebige Raster und Vektor-Datenebnenen kombiniert werden
# Offensichtlich ist der Clusteralgorithmus überfordert mit zwei Klassen ein sinnvolles Ergebnis zu erzeugen
mapview(mask*prediction_kmeans_2019$map,
        col.regions = mapviewPalette("mapviewRasterColors"),
        at = seq(0, nclasses, 1),
        legend = TRUE,
        alpha.regions = 1,
        maxpixels =  1693870) +  # neue Ebene
  viewRGB(mask * pred_stack_2019,
          r=4,g=5,b=6,
          maxpixels =  1693870) +  # neue Ebene
  mapview(mask*prediction_kmeans_2020$map,
          col.regions = mapviewPalette("mapviewRasterColors"),
          at = seq(0, nclasses, 1),
          legend = FALSE,
          alpha.regions = 1,
          maxpixels =  1693870) +  # neue Ebene
  viewRGB(mask * pred_stack_2020,
          r=4,g=5,b=6,
          maxpixels =  1693870)

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
#names(train_areas_2019) = c("class","nir","red_1","green_1","red_2","green_2","blue","EVI","MSAVI2","NDVI","SAVI","X","Y","cell","coverage_fraction")
} else {
  train_areas_2019 = readRDS(paste0(envrmt$path_data,"train_areas_2019.rds"))
  train_areas_2020 = readRDS(paste0(envrmt$path_data,"train_areas_2020.rds"))
}



## ---- Überwachte  mit Klassifikation Random Forest ----
## Hier wird der Random Forest über das Utility Paket caret aufgerufen

# Setzen eines "seed" ermöglicht reproduzierbaren Zufall
set.seed(123)

# Zufälliges Ziehen von 25% der Daten
trainDat_2019 =  train_areas_2019[createDataPartition(train_areas_2019$class,list = FALSE,p = 0.25),]
trainDat_2020 =  train_areas_2020[createDataPartition(train_areas_2020$class,list = FALSE,p = 0.25),]
#  Response-Variable (=Spalte "class") muss den Datentyp "factor" haben
trainDat_2019$class <- as.factor(trainDat_2019$class)
trainDat_2020$class <- as.factor(trainDat_2020$class)

# Einstellungen Modelltraining: cross-validation, 10 Wiederholungen
ctrlh = trainControl(method = "cv",
                     number = 10,
                     savePredictions = TRUE)

#--- random forest model training
cv_model_2019 = train(trainDat_2019[,2:11], # in den Spalten 2 bis 20 stehen die Trainingsdaten (Prediktoren genannt)
                 trainDat_2019[,1],         # in der Spalte 1 steht die zu klassizierende Variable (Response genannt)
                 method = "rf",             # Methode hier rf für random forest
                 metric = "Kappa",          # Qualitäts/Performanzmaß KAppa
                 trControl = ctrlh,         # obig erzeugte Trainingssteuerung soll eingelsen werden
                 importance = TRUE)         # Die Bedeung der Variablen wird mit abgespeichert
# Klassifikation (auch Vorhersage genannt)
prediction_rf_2019  = raster::predict(pred_stack_2019 ,cv_model_2019, progress = "text")

#--- 2020
# Einstellungen wie 2019
cv_model_2020 = train(trainDat_2020[,2:11],
                      trainDat_2020[,1],
                      method = "rf",
                      metric = "Kappa",
                      trControl = ctrlh,
                      importance = TRUE)
prediction_rf_2020  = raster::predict(pred_stack_2020 ,cv_model_2020, progress = "text")


#---- Visualisierung mit tmap qtm (interaktive = "view", statisch = "plot")----
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


## Maximum Likelihood Classification
#  superClass() Funktion aus dem Paket RSToolbox

# umwandeln der Tabelle in das geforderte SpatialdataPoint Objekt
# Identifikation der Koordinaten
sp::coordinates(trainDat_2019) = ~X+Y
sp::coordinates(trainDat_2020) = ~X+Y
# Zuweisen der korrekten Projektion (hier aus dem zugehörigen Raster-Stack)
crs(trainDat_2019) = crs(pred_stack_2019)
crs(trainDat_2020) = crs(pred_stack_2020)

# superClass method "mlc" berechnet die Statisik und klassifiziert in einem aufruf
prediction_mlc_2019       <- superClass(pred_stack_2019, trainData = trainDat_2019[,1:11], responseCol = "class",
                                  model = "mlc", tuneLength = 1, trainPartition = 0.5)
prediction_mlc_2020       <- superClass(pred_stack_2020, trainData = trainDat_2020[,1:11], responseCol = "class",
                                        model = "mlc", tuneLength = 1, trainPartition = 0.5)
# Ergebnisse der rf Klassifikation
prediction_mlc_2019
prediction_mlc_2020

## Interaktive Karte
mapview(mask*prediction_rf_2019 , alpha.regions = 1, maxpixels =  1693870,
        col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = TRUE) +
  mapview(mask*prediction_rf_2020, alpha.regions = 1, maxpixels =  1693870,
          col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = FALSE) +
  mapview(mask*prediction_mlc_2019$map,alpha.regions = 1, maxpixels =  1693870,
          col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = FALSE) +
  mapview(mask*prediction_mlc_2020$map,alpha.regions = 1, maxpixels =  1693870,
          col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = FALSE)


# Differenzbild random forest
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
# Berechnung der diversen Kappa Werte
# https://giswerk.org/doku.php?id=r:r-tutorials:calculatekappa
kstat(prediction_mlc_2019$map,prediction_rf_2019, perCategory = FALSE)
kstat(prediction_mlc_2020$map,prediction_rf_2020,perCategory = FALSE)


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
# Interaktive und synchronisierte Karten
tmap::tmap_arrange(m,sync = TRUE)

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
mrf_compare_2020_2019 = lsp_compare(st_as_stars(prediction_rf_2019), st_as_stars(prediction_rf_2020),
                        type = "cove", dist_fun = "jensen-shannon",
                        window = win_size, threshold = 0.9)


# Visualisierung der Gesamtwahrscheinlichkeiten
# logarithmische Skala
tmap_mode("plot")
tm_compare2 = tm_shape(mrf_compare_2020_2019) +
  tm_raster("dist",
            breaks = c(0, 0.001, 0.01, 0.1, 1.01),
            legend.is.portrait = F,
            palette = "viridis",
            title =  "Distance (JSD)") +
  tm_layout(legend.show = TRUE,
            legend.text.size = 0.3,
            legend.position = c(.0, 0.15),
            legend.outside.position = "bottom",
            legend.outside = TRUE)
tm_compare2

# Identifikation der Gebiete (Referenz ist win_size) die maximale Veränderungen aufweisen
# hier der logarithmisch skalierte Paramter "dist" soll größer 0.1 sein
lc_am_compare_sel = st_as_sf(mrf_compare_2020_2019) %>%
  subset(dist > 0.4)
# Sortierung des Ergebnis nach Größe
lc_am_compare_sel = lc_am_compare_sel[order(lc_am_compare_sel$dist,
                                            decreasing = TRUE), ]
# Extraktion der Top 3 Gebiete
lc=list()
for (i in 1:length(lc_am_compare_sel$id)){
lc[[i]] = lsp_extract(c(st_as_stars(prediction_rf_2019), st_as_stars(prediction_rf_2020)), window = win_size, id = lc_am_compare_sel$id[i])
# lc_2 = lsp_extract(c(st_as_stars(prediction_rf_2019), st_as_stars(prediction_rf_2020)), window = win_size, id = lc_am_compare_sel$id[2])
# lc_3 = lsp_extract(c(st_as_stars(prediction_rf_2019), st_as_stars(prediction_rf_2020)), window = win_size, id = lc_am_compare_sel$id[3])
}
# Interaktive Visualisierung
tm_plot = function(x,vt="plot",ls=FALSE,lo=FALSE){
  tmap_mode(vt)
    tm_shape(x) +
    tm_raster(style = "cat",
              title = "Land cover:",
              legend.is.portrait = FALSE) +
    tm_facets(ncol = 2) +
    tm_layout(legend.show = ls,
              legend.only = lo,
              legend.position = c(.0, 0.15),
              legend.outside.position = "bottom",
              legend.outside = T,
              panel.label.height=0.6,
              panel.label.size=0.6,
              panel.labels = c(2019, 2020))
}

tm = list()
for (i in 1:length(lc)){
tm[[i]] = tm_plot(lc[[i]])
if (i == length(lc))
  tm[[i]] = tm_plot(lc[[i]],ls=T)
}

tmap_arrange(lapply(1:length(tm),function(i){tm[[i]]}),ncol=2,nrow = ceiling(length(tm)/2))

#####
