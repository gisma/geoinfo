#
# Type: script
# Name: change_detection_workflow.R
# Description:  basic reproducible Change Detection workflow for Sentinel data
#               data download -> digitizing training areas ->
#               extracting and preparing training data ->
#               example classifications: cluster analysis, random forest, MaxLike
#               quality assesment and comparison with basic tools, kappa, and motif
#               basic visualisation examples with tmap and mapview
# Dependencies:
# Output: original sentinel tile
#         AOI window of this tile (research_area)
#         training areas
# Copyright: GPL (>= 3), Chris Reudenbach, creuden@gmail.com
#

#--- laden der notwendigen Bibliotheken
# Achtung Pakete müssen evtl. manuell installiert werden
library(envimaR)
library(tmaptools)
library(motif)
library(kableExtra)

#--- Schalter für den Download der sentinel daten
get_sen = FALSE

#--- schalter ob digitalisiert werden muss falls er auf FALSE gesetzt ist werden die
# (zuvor erstellten und gesicherten Daten ) im else Teil der Verzweigung eingelesen
digitize = FALSE

## setzen des aktuellen Projektverzeichnisses (erstellt mit envimaR) als root_folder
#root_folder = find_rstudio_root_file()
root_folder = "~/edu/geoinfo/"

# einlesen des zuvor erstellten Setup-Skripts
source(file.path(root_folder, "src/functions/000_setup.R"))

# weitere Parameter

# größe der Aggregationsfesnter für die CD Auswertung Distanz-Schwellenwert
# "dist" ist abhängig von wins_size je kleiner win_size desto größer der
# Schwellenwert
win_size = 25

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
# Daher die Daten manuell herunterladen und in das Verzeichnis kopieren und entpacken
# Dann auskommentierten Tail ausführen
# corine_eu = raster(file.path(envrmt$path_data_lev0,"u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"))
# tmp = projectRaster(pred_stack_2019[[1]],crs = crs(corine_eu))
# corine_crop = raster::crop(corine_eu,tmp)
# corine_utm = projectRaster(corine_crop,crs = crs(pred_stack_2019))
# corine = resample(corine_utm,pred_stack_2019[[1]])
# raster::writeRaster(corine,file.path(envrmt$path_data_lev0,"/corine.tif"),overwrite=TRUE)

# ansonsten den Beipieldatensatz laden
corine = raster::raster(file.path(envrmt$path_data_lev0,"corine.tif"))

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

# Visualisierung mit Mapview: mit dem + Zeichen können beliebige Raster und
# Vektor-Datenebnenen kombiniert werden Bei der Ansicht fällt auf dass
# offensichtlich der Clusteralgorithmus keine zielführenden Klassen der
# Einstellung nClasses=2 erzeugt
mapview(mask*prediction_kmeans_2019$map,
        at = seq(0, nclasses, 1), # Anzahl der Legendenklassen
        legend = TRUE,            # Legende zeigen
        alpha.regions = 1,        # layer undurchsichtig
        maxpixels =  1693870) +  #  volle auflösung
  viewRGB(mask * pred_stack_2019, # funktion viewRGB zeigt dreikanalige Raster Bilder
          r=4,g=5,b=6,            # Kanäle aus dem angegebenen Rasterstack
          maxpixels =  1693870) +
  mapview(mask*prediction_kmeans_2020$map, # RSToolbox Rasterobjekte werden unter $map gespeichert
          at = seq(0, nclasses, 1),
          legend = FALSE,
          alpha.regions = 1,
          maxpixels =  1693870) +  # neue Ebene
  viewRGB(mask * pred_stack_2020,
          r=4,g=5,b=6,
          maxpixels =  1693870)

#---- Digitalisierung der Trainingsdaten ----
# Muß natürlich nur einmal gemacht werden
if (digitize) {
  source(file.path(root_folder, "src/functions/digitize.R"))
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


#---- Visualisierung mit tmap qtm() ----

# (interaktiv = "view", statisch = "plot")
tmap_mode("plot")
w1 = qtm(prediction_rf_2019,projection = 32632,inner.margins=0.01) +
  tm_legend(scale= 0.5,
            legend.outside=T,
            title = "RF 2019",
            title.size = 1.0) +
  tm_grid()

w2 = qtm(prediction_kmeans_2019$map,projection = 32632,inner.margins=0.01) +
  tm_legend(scale= 0.5,
            legend.outside=T,
            title = "kmeans 2019",
            title.size = 1.0) +
  tm_grid()

w3 = qtm(prediction_rf_2020,projection = 32632,inner.margins=0.01) +
  tm_legend(scale= 0.5,
            legend.outside=T,
            title = "RF 2020",
            title.size = 1.0) +
  tm_grid()

w4 = qtm(prediction_kmeans_2020$map,projection = 32632,inner.margins=0.01) +
  tm_legend(scale= 0.5,
            legend.outside=T,
            title = "kmeans 2020",
            title.size = 1.0) +
  tm_grid()
# ohne jede Anordnung
w2;w4;w1;w3

# mit der bordeigenen tmap fnktion tmap_arrrange()
tmap_arrange(w1, w3, w2, w4,asp = NA)

# mit dem paket gid (universell nutzbar um graphische ausgabe zu gestalten)
# hier die aufteilung wi tmap_arrange()
library(grid)
grid.newpage()
print(w2, vp=viewport(0.25, 0.25, width = 0.5, height = 0.5))
print(w1, vp=viewport(0.25, 0.75, width = 0.5, height = 0.5))
print(w4, vp=viewport(0.75, 0.25, width = 0.5, height = 0.5))
print(w3, vp=viewport(0.75, 0.75, width = 0.5, height = 0.5))

# mit dem paket cowplot
library(cowplot)
plot_grid(tmap_grob(w1),tmap_grob(w3),tmap_grob(w2),tmap_grob(w4),
          labels = c('FIG A', 'FIG B', 'FIG C', 'FIG D'),
          label_size= 10, vjust = 5.)


# ---- Maximum Likelihood Classification ----
# superClass() Funktion aus dem Paket RSToolbox umwandeln der Tabelle in das
# geforderte SpatialdataPoint Objekt Identifikation der Koordinaten
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

## ---- Visualisierung mit mapview ----
mapview::viewRGB(mask*pred_stack_2020, r = 4, g =5, b = 6,maxpixels =  1693870)+
  mapview(mask*prediction_rf_2019 , alpha.regions = 0.5, maxpixels =  1693870,
          col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = TRUE) +
  mapview(mask*prediction_rf_2020, alpha.regions = 0.5, maxpixels =  1693870,
          col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = FALSE) +
  mapview(mask*prediction_mlc_2019$map,alpha.regions = 0.5, maxpixels =  1693870,
          col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = FALSE) +
  mapview(mask*prediction_mlc_2020$map,alpha.regions = 0.5, maxpixels =  1693870,
          col.regions = mapviewPalette("mapviewRasterColors"),at = seq(0, nclasses, 1), legend = FALSE)


## ---- Change Detection Auswertung ----

# Differenzbild random forest
qtm_min(prediction_rf_2020 - prediction_rf_2019,title = "kmeans 2019")

# ---- CDA Extraktion der Klassennamen ----
categories = levels(prediction_rf_2019)[[1]]$value
categories

# ---- Berechnung der Vierfeld Tabelle mit raster Basisfunktion ----
ct = raster::crosstab(prediction_rf_2019,prediction_rf_2020)
rownames(ct) = categories
colnames(ct) = categories
ct %>%
  kbl(caption = "Crosstab 2019 vs 2020",) %>%
  kable_classic(full_width = F)

# add_header_above(c("layer 1 " = 1, "layer 2" = 2))

# ---- kappa ----
# Vergleich der Übereinstimmung unterschiedlicher Klassifikationen (hier MaxLike
# und RF) mit Hilfe diverser Kappa Werte
# https://giswerk.org/doku.php?id=r:r-tutorials:calculatekappa
# 2019
kappa_2019 = kstat(prediction_mlc_2019$map,prediction_rf_2019, perCategory = FALSE)
kappa_2019
kappa_2019  %>%
  kbl(caption = "Kappa 2019 MLC vs RF",) %>%
  kable_classic(full_width = F)

# 2020
kappa_2020 = kstat(prediction_mlc_2020$map,prediction_rf_2020,perCategory = FALSE)
kappa_2020
kappa_2020  %>%
  kbl(caption = "Kappa 2019 MLC vs RF",) %>%
  kable_classic(full_width = F) %>%
  column_spec(1,color = "red", link = "https://giswerk.org/doku.php?id=r:r-tutorials:calculatekappa")


# Erzeugen aller Vergleichsraster der Kontingenztabelle
# Die lapply funktionen sind integrierte FOR Schleifen die über die Liste der
# Kategorien die Funktion changefrom() für die Kreuztabelle anwenden
r = lapply(1:length(categories),
           function(i){lapply(1:length(categories),
                              function(j){changefrom(prediction_rf_2019, prediction_rf_2020, i,j)})})
r

# ---- Visualsierung der Kreuztabellierten Von-Zu-Raster ----
# Plotten der Raster hierzu werden erneut alle Kategorien einzeln geplottet i
# und j sind hilfsvariablen um die korrekten Raster Layer ansprechen zu können.
# t ist eine Hilfvariable um eine Liste für die Ergebnisbilder hochzählen zu
# können
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


# ---- Analyse von Veränderungen mit dem paket motif ----
# lsp_compare vergleicht zwei (oder mehr) kategoriale Karten (change detection
# klassifikationen etc.) miteinander und nutzt für die Ausgabe der
# Wahrscheinlichkeiten verschiedener räumlicher Aggregierungsstufen und
# Merkmalsraum-Distanzen um Veränderungswahrscheinlichkeiten zu ermitteln
mrf_compare_2020_2019 = lsp_compare(st_as_stars(mask*prediction_rf_2019),
                                    st_as_stars(mask*prediction_rf_2020),
                                    type = "cove",
                                    dist_fun = "jensen-shannon",
                                    window = win_size,
                                    threshold = 0.9)

# Visualisierung der Gesamtwahrscheinlichkeiten Achtung logarithmische Skala

tmap_mode("plot")
tm_compare2 = tm_shape(mrf_compare_2020_2019) +
  tm_raster("dist",
            breaks = c(0, 0.001, 0.01, 0.1, 1.01),
            legend.is.portrait = F,
            palette = "viridis",
            title =  "Distance (JSD)") +
  tm_layout(legend.show = TRUE,
            legend.text.size = 0.3,
            legend.outside = TRUE) +
  tm_legend(scale= 0.5,
            legend.outside=T,
            title = paste("RF 2019 vs 2020 ",win_size*10,"m**2" ),
            title.size = 1.0) +
  tm_grid()

tm_compare2

# ---- Detail-Analyse Teil 1 ----
# Identifikation der Gebiete (Referenz ist win_size) die maximale Veränderungen
# aufweisen im Beispiel soll  "dist" soll größer 0.001 sein (logarithmische
# Skalierung!)
lc_am_compare_sel = st_as_sf(mrf_compare_2020_2019) %>%
  subset(dist > 0.01)

# Sortierung nach Größe
lc_am_compare_sel = lc_am_compare_sel[order(lc_am_compare_sel$dist,decreasing = TRUE), ]
lc_am_compare_sel

##- Visualsierung
# Plotten der top ten gebiete wir nutzen die selbst geschriebene Funktion tm_plot()
# siehe src/functions/fun_tmap.R
tm_plot(sel=lc_am_compare_sel[1:nrow(lc_am_compare_sel),],ov = T)
# als interaktive karte
tm_plot(sel=lc_am_compare_sel[1:nrow(lc_am_compare_sel),],vt = "view",ov = T)


# ---- Visualsierng Detail-Analyse ----
# Extraktion nach der sortierten Liste lc_comapare_sel  werden die top 10 change
# detection hotspots Daten extrahiert und in die liste lc geschrieben
lc = lapply(seq(1:10),function(i){
              lsp_extract(c(st_as_stars(prediction_rf_2019),
                            st_as_stars(prediction_rf_2020)),
                          window = win_size,
                          id = lc_am_compare_sel$id[i])
            })

# Erzeugen der top 3 Change Detection Hotspots Karten-Ausschnitte
tm_lc = lapply(seq(1:3),function(i){
               tm_plot(lc[[i]])
             })

# plotten der erzeugten Karten
tmap_arrange(tm_lc)


#####

