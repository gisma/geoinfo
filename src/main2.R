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

# Stack-Loop über die Daten
for (pat in c("RGB432B","EVI","MSAVI2","NDVI","SAVI")){
  pred_stack_2019 = raster::stack(pred_stack_2019,raster::stack(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20190619",full.names = TRUE)))
  pred_stack_2020 = raster::stack(pred_stack_2020,raster::stack(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20200623",full.names = TRUE)))
}

# Zuweisen von leserlichen Namen auf die Datenebenen
names(pred_stack_2019) = c("nir","red","green","red","green","blue","EVI","MSAVI2","NDVI","SAVI")
names(pred_stack_2020) = c("nir","red","green","red","green","blue","EVI","MSAVI2","NDVI","SAVI")


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
# train_area <- mapview::viewRGB(pred_stack_2020, r = 1, g = 2, b = 3,maxpixels =  1693870) %>% mapedit::editMap()
# clearcut <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "clearcut", id = 1)
# train_area <- mapview::viewRGB(pred_stack_2020, r = 1, g = 2, b = 3) %>% mapedit::editMap()
# other <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "other", id = 2)
# train_areas_2020 <- rbind(clearcut, other)
# train_areas_2020 = sf::st_transform(train_areas_2020,crs = sf::st_crs(pred_stack_2020))
# tDF_2020 = exactextractr::exact_extract(pred_stack_2020, train_areas_2020,  force_df = TRUE,
#                                         include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "class")
# tDF_2020 = dplyr::bind_rows(tDF_2020)
# tDF_2020 = tDF_2020[  rowSums(is.na(tDF_2020)) == 0,]
# saveRDS(tDF_2020, paste0(envrmt$path_data,"train_areas_2020.rds"))

} else {
  train_areas_2019 = readRDS(paste0(envrmt$path_data,"train_areas_2019.rds"))
#  train_areas_2020 = readRDS(paste0(envrmt$path_data,"train_areas_2020.rds"))

}

## ---- Überwachte  mit Klassifikation Random Forest ----

## Hier wird der Entscheidungsbaum Algorithmus Random Forest über das Utility Paket caret aufgerufen
## ACHTUNG Das Beispiel ist nur für 2019

# das setzen eines seed ermöglicht reproduzierbaren zufall
set.seed(123)

# Aufsplitten  der Daten in training und test , zufällige extraktion von 25% der Daten
trainDat_2019 =  train_areas_2019[createDataPartition(train_areas_2019$class,list = FALSE,p = 0.25),]
# die response variable muss auf den Datentyp "factor" gesetzt werden
trainDat_2019$class <- as.factor(trainDat_2019$class)
# Training Steuerung mit  cross-validation, 10 wiederholungen
ctrlh = trainControl(method = "cv",
                     number = 10,
                     savePredictions = TRUE)
#--- random forest model training

# Hiermit wird das Modell berechnet "trainiert"
cv_model_2019 = train(trainDat_2019[,2:11], # in den Spalten 2 bis 20 stehen die Trainingsdaten (Prediktoren genannt)
                 trainDat_2019[,1],         # in der Spalte 1 stet die zu Klassizierende Variable (Response genannt)
                 method = "rf",             # Methode hier rf für random forest
                 metric = "Kappa",          # Qualitäts/Performanzmaß KAppa
                 trControl = ctrlh,         # obig erzeugte Trainingssteuerung soll eingelsen werden
                 importance = TRUE)         # Die Bedeung der Variablen wird mit abgespeichert


# Klassifikation wird häufig auch Vorhersage genannt.
prediction_rf_2019  = predict(pred_stack_2019 ,cv_model_2019, progress = "text")
mapview(prediction_rf_2019,col.regions = mapviewPalette("mapviewSpectralColors"), at = seq(0, 2, 1), legend = TRUE,alpha.regions = 0.5)



