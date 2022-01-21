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
library(rprojroot)
library(sen2r)
#--- Schalter für Download
get_sen = FALSE


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
pred_stack_2019 = raster::stack(list.files(file.path(envrmt$path_data_lev1,"RGB432B"),pattern = "20190619",full.names = TRUE))
pred_stack_2020 = raster::stack(file.path(envrmt$path_data_lev1,"RGB843B",basename(list.files(file.path(envrmt$path_data_lev1,"RGB432B"),pattern = "20200730"))))

# Stack-Loop über die Daten
for (pat in c("RGB432B","EVI","MSAVI2","NDVI","SAVI")){
  pred_stack_2019 = raster::stack(pred_stack_2019,file.path(envrmt$path_data_lev1,pat,basename(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20190619"))))
  #pred_stack_2020 = raster::stack(pred_stack_2020,file.path(envrmt$path_data_lev1,pat,basename(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20200730"))))
}
# Zuweisen von leserlichen Namen auf die Datenebenen

names(pred_stack_2019) = c("nir","red","green","red","green","blue","EVI","MSAVI2","NDVI","SAVI")
names(pred_stack_2020) = c("nir","red","green","red","green","blue","EVI","MSAVI2","NDVI","SAVI")

#--- Digitalisierung der Trainingsdaten
#--- 2019
# Kahlschlag
train_area <- mapview::viewRGB(pred_stack_2019, r = 1, g = 2, b = 3) %>% mapedit::editMap()
# Hinzufügen der Attribute class (text) und id (integer)
clearcut <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "clearcut", id = 1)
# kein Wald
train_area <- mapview::viewRGB(pred_stack_2019, r = 1, g = 2, b = 3) %>% mapedit::editMap()
other <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "other", id = 2)
# merge in eine Datei
train_areas_2019 <- rbind(clearcut, other)
# sichern
saveRDS(train_areas_2019, paste0(envrmt$path_data,"train_areas_2019.rds"))
#--- 2020
# Kahlschlag
train_area <- mapview::viewRGB(pred_stack_2020, r = 1, g = 2, b = 3) %>% mapedit::editMap()
# Hinzufügen der Attribute class (text) und id (integer)
clearcut <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "clearcut", id = 1)
# kein Wald
train_area <- mapview::viewRGB(pred_stack_2020, r = 1, g = 2, b = 3) %>% mapedit::editMap()
other <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "other", id = 2)
# merge zu einer Datei
train_areas_2020 <- rbind(clearcut, other)
# sichern
saveRDS(train_areas_2020, paste0(envrmt$path_data,"train_areas_2020.rds"))

#--- 2019
# Projektion
tp_2019 = sf::st_transform(train_areas_2019,crs = sf::st_crs(pred_stack))
# Extraktion
tDF_2019 = exactextractr::exact_extract(pred_stack_2019, train_areas_2019,  force_df = TRUE,
                                   include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "class")
# merge in eine Datei
tDF_2019 = dplyr::bind_rows(tDF_2019)
# Löschen von NA Zeilen
tDF = tDF[  rowSums(is.na(tDF)) == 0,]

#--- 2020
tp_2020 = sf::st_transform(train_areas_2020,crs = sf::st_crs(pred_stack))
tDF_2020 = exactextractr::exact_extract(pred_stack_2020, train_areas_2020,  force_df = TRUE,
                                        include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "class")
tDF_2020 = dplyr::bind_rows(tDF_2020)
tDF = tDF[  rowSums(is.na(tDF)) == 0,]


## k-means über RStoolbox
# Modell
prediction_kmeans_2019 = unsuperClass(pred_stack_2019, nClasses = 2,norm = TRUE, algorithm = "MacQueen")
# Klassifikation
mapview(prediction_kmeans_2019$map, col = c('darkgreen', 'burlywood', 'green'))
prediction_kmeans_2020 = unsuperClass(pred_stack_2020, nClasses = 2,norm = TRUE, algorithm = "MacQueen")
mapview(prediction_kmeans2020$map, col = c('darkgreen', 'burlywood', 'green'))


## random forest via caret
# seed ermglicht reproduzierbaren zufall
set.seed(123)
# auteilung der daten in training und test , zufällige extraktion von 25% der Daten
trainDat_2019 =  tDF_2019[createDataPartition(tDF_2019$class,list = FALSE,p = 0.25),]
# Training Steuerung mit  cross-validation, 10 wiederholungen
ctrlh = trainControl(method = "cv",
                     number = 10,
                     savePredictions = TRUE)
#--- random forest model training
# Paralleisierung
cl = parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
set.seed(123)
# Modelltraining
cv_model_2019 = train(trainDat_2019[,2:20],
                 trainDat_2019[,1],
                 method = "rf",
                 metric = "Kappa",
                 trControl = ctrlh,
                 importance = TRUE)
stopCluster(cl)

# Klassifikation
prediction_rf_2019  = predict(pred_stack_2019 ,cv_model_2019, progress = "text")
mapview(prediction_rf_2019,col.regions = c('darkgreen', 'burlywood', 'green'))



