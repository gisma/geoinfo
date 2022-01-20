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
if (get_sen)
  out_paths_3 <- sen2r(
    gui = TRUE,
    param_list = "~/edu/geoinfo/data/harz.json",
    tmpdir = envrmt$path_tmp,
  )

#--- Einlesen der Daten aus den Verzeichnissen
# RGB stack der beiden Jahre
rgb_2019 = raster::stack(file.path(envrmt$path_data_lev1,"RGB432B",basename(list.files(file.path(envrmt$path_data_lev1,"RGB432B"),pattern = "20190724"))))
rgb_2020 = raster::stack(file.path(envrmt$path_data_lev1,"RGB432B",basename(list.files(file.path(envrmt$path_data_lev1,"RGB432B"),pattern = "20200730"))))
names(rgb_2019) = c("red","green","blue")
names(rgb_2020) = c("red","green","blue")
# Übergabe des RGB stacks an den Prädiktoren-Stack
pred_stack_2019 = rgb_2019
pred_stack_2020 = rgb_2020

# Stack-Loop über die Daten
for (pat in c("EVI","MSAVI2","NDVI","SAVI")){
  pred_stack_2019 = raster::stack(pred_stack_2019,file.path(envrmt$path_data_lev1,pat,basename(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20190724"))))
  pred_stack_2020 = raster::stack(pred_stack_2020,file.path(envrmt$path_data_lev1,pat,basename(list.files(file.path(envrmt$path_data_lev1,pat),pattern = "20200730"))))
}
# Zuweisen von leserlichen Namen auf die Datenebenen

names(pred_stack_2019) = c("red","green","blue","EVI","CSI","MSAVI2","NDVI","SAVI")
names(pred_stack_2020) = c("red","green","blue","EVI","NDVI","MSAVI2","NDVI","SAVI")

#--- Digitalisierung der Trainingsdaten
# Wald
train_area <- mapview::viewRGB(rgb, r = 1, g = 2, b = 3) %>% mapedit::editMap()
# Hinzufügen der Attribute class (text) und id (integer)
forest <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "forest", id = 1)

# kein wald
train_area <- mapview::viewRGB(rgb, r = 1, g = 2, b = 3) %>% mapedit::editMap()
no_forest <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "no_forest", id = 2)

# Kahlschlag
train_area <- mapview::viewRGB(rgb, r = 4, g = 5, b = 6) %>% mapedit::editMap()
clearcut <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "clearcut", id = 3)

# bind it together to one file
train_areas <- rbind(forest, no_forest, clearcut)

# save results
saveRDS(train_areas, paste0(envrmt$path_data,"train_areas.rds"))


# first we have to project the data into the correct crs
tp = sf::st_transform(train_areas,crs = sf::st_crs(pred_stack))
DF <- raster::extract(stack, tp, df=TRUE)
tDF = exactextractr::exact_extract(pred_stack, train_areas,  force_df = TRUE,
                                   include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "class")
tDF = dplyr::bind_rows(tDF)

# brute force approach to get rid of NA
tDF = tDF[  rowSums(is.na(tDF)) == 0,]

# # conversion to a spatial sf object
# tDF_ex_sf = sf::st_as_sf(tDF_ex ,
#                          coords = c("x", "y"),
#                          crs = projection(predStack),
#                          agr = "constant")



## k-means über RStoolbox
prediction_kmeans = unsuperClass(pred_stack, nSamples = 25000, nClasses = 3, nStarts = 25,
                                 nIter = 250, norm = TRUE, clusterMap = TRUE,
                                 Algorithmus = "MacQueen")
mapview(prediction_kmeans$map, col = c('darkgreen', 'burlywood', 'green'))



## random forest via caret
set.seed(123)
# split data into train and test data and take only a fraction of them
trainDat =  tDF[createDataPartition(tDF$class,list = FALSE,p = 0.25),]
# define a training control object for caret with cross-validation, 10 repeats
ctrlh = trainControl(method = "cv",
                     number = 10,
                     savePredictions = TRUE)
# train random forest via caret model
cl = parallel::makeCluster(24)
doParallel::registerDoParallel(cl)
set.seed(seed)
cv_model = train(trainDat[,2:20],
                 trainDat[,1],
                 method = "rf",
                 metric = "Kappa",
                 trControl = ctrlh,
                 importance = TRUE)
stopCluster(cl)

prediction_rf  = predict(pred_stack ,cv_model, progress = "text")
mapview(prediction_rf,col.regions = c('darkgreen', 'burlywood', 'green'))



