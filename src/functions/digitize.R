# Für die überwachte Klassifikation benötigen wir Trainingsgebiete.
# Sie können Sie wie nachfolgend digitalisieren oder alternativ z.B. QGis verwenden
# Dieses Skript wird in das change_detection_workflow.R Skript eingelesen
# falls digitize = TRUE definiert ist


#--- 2019
# Kahlschlag
# Es wird das Falschfarbenkomposit in originaler Auflösung genutzt (maxpixels =  1693870)
# Bitte beachten Sie dass es in voller Auflösung (1) deutlich länger lädt und
# (2) Vegetation in Rot dargestellt wird. Sollten sie Echtfarben bevorzugen
# nutzen sie die Kanäle r=4, g=5, b=6

# Kahlschläge/Trockenschäden etc.
train_area <- mapview::viewRGB(pred_stack_2019, r = 1, g = 2, b = 3, maxpixels =  1693870) %>% mapedit::editMap()
# Hinzufügen der Attribute class (text) und id (integer)
clearcut <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "clearcut", id = 1)

# other: hier gilt es möglichst gut verteilt übers Bild
# alle nicht zu Kahlschlag  gehörenden Flächen zu erfassen
train_area <- mapview::viewRGB(pred_stack_2019, r = 1, g = 2, b = 3) %>% mapedit::editMap()
other <- train_area$finished$geometry %>% st_sf() %>% mutate(class = "other", id = 2)

# rbind  kopiert die beiden obigen Vektorobjekte in eine Datei
train_areas_2019 <- rbind(clearcut, other)

# Umprojizieren auf die Raster Datei
train_areas_2019 = sf::st_transform(train_areas_2019,crs = sf::st_crs(pred_stack_2019))


# Extraktion der Trainingsdaten für die digitalisierten Flächen
tDF_2019 = exactextractr::exact_extract(pred_stack_2019, train_areas_2019,  force_df = TRUE,
                                        include_cell = TRUE,include_xy = TRUE,full_colnames = TRUE,include_cols = "class")
# auch hier wieder merge in eine Datei
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


# ---- Das gleiche muss für 2020 wiederholt werden
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

