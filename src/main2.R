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

# laden der notwendigen Bibliotheken
## Achtung sie müssen evtl. installiert werden
library(envimaR)
library(rprojroot)
## setzen des aktuellen Projektverzeichnisses (erstellt mit envimaR) als rootDIR
#rootDIR = find_rstudio_root_file()
root_folder = "~/edu/courses/src_courses/geoinfo/"
library(sen2r)

# einlesen des zuvor erstellten setup Srkiptes
source(file.path(root_folder, "src/functions/000_setup.R"))
tmpDir()=envrmt$path_tmp

#safe_is_online("/home/creu/.sen2r/lta_orders/lta_20211213_163011.json")
json_path <- "~/edu/courses/src_courses/geoinfo/data/harz.json"

out_paths_3 <- sen2r(
  gui = TRUE,
  param_list = json_path,
)



# st = myextent,
# extent_name = "harz",
# timewindow = c(as.Date("2019-06-01"), as.Date("2019-07-01")),
# list_prods stack=raster::stack(paste0(envrmt$path_research_area,"/BOA/",fn))

# subsetting the filename(s) of the interesting file(s)
fn_noext=xfun::sans_ext(basename(list.files(paste0(envrmt$path_data_lev1,"/BOA/"),pattern = "S2B2A")))
fn = basename(list.files(paste0(envrmt$envrmt$path_data_lev1,"/BOA/"),pattern = "S2B2A"))

# creating a raster stack
stack=raster::stack(paste0(envrmt$path_research_area,"/BOA/",fn))


# first we have to project the data into the correct crs
tp = sf::st_transform(train_areas,crs = sf::st_crs(stack))
## next we extract the values from every band of the raster stack
# we force the values to be returned as an data frame
# because extracting the raster way is very slow
DF <- raster::extract(stack, tp, df=TRUE)
# we rename the layers for simplicity
# now we add the "class" category which we need later on for training
# it was dropped during extraction
DF_sf =st_as_sf(inner_join(DF,tp))
# finally we produce a simple data frame without geometry
DF2 = DF_sf
st_geometry(DF2)=NULL

# berechnen der Pberflächen Albedo einer angepassten Regression aus dem Paket 'agriwater'
# Einlesen der notwendigen Knäle aus dem Sentinel Stack  in einzelne raster ebenen
# /10000 ist notwendig um die Originalwerte korrekt zu skalieren
b2 <- stack[[2]]/10000
b3 <- stack[[3]]/10000
b4 <- stack[[4]]/10000
b8 <- stack[[8]]/10000

# lineare regression albedo top of Atmosphere
alb_top = b2 * 0.32 + b3 * 0.26 + b4 * 0.25 + b8 * 0.17
# lineare regression oberflächenalbedo
alb_surface = 0.6054 * alb_top + 0.0797

# Visualierung
plot(alb_top)
plot(alb_surface)
