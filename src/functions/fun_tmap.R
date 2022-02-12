# Interaktive Visualisierung
tm_plot = function(x,sel=lc_am_compare_sel,vt="plot",ls=FALSE,lo=FALSE,ov=F){
  tmap_mode(vt)
  if (ov){
    tm_shape(sel)+
      tm_polygons("dist",alpha = 0.3,group = "sampling areas",
                  legend.show=T,convert2density=F,legend.is.portrait = F,
                  breaks = c(0, 0.001, 0.01, 0.1, 1.01),
                  legend.is.portrait = F,
                  palette = "viridis",
                  title =  "Distance (JSD)") +
      tm_layout(legend.position = c(.0, 0.15),
                legend.outside.position = "bottom",
                legend.outside = T,
                panel.label.height=0.6,
                panel.label.size=0.6,
                panel.labels = c("Aggregated Areas"))
  } else {
    tm_shape(x) +
      tm_raster(style = "cat",
                title = "Land cover:",
                legend.is.portrait = F) +
      tm_facets(ncol = 2) +
      tm_layout(legend.show = ls,
                legend.only = lo,
                legend.position = c(.0, 0.15),
                legend.outside.position = "bottom",
                legend.outside = T,
                panel.label.height=0.6,
                panel.label.size=0.6,
                panel.labels = c(2019, 2020))}
}

