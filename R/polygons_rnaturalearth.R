
polygons_rnaturalearth = function(  countries= c("United States of America", "Canada"), xlim=c(-80,-40), ylim=c(38, 60) ) {

  if (!require("rnaturalearth"))  install.packages("rnaturalearth")
  if (!require("rnaturalearthdata"))  install.packages("rnaturalearthdata")
    if (0) {
      # or install directly
      devtools::install_github("ropenscilabs/rnaturalearth")
      devtools::install_github("ropenscilabs/rnaturalearthdata")
      install.packages("rnaturalearthhires",
            repos = "http://packages.ropensci.org",
            type = "source")
    }

    require("rnaturalearth")
    require("rnaturalearthhires")
    require("rnaturalearthdata")
    require(sf)
    
    if (is.null(countries)) {
      out = ne_countries(type = 'countries', returnclass="sf")
    } else {
      out = ne_states( countries, returnclass="sf" )
    }

    st_crs(out) = st_crs( "epsg:4326" )
    xy = as.matrix( expand.grid( xlim, ylim) )
    mp = st_multipoint( xy )
    mp = st_bbox (mp )
    bbox =  st_as_sfc( st_bbox( mp ) )
    st_crs(bbox) = st_crs( "epsg:4326" )
    out = st_difference( bbox, out )

  return(out)
}


