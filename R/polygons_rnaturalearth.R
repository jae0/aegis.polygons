
polygons_rnaturalearth = function(  
  countries=NULL, 
  xlim=c(-80,-40), 
  ylim=c(38, 60), 
  scale="large", returnclass="sf" ) {


if (0){
  loadfunctions("aegis.polygons")

  # to debug
  countries=NULL
  xlim=c(-80,-40)
  ylim=c(38, 60)
  scale="large"
  returnclass="sf" 

  # to install directly
  devtools::install_github("ropenscilabs/rnaturalearth")
  devtools::install_github("ropenscilabs/rnaturalearthdata")
  install.packages("rnaturalearthhires",
        repos = "http://packages.ropensci.org",
        type = "source")

}

  if (!require("rnaturalearth"))  install.packages("rnaturalearth")
  if (!require("rnaturalearthdata"))  install.packages("rnaturalearthdata")
  if (!require("rnaturalearthhires"))  install.packages("rnaturalearthhires")

  require("rnaturalearth")
  require("rnaturalearthhires")
  require("rnaturalearthdata")
  require(sf)
  
  if (!is.null(countries)) {

     countries = ne_countries(type = "countries", scale=scale, returnclass=returnclass)

     out = countries[ countries$name %in% countries, "name" ] 

  } else {

    ne_download( scale = 110L, type = "states", category ="cultural" )
    countries = ne_countries(type = "countries", scale=scale, returnclass=returnclass)
    
    spmg = countries[ countries$name %in% 
      c(  "St. Pierre and Miquelon", 
          "Greenland" ), 
      "name" ] 

    subareas = ne_states( 
      c( "United States of America", 
        "Canada" 
      ), 
      returnclass=returnclass 
    )
      
    subareas = subareas[ subareas$name %in% c(
        "Québec", 
        "New Brunswick", 
        "Newfoundland and Labrador",
        "Nova Scotia",
        "Prince Edward Island",
        "Maine"
        ) , "name" ]

    out = rbind( subareas, spmg )
  }

  st_crs(out) = st_crs( "epsg:4326" )
  xy = as.matrix( expand.grid( xlim, ylim) )
  mp = st_multipoint( xy )
  mp = st_bbox (mp )
  bbox =  st_as_sfc( st_bbox( mp ) )
  st_crs(bbox) = st_crs( "epsg:4326" )
 
  out = st_intersection( out, bbox )
  
  return(out)
}


