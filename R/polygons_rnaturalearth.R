
polygons_rnaturalearth = function( ) {

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
  world = ne_countries(type = 'countries', scale = 'large', returnclass="sf")

  return(world)
}