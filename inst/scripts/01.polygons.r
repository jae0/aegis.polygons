
# Roll your own polygon via mouse interaction

# Example

  # load helper functions
  project.library( "aegis" )
  Rlibrary( "rgdal", "sp", "raster", "maps", "mapdata")

  p = list() # start parameter list

  # define "extents" or "bounding box" in lon/lat
  p$xlim = c(-72, -56 ) # longitudes
  p$ylim = c(42, 49) # latitudes
  p$regions = c("Canada", "USA") # library "map" polygon designations
  p$output.directory = project.datadirectory("aegis", "data","Science")
  p$output.filename = "test.dat"   # must have an "extention" like dat or txt etc,

  # plot background map and using mouse interaction define region:
  # left mouse click to register, right mouse click to finish (or [Esc] is using Rstudio)
  polygon.db( DS="create", p=p )


  # Access method 1: low level
  # read the data with read.table or read.csv
  polygon.db( DS="map.background", p=p )
  scotianshelf = read.table( aegis::polygon_file( "test"  ) )
  lines( scotianshelf, col="green" )  # plot to confirm we have the right data

  # Access method 2: medium
  polygon.db ( DS="map.background", p=p )
  scotianshelf = polygon.db( id="test" )
  lines( scotianshelf, col="orange" )

  # Access method 3: one step .. just the data (projected or not)
  scotianshelf = polygon.db( id="test", crs="+proj=utm +ellps=WGS84 +zone=20 +units=km" )

  # Access method 4: one step data and plot it too
  scotianshelf = polygon.db( id="test", crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", p=p, plotmap=TRUE ) # p contains the extent

