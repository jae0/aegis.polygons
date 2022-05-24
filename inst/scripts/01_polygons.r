



# ------------------------------------------------
# Create unerlying areal units for modelling of bathymetry, etc using lattice methods via carstm

require(aegis.bathymetry)

for ( areal_units_resolution_km in c(5, 10, 20, 25) ) {
  # for ( spatial_domain in c("snowcrab", "SSE")) {
   for ( spatial_domain in c("snowcrab", "SSE")) {
    areal_units_overlay = "snowcrab_managementareas"
    if ( spatial_domain=="SSE") areal_units_overlay = "none"
    p = bathymetry_parameters(
      project_class = "carstm", # defines which parameter set to load
      inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      spatial_domain = spatial_domain,  # defines spatial area
      areal_units_resolution_km = areal_units_resolution_km, # km dim of lattice
      areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km",  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
      areal_units_type = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_overlay = areal_units_overlay, # additional polygon layers for subsequent analysis such as management area: "snowcrab" or "groundfish"  # for now ..
      libs = RLibrary ( "sp", "rgeos", "INLA", "raster", "aegis",  "aegis.polygons", "aegis.bathymetry", "carstm" )
    )

    sppoly = areal_units( p=p, redo=TRUE )

  }

}



# Roll your own polygon via mouse interaction

# Example

  # load helper functions
  project.library( "aegis", "aegis.polygons" )


  RLibrary( "rgdal", "sp", "raster", "maps", "mapdata")

  p = list() # start parameter list

  # define "extents" or "bounding box" in lon/lat
  p$xlim = c(-72, -56 ) # longitudes
  p$ylim = c(42, 49) # latitudes
  p$regions = c("Canada", "USA") # library "map" polygon designations
  p$output.directory = project.datadirectory("aegis", "data","Science")
  p$output.filename = "test.dat"   # must have an "extention" like dat or txt etc,

  # plot background map and using mouse interaction define region:
  # left mouse click to register, right mouse click to finish (or [Esc] is using Rstudio)
  polygon_db( DS="create", p=p )


  # Access method 1: low level
  # read the data with read.table or read.csv
  polygon_db( DS="map.background", p=p )
  scotianshelf = read.table( aegis.polygons::polygon_file( "test"  ) )
  lines( scotianshelf, col="green" )  # plot to confirm we have the right data

  # Access method 2: medium
  polygon_db ( DS="map.background", p=p )
  scotianshelf = polygon_db( polyid="test" )
  plot( scotianshelf, col="orange" )

  # Access method 3: one step .. just the data (projected or not)
  scotianshelf = polygon_db( polyid="test", project_to="+proj=utm +ellps=WGS84 +zone=20 +units=km" )

  # Access method 4: one step data and plot it too
  scotianshelf = polygon_db( polyid="test", project_to="+proj=utm +ellps=WGS84 +zone=20 +units=km", p=p ) # p contains the extent

    u = maps::map( database="worldHires", regions=p$regions, xlim=p$xlim, ylim=p$ylim, fill=FALSE, plot=FALSE )
    v = data.frame( cbind( u$x, u$y) )
    w = sf::sf_project( pts=as.matrix(v), from=sf::st_crs("EPSG:4326"), to=project_to )
    plot( w, pch=".", col="gray", xlab="Easting", ylab="Northing")
    plot( scotianshelf, col="green", add=TRUE)

