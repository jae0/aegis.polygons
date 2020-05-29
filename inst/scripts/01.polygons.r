



# ------------------------------------------------
# Create unerlying areal units for modelling of bathymetry, etc using lattice methods via carstm

for ( areal_units_resolution_km in c(10, 20, 25) ) {
  # for ( spatial_domain in c("snowcrab", "SSE")) {
   for ( spatial_domain in c("snowcrab", "SSE")) {
    areal_units_overlay = "snowcrab_managementareas"
    if ( spatial_domain=="SSE") areal_units_overlay = "groundfish_strata"
    p = aegis.carstm::bathymetry_carstm(
      DS = "parameters",
      project_class = "carstm", # defines which parameter set to load
      inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
      spatial_domain = spatial_domain,  # defines spatial area
      areal_units_resolution_km = areal_units_resolution_km, # km dim of lattice
      areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
      # areal_units_proj4string_planar_km = "+proj=omerc +lat_0=44.0 +lonc=-63.0 +gamma=0.0 +k=1 +alpha=325 +x_0=0 +y_0=0 +ellps=WGS84 +units=km",  # oblique mercator, centred on Scotian Shelf rotated by 325 degrees
      areal_units_source = "lattice", # "stmv_fields" to use ageis fields instead of carstm fields ... note variables are not the same
      areal_units_overlay = areal_units_overlay, # additional polygon layers for subsequent analysis such as management area: "snowcrab" or "groundfish"  # for now ..
      libs = RLibrary ( "sp", "rgeos", "INLA", "raster", "aegis",  "aegis.polygons", "aegis.bathymetry", "carstm" )
    )

    sppoly = areal_units( p=p, redo=TRUE )

  }

}



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
  scotianshelf = read.table( aegis.polygons::polygon_file( "test"  ) )
  lines( scotianshelf, col="green" )  # plot to confirm we have the right data

  # Access method 2: medium
  polygon.db ( DS="map.background", p=p )
  scotianshelf = polygon.db( polyid="test" )
  lines( scotianshelf, col="orange" )

  # Access method 3: one step .. just the data (projected or not)
  scotianshelf = polygon.db( polyid="test", crs="+proj=utm +ellps=WGS84 +zone=20 +units=km" )

  # Access method 4: one step data and plot it too
  scotianshelf = polygon.db( polyid="test", crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", p=p, plotmap=TRUE ) # p contains the extent
