# aegis.polygons

This R-library serves three functions:

- A front end to save/load polygons stored as Rdata or RDS files. 

- To create polygons (areal units) for use in areal-unit modelling
    
    - using tesselation of Voronoi triangles (see carstm examples in their own projects)

    - using simple lattice structure defined by "areal_units_resolution_km" (see below)

- To create polygons interactively with mouse clicks.

---

The following snippets show examples of these use cases.

```r

# An example creating underlying areal units for modelling of bathymetry via carstm
# In this case it is a simple lattice

require(raster)  
require(fasterize)  
require(sf)  

require(aegis)
require(aegis.polygons)
require(aegis.bathymetry)


spatial_domain  = "SSE"  # target spatial domain

areal_units_resolution_km = 10 # resoltuion of spatial lattice

# set up required parameters
p = bathymetry_parameters(
  project_class = "carstm", # defines which parameter set to load (areal units)
  inputdata_spatial_discretization_planar_km = 1,  # 1 km .. requires 32 GB RAM and limit of speed -- controls resolution of data prior to modelling to reduce data set and speed up modelling
  spatial_domain = spatial_domain,  # defines spatial area
  areal_units_resolution_km = areal_units_resolution_km, # km dim of lattice
  areal_units_proj4string_planar_km = projection_proj4string("utm20"),  # coord system to use for areal estimation and gridding for carstm
  areal_units_type = "lattice",  
  areal_units_overlay =  "none"  # additional polygon layers for subsequent analysis  
)

sppoly = areal_units( p=p, redo=TRUE )

plot(sppoly["AUID"])

```


To roll your own polygon via mouse interaction is easy:

```r

# load helper functions (if not installed as packages)
loadfunctions( "aegis")
loadfunctions( "aegis.polygons" )

message("FIXE ME::: deprecated libs, use sf/stars")

RLibrary( "raster" )

temp_directory = tempdir()
dir.create( dir.local, recursive=TRUE, showWarnings=FALSE )

p = list() # start parameter list

# define "extents" or "bounding box" in lon/lat
p$xlim = c(-72, -56 ) # longitudes
p$ylim = c(42, 49) # latitudes
p$regions = c("Canada", "USA") # library "map" polygon designations
p$output.directory = temp_directory
p$output.filename = "test.dat"   # must have an "extention" like dat or txt etc,


# plot background map and using mouse interaction define region:
require(GADMTools)

message( "Downloading .. \n" )
message( "Warnings about 'old-style crs object detected ...' ' can be ignored .. the maintainer needs to update their files" )

maritimes = GADMTools::gadm_subset(GADMTools::gadm_sf.loadCountries( fileNames="CAN", level=1, basefile=dir.local   ),
  level=1, regions=c("Nova Scotia", "Prince Edward Island", "Newfoundland and Labrador", "Québec","New Brunswick"  )  )$sf

plot(maritimes["NAME_1"])

# left mouse click to register, right mouse click to finish (or [Esc] is using Rstudio)
res = polygon_db( DS="create", p=p )

str(res)

```

