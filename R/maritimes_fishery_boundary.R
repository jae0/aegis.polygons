maritimes_fishery_boundary = function( DS="maritimes", areal_units_timeperiod="pre2014", internal_resolution_km=0.1, crs_km="+proj=utm +ellps=WGS84 +zone=20 +units=km", redo=FALSE ) {

  polydir = project.datadirectory("aegis", "polygons")

  if (DS == "maritimes") {
    outdir = file.path(  )
    fn = file.path( polydir, "Science", "Management_Areas", "Fisheries",
      paste( "maritimes_boundary", "rdz", sep="." ) )
  }

  if (DS == "snowcrab" ) {
    fn = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Snowcrab",
      paste( "maritimes_snowcrab_boundary", "rdz", sep="." ) )
  }

  if (DS=="groundfish") {
    fn = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Groundfish",
    paste( "maritimes_groundfish_boundary", "rdz", sep="." ) )
  }

  boundary = NULL
  if (!redo) {
    if (file.exists(fn)) {
      boundary = read_write_fast(fn)
      return(boundary)
    }
  }

  if (DS %in% c("maritimes") ) {
    boundary = (
      polygon_managementareas( species="maritimes" )
      %>% st_transform( st_crs(crs_km) )
      %>% st_simplify()
      %>% st_buffer(0.1)
      %>% st_union()
      %>% st_simplify()
    )
  }

  if (DS=="snowcrab") {
    boundary = (
      polygon_managementareas( species="snowcrab" )
      %>% st_transform( st_crs(crs_km) )
      %>% st_simplify()
      %>% st_buffer(0.1)
      %>% st_union()
      %>% st_simplify()
    )
  }

  if (DS=="groundfish") {
    boundary = maritimes_groundfish_strata( areal_units_timeperiod=areal_units_timeperiod )
    internal_resolution_km = 1 # anything smaler gives incomplete polygons
  }

  boundary = (
    st_transform(boundary, st_crs(crs_km) )
    %>% st_buffer( internal_resolution_km  )
    %>% st_union()
    %>% st_simplify()
    %>% st_make_valid()
    %>% st_transform( st_crs(projection_proj4string("lonlat_wgs84")) )
  )

  if (!file.exists( dirname(fn) )) dir.create( dirname(fn), recursive=TRUE, showWarnings=FALSE )
  read_write_fast( boundary, file=fn )
  return( boundary)
}
