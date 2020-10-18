maritimes_fishery_boundary = function( DS="maritimes", areal_units_timeperiod="pre2014", internal_resolution_km=0.1, crs_km="", redo=FALSE ) {

  polydir = project.datadirectory("aegis", "polygons")

  if (DS == "maritimes") {
    outdir = file.path(  )
    fn = file.path( polydir, "Science", "Management_Areas", "Fisheries",
      paste( "maritimes_boundary", "rdata", sep="." ) )
  }

  if (DS == "snowcrab" ) {
    fn = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Snowcrab",
      paste( "maritimes_snowcrab_boundary", "rdata", sep="." ) )
  }

  if (DS=="groundfish") {
    fn = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Groundfish",
    paste( "maritimes_groundfish_boundary", "rdata", sep="." ) )
  }

  boundary = NULL
  if (!redo) {
    if (file.exists(fn)) {
      load(fn)
      return(boundary)
    }
  }

  if (DS %in% c("maritimes", "snowcrab") ) {
    boundary = (
      as( polygon_managementareas( species="maritimes" ) , "sf" )
      %>% st_transform( st_crs(sppoly) )
      %>% st_simplify()
      %>% st_buffer(0.1)
      %>% st_union()
      %>% st_simplify()
    )
  }

  if (DS=="groundfish") {
    if ( areal_units_timeperiod=="default") areal_units_timeperiod = "pre2014"
    boundary = as( maritimes_groundfish_strata( areal_units_timeperiod=areal_units_timeperiod, returntype="polygons" ), "sf")
  }

  boundary = (
    st_transform(boundary, crs_km )
    %>% st_buffer( internal_resolution_km  )
    %>% st_union()
    %>% st_simplify()
    %>% st_make_valid()
    %>% st_transform( st_crs(projection_proj4string("lonlat_wgs84")) )
  )

  if (!file.exists( dirname(fn) )) dir.create( dirname(fn), recursive=TRUE, showWarnings=FALSE )
  save( boundary, file=fn, compress=TRUE )
  return( boundary)
}
