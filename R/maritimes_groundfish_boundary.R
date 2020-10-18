maritimes_groundfish_boundary = function( areal_units_timeperiod="pre2014", internal_resolution_km=0.1, crs_km="", redo=FALSE ) {

    polydir = project.datadirectory("aegis", "polygons")
    outdir = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Groundfish" )
    if (!file.exists( outdir )) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )
    fn = file.path( outdir, paste( "maritimes_groundfish_boundary", "rdata", sep="." ) )
    boundary = NULL
    if (!redo) {
      if (file.exists(fn)) {
        load(fn)
        return(boundary)
      }
    }
    if ( areal_units_timeperiod=="default") areal_units_timeperiod = "pre2014"

    boundary = (
      as( maritimes_groundfish_strata( areal_units_timeperiod=areal_units_timeperiod, returntype="polygons" ), "sf")
      %>% st_transform(crs_km )
      %>% st_buffer( internal_resolution_km  )
      %>% st_union()
      %>% st_simplify()
      %>% st_make_valid()
      %>% st_transform( st_crs(projection_proj4string("lonlat_wgs84")) )
     )

    save( boundary, file=fn, compress=TRUE )
    return( boundary)

}
