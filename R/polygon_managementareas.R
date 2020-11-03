
polygon_managementareas = function( species="maritimes", area="cfaall", redo=FALSE, project_to="+proj=utm +ellps=WGS84 +zone=20 +units=km", returntype="sf" ) {

  if (species %in% c("maritimes", "snowcrab") ) {
    polydir = project.datadirectory("aegis", "polygons")
    outdir = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Snowcrab" )
    if (!file.exists( outdir )) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )
    fn = file.path( outdir, paste( area, "rdata", sep="." ) )

    shp = NULL

    if (!redo) {
      if (file.exists(fn)) load(fn)
      if( !is.null(shp) ) {
        if (returntype=="sp") as( shp, "Spatial")
        return(shp)
      }
    }

    w = NULL
    for (su in aegis.polygons::polygon_internal_code(area)) {
      v = st_as_sf( polygon_db( polyid=su, returntype="sf"  ) )
      v[, "ID"] = su
      w = rbind( w, v )
    }

    wsp = st_transform( w, st_crs(project_to) )
    wsp = st_simplify(wsp)
    wsp = st_buffer(wsp, 0.1)

    coast = coastline_db( p=p, DS="eastcoast_gadm" )
    coast = st_transform( coast, st_crs(project_to) )
    coast = st_simplify(coast)
    coast = st_buffer(coast, 0.1)

    shp = st_difference( st_buffer( st_union( wsp), 0.2), st_union(coast) )
    shp = st_transform( shp, st_crs(projection_proj4string("lonlat_wgs84")) )

    save(shp, file=fn, compress=TRUE)

    if (returntype=="sp") as( shp, "Spatial")
    return(shp)

    plot(shp)
  }

  ## add more as required

}
