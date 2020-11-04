
polygon_managementareas = function( species="maritimes", area="cfaall", redo=FALSE, project_to=projection_proj4string("lonlat_wgs84"), project_km="+proj=utm +ellps=WGS84 +zone=20 +units=km" ) {

  if (species %in% c("maritimes", "snowcrab") ) {
    polydir = project.datadirectory("aegis", "polygons")
    outdir = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Snowcrab" )
    if (!file.exists( outdir )) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )
    fn = file.path( outdir, paste( area, "rdata", sep="." ) )

    shp = NULL

    if (!redo) {
      if (file.exists(fn)) load(fn)
      if( !is.null(shp) ) {
        return(shp)
      }
    }

    w = NULL
    for (su in aegis.polygons::polygon_internal_code(area)) {
      v = st_as_sf( polygon_db( polyid=su ) )
      v[, "ID"] = su
      w = rbind( w, v )
    }

    wsp = st_transform( w, st_crs(project_km) )
    wsp = st_simplify(wsp)
    wsp = st_buffer(wsp, 0.1)

    coast = aegis.coastline::coastline_db()
    coast = st_transform( coast, st_crs(wsp) )
    coast = st_simplify(coast)
    coast = st_buffer(coast, 0.1)

    shp = st_difference( st_buffer( st_union( wsp), 0.2), st_union(coast) )

    shp = st_transform( shp, st_crs(project_to) )

    save(shp, file=fn, compress=TRUE)
    return(shp)

    plot(shp)
  }

  ## add more as required

}
