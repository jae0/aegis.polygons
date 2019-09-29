
polygons_managementarea = function( species="snowcrab", area="cfall", redo=FALSE, crs="+proj=utm +ellps=WGS84 +zone=20 +units=km" ) {

  if (species=="snowcrab") {
    polydir = project.datadirectory("aegis", "polygons")
    fn = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Snowcrab", paste( area, "rdata", sep="." ) )
    shp = NULL

    if (!redo) {
      if (file.exists(fn)) load(fn)
      if( !is.null(shp) ) return(shp)
    }

    w = list()
    for (su in aegis.polygons::polygon_internal_code(area)) {
      v = aegis.polygons::polygon.db( polyid=su  )
      vp = Polygon( as.matrix( v ) )
      w = c( w, list( Polygons(list(vp), ID=su )) )
    }

    wsp = SpatialPolygons( w, proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") )  # expect lon/lat
    wsp = spTransform( wsp, sp::CRS(crs) )
    wsp = as(wsp, "sf")
    wsp = st_simplify(wsp)
    wsp = st_buffer(wsp, 0.1)

    coast = coastline.db( p=p, DS="eastcoast_gadm" )
    coast = spTransform( coast, sp::CRS(crs) )
    coast = as( coast, "sf")
    coast = st_simplify(coast)
    coast = st_buffer(coast, 0.1)

    shp = st_difference( st_buffer( st_union( wsp), 0.2), st_union(coast) )
    shp = spTransform( as( shp, "Spatial" ), sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") )
    save(shp, file=fn, compress=TRUE)
    return(shp)

    plot(shp)
  }

  ## add more as required

}
