
polygon_managementareas = function( species="maritimes", area="cfaall", redo=FALSE, project_to=projection_proj4string("lonlat_wgs84"), project_km="+proj=utm +ellps=WGS84 +zone=20 +units=km" ) {

  if (species %in% c("maritimes" ) ) {
    polydir = project.datadirectory("aegis", "polygons")
    outdir = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Snowcrab" )
    if (!file.exists( outdir )) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )
    fn = file.path( outdir, paste( area, "rdz", sep="." ) )

    shp = NULL

    if (!redo) {
      if (file.exists(fn)) shp = read_write_fast(fn)
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

    read_write_fast(shp, file=fn)
    return(shp)

    plot(shp)
  }


  if (species %in% c( "snowcrab") ) {
    polydir = project.datadirectory("aegis", "polygons")
    outdir = file.path( polydir, "Science", "Management_Areas", "Fisheries", "Snowcrab" )
    if (!file.exists( outdir )) dir.create( outdir, recursive=TRUE, showWarnings=FALSE )
    fn = file.path( outdir, paste( "snowcrab", "rdz", sep="." ) )

    shp = NULL

    if (!redo) {
      if (file.exists(fn)) shp = read_write_fast(fn)
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

    pB = bathymetry_parameters( spatial_domain="snowcrab", project_class="stmv" )
    Z = bathymetry_db( p=pB, DS="complete" )
    Z = Z[ geo_subset( spatial_domain="snowcrab", Z=Z ), ]
    Z = st_as_sf( Z[, c("plon", "plat", "z")], coords= c("plon", "plat"), crs=st_crs(shp)  )
    inside = st_points_in_polygons( Z, shp, method="sp::point.in.polygon" )
    Z = Z[which(inside), ]

    ZG = st_as_sf( st_make_grid( Z, cellsize=1, what="polygons", square=TRUE ) )
    ZG = st_join(ZG, Z )
    ZG = ZG [ which(ZG$z > 5 & ZG$z < 350), ]
    
    shp = st_union(ZG)
    shp = st_simplify( shp)
    shp = st_make_valid( shp)

    shp = st_transform( shp, st_crs(project_to) )

    read_write_fast(shp, file=fn)
    return(shp)

    plot(shp)
  }
    


  ## add more as required

}
