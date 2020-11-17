
polygon_db = function( DS="load", p=NULL, polyid=NULL, project_to=projection_proj4string("lonlat_wgs84"), plotmap=FALSE ) {
  #\\ create/extract polygons and/or return on a map
  #\\ if project_to is passed, default storage/load CRS is assumed lonlat
  #\\ default return value is lon/lat in data frame, also possible to return as a polygon

  if (DS=="map.background") {
    # load libraries to quickly map coastlines
    maps::map( database="worldHires", regions=p$regions, xlim=p$xlim, ylim=p$ylim, fill=FALSE, plot=TRUE )
    maps::map.scale()
    box()
  }

  if (DS=="load") {
    fn = NULL
    fn = try( aegis.polygons::polygon_file( polyid ) )
    if ( class(fn) %in% "try-error") {
      print( "Something went wrong. See error message below:" )
      print( fn)
    }
    X = read.table (fn)
    colnames( X ) = c("lon", "lat" )
    X = (
      as.matrix( X )
      %>% st_multipoint()
      %>% st_sfc( crs=st_crs(project_to) )
      %>% st_cast("POLYGON" )
      %>% st_make_valid()
    )
    return( X )
  }

  if (DS=="create") {
    print( "Left mouse click and then to end right-mouse click (or [Esc] in Windows ) " )
    polygon_db ( DS="map.background", p=p )
    X = locator(type="o" )
    X = as.data.frame( X)
    colnames(X) = c("lon", "lat")
    X = rbind(X, X[1,])
    lines (X)
    u = readline("If it looks good, type 'good' and [enter], otherwise interrrupt and start over.. " )
    save.filename = file.path( p$output.directory, p$output.filename )
    print( save.filename )
    if (file.exists(save.filename) ) {
      u = readline( "Filename (above) exists. Interrupt now otherwise it will be overwritten" )
    }
    write.table( X, file=save.filename )
    return ("save completed")
  }

}
