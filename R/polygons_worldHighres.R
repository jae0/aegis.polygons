
polygons_worldHighres = function( regions="Canada", xlim=c(-80,-40), ylim=c(38, 60) ) {

  require(mapdata)
  require(map)
  require(sf)

  o = sf::st_as_sf(maps::map("world", regions=regions, plot = FALSE, fill = TRUE))
  st_crs(o) = st_crs( "epsg:4326" )

  xy = as.matrix( expand.grid( xlim, ylim) )
  mp = st_multipoint( xy )
  mp = st_bbox (mp )
  bbox =  st_as_sfc( st_bbox( mp ) )
  st_crs(bbox) = st_crs( "epsg:4326" )
  o = st_difference( bbox, o )

  return(o)
}