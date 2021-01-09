
areal_units_inla_mesh = function(locs, areal_units_resolution_km, areal_units_proj4string_planar_km ) {

    require(INLA)
    sppoly = inla.mesh.2d (
      loc=locs,
      max.edge = c( 0.5, 5 ) * areal_units_resolution_km,  #   # max size of a triange (in, out)
      offset = c( 0.1, 1 ) * areal_units_resolution_km , # how much to extend inside and outside of boundary,
      cutoff = c( 0.5, 5 ) * areal_units_resolution_km # min distance allowed between points #,
      # boundary =  inla.mesh.segment(st_coordinates( as(boundary, "sf") )[,c(1,2)])
    )

    # convert to sp*
    sppoly = SpatialPolygonsDataFrame(
      Sr = SpatialPolygons( lapply(
        1:nrow(sppoly$graph$tv),
        function(x) {
          tv = sppoly$graph$tv[x, , drop = TRUE]
          Polygons(list(Polygon(
            sppoly$loc[tv[c(1, 3, 2, 1)], 1:2, drop = FALSE])),
          ID = x)
        }
        ),
      proj4string = sp::CRS( areal_units_proj4string_planar_km )
      ),
      data = as.data.frame(sppoly$graph$tv[, c(1, 3, 2), drop = FALSE]),
      match.ID = FALSE
    )
 
  return(sppoly)

}
