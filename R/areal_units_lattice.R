
areal_units_lattice = function(spatial_domain, areal_units_resolution_km, areal_units_proj4string_planar_km, use_stmv_solution=TRUE, rastermethod="sf" ) {

     # res based on grids ... rather than arbitrary polygons
    # static features only so far
    if (use_stmv_solution) {
      pB = bathymetry_parameters( spatial_domain=spatial_domain, project_class="stmv" )
      Z = bathymetry_db( p=pB, DS="complete" )

    } else {
      pB = bathymetry_parameters( spatial_domain=spatial_domain, project_class="core"  )  # default is the "best" performing method
      # as points/grids
      Z = bathymetry_db( p=pB, DS="aggregated_data" )
      names(Z)[which(names(Z)=="z.mean" )] = "z"
      # Z = lonlat2planar(Z, aegis_proj4string_planar_km)  # should not be required but to make sure
    }
    Z = Z[ filter_by_spatial_domain( spatial_domain=spatial_domain, Z=Z ), ]
    Z = st_as_sf( Z[, c("plon", "plat", "z")], coords= c("plon", "plat"), crs=st_crs( pB$aegis_proj4string_planar_km ) )
    Z = st_transform(Z, st_crs(areal_units_proj4string_planar_km))

    if (rastermethod=="sf") {
      sppoly = st_as_sf( st_make_grid( Z, cellsize=areal_units_resolution_km,  what="polygons", square=TRUE ) )
      sppoly$AUID = as.character( 1:nrow(sppoly) )  # row index
      row.names(sppoly) = sppoly$AUID
    } else if (rastermethod=="raster") {
      require(raster)
      raster_template = raster::raster(extent(Z), res=areal_units_resolution_km, crs=projection(Z) ) # +1 to increase the area
      sppoly = raster::rasterize( Z, raster_template, field=Z$z )
      sppoly = as(sppoly, "SpatialPixelsDataFrame")
      sppoly = as( as(sppoly, "SpatialPolygonsDataFrame"), "sf")
      raster_template = NULL
    }
    Z = NULL

  return(sppoly)

}
