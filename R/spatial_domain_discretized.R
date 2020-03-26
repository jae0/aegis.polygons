
spatial_domain_discretized = function( spatial_domain ) {
  # as points/grids
  pn = spatial_parameters( spatial_domain=spatial_domain )
  Z = bathymetry.db( p=pn, DS="aggregated_data" )
  names(Z)[which(names(Z)=="z.mean" )] = "z"
  Z = lonlat2planar(Z, pn$aegis_proj4string_planar_km)  # should not be required but to make sure
  Z = geo_subset( spatial_domain=spatial_domain, Z=Z )
  spdf0 = SpatialPointsDataFrame(Z[, c("plon", "plat")], data=Z, proj4string=sp::CRS(pn$aegis_proj4string_planar_km) )
  return(spdf0)
}
