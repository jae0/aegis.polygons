
areal_units_overlay_filter = function( sppoly, areal_units_overlay, inputdata_spatial_discretization_planar_km, areal_units_resolution_km, areal_units_proj4string_planar_km, areal_units_timeperiod="pre2014" ) {


  if ( grepl("groundfish_strata", areal_units_overlay) ) {

    ### next section is the same as "stratanal polygons" below .. copied here to avoid recursion
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers

    gf = maritimes_groundfish_strata( areal_units_timeperiod=areal_units_timeperiod )
    gf = st_make_valid(gf)
    gf = st_transform(gf, st_crs(sppoly) )
    gf$gfUID = as.character(gf$AUID)
    gf$AUID = NULL
    row.names(gf) = gf$gfUID


    boundary = maritimes_groundfish_boundary( areal_units_timeperiod=areal_units_timeperiod, internal_resolution_km=inputdata_spatial_discretization_planar_km, crs_km=st_crs(sppoly) )
    boundary = st_transform(boundary, st_crs(sppoly) )

    sp0 = sppoly
    sp0$AUID = as.character(1:length(sp0))

    sppoly = (
      st_intersection( st_intersection( boundary, sppoly ), gf )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )
    boundary = NULL

    sppoly = (
      st_join( st_sf(sppoly), gf, join=st_intersects,  largest=TRUE )
    )
    gf = NULL

    # recover the data from sp0
    sppoly = (
      st_join( st_sf(sppoly), sp0, join=st_intersects,  largest=TRUE )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )
    sp0 = NULL
    sppoly = st_make_valid(sppoly)
    sppoly$index = sppoly$AUID
    sppoly$AUID = paste( sppoly$AUID , sppoly$gfUID, sep="_")

  }

  # --------------------
  if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {
    # main domain/boundary .. lactual management areas applied elsewhere (later after filtering to compute SA's)
    sppoly = (
      st_intersection(
        sppoly,
        st_transform( maritimes_fishery_boundary( DS="maritimes", internal_resolution_km=areal_units_resolution_km / 10, crs_km=st_crs(areal_units_proj4string_planar_km) ), st_crs(sppoly) ) )
      %>% st_simplify()
    )
  }
  return(sppoly)

}
