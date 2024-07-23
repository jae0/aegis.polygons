
areal_units_filename = function( 
    p = NULL,
    spatial_domain=NULL,
    areal_units_type=NULL,
    areal_units_resolution_km=NULL,
    areal_units_constraint="none",
    areal_units_constraint_ntarget=NULL,
    areal_units_constraint_nmin=NULL,
    areal_units_timeperiod=NULL,
    areal_units_overlay=NULL ) {
    
    # though simple, make a function to have only one unique method of construction
    if (!is.null(p)) {
      # assume all params are present in p
      areal_units_fn = paste(
        p$spatial_domain,
        p$areal_units_type,
        p$areal_units_resolution_km,
        p$areal_units_constraint,
        p$areal_units_constraint_ntarget,
        p$areal_units_constraint_nmin,
        p$areal_units_timeperiod,
        paste0(p$areal_units_overlay, collapse="~"),
        sep="~"
      )

    } else {
      areal_units_fn = paste(
        spatial_domain,
        areal_units_type,
        areal_units_resolution_km,
        areal_units_constraint,
        areal_units_constraint_ntarget,
        areal_units_constraint_nmin,
        areal_units_timeperiod,
        paste0(areal_units_overlay, collapse="~"),
        sep="~"
      )
    }
    return( areal_units_fn )
}
