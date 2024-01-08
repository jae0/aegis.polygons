
areal_units_filename = function( 
    spatial_domain,
    areal_units_type,
    areal_units_resolution_km,
    areal_units_constraint,
    areal_units_constraint_ntarget,
    areal_units_constraint_nmin,
    areal_units_timeperiod,
    areal_units_overlay ) {
    
    # though simple, make a function to have only one unique method of construction

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
    return( areal_units_fn )
}
