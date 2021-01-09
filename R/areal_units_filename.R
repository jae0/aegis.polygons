
areal_units_filename = function(p, areal_units_directory = project.datadirectory("aegis", "polygons", "areal_units" ), ...) {
  #// construct a unique filename from parameters

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args

  areal_units_fn = paste(
    p$project_name,
    p$spatial_domain,
    p$areal_units_type,
    p$areal_units_resolution_km,
    p$areal_units_constraint,
    p$areal_units_constraint_nmin,
    paste0(p$areal_units_overlay, collapse="~"),
    sep="|"
  )
  
  areal_units_fn_full = file.path( areal_units_directory, paste(areal_units_fn, "rdata", sep="." ) )

  return( areal_units_fn_full )
}
