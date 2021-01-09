
areal_units_filename = function(p, areal_units_directory = project.datadirectory("aegis", "polygons", "areal_units" ), ...) {
  #// construct a unique filename from parameters

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args

  # areal_units_type:  "tessilation", "lattice", "stratanal_polygons", "groundfish_strata",  "inla_mesh" 
  # areal_units_overlay: "groundfish_strata", "snowcrab_managementareas", "none"

  areal_units_resolution_km =  ifelse (exists("areal_units_resolution_km", p), p$areal_units_resolution_km, 25 )
  aegis_internal_resolution_km = ifelse (exists("aegis_internal_resolution_km", p), p$aegis_internal_resolution_km, min(1, areal_units_resolution_km) )
  areal_units_proj4string_planar_km =  ifelse (exists("areal_units_proj4string_planar_km", p), p$areal_units_proj4string_planar_km, p$aegis_proj4string_planar_km )

  # these are required:
  project_name =  ifelse (exists("project_name", p), p$project_name, "default" )
  spatial_domain =  ifelse (exists("spatial_domain", p), p$spatial_domain, "SSE" )
  areal_units_type =  ifelse (exists("areal_units_type", p), p$areal_units_type, "tessilation" )
  areal_units_overlay =  ifelse (exists("areal_units_overlay", p), p$areal_units_overlay, "none" )
  areal_units_constraint =  ifelse (exists("areal_units_constraint", p), p$areal_units_constraint, "none" )
  areal_units_constraint_nmin =  ifelse (exists("areal_units_constraint_nmin", p), p$areal_units_constraint_nmin, 0)

  areal_units_fn = paste(
    project_name,
    spatial_domain,
    areal_units_type,
    areal_units_resolution_km,
    areal_units_constraint,
    areal_units_constraint_nmin,
    paste0(areal_units_overlay, collapse="_"),
    sep="_"
  )
  
  areal_units_directory = project.datadirectory("aegis", "polygons", "areal_units" )
  areal_units_fn_full = file.path( areal_units_directory, paste(areal_units_fn, "rdata", sep="." ) )

  return( areal_units_fn_full )
}
