

areal_units = function( p=NULL, areal_units_fn_full=NULL, plotit=FALSE, sa_threshold_km2=0, redo=FALSE, use_stmv_solution=FALSE, rastermethod="sf",  xydata=NULL, constraintdata=NULL, spbuffer=5, hull_multiplier = 6, ... ) {

  if (0) {
    plotit=FALSE
    sa_threshold_km2=0
    redo=FALSE
    use_stmv_solution=FALSE
    spbuffer=5
    hull_multiplier = 6
  }

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args


  # areal_units_type:  "tessilation", "lattice", "stratanal_polygons", "groundfish_strata",  "inla_mesh" 
  # areal_units_overlay: "groundfish_strata", "snowcrab_managementareas", "none"

  areal_units_resolution_km =  ifelse (exists("areal_units_resolution_km", p), p$areal_units_resolution_km, 25 )
  inputdata_spatial_discretization_planar_km = ifelse (exists("inputdata_spatial_discretization_planar_km", p), p$inputdata_spatial_discretization_planar_km, min(1, areal_units_resolution_km) )
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
    paste0(areal_units_overlay, collapse="~"),
    sep="|"
  )
  
  areal_units_directory = project.datadirectory("aegis", "polygons", "areal_units" )
  areal_units_fn_full = file.path( areal_units_directory, paste(areal_units_fn, "rdata", sep="." ) )
    

  sppoly = NULL
  boundary = NULL

  if (!redo) {
    if ( file.exists(areal_units_fn_full) ) load(areal_units_fn_full)
    if ( !is.null(sppoly) ) return(sppoly)
  }

  message( "Creating areal units:", areal_units_fn_full)

  if ( areal_units_type == "lattice" ) {
    sppoly = aegis.polygons::areal_units_lattice(
      spatial_domain = spatial_domain, 
      areal_units_resolution_km=areal_units_resolution_km, 
      areal_units_proj4string_planar_km=areal_units_proj4string_planar_km ,
      rastermethod = rastermethod,
      use_stmv_solution=FALSE
    )
  }
   
  # ------------------------------------------------

  if (areal_units_type %in% c( "stratanal_polygons_pre2014", "stratanal_polygons", "groundfish_strata")  ) {
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    areal_units_timeperiod = "pre2014"
    sppoly = maritimes_groundfish_strata( areal_units_timeperiod = areal_units_timeperiod  )
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
    sppoly$AUID = as.character(sppoly$AUID)
    row.names(sppoly) = sppoly$AUID
  }


  if (areal_units_type %in% c( "stratanal_polygons_post2014")  ) {
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    areal_units_timeperiod = "post2014"
    sppoly = maritimes_groundfish_strata( areal_units_timeperiod = areal_units_timeperiod  )
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
    sppoly$AUID = as.character(sppoly$AUID)
    row.names(sppoly) = sppoly$AUID
  }



  # ------------------------------------------------
  if (is.null(xydata)) {
    # if (project_name == "aegis") {
    #   xydata = survey_db( DS="filter", p=p )
    #   xydata = lonlat2planar(xydata, areal_units_proj4string_planar_km)  # should not be required but to make sure
    #   xydata = xydata[ geo_subset( spatial_domain=spatial_domain, Z=xydata ), ]
    #   xydata$AUID = xydata$id
    #   xydata = st_as_sf ( xydata, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
    #   xydata = st_transform( xydata, st_crs( areal_units_proj4string_planar_km ))
    # }

    if (project_name == "temperature") {
      xydata = temperature.db( p=p, DS="bottom.all"  )  #
      xydata = xydata[ , c("lon", "lat", "yr" )]
      xydata = lonlat2planar(xydata, areal_units_proj4string_planar_km)  # should not be required but to make sure
      xydata = st_as_sf ( xydata, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
      xydata = st_transform( xydata, st_crs( areal_units_proj4string_planar_km ))

      locs = st_coordinates( xydata )
      locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

      boundary = st_sfc( st_multipoint( non_convex_hull( locs, alpha=spbuffer*hull_multiplier  ) ), crs=st_crs(areal_units_proj4string_planar_km) )

      # aegis_mesh tweaks
      tus="yr"
      fraction_cv = 0.5 
      fraction_good_bad = 0.9 
      nAU_min = 5  

      areal_units_timeperiod = "none"
  
    }
  


    if (project_name == "survey") {
      xydata = survey_db( DS="set.base", p=p )
      xydata = lonlat2planar(xydata, areal_units_proj4string_planar_km)  # should not be required but to make sure
      xydata$AUID = xydata$id
      xydata = st_as_sf ( xydata, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
      xydata = st_transform( xydata, st_crs( areal_units_proj4string_planar_km ))
      
      locs = st_coordinates( xydata )
      locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

      boundary = 
        maritimes_fishery_boundary( DS="groundfish", internal_resolution_km=1, crs_km=st_crs(areal_units_proj4string_planar_km) ) # post 2014 is larger

      # aegis_mesh defaults
      tus="yr"
      fraction_cv=1.0
      fraction_good_bad=0.8
      nAU_min=5

      areal_units_timeperiod = "pre2014"

    }


    if (project_name == "bio.snowcrab") {
      xydata = snowcrab.db( p=p, DS="set.clean"  )  #
      xydata = xydata[ , c("lon", "lat", "yr"  )]
      xydata = lonlat2planar(xydata, areal_units_proj4string_planar_km)  # should not be required but to make sure
      xydata = st_as_sf ( xydata, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
      xydata = st_transform( xydata, st_crs( areal_units_proj4string_planar_km ))

      locs = st_coordinates( xydata )
      locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

      boundary = st_sfc( st_multipoint( non_convex_hull( locs, alpha=spbuffer*hull_multiplier  ) ), crs=st_crs(areal_units_proj4string_planar_km) )

      # aegis_mesh tweaks
      tus="yr"
      fraction_cv = 0.5 
      fraction_good_bad = 0.9 
      nAU_min = 5  

      areal_units_timeperiod = "none"
  
    }
  

    # if (project_name == "snowcrab_biological_data") {
    #   xydata = snowcrab.db( p=p, DS="biological_data"  )  #
    #   xydata = xydata[ , c("lon", "lat", "yr" )]
    #   xydata = lonlat2planar(xydata, areal_units_proj4string_planar_km)  # should not be required but to make sure
    #   xydata = xydata[ geo_subset( spatial_domain=spatial_domain, Z=xydata ), ]
    #   xydata = st_as_sf ( xydata, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
    #   xydata = st_transform( xydata, st_crs( areal_units_proj4string_planar_km ))
    # }

  }  ## end xydata

  if (!is.null(boundary)) {
      boundary = (
        boundary
        %>% st_cast("POLYGON" )
        %>% st_make_valid()
        %>% st_buffer( areal_units_resolution_km )
        %>% st_union()
        %>% st_cast("POLYGON" )
        %>% st_make_valid()
      )
  }


    if ( areal_units_type == "inla_mesh" ) {
      sppoly = aegis.polygons::areal_units_inla_mesh(
        locs=locs, 
        areal_units_resolution_km=areal_units_resolution_km, 
        areal_units_proj4string_planar_km=areal_units_proj4string_planar_km 
      )
    }


    if ( areal_units_type == "tesselation" ) {
      sppoly = aegis.polygons::aegis_mesh( 
        pts=xydata, 
        boundary=boundary, 
        resolution=areal_units_resolution_km, 
        spbuffer=areal_units_resolution_km, 
        areal_units_constraint_nmin=areal_units_constraint_nmin, 
        tus=tus,
        fraction_cv = fraction_cv, 
        fraction_good_bad = fraction_good_bad, 
        nAU_min = nAU_min  
      )  # voroni tesslation and delaunay triagulation
    }
    
      
  

  if (is.null(sppoly)) stop()


  if (!is.null(boundary)) {
    # must be done separately (after all areal)unit_types has been processed)
    sppoly = (
      st_intersection( st_sf(sppoly), st_transform( boundary, st_crs(sppoly )) )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )
    boundary = NULL
  }  

    #  remove.coastline
    require(aegis.coastline)
    coast = (
        coastline_db( p=p, DS="eastcoast_gadm" )
        %>% st_transform( st_crs( areal_units_proj4string_planar_km ))
        %>% st_simplify()
        %>% st_buffer(inputdata_spatial_discretization_planar_km )
        %>% st_union()
    )
    sppoly = st_difference( sppoly, coast)
    coast = NULL


  # --------------------
  # Overlays
    sppoly = areal_units_overlay(
      sppoly = sppoly, 
      areal_units_overlay = areal_units_overlay, 
      areal_units_resolution_km = areal_units_resolution_km, 
      areal_units_proj4string_planar_km = areal_units_proj4string_planar_km, 
      inputdata_spatial_discretization_planar_km = inputdata_spatial_discretization_planar_km ,
      areal_units_timeperiod = areal_units_timeperiod   # only useful for groundfish
    ) 

  # --------------------
  # Constraints

  if (areal_units_constraint == "groundfish")  constraintdata = survey_db( DS="set.base", p=p )[, c("lon", "lat")]  #
  if (areal_units_constraint == "snowcrab")    constraintdata = snowcrab.db( p=p, DS="set.clean" )[, c("lon", "lat")]  #

  if ( sa_threshold_km2  == 0 ) {
    if ( exists("sa_threshold_km2", p)) sa_threshold_km2 = p$sa_threshold_km2
  }

  sppoly = areal_units_constraint_filter( 
    sppoly=sppoly, 
    areal_units_type = areal_units_type,
    areal_units_constraint_nmin = areal_units_constraint_nmin,
    areal_units_proj4string_planar_km=areal_units_proj4string_planar_km,
    sa_threshold_km2 = sa_threshold_km2, 
    constraintdata = constraintdata 
  ) 


  if (0) {
    W.nb = poly2nb(sppoly, row.names=sppoly$internal_id, queen=TRUE)  # slow .. ~1hr?
    plot(W.nb, st_geometry(sppoly))
    edit.nb(W.nb, polys=as(sppoly, "Spatial"), use_region.id=TRUE)
  }

  # --------------------
  # completed mostly, final filters where required

  # ------------------------------------------------
  # overalys of other areal units not directly related to AU's ( eg. large management areas /zones
  # ... but they require SA estimates after any filters (above))
    if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {
      # as a last pass, calculate surface areas of each subregion .. could be done earlier but it is safer done here due to deletions above
      message("Computing surface areas for each subarea ... can be slow if this is your first time")
      for (subarea in c("cfanorth", "cfasouth", "cfa23", "cfa24", "cfa4x" ) ) {
        print(subarea)
        csa = polygon_managementareas( species="snowcrab", area=subarea )
        csa = st_transform( csa, st_crs( areal_units_proj4string_planar_km ) )
        ooo = st_intersection( csa, sppoly )
        ooo$surfacearea = st_area( ooo )
        vn = paste(subarea, "surfacearea", sep="_")
        sppoly[[ vn ]] = 0
        j = match( ooo$AUID, sppoly$AUID )
        if (length(j) > 0)  sppoly[[ vn ]][j] = ooo$surfacearea
      }
    }


    # ------------------------------------------------


  sppoly = st_make_valid(sppoly)
  row.names(sppoly) = sppoly$AUID
  W.nb = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE, snap=areal_units_resolution_km )  # slow .. ~1hr?
  W.remove = which(card(W.nb) == 0)


  if (0) {
    # https://cran.r-project.org/web/packages/spdep/vignettes/nb_sf.html
     # sf::st_relate takes about 136 s. for a total of 139 s. to generate a queen neighbour object. The contiguity neighbour objects using st_queen
    edit(W.nb, polys=sppoly, use_region.id=TRUE)
  }


  if ( length(W.remove) > 0 ) {
    # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
    W.keep = which(card(W.nb) > 0)
    W.nb = nb_remove( W.nb, W.remove )
    sppoly = sppoly[W.keep,]

    row.names(sppoly) = sppoly$AUID
    sppoly = sppoly[order(sppoly$AUID),]
  }

  attr(sppoly, "nb") = W.nb  # adding neighbourhood as an attribute to sppoly
  attr(sppoly, "project_name") = project_name
  attr(sppoly, "spatial_domain") = spatial_domain
  attr(sppoly, "inputdata_spatial_discretization_planar_km") = inputdata_spatial_discretization_planar_km
  attr(sppoly, "areal_units_directory") = dirname(areal_units_fn_full)
  attr(sppoly, "areal_units_fn") = basename(areal_units_fn_full)
  attr(sppoly, "areal_units_fn_full") = areal_units_fn_full
  attr(sppoly, "areal_units_overlay") = areal_units_overlay
  attr(sppoly, "areal_units_resolution_km") = areal_units_resolution_km
  attr(sppoly, "areal_units_type") = areal_units_type
  attr(sppoly, "areal_units_constraint") = areal_units_constraint 
  attr(sppoly, "areal_units_constraint_nmin") = areal_units_constraint_nmin

  save(sppoly, file=areal_units_fn_full, compress=TRUE)

  if (plotit) plot(sppoly)

  return( sppoly )

}
