

areal_units = function( p=NULL, areal_units_fn_full=NULL, areal_units_directory=NULL, plotit=FALSE, sa_threshold_km2=0, redo=FALSE,
  use_stmv_solution=TRUE, rastermethod="sf",  xydata=NULL, constraintdata=NULL, spbuffer=5, hull_alpha =15, duplications_action="union",  areal_units_timeperiod=NULL, verbose=FALSE, return_crs=NULL, 
      count_time=TRUE, respect_spatial_domain=TRUE, ... ) {

  if (0) {
    plotit=FALSE
    sa_threshold_km2=0
    redo=FALSE
    use_stmv_solution=TRUE
    spbuffer=5
    hull_alpha = 15
    rastermethod="sf"
    xydata=NULL
    constraintdata=NULL
    duplications_action="union"
    areal_units_timeperiod=NULL
    verbose=TRUE
  }

  require(spdep)

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args

  if (exists("hull_alpha", p)) hull_alpha = p$hull_alpha

  # areal_units_type:  "tesselation", "lattice", "stratanal_polygons", "groundfish_strata",  "inla_mesh"
  # areal_units_overlay: "groundfish_strata", "snowcrab_managementareas", "none"

  areal_units_resolution_km =  ifelse (exists("areal_units_resolution_km", p), p$areal_units_resolution_km, 25 )
  inputdata_spatial_discretization_planar_km = ifelse (exists("inputdata_spatial_discretization_planar_km", p), p$inputdata_spatial_discretization_planar_km, min(1, areal_units_resolution_km) )
  areal_units_proj4string_planar_km =  ifelse (exists("areal_units_proj4string_planar_km", p), p$areal_units_proj4string_planar_km, p$aegis_proj4string_planar_km )

  # these are required:
  project_name =  ifelse (exists("project_name", p), p$project_name, "default" )
  spatial_domain =  ifelse (exists("spatial_domain", p), p$spatial_domain, "SSE" )
  areal_units_type =  ifelse (exists("areal_units_type", p), p$areal_units_type, "tesselation" )
  areal_units_overlay =  ifelse (exists("areal_units_overlay", p), p$areal_units_overlay, "none" )
  areal_units_constraint =  ifelse (exists("areal_units_constraint", p), p$areal_units_constraint, "none" )
  areal_units_timeperiod =  ifelse (exists("areal_units_timeperiod", p), p$areal_units_timeperiod, "none" )
  
  areal_units_constraint_nmin =  ifelse (exists("areal_units_constraint_nmin", p), p$areal_units_constraint_nmin, 0)
  areal_units_constraint_ntarget =  ifelse (exists("areal_units_constraint_ntarget", p), p$areal_units_constraint_ntarget, 0)


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

  if ( !is.null(areal_units_fn_full) ) areal_units_directory =  dirname(areal_units_fn_full)

  if ( is.null(areal_units_fn_full) )  {
    if ( is.null(areal_units_directory) )  areal_units_directory = file.path( p$datadir, "areal_units" )
    areal_units_fn_full = file.path( areal_units_directory, paste(areal_units_fn, "rdata", sep="." ) )
  }

  dir.create( areal_units_directory, showWarnings = FALSE, recursive = TRUE )
  
  sppoly = NULL
  boundary = NULL

  if (!redo) {
    if ( file.exists(areal_units_fn_full) ) {
      load(areal_units_fn_full)
      if (!is.null(sppoly)) {
        if ( !is.null( return_crs )) sppoly=st_transform(sppoly, crs=st_crs(return_crs) )
        if ( !is.null( sppoly ) ) return(sppoly)
      }
    }
  }

  message( "Creating areal units: ", areal_units_fn_full)

  message( "Areal units base structure (type): ", areal_units_type)

  if ( areal_units_type == "lattice" ) {
    sppoly = areal_units_lattice(
      spatial_domain = spatial_domain,
      areal_units_resolution_km=areal_units_resolution_km,
      areal_units_proj4string_planar_km=areal_units_proj4string_planar_km ,
      rastermethod = rastermethod,
      use_stmv_solution=use_stmv_solution
    )
  }

  # ------------------------------------------------

  if (areal_units_type %in% c( "stratanal_polygons_pre2014", "stratanal_polygons", "groundfish_strata")  ) {
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    areal_units_timeperiod = "pre2014"
    sppoly = maritimes_groundfish_strata( areal_units_timeperiod = areal_units_timeperiod  )
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers`
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
    if (exists("areal_units_xydata", p)) {
       assign("xydata", eval(parse(text=p$areal_units_xydata) ) )
    } else {
       message( "To create areal units, xydata is required.")
       stop()
    }
  }
  if (!exists("lon", xydata)) {
    if (! "sf" %in% class(xydata) ) {
      stop( "areal_units requires 'lon', 'lat', and possibly 'yr' or that it be an 'sf' object")
    }
  }

  if (! "sf" %in% class(xydata) ) {
    xydata = st_as_sf ( xydata, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
  }
  xydata = st_transform( xydata, st_crs( areal_units_proj4string_planar_km ))

  if (respect_spatial_domain) {
    xydata = xydata[ geo_subset( spatial_domain=spatial_domain, Z=xydata ) , ]
  }

  if ( areal_units_type %in% c( "stratanal_polygons_pre2014", "stratanal_polygons", "groundfish_strata", "stratanal_polygons_post2014") ) {
    message( "Determining areal unit domain boundary from groundfish survey")
    boundary = maritimes_fishery_boundary( DS="groundfish", internal_resolution_km=1, crs_km=st_crs(areal_units_proj4string_planar_km) ) # post 2014 is larger (crm is for incoming data)
    boundary = st_transform( boundary, st_crs(areal_units_proj4string_planar_km) ) # output of above is lonlat
  
  } else {
    
    if ( project_name == "bio.snowcrab") {
    
      message( "Determining areal unit domain boundary from snowcrab survey")
      boundary = polygon_managementareas( species="snowcrab" )
      boundary = st_transform( boundary, st_crs(areal_units_proj4string_planar_km) )
      boundary = st_buffer(boundary, 0)
      boundary = st_make_valid(boundary)

      data_boundary = st_sfc( st_multipoint( non_convex_hull(
        st_coordinates( xydata ) + runif( nrow(xydata)*2, min=-1e-3, max=1e-3 ) ,  # noise increases complexity of edges -> better discrim of polys
        alpha=hull_alpha
        ) ), crs=st_crs(areal_units_proj4string_planar_km) )
      data_boundary = (
          data_boundary
          %>% st_cast("POLYGON" )
          %>% st_make_valid()
          %>% st_buffer( areal_units_resolution_km/10 )
          %>% st_union()
          %>% st_cast("POLYGON" )
          %>% st_make_valid()
        )
      boundary = st_intersection(data_boundary, boundary)

    }  else {

      message( "Determining areal unit domain boundary from input: xydata")
      message( "If you get strange boundary outlines (due to data sparsity),")
      message(" try changing (increasing) the value of hull_alpha, the granularity of boundary tracing, in km.")
      message( "The current hull_alpha is: ", hull_alpha )
      boundary = st_sfc( st_multipoint( non_convex_hull(
        st_coordinates( xydata ) + runif( nrow(xydata)*2, min=-1e-3, max=1e-3 ) ,  # noise increases complexity of edges -> better discrim of polys
        alpha=hull_alpha
      ) ), crs=st_crs(areal_units_proj4string_planar_km) )
  
    }
  
  }

  if (0){
    plot( xydata, reset=FALSE )
    plot( boundary, add=TRUE )

  }

  if (!is.null(boundary)) {
      boundary = (
        boundary
        %>% st_cast("POLYGON" )
        %>% st_make_valid()
        %>% st_simplify()
        %>% st_buffer( areal_units_resolution_km/10 )
        %>% st_union()
        %>% st_cast("POLYGON" )
        %>% st_make_valid()
      )
  }



    if ( areal_units_type == "inla_mesh" ) {
      message( "Determining areal units as an INLA mesh")
      sppoly = areal_units_inla_mesh(
        locs=st_coordinates( xydata ) + runif( nrow(xydata)*2, min=-1e-3, max=1e-3 ) , # add  noise  to prevent a race condition
        areal_units_resolution_km=areal_units_resolution_km,
        areal_units_proj4string_planar_km=areal_units_proj4string_planar_km
      )
    }

    if ( areal_units_type == "tesselation" ) {
      message( "Determining areal units via iterative Voroni tesselation of AU centroids and dissolution of AUs")
      sppoly = aegis_mesh(
        pts=xydata,
        boundary=boundary,
        resolution=areal_units_resolution_km,
        spbuffer=areal_units_resolution_km,
        areal_units_constraint_ntarget=areal_units_constraint_ntarget,
        areal_units_constraint_nmin=areal_units_constraint_nmin,
        tus=p$tus,
        fraction_todrop = p$fraction_todrop,
        fraction_cv = p$fraction_cv,  # stopping criterion: when cv drops below this value
        nAU_min = p$nAU_min,   # stoppping criterion: allow no less than this number of areal units
        count_time = count_time,
        verbose = verbose
      )  # voroni tesslation and delaunay triagulation
    }

    xydata = NULL


  if (is.null(sppoly)) stop("Error in areal units: none found")


  if (!is.null(boundary)) {
    # must be done separately (after all areal)unit_types has been processed)
    sppoly = (
      st_intersection( st_sf(sppoly), st_transform( boundary, st_crs(sppoly )) )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )
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

    sppoly = st_transform( sppoly, st_crs(areal_units_proj4string_planar_km) )
    sppoly = st_difference( sppoly, coast)
    coast = NULL

 
    # --------------------
    # Overlays
    message( "Filtering areal units on overlays")

    sppoly = areal_units_overlay_filter(
      sppoly = sppoly,
      areal_units_overlay = areal_units_overlay,
      areal_units_resolution_km = areal_units_resolution_km,
      areal_units_proj4string_planar_km = areal_units_proj4string_planar_km,
      inputdata_spatial_discretization_planar_km = inputdata_spatial_discretization_planar_km ,
      areal_units_timeperiod = areal_units_timeperiod   # only useful for groundfish
    )



  if ( !is.null(constraintdata)) {
    # --------------------
    # Additional Constraints from other data

    if (areal_units_constraint == "groundfish")  {
      message( "Constrain areal units to required number of groundfish survey locations by merging into adjacent AUs")
      constraintdata = groundfish_survey_db(DS="set.base", yrs=p$yrs )[, c("lon", "lat")]  #
    }
    
    if (areal_units_constraint == "snowcrab")    {
      message( "Constrain areal units to required number of snowcrab survey locations by merging into adjacent AUs")
      constraintdata = snowcrab.db( p=p, DS="set.clean" )[, c("lon", "lat")]  #
    }

    if (areal_units_constraint == "survey")    {
      message( "Constrain areal units to required number of aegis.survey locations by merging into adjacent AUs")
      constraintdata = survey_db( p=p, DS="set.base" )[, c("lon", "lat")]  #
    }


    sppoly$npts  = 0
    sppoly$internal_id = 1:nrow(sppoly)
    row.names( sppoly) = sppoly$internal_id

    sppoly = st_transform( sppoly, st_crs( areal_units_proj4string_planar_km ))
    sppoly$au_sa_km2 = st_area(sppoly)
    attributes( sppoly$au_sa_km2 ) = NULL

    constraintdata = sf::st_as_sf( constraintdata, coords = c("lon","lat"), crs=st_crs(projection_proj4string("lonlat_wgs84")) )
    constraintdata = st_transform( constraintdata, st_crs(areal_units_proj4string_planar_km ))
    # constraintdata = st_join( constraintdata, sppoly, join=st_within )
    constraintdata$internal_id = st_points_in_polygons( constraintdata, sppoly, varname="internal_id" )
    ww = tapply( rep(1, nrow(constraintdata)), constraintdata$internal_id, sum, na.rm=T )
    sppoly$npts[ match( names(ww), as.character(sppoly$internal_id) )] = ww
 
    zeros = which( sppoly$npts == 0 )
    if ( length(zeros) > 0 ) sppoly = sppoly[-zeros,]
    

    # This section needs to be reworked :: 
    # NOTE:  areal_units_constraint_nmin already operates in the tesselation
    # dropping causes too many data points from being dropped too (esp near boudndaries .,,  must make it respect boundary)

    # todrop = which( sppoly$npts < areal_units_constraint_nmin )

    # if (length(todrop) > 0 ) {

    #   if (areal_units_type == "lattice" ) {
    #     # laatice structure is required, simply drop where there is no data
    #     sppoly = sppoly[ -todrop , ]

    #   } else {
    #     # already done in tesselation method so not really needed but
    #     # try to join to adjacent au's
    #     # for other methods, this ensures a min no samples in each AU
    #     sppoly$dropflag = FALSE
    #     sppoly$nok = TRUE
    #     sppoly$nok[todrop] = FALSE
    #     W.nb = poly2nb(sppoly, row.names=sppoly$internal_id, queen=TRUE)  
    #     for (i in order(sppoly$npts) ) {
    #       if ( sppoly$nok[i]) next()
    #       lnb = W.nb[[ i ]]
    #       if (lnb < 1) next()
    #       local_finished = FALSE
    #       for (f in 1:length(lnb) ) {
    #         if (local_finished) break()
    #         v = setdiff( 
    #           intersect( lnb, which(sppoly$nok) ), ## AU neighbours that are OK and so can consider dropping
    #           which(sppoly$dropflag)                       ## AUs confirmed already to drop
    #         )
    #         if (length(v) > 0) {
    #           j = v[ which.min( sppoly$npts[v] )]
    #           g_ij = try( st_union( st_geometry(sppoly)[j] , st_geometry(sppoly)[i] ) )
    #           if ( !inherits(g_ij, "try-error" )) {
    #             st_geometry(sppoly)[j] = g_ij
    #             sppoly$npts[j] = sppoly$npts[j] + sppoly$npts[i]
    #             sppoly$nok[i] = FALSE
    #             sppoly$dropflag[i] = TRUE
    #             if ( sppoly$npts[j] >= areal_units_constraint_nmin) {
    #               local_finished=TRUE
    #               sppoly$nok[j] = TRUE
    #             }
    #           }
    #         }
            
    #       }

    #     }
    #     # final check
    #     toofew = which( sppoly$npts < areal_units_constraint_nmin )
    #     if (length(toofew) > 0) sppoly$dropflag[toofew] = TRUE

    #     sppoly = sppoly[ - which( sppoly$dropflag ), ]
    #     sppoly$nok =NULL
    #   }

    #   # update counts
    #   ww = tapply( rep(1, nrow(constraintdata)), constraintdata$internal_id, sum, na.rm=T )
    #   sppoly$npts = 0
    #   oo = match( names(ww), as.character(sppoly$internal_id) )
    #   ii = which(is.finite(oo))
    #   sppoly$npts[ oo[ii] ] = ww[ii]
      
    #   sppoly$internal_id = NULL
    #   sppoly = st_make_valid(sppoly)

    #   message( "Dropping due to areal_units_constraint_nmin, now there are : ", nrow(sppoly), " areal units." )
    # }

  # message( "Applying constraints leaves: ",  nrow(sppoly), " areal units." )

  # return(sppoly)

  
  if ( nrow( sppoly) == 0 ) {
    message ( "No polygons meet the specified criteria (check sa_threshold_km2 ?)." )
    browser()
  }

}


  # SA check
  if ( sa_threshold_km2  == 0 ) {
    if ( exists("sa_threshold_km2", p)) sa_threshold_km2 = p$sa_threshold_km2
  }
  sppoly = st_transform( sppoly, st_crs( areal_units_proj4string_planar_km ))
  sppoly$au_sa_km2 = st_area(sppoly)
  attributes( sppoly$au_sa_km2 ) = NULL
  toremove = which( sppoly$au_sa_km2 < sa_threshold_km2 )
  if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?


  # --------------------
  # completed mostly, final filters where required

  sppoly = st_make_valid(sppoly)
  if (!exists( "AUID", sppoly)) sppoly[, "AUID"]  = as.character( 1:nrow(sppoly) )
  oo = which(duplicated(sppoly$AUID) )
  if (length(oo)>0) {
    if ( duplications_action=="union" ) {
      # adding features (islands, coastlines) can break areal units that might be best left together
      todrop = NULL
      for (o in oo) {
        uu = which(sppoly$AUID == sppoly$AUID[o])
        vv = st_union( sppoly[uu,])
        st_geometry( sppoly[uu[1],]) = st_geometry(vv)
        todrop = c(todrop, setdiff(uu, uu[1]))
      }
      sppoly = sppoly[- todrop, ]
    }

    if ( duplications_action=="separate" ) {
      # adding features (islands, coastlines) can break areal units that might be best left together
      for (o in oo) {
        uu = which(sppoly$AUID == sppoly$AUID[o])
        sppoly$AUID[ uu] = paste( sppoly$AUID[uu], 1:length(uu), sep="_")
      }
    }
    sppoly = st_make_valid(sppoly)
  }
   
  row.names(sppoly) = sppoly$AUID

  require(spdep)
  
  W.nb = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE, snap=areal_units_resolution_km )  # slow .. ~1hr?
  W.remove = which(card(W.nb) == 0)

  if ( length(W.remove) > 0 ) {
    # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
    W.keep = which(card(W.nb) > 0)
    W.nb = nb_remove( W.nb, W.remove )
    sppoly = sppoly[W.keep,]
    row.names(sppoly) = sppoly$AUID
    sppoly = sppoly[order(sppoly$AUID),]
    sppoly = st_make_valid(sppoly)
  }



  if (0) {
    W.nb = poly2nb(sppoly, row.names=sppoly$internal_id, queen=TRUE)  # slow .. ~1hr?
    plot(W.nb, st_geometry(sppoly))
    edit.nb(W.nb, polys=as(sppoly, "Spatial"), use_region.id=TRUE)

    # https://cran.r-project.org/web/packages/spdep/vignettes/nb_sf.html
     # sf::st_relate takes about 136 s. for a total of 139 s. to generate a queen neighbour object. The contiguity neighbour objects using st_queen
    edit(W.nb, polys=sppoly, use_region.id=TRUE)
  }


# ------------------------------------------------
# overalys of other areal units not directly related to AU's ( eg. large management areas /zones
# ... but they require SA estimates after any filters (above))
  if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {
 
    # as a last pass, calculate surface areas of each subregion .. could be done earlier but it is safer done here due to deletions above
    if (verbose) message("Computing surface areas for each subarea ... can be slow if this is your first time")
    for (subarea in c("cfanorth", "cfasouth", "cfa23", "cfa24", "cfa4x" ) ) {
      print(subarea)
      csa = polygon_managementareas( species="maritimes", area=subarea )
      csa = st_transform( st_as_sf(csa), st_crs( areal_units_proj4string_planar_km ) )
      ooo = st_intersection( csa, sppoly )
      ooo$surfacearea = st_area( ooo )
      vn = paste(subarea, "surfacearea", sep="_")
      sppoly[[ vn ]] = 0
      j = match( ooo$AUID, sppoly$AUID )      
      if (length(j) > 0)  sppoly[[ vn ]][j] = ooo$surfacearea
    }
  }


    # ------------------------------------------------
  sppoly$au_sa_km2 = st_area( sppoly )

  if (!exists("strata_to_keep", sppoly) ) sppoly$strata_to_keep = TRUE  # flag for aggregations

  nb = INLA::inla.read.graph( spdep::nb2mat( W.nb ))
  attr(nb, "region.id") = sppoly$AUID

  attr(sppoly, "nb") = nb  # adding neighbourhood as an attribute to sppoly
  attr(sppoly, "W.nb") = W.nb  # adding neighbourhood as an attribute to sppoly
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

  if (!is.null( return_crs )) sppoly=st_transform(sppoly, crs=st_crs(return_crs) )
  save(sppoly, file=areal_units_fn_full, compress=TRUE)
  message( "Saved polygons as: ", areal_units_fn_full )

  if (plotit) plot(sppoly["au_sa_km2"])
   
  return( sppoly )

}
