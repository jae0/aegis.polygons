

areal_units = function( 
  p=NULL, 
  areal_units_fn_full=NULL, 
  areal_units_directory=NULL, 
  plotit=FALSE, 
  sa_threshold_km2=0, 
  redo=FALSE,
  use_stmv_solution=TRUE, 
  rastermethod="sf",  
  xydata=NULL, 
  spbuffer=5, 
  n_iter_drop=1, 
  hull_noise=1e-4, 
  lenprob=0.9,
  duplications_action="union",  
  areal_units_timeperiod=NULL, 
  verbose=FALSE, 
  return_crs=NULL, 
  count_time=TRUE, 
  respect_spatial_domain=TRUE, 
  ... ) {

  if (0) {
    plotit=FALSE
    sa_threshold_km2=0
    redo=FALSE
    use_stmv_solution=TRUE
    spbuffer=5
    rastermethod="sf"
    xydata=NULL
    lenprob = 0.9
    duplications_action="union"
    areal_units_timeperiod=NULL
    areal_units_fn_full = NULL
    areal_units_directory=NULL
    verbose=TRUE
    n_iter_drop=1 
    hull_noise=1e-4 
    return_crs=NULL
    count_time=TRUE 
    respect_spatial_domain=TRUE
  }

  require(spdep)

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args

  # hull (boundary) related:
  if (exists("spbuffer", p)) spbuffer = p$spbuffer
  if (exists("lenprob", p)) lenprob = p$lenprob
  if (exists("hull_noise", p)) hull_noise = p$hull_noise
  if (exists("n_iter_drop", p)) n_iter_drop = p$n_iter_drop
  if (exists("sa_threshold_km2", p)) sa_threshold_km2 = p$sa_threshold_km2

  if (exists("rastermethod", p)) rastermethod = p$rastermethod
  if (exists("use_stmv_solution", p)) use_stmv_solution = p$use_stmv_solution
  if (exists("areal_units_directory", p)) {
    if (is.null(areal_units_directory)) {
      areal_units_directory = p$areal_units_directory  # priority to arg passed to this functions
    }
  }
  
  if (exists("areal_units_fn_full", p)) areal_units_fn_full = p$areal_units_fn_full
  
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
  areal_units_constraint_ntarget =  ifelse (exists("areal_units_constraint_ntarget", p), round(p$areal_units_constraint_ntarget), 0)

  areal_units_fn = areal_units_filename(
    spatial_domain = spatial_domain,
    areal_units_type = areal_units_type,
    areal_units_resolution_km = areal_units_resolution_km,
    areal_units_constraint = areal_units_constraint,
    areal_units_constraint_ntarget = areal_units_constraint_ntarget,
    areal_units_constraint_nmin = areal_units_constraint_nmin,
    areal_units_timeperiod = areal_units_timeperiod,
    areal_units_overlay = areal_units_overlay )
    

  if ( !is.null(areal_units_fn_full) ) {
    if (is.null(areal_units_directory)) {
      areal_units_directory =  dirname(areal_units_fn_full)
    }
  }

  if ( is.null(areal_units_fn_full) )  {
    if ( is.null(areal_units_directory) )  areal_units_directory = file.path( p$datadir, "areal_units" )
    areal_units_fn_full = file.path( areal_units_directory, paste(areal_units_fn, "rdz", sep="." ) )
  }

  dir.create( areal_units_directory, showWarnings = FALSE, recursive = TRUE )
  
  sppoly = NULL

  if (!redo) {
    if ( file.exists(areal_units_fn_full) ) {
      sppoly = read_write_fast(areal_units_fn_full)
      if (!is.null(sppoly)) {
        if ( !is.null( return_crs )) sppoly=st_transform(sppoly, crs=st_crs(return_crs) )
        if ( !is.null( sppoly ) ) return(sppoly)
      }
    }
  }

  message("***\n
      If any new parameter settings are used for sppoly creation, \n
      then move them into *_parameters.R as the lookup mechanism \n 
      uses these parameter settings for file lookup \n***\n")

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
    xydata = st_as_sf ( xydata, coords= c('lon', 'lat'))
    st_crs(xydata) = st_crs(projection_proj4string("lonlat_wgs84")) 
  }
  xydata = st_transform( xydata, st_crs( areal_units_proj4string_planar_km ))


  boundary = NULL

  if ( areal_units_type %in% c( "stratanal_polygons_pre2014", "stratanal_polygons", "groundfish_strata", "stratanal_polygons_post2014") ) {
    message( "Determining areal unit domain boundary from groundfish survey")
    boundary = maritimes_fishery_boundary( DS="groundfish", internal_resolution_km=1, crs_km=st_crs(areal_units_proj4string_planar_km) ) # post 2014 is larger (crm is for incoming data)
    boundary = st_transform( boundary, st_crs(areal_units_proj4string_planar_km) ) # output of above is lonlat
    boundary = (
        boundary
        %>% st_cast("POLYGON" )
        %>% st_make_valid()
        %>% st_simplify()
        %>% st_buffer( areal_units_resolution_km )
        %>% st_union()
        %>% st_cast("POLYGON" )
        %>% st_make_valid()
      )

  } 
    
  if ( project_name == "bio.snowcrab") {
 
      message( "Determining areal unit domain boundary from snowcrab survey")
      boundary = polygon_managementareas( species="snowcrab" )
      boundary = st_transform( boundary, st_crs(areal_units_proj4string_planar_km) )
      boundary = st_buffer(boundary, 0)  
      boundary =  st_simplify(boundary, TRUE, 1 )
      boundary = st_cast(boundary, "POLYGON" )
      boundary = st_make_valid(boundary)
 
      data_boundary = (
        st_combine(xydata)
        %>% st_concave_hull( ratio=0.01, allow_holes=FALSE )  #      xyd = non_convex_hull( xydata, lengthscale=spbuffer, lenprob=lenprob )  
        %>% st_sfc(crs=st_crs(areal_units_proj4string_planar_km))
        %>% st_cast("POLYGON" )
        %>% st_simplify(dTolerance=inputdata_spatial_discretization_planar_km)
        %>% st_union()
        %>% st_make_valid()
      )
 
      boundary = st_intersection(data_boundary, boundary)
  
  }  
 
  
  if (is.null(boundary)) {
      message( "Determining areal unit domain boundary from input: xydata") 
      boundary = (
        st_combine(xydata)
        %>% st_concave_hull( ratio=0.01, allow_holes=FALSE )  #      xyd = non_convex_hull( xydata, lengthscale=spbuffer, lenprob=lenprob )  
        %>% st_sfc(crs=st_crs(areal_units_proj4string_planar_km))
        %>% st_cast("POLYGON" )
        %>% st_simplify(dTolerance=inputdata_spatial_discretization_planar_km)
        %>% st_union()
        %>% st_make_valid()
      )
 
  }
    #  remove.coastline
    require(aegis.coastline)
    coast = (
        coastline_db( p=p, DS="eastcoast_gadm" )
        %>% st_transform( st_crs( areal_units_proj4string_planar_km ))
        %>% st_simplify(dTolerance=inputdata_spatial_discretization_planar_km)
        %>% st_union()
        %>% st_cast("MULTIPOLYGON")
    )

    boundary = st_difference( boundary, coast )
    coast = NULL

    
    if (0){
      plot( xydata, reset=FALSE )
      plot( boundary, add=TRUE )
    }


    if ( areal_units_type == "inla_mesh" ) {
      message( "Determining areal units as an INLA mesh")
      sppoly = areal_units_inla_mesh(
        locs=st_coordinates( xydata ) [,c("X", "Y")] + runif( nrow(xydata)*2, min=-hull_noise, max=hull_noise ) , # add  noise  to prevent a race condition
        areal_units_resolution_km=areal_units_resolution_km,
        areal_units_proj4string_planar_km=areal_units_proj4string_planar_km
      )
    }
    
    if ( areal_units_type == "tesselation" ) {
      message( "Determining areal units via iterative Voronoi tesselation of AU centroids and dissolution of AUs")
      # boundary precomputed specific to project here if possible to give more control

      sppoly = aegis_mesh(
        pts=xydata,
        boundary=boundary,
        output_type="polygons",
        resolution=areal_units_resolution_km,
        spbuffer=spbuffer,
        hull_lengthscale=spbuffer,  # for rasterization .. not used if boundary is provided
        areal_units_constraint_ntarget=areal_units_constraint_ntarget,
        areal_units_constraint_nmin=areal_units_constraint_nmin,
        tus=p$tus,
        fraction_todrop = p$fraction_todrop,
        fraction_cv = p$fraction_cv,  # stopping criterion: when cv drops below this value
        nAU_min = p$nAU_min,   # stoppping criterion: allow no less than this number of areal units
        count_time = count_time,
        verbose = verbose
      )  # Voronoi tesslation and delaunay triagulation
    }


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
 

    # --------------------
    # Additional Constraints from other data

    # # SA check
    if ( sa_threshold_km2  == 0 ) {
      if ( exists("sa_threshold_km2", p)) sa_threshold_km2 = p$sa_threshold_km2
    }


    sppoly = st_transform( sppoly, st_crs( areal_units_proj4string_planar_km ))
    sppoly$au_sa_km2 = st_area(sppoly)
    attributes( sppoly$au_sa_km2 ) = NULL
    toremove = which( sppoly$au_sa_km2 < sa_threshold_km2 )
    if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?


    constraintdata = NULL

    if (areal_units_constraint == "groundfish")  {
      message( "Constrain areal units to required number of groundfish survey locations by merging into adjacent AUs")
      constraintdata = groundfish_survey_db(DS="set.base", yrs=p$yrs )[, c("lon", "lat")]  #
    }
 
    if (areal_units_constraint == "survey")    {
      message( "Constrain areal units to required number of aegis.survey locations by merging into adjacent AUs")
      constraintdata = survey_db( p=p, DS="set.base" )[, c("lon", "lat")]  #
    }

    if (!is.null(constraintdata)) {
      constraintdata = sf::st_as_sf( constraintdata, coords = c("lon","lat"), crs=st_crs(projection_proj4string("lonlat_wgs84")) )
      constraintdata = st_transform( constraintdata, st_crs(areal_units_proj4string_planar_km ))
    }

    if (is.null(constraintdata)) constraintdata = xydata
    xydata = NULL
  
    # count
    sppoly$internal_id = 1:nrow(sppoly)
    row.names( sppoly) = sppoly$internal_id
    cdd = setDT( st_drop_geometry( constraintdata ) )
    cdd$internal_id = st_points_in_polygons( constraintdata, sppoly, varname="internal_id" )
    cdd = cdd[,.(npts=.N), .(internal_id) ]
    sppoly$npts = cdd$npts[ match( as.character(sppoly$internal_id),  as.character(cdd$internal_id) )]

    uu = which( !is.finite(sppoly$npts) )
    if ( length(uu) > 0 ) sppoly$npts[uu] = 0
 
    if (areal_units_type == "lattice" ) {
      # lattice structure is required, simply drop where there is no data
        todrop = which( sppoly$npts < areal_units_constraint_nmin )
      if (length(todrop) > 0 ) {
        sppoly = sppoly[ -todrop , ]
      }

    } else if (areal_units_type == "tesselation" ) {

      # try to join to adjacent au's
      todrop = which( sppoly$npts < areal_units_constraint_nmin )
      if (length(todrop) == 0 ) break()

      sppoly$already_dropped = FALSE
      sppoly$count_is_ok = TRUE
      sppoly$count_is_ok[todrop] = FALSE

      NB_graph = poly2nb(sppoly, row.names=sppoly$internal_id, queen=TRUE)  
      message( "Joining adjacent areal units to obtain target minimum number of data points in a cell: ", areal_units_constraint_nmin )
      
      for (i in order(sppoly$npts) ) {
        if ( sppoly$count_is_ok[i]) next()
        lnb = NB_graph[[ i ]]
        if (length(lnb) < 1) next()
        local_finished = FALSE
        for (f in 1:length(lnb) ) {
          if (!local_finished) {
            v = setdiff( 
              intersect( lnb, which(sppoly$count_is_ok) ), ## AU neighbours that are OK and so can consider dropping
              which(sppoly$already_dropped)                       ## AUs confirmed already to drop
            )
            if (length(v) == 0 ) break()

            j = v[ which.min( sppoly$npts[v] )]
            if (length(j)>0) {
              g_ij = try( st_union( st_geometry(sppoly)[j] , st_geometry(sppoly)[i] ) )
              if ( !inherits(g_ij, "try-error" )) {
                st_geometry(sppoly)[j] = g_ij
                sppoly$npts[j] = sppoly$npts[j] + sppoly$npts[i]
                sppoly$count_is_ok[i] = FALSE
                sppoly$already_dropped[i] = TRUE
                if ( sppoly$npts[j] >= areal_units_constraint_nmin) {
                  local_finished=TRUE
                  sppoly$count_is_ok[j] = TRUE
                }
              }
            }

          }
        }
      }
      # final check
      toofew = which( sppoly$npts < areal_units_constraint_nmin )
      if (length(toofew) > 0) sppoly$already_dropped[toofew] = TRUE
      sppoly = sppoly[ - which( sppoly$already_dropped ), ]
      sppoly$count_is_ok =NULL
      message( "After merge, there are:  ", nrow(sppoly), " areal units." )
    }
        
    message( "Last pass: Force dropping of low count locations:" )
    # update counts: iterativ, so do it a few times
    for (i in 1:n_iter_drop) {

      uu = jitter( st_coordinates(st_centroid(sppoly))[,c("X", "Y")] )  # jitter noise to keep from getting stuck
      sppoly = tessellate( uu, outformat="sf", crs=st_crs( sppoly )) # centroids via voronoi
      sppoly = st_sf( st_intersection( sppoly, boundary ) ) # crop
 
      sppoly$internal_id = 1:nrow(sppoly)
      row.names( sppoly) = sppoly$internal_id
      cdd = setDT( st_drop_geometry( constraintdata ) )
      cdd$internal_id = st_points_in_polygons( constraintdata, sppoly, varname="internal_id" )
      cdd = cdd[,.(npts=.N), .(internal_id) ]
      sppoly$npts = cdd$npts[ match( as.character(sppoly$internal_id),  as.character(cdd$internal_id) )]
      sppoly = st_make_valid(sppoly)

      # drop 
      o =  which( !is.finite(sppoly$npts ) | sppoly$npts < areal_units_constraint_nmin )
      if (length(o)==0) break()

      sppoly = sppoly[ -o, ]

    } 

    # final update      

    uu =  st_coordinates(st_centroid(sppoly))[,c("X", "Y")] 

    sppoly = tessellate( uu, outformat="sf", crs=st_crs( sppoly )) # centroids via voronoi
    sppoly = st_sf( st_intersection( sppoly, boundary ) ) # crop
   
    sppoly$internal_id = 1:nrow(sppoly)
    row.names( sppoly) = sppoly$internal_id
    cdd = setDT( st_drop_geometry( constraintdata ) )
    cdd$internal_id = st_points_in_polygons( constraintdata, sppoly, varname="internal_id" )
    cdd = cdd[,.(npts=.N), .(internal_id) ]
    sppoly$npts = cdd$npts[ match( as.character(sppoly$internal_id),  as.character(cdd$internal_id) )]
    sppoly$internal_id = NULL
    sppoly = st_make_valid(sppoly)
 
    message( "Applying additional constraints leaves: ",  nrow(sppoly), " areal units." )

    if ( nrow( sppoly) == 0 ) {
      message ( "No polygons meet the specified criteria (check sa_threshold_km2 ?)." )
      browser()
    }

    
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
  
  NB_graph = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE, snap=areal_units_resolution_km )  # slow .. ~1hr?
  NB_remove = which(card(NB_graph) == 0)

  if ( length(NB_remove) > 0 ) {
    # remove isolated locations and recreate sppoly .. alternatively add links to NB_graph
    NB_keep = which(card(NB_graph) > 0)
    NB_graph = nb_remove( NB_graph, NB_remove )
    sppoly = sppoly[NB_keep,]
    row.names(sppoly) = sppoly$AUID
    sppoly = sppoly[order(sppoly$AUID),]
    sppoly = st_make_valid(sppoly)
  }



  if (0) {
    NB_graph = poly2nb(sppoly, row.names=sppoly$internal_id, queen=TRUE)  # slow .. ~1hr?
    plot(NB_graph, st_geometry(sppoly))
    edit.nb(NB_graph, polys=as(sppoly, "Spatial"), use_region.id=TRUE)

    # https://cran.r-project.org/web/packages/spdep/vignettes/nb_sf.html
     # sf::st_relate takes about 136 s. for a total of 139 s. to generate a queen neighbour object. The contiguity neighbour objects using st_queen
    edit(NB_graph, polys=sppoly, use_region.id=TRUE)
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

  # update counts
  sppoly = st_make_valid(sppoly)
  
  sppoly$internal_id = 1:nrow(sppoly)
  row.names( sppoly) = sppoly$internal_id
  cdd = setDT( st_drop_geometry( constraintdata ) )
  cdd$internal_id = st_points_in_polygons( constraintdata, sppoly, varname="internal_id" )
  cdd = cdd[,.(npts=.N), .(internal_id) ]
  sppoly$npts = cdd$npts[ match( as.character(sppoly$internal_id),  as.character(cdd$internal_id) )]
    
  sppoly$internal_id = NULL

  # ------------------------------------------------

  sppoly = st_make_valid(sppoly)  # once more, just in case
  sppoly$au_sa_km2 = st_area( sppoly )

  if (!exists("strata_to_keep", sppoly) ) sppoly$strata_to_keep = TRUE  # flag for aggregations

  attr(NB_graph, "region.id") = sppoly$AUID

  attr(sppoly, "NB_graph") = NB_graph  # adding neighbourhood as an attribute to sppoly
  attr(sppoly, "nb") = INLA::inla.read.graph( spdep::nb2mat( NB_graph ))  # adding neighbourhood as an attribute to sppoly
  
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
  read_write_fast(sppoly, file=areal_units_fn_full)
  message( "Saved polygons as: ", areal_units_fn_full )

  if (plotit) plot(sppoly["au_sa_km2"])
   
  return( sppoly )

}


