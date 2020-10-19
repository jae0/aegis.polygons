

areal_units = function( p=NULL,  plotit=FALSE, sa_threshold_km2=0, redo=FALSE, use_stmv_solution=FALSE, areal_units_constraint="none", ... ) {

  if (0) {
    areal_units_timeperiod="default"
    plotit=FALSE
    # sa_threshold_km2=0
    redo=FALSE
    use_stmv_solution=FALSE
  }

  p = parameters_add(p, list(...) ) # add passed args to parameter list, priority to args

  aegis_internal_resolution_km = ifelse (exists("aegis_internal_resolution_km", p), p$aegis_internal_resolution_km, 1 )
  areal_units_proj4string_planar_km =  ifelse (exists("areal_units_proj4string_planar_km", p), p$areal_units_proj4string_planar_km, "+proj=utm +ellps=WGS84 +zone=20 +units=km" )


  # these are required:
  spatial_domain =  ifelse (exists("spatial_domain", p), p$spatial_domain, "SSE" )
  areal_units_overlay =  ifelse (exists("areal_units_overlay", p), p$areal_units_overlay, "none" )
  areal_units_resolution_km =  ifelse (exists("areal_units_resolution_km", p), p$areal_units_resolution_km, 25 )
  areal_units_source =  ifelse (exists("areal_units_source", p), p$areal_units_source, "lattice" )
  areal_units_timeperiod =  ifelse (exists("areal_units_timeperiod", p), p$areal_units_timeperiod, "default" )
  areal_units_constraint =  ifelse (exists("areal_units_constraint", p), p$areal_units_constraint, "none" )
  areal_units_constraint_nmin =  ifelse (exists("areal_units_constraint_nmin", p), p$areal_units_constraint_nmin, 0)

  areal_units_fn = paste(
    spatial_domain,
    paste0(areal_units_overlay, collapse="_"),
    areal_units_resolution_km,
    areal_units_source,
    areal_units_timeperiod,
    areal_units_constraint,
    areal_units_constraint_nmin,
    sep="_"
  )

  areal_units_directory = project.datadirectory("aegis", "polygons", "areal_units" )

  areal_units_fn_full = file.path( areal_units_directory, paste(areal_units_fn, "rdata", sep="." ) )

  sppoly = NULL
  boundary = NULL

  if (!redo) {
    if (file.exists(areal_units_fn_full)) {
      load(areal_units_fn_full)
      # message( "Using areal units specified in ", areal_units_fn_full)
    }
    if( !is.null(sppoly) ) return(sppoly)
  }

  message( "Creating/over-writing areal units specified in ", areal_units_fn_full)

  if (areal_units_source == "lattice" ) {
    # res based on grids ... rather than arbitrary polygons
    # static features only so far
    if (use_stmv_solution) {
      sppoly = aegis_db_spatial_object(
        spatial_domain=spatial_domain,
        proj4string=areal_units_proj4string_planar_km,
        areal_units_resolution_km=areal_units_resolution_km,
        returntype="SpatialPolygonsDataFrame"
      )

    } else {

      # as points/grids
      pn = spatial_parameters( spatial_domain=spatial_domain )
      Z = bathymetry_db( p=pn, DS="aggregated_data" )
      names(Z)[which(names(Z)=="z.mean" )] = "z"
      Z = lonlat2planar(Z, pn$aegis_proj4string_planar_km)  # should not be required but to make sure
      Z = geo_subset( spatial_domain=spatial_domain, Z=Z )  # position and depth filters
      spdf0 = SpatialPointsDataFrame(Z[, c("plon", "plat")], data=Z, proj4string=sp::CRS(pn$aegis_proj4string_planar_km) )

      require(raster)
      raster_template = raster(extent(spdf0)) # +1 to increase the area
      res(raster_template) = areal_units_resolution_km  # in units of crs (which should be in  km)
      crs(raster_template) = projection(spdf0) # transfer the coordinate system to the raster

      sppoly = rasterize( spdf0[, c("plon", "plat")], raster_template, field=spdf0$z )
      sppoly = as(sppoly, "SpatialPolygonsDataFrame")
    }

    sppoly = as( sppoly, "sf")

  }


  # ------------------------------------------------

  if (areal_units_source %in% c( "stratanal_polygons", "groundfish_strata")  ) {
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    sppoly = maritimes_groundfish_strata( areal_units_timeperiod=ifelse( areal_units_timeperiod=="default", "pre2014", areal_units_timeperiod ), returntype="polygons" )
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
    sppoly = as(sppoly, "sf")
    sppoly$AUID = as.character(sppoly$AUID)
    row.names(sppoly) = sppoly$AUID
  }



  if (areal_units_source %in% c("groundfish_polygons_inla_mesh",  "groundfish_polygons_tesselation") ) {

    ## method 1: Voronoi triangulation
    gfset = survey_db( DS="set.base", p=p )
    gfset = lonlat2planar(gfset, areal_units_proj4string_planar_km)  # should not be required but to make sure
    gfset = geo_subset( spatial_domain=spatial_domain, Z=gfset )
    gfset$AUID = gfset$id

    gfset =  st_as_sf ( gfset, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
    gfset = st_transform( gfset, st_crs( areal_units_proj4string_planar_km ))

    locs = st_coordinates( gfset )
    locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

    boundary = maritimes_fishery_boundary( DS="groundfish", internal_resolution_km=1, crs_km=st_crs(sppoly) ) # post 2014 is larger
    boundary = st_transform(boundary, st_crs(sppoly) )

    if (0) {
      #altenate way of determining boundary based upon data .. slower so turned off
      spbuffer = 5
      hull_multiplier = 6
      boundary = non_convex_hull( locs, alpha=spbuffer*hull_multiplier  )
      boundary = list( Polygons(list( Polygon( as.matrix( boundary ) )), ID="boundary" ))
      boundary = SpatialPolygons( boundary, proj4string=sp::CRS(projection_proj4string("utm20")) )
      boundary = gBuffer( gUnaryUnion( gBuffer( boundary, width=spbuffer, byid=TRUE) ), width=spbuffer)
      # plot(boundary)
    }

    if (areal_units_source == "groundfish_polygons_tesselation" ) {

     spmesh = aegis_mesh( SPDF=as(gfset, "Spatial"), SPDF_boundary=as(boundary, "Spatial"), resolution=areal_units_resolution_km,
       spbuffer=areal_units_resolution_km, areal_units_constraint_nmin=areal_units_constraint_nmin, tus="yr", output_type="polygons" )  # voroni tesslation and delaunay triagulation
    }

    if (areal_units_source == "groundfish_polygons_inla_mesh" ) {
      require(INLA)
      spmesh = inla.mesh.2d (
        loc=locs,
        max.edge = c( 0.5, 5 ) * areal_units_resolution_km,  #   # max size of a triange (in, out)
        offset = c( 0.1, 1 ) * areal_units_resolution_km , # how much to extend inside and outside of boundary,
        cutoff = c( 0.5, 5 ) * areal_units_resolution_km # min distance allowed between points #,
        # boundary =  inla.mesh.segment(st_coordinates( as(boundary, "sf") )[,c(1,2)])
      )

      # convert to sp*
      spmesh = SpatialPolygonsDataFrame(
        Sr = SpatialPolygons( lapply(
          1:nrow(spmesh$graph$tv),
          function(x) {
            tv = spmesh$graph$tv[x, , drop = TRUE]
            Polygons(list(Polygon(
              spmesh$loc[tv[c(1, 3, 2, 1)], 1:2, drop = FALSE])),
            ID = x)
          }
          ),
        proj4string = sp::CRS( areal_units_proj4string_planar_km )
        ),
        data = as.data.frame(spmesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
        match.ID = FALSE
      )
    }


    # must be done separately
    sppoly = (
      st_intersection( as( spmesh, "sf"), boundary )
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )

  #  remove.coastline
    require(aegis.coastline)
    coast = (
        as( coastline_db( p=p, DS="eastcoast_gadm" ), "sf")
        %>% st_transform( sp::CRS( areal_units_proj4string_planar_km ))
        %>% st_simplify()
        %>% st_buffer(0.1)
        %>% st_union()
    )
    sppoly = st_difference( sppoly, coast)

  }


  if (areal_units_source %in% c("snowcrab_polygons_inla_mesh",  "snowcrab_polygons_tesselation") ) {

    snset = snowcrab.db( p=p, DS="set.clean"  )  #
    snset = lonlat2planar(snset, areal_units_proj4string_planar_km)  # should not be required but to make sure
    snset = geo_subset( spatial_domain=spatial_domain, Z=snset )

    snset =  st_as_sf ( snset, coords= c('lon', 'lat'), crs = st_crs(projection_proj4string("lonlat_wgs84")) )
    snset = st_transform( snset, st_crs( areal_units_proj4string_planar_km ))
    # snset$AUID = snset$id
    # rownames( snset ) = snset$AUID
    # sp::proj4string(snset) = projection_proj4string("lonlat_wgs84")
    # snset = spTransform( snset, sp::CRS( areal_units_proj4string_planar_km ) )  # in km

    locs = st_coordinates( snset )
    locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

    spbuffer = 5
    hull_multiplier = 6

    boundary = non_convex_hull( locs, alpha=spbuffer*hull_multiplier  )
    boundary = st_sfc( st_multipoint( as.matrix(boundary) ), crs=st_crs(areal_units_proj4string_planar_km)  )

    if (0) {
      # using sp, defunct
      boundary = non_convex_hull( locs, alpha=spbuffer*hull_multiplier  )
      boundary = list( Polygons(list( Polygon( as.matrix( boundary ) )), ID="boundary" ))
      boundary = SpatialPolygons( boundary, proj4string=sp::CRS(projection_proj4string("utm20")) )
      boundary = gBuffer( gUnaryUnion( gBuffer( boundary, width=spbuffer, byid=TRUE) ), width=spbuffer)
            # plot(boundary)
    }


    if (areal_units_source == "snowcrab_polygons_tesselation" ) {

      spmesh = aegis_mesh( SPDF=as(snset, "Spatial"), SPDF_boundary=as(boundary, "Spatial"), resolution=areal_units_resolution_km, spbuffer=areal_units_resolution_km, areal_units_constraint_nmin=areal_units_constraint_nmin, tus="yr", output_type="polygons", boundary=boundary )  # voroni tesslation and delaunay triagulation
    }

    if (areal_units_source == "snowcrab_polygons_inla_mesh" ) {

      require(INLA)
      spmesh = inla.mesh.2d (
        loc=locs,
        max.edge = c( 0.5, 5 ) * areal_units_resolution_km,  #   # max size of a triange (in, out)
        offset = c( 0.1, 1 ) * areal_units_resolution_km , # how much to extend inside and outside of boundary,
        cutoff = c( 0.5, 5 ) * areal_units_resolution_km # min distance allowed between points #,
        # boundary =  inla.mesh.segment(st_coordinates( as(boundary, "sf") )[,c(1,2)])
      )

      # convert to sp*
      spmesh = SpatialPolygonsDataFrame(
        Sr = SpatialPolygons( lapply(
          1:nrow(spmesh$graph$tv),
          function(x) {
            tv = spmesh$graph$tv[x, , drop = TRUE]
            Polygons(list(Polygon(
              spmesh$loc[tv[c(1, 3, 2, 1)], 1:2, drop = FALSE])),
            ID = x)
          }
          ),
        proj4string = sp::CRS( areal_units_proj4string_planar_km )
        ),
        data = as.data.frame(spmesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
        match.ID = FALSE
      )
    }

    # must be done separately
    sppoly = (
      st_intersection( as( spmesh, "sf"), st_transform( st_cast(boundary, "POLYGON" ), st_crs(spmesh )) )
      %>% st_simplify()
      %>% st_cast( "POLYGON" )
    )

  }

  if (is.null(sppoly)) stop( "areal_units_source was not recognized!")

  # --------------------
  # Overlays

  if ( grepl("groundfish_strata", areal_units_overlay) ) {

    areal_units_timeperiod=ifelse( areal_units_timeperiod=="default", "pre2014", areal_units_timeperiod )

    ### next section is the same as "stratanal polygons" below .. copied here to avoid recursion
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers

    gf = as( maritimes_groundfish_strata( areal_units_timeperiod=areal_units_timeperiod, returntype="polygons" ), "sf")
    gf = st_make_valid(gf)
    gf = st_transform(gf, st_crs(sppoly) )
    gf$gfUID = as.character(gf$AUID)
    gf$AUID = NULL
    row.names(gf) = gf$gfUID


    boundary = maritimes_groundfish_boundary( areal_units_timeperiod=areal_units_timeperiod, internal_resolution_km=aegis_internal_resolution_km, crs_km=st_crs(sppoly) )
    boundary = st_transform(boundary, st_crs(sppoly) )

    sp0 = sppoly
    sp0$AUID = as.character(1:length(sp0))

    sppoly = (
      st_intersection( st_intersection( boundary, sppoly ), gf )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )

    sppoly = (
      st_join( st_sf(sppoly), gf, join=st_intersects,  largest=TRUE )
    )

    # recover the data from sp0
    sppoly = (
      st_join( st_sf(sppoly), sp0, join=st_intersects,  largest=TRUE )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )
    sppoly = st_make_valid(sppoly)
    sppoly$index = sppoly$AUID
    sppoly$AUID = paste( sppoly$AUID , sppoly$gfUID, sep="_")

  }

  # --------------------

  if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {
    # main domain/boundary
    sppoly = (
      st_intersection(
        sppoly,
        st_transform(maritimes_fishery_boundary( DS="maritimes", internal_resolution_km=0.1, crs_km=st_crs(areal_units_proj4string_planar_km) ), st_crs(sppoly) ) )
      %>% st_simplify()
    )
  }


  # --------------------


  if (class( areal_units_constraint ) == "character") {
    if (areal_units_constraint == "snowcrab") {
      areal_units_constraint = snowcrab.db( p=p, DS="set.clean" )[, c("lon", "lat")]  #
    }
    if (any( areal_units_constraint %in% c("groundfish", "aegis.survey" ) )) {
      areal_units_constraint = gfset = survey_db( DS="set.base", p=p )[, c("lon", "lat")]  #
    }
  }

  if (inherits( areal_units_constraint, c("data.frame", "matrix") ) ) {
    # this part done as a "Spatial" object
    # already done in tessilation method so not really needed but
    # for other methods, this ensures a min no samples in each AU
    # todo:  convert to sf:: st_join

    cst = st_sfc( st_multipoint( as.matrix(areal_units_constraint) ), crs=st_crs(projection_proj4string("lonlat_wgs84"))  )
    cst = st_transform( cst, st_crs(sppoly ))
    sppoly$internal_id = 1:nrow(sppoly)
    sppolywithdata = st_intersection(  sppoly, cst )$internal_id  # all with at least 1
    sppoly = sppoly[ match( sppolywithdata, sppoly$internal_id), ]
    if (areal_units_constraint_nmin > 0 ) {
      sppoly$npts  = lengths( st_intersection( cst, sppoly ) )
      todrop = which( sppoly$npts < areal_units_constraint_nmin )
      if (length(todrop) > 0 ) sppoly = sppoly[ - todrop, ]
    }
  }

  # --------------------
  # completed mostly, final filters where required
  if (!exists("AUID", sppoly)) {
    sppoly[, "AUID"]  = as.character( 1:nrow(sppoly) )
    # rownames(sppoly) = sppoly$AUID
  }

  sppoly = st_transform( sppoly, sp::CRS( areal_units_proj4string_planar_km ))
  sppoly$au_sa_km2 = st_area(sppoly)

  # plot(st_geometry(sppoly))
  # plot(sppoly[,"au_sa_km2"])

  attributes( sppoly$au_sa_km2 ) = NULL
  if ( sa_threshold_km2 == 0 ) {
    if ( exists("sa_threshold_km2", p)) sa_threshold_km2 = p$sa_threshold_km2
  }
  toremove = which( sppoly$au_sa_km2 < sa_threshold_km2 )
  if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?

  if ( nrow( sppoly) == 0 ) {
    message ( "No polygons meet the specified criteria (check sa_threshold_km2 ?)." )
    return (sppoly)
  }

  if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {
    # as a last pass, calculate surface areas of each subregion .. could be done earlier but it is safer done here due to deletions above
    message("Computing surface areas for each subarea ... can be slow if this is your first time")
    for (subarea in c("cfanorth", "cfasouth", "cfa23", "cfa24", "cfa4x" ) ) {
      print(subarea)
      csa = polygon_managementareas( species="snowcrab", area=subarea )
      csa = as(csa, "sf")
      csa = st_transform( csa, sp::CRS( areal_units_proj4string_planar_km ) )
      ooo = st_intersection( csa, sppoly )
      ooo$surfacearea = st_area( ooo )
      vn = paste(subarea, "surfacearea", sep="_")
      sppoly[[ vn ]] = 0
      j = match( ooo$AUID, sppoly$AUID )
      if (length(j) > 0)  sppoly[[ vn ]][j] = ooo$surfacearea
    }
  }


  # ------------------------------------------------


  # force saves as Spatial* data (for poly2nb, really)
  sppoly = st_make_valid(sppoly)
  sppoly = as(sppoly, "Spatial")

  # poly* function operate on Spatial* data
  W.nb = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE)  # slow .. ~1hr?
  W.remove = which(card(W.nb) == 0)

  if ( length(W.remove) > 0 ) {
    # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
    W.keep = which(card(W.nb) > 0)
    W.nb = nb_remove( W.nb, W.remove )
    sppoly = sppoly[W.keep,]

    row.names(sppoly) = sppoly$AUID
    sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
    sppoly = sppoly[order(sppoly$AUID),]
  }

  attr(sppoly, "nb") = W.nb  # adding neighbourhood as an attribute to sppoly
  attr(sppoly, "spatial_domain") = spatial_domain
  attr(sppoly, "areal_units_fn") = areal_units_fn
  attr(sppoly, "areal_units_fn_full") = areal_units_fn_full
  attr(sppoly, "areal_units_directory") = areal_units_directory
  attr(sppoly, "areal_units_overlay") = areal_units_overlay
  attr(sppoly, "areal_units_resolution_km") = areal_units_resolution_km
  attr(sppoly, "areal_units_source") = areal_units_source
  attr(sppoly, "areal_units_timeperiod") = areal_units_timeperiod
  attr(sppoly, "areal_units_constraint_nmin") = areal_units_constraint_nmin

  save(sppoly, file=areal_units_fn_full, compress=TRUE)

  if (plotit) plot(sppoly)

  return( sppoly )

}
