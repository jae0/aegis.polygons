
areal_units = function( p=NULL, timeperiod="default", plotit=FALSE, sa_threshold_km2=0, redo=FALSE, use_stmv_solution=FALSE, areal_units_constraint="none", ... ) {

  if (0) {
    timeperiod="default"
    plotit=FALSE
    sa_threshold_km2=0
    redo=FALSE
    use_stmv_solution=FALSE
  }

  p = parameters_control(p, list(...), control="add") # add passed args to parameter list, priority to args

  spatial_domain =  ifelse (exists("spatial_domain", p), p$spatial_domain, "SSE" )

  aegis_internal_resolution_km = ifelse (exists("aegis_internal_resolution_km", p), p$aegis_internal_resolution_km, 1 )
  areal_units_proj4string_planar_km =  ifelse (exists("areal_units_proj4string_planar_km", p), p$areal_units_proj4string_planar_km, "+proj=utm +ellps=WGS84 +zone=20 +units=km" )
  areal_units_source =  ifelse (exists("areal_units_source", p), p$areal_units_source, "lattice" )
  areal_units_resolution_km =  ifelse (exists("areal_units_resolution_km", p), p$areal_units_resolution_km, 20 )
  areal_units_overlay =  ifelse (exists("areal_units_overlay", p), p$areal_units_overlay, "none" )

  areal_units_fn = ifelse ( exists( "areal_units_fn", p), p$areal_units_fn,
    paste(
      spatial_domain,
      paste0(areal_units_overlay, collapse="_"),
      areal_units_resolution_km,
      areal_units_source,
      timeperiod,
      sep="_"
    )
  )

  fn = file.path( project.datadirectory("aegis", "polygons", "areal_units" ), paste(areal_units_fn, "rdata", sep="." ) )
  print(fn)

  sppoly = NULL

  if (!redo) {
    if (file.exists(fn)) load(fn)
    if( !is.null(sppoly) ) return(sppoly)
  }


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

      spdf0 = spatial_domain_discretized(spatial_domain)

      require(raster)
      raster_template = raster(extent(spdf0)) # +1 to increase the area
      res(raster_template) = areal_units_resolution_km  # in units of crs (which should be in  km)
      crs(raster_template) = projection(spdf0) # transfer the coordinate system to the raster

      sppoly = rasterize( spdf0[, c("plon", "plat")], raster_template, field=spdf0$z )
    }

    sppoly = as( as(sppoly, "SpatialPolygonsDataFrame"), "sf")

  }


  # ------------------------------------------------

  if (areal_units_source == "stratanal_polygons") {
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    sppoly = maritimes_groundfish_strata( timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod ), returntype="polygons" )
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
    sppoly = as(sppoly, "sf")
    sppoly$AUID = as.character(sppoly$AUID)
    row.names(sppoly) = sppoly$AUID
  }


  if (areal_units_source %in% c("groundfish_polygons_inla_mesh",  "groundfish_polygons_tesselation") ) {

    ## method 1: Voronoi triangulation
    gfset = survey.db( DS="set.base", p=p )
    gfset = lonlat2planar(gfset, areal_units_proj4string_planar_km)  # should not be required but to make sure
    gfset = geo_subset( spatial_domain=spatial_domain, Z=gfset )

    gfset$AUID = gfset$id
    coordinates(gfset) = ~ lon+lat
    rownames(gfset@data) = gfset$AUID
    sp::proj4string(gfset) = projection_proj4string("lonlat_wgs84")

    gfset = spTransform( gfset, sp::CRS( areal_units_proj4string_planar_km ) )  # in km
    locs = coordinates( gfset )
    locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

    spbuffer = 5
    hull_multiplier = 6
    data_boundary = non_convex_hull( locs, alpha=spbuffer*hull_multiplier  )
    data_boundary = list( Polygons(list( Polygon( as.matrix( data_boundary ) )), ID="boundary" ))
    data_boundary = SpatialPolygons( data_boundary, proj4string=sp::CRS(projection_proj4string("utm20")) )
    data_boundary = gBuffer( gUnaryUnion( gBuffer( data_boundary, width=spbuffer, byid=TRUE) ), width=spbuffer)
          # plot(data_boundary)

    if (areal_units_source == "groundfish_polygons_tesselation" ) {
      areal_units_tessilation_factor = ifelse( exists("areal_units_tessilation_factor", p), p$areal_units_tessilation_factor, 5 )  # reduction for initial pass
      spmesh = aegis_mesh( SPDF=gfset,  resolution=areal_units_resolution_km/ areal_units_tessilation_factor, output_type="grid.count" )  # coarse grid representation
      if (exists("areal_units_tessilation_nmin", p)) {
        spmesh = spmesh[ spmesh$layer > p$areal_units_tessilation_nmin, ]
        spmesh$AUID = as.character( 1:length(spmesh) )
        row.names(spmesh) = spmesh$AUID
      }
      spmesh = aegis_mesh( SPDF=spmesh, resolution=areal_units_resolution_km, spbuffer=areal_units_resolution_km, output_type="polygons" )  # voroni tesslation and delaunay triagulation
    }

    if (areal_units_source == "groundfish_polygons_inla_mesh" ) {

      require(INLA)
      spmesh = inla.mesh.2d (
        loc=locs,
        max.edge = c( 0.5, 5 ) * areal_units_resolution_km,  #   # max size of a triange (in, out)
        offset = c( 0.1, 1 ) * areal_units_resolution_km , # how much to extend inside and outside of boundary,
        cutoff = c( 0.5, 5 ) * areal_units_resolution_km # min distance allowed between points #,
        # boundary =  inla.mesh.segment(st_coordinates( as(data_boundary, "sf") )[,c(1,2)])
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

    require(aegis.coastline)
    coast = (
        as( coastline.db( p=p, DS="eastcoast_gadm" ), "sf")
        %>% st_transform( sp::CRS( areal_units_proj4string_planar_km ))
        %>% st_simplify()
        %>% st_buffer(0.1)
        %>% st_union()
    )

    sppoly = (
      as(data_boundary, "sf")
      %>% st_union()
      %>% st_buffer(0.1)
      %>% st_difference( coast)
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )

    # must be done separately
    sppoly = (
      st_intersection( as( spmesh, "sf"), sppoly )
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )

  }


  if (areal_units_source %in% c("snowcrab_polygons_inla_mesh",  "snowcrab_polygons_tesselation") ) {

    snset = snowcrab.db( p=p, DS="set.clean"  )  #
    snset = lonlat2planar(snset, areal_units_proj4string_planar_km)  # should not be required but to make sure
    snset = geo_subset( spatial_domain=spatial_domain, Z=snset )
    snset$AUID = snset$id
    coordinates(snset) = ~ lon+lat
    rownames(snset@data) = snset$AUID
    sp::proj4string(snset) = projection_proj4string("lonlat_wgs84")

    snset = spTransform( snset, sp::CRS( areal_units_proj4string_planar_km ) )  # in km
    locs = coordinates( snset )
    locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

    spbuffer = 5
    hull_multiplier = 6
    data_boundary = non_convex_hull( locs, alpha=spbuffer*hull_multiplier  )
    data_boundary = list( Polygons(list( Polygon( as.matrix( data_boundary ) )), ID="boundary" ))
    data_boundary = SpatialPolygons( data_boundary, proj4string=sp::CRS(projection_proj4string("utm20")) )
    data_boundary = gBuffer( gUnaryUnion( gBuffer( data_boundary, width=spbuffer, byid=TRUE) ), width=spbuffer)
          # plot(data_boundary)

    if (areal_units_source == "snowcrab_polygons_tesselation" ) {
      areal_units_tessilation_factor = ifelse( exists("areal_units_tessilation_factor", p), p$areal_units_tessilation_factor, 5 )  # reduction for initial pass
      spmesh = aegis_mesh( SPDF=snset,  resolution=areal_units_resolution_km/ areal_units_tessilation_factor, output_type="grid.count" )  # coarse grid representation
      if (exists("areal_units_tessilation_nmin", p)) {
        spmesh = spmesh[ spmesh$layer > p$areal_units_tessilation_nmin, ]
        spmesh$AUID = as.character( 1:length(spmesh) )
        row.names(spmesh) = spmesh$AUID
      }
      spmesh = aegis_mesh( SPDF=spmesh, resolution=areal_units_resolution_km, spbuffer=areal_units_resolution_km, output_type="polygons" )  # voroni tesslation and delaunay triagulation
    }

    if (areal_units_source == "snowcrab_polygons_inla_mesh" ) {

      require(INLA)
      spmesh = inla.mesh.2d (
        loc=locs,
        max.edge = c( 0.5, 5 ) * areal_units_resolution_km,  #   # max size of a triange (in, out)
        offset = c( 0.1, 1 ) * areal_units_resolution_km , # how much to extend inside and outside of boundary,
        cutoff = c( 0.5, 5 ) * areal_units_resolution_km # min distance allowed between points #,
        # boundary =  inla.mesh.segment(st_coordinates( as(data_boundary, "sf") )[,c(1,2)])
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
      st_intersection( as( spmesh, "sf"), as(data_boundary, "sf") )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )
  }


  # --------------------
  # Overlays

  if ( grepl("groundfish_strata", areal_units_overlay) ) {

    timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod )

    ### next section is the same as "stratanal polygons" below .. copied here to avoid recursion
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers

    gf = as( maritimes_groundfish_strata( timeperiod=timeperiod, returntype="polygons" ), "sf")
    gf$AUID = as.character(gf$AUID)
    row.names(gf) = gf$AUID

    o = over( as(sppoly, "Spatial"), as(gf, "Spatial") ) # match each datum to an area
    sppoly = sppoly[ which(!is.na(o$AUID)), ]
    sppoly = as( sppoly, "sf")

    gf = (
      as( gf, "sf" )
      %>% st_simplify()
      %>% st_buffer(aegis_internal_resolution_km)
      %>% st_union()
      %>% st_simplify()
    )

    sppoly = (
      st_intersection( gf, sppoly )
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )

  }

  # --------------------

  if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {

    spdf0 = spatial_domain_discretized(spatial_domain)

    require(raster)

    raster_resolution = p$pres

    # higher resolution raster
    raster_template = raster(extent(spdf0)) # +1 to increase the area
    res(raster_template) = raster_resolution  # in units of crs (which should be in 1 km .. so 100m)
    crs(raster_template) = projection(spdf0) # transfer the coordinate system to the raster


    rasterized_depths = rasterize( spdf0[, c("plon", "plat")], raster_template, field=spdf0$z )

    rasterized_depths = (
      as(as(rasterized_depths, "SpatialPolygons"), "sf")
      %>% st_buffer( raster_resolution )
      %>% st_simplify()
      %>% st_union( )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )

    sppoly = (
      st_intersection( sppoly, rasterized_depths )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )


    domain = (
      as( polygons_managementarea( species="snowcrab", area="cfaall") , "sf" )
      %>% st_transform( st_crs(sppoly) )
      %>% st_simplify()
      %>% st_buffer(0.1)
      %>% st_union()
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )

    sppoly = (
      st_intersection( sppoly, domain )
      %>% st_simplify()
      %>% st_collection_extract("POLYGON")
      %>% st_cast( "POLYGON" )
    )
 }


  # --------------------


  if (class( areal_units_constraint ) %in% c("data.frame", "matrix") ) {
    # this part done as a "Spatial" object

    sppoly = as(sppoly, "Spatial")
    sppoly$uid_internal = as.character( 1:nrow(sppoly) )
    cst = SpatialPoints( coords=areal_units_constraint, CRS(projection_proj4string("lonlat_wgs84")) )
    cst = spTransform( cst, sp::proj4string(sppoly) )
    oo = over( sppoly, cst )
    sppoly = sppoly[ which(!is.na(oo) ), ]
    if (exists("areal_units_constraint_nmin", p)) {
      vv = over( cst, sppoly  )
      ww = tapply( rep(1, nrow(vv)), vv$uid_internal, sum, na.rm=T )
      good = which(ww > p$areal_units_constraint_nmin )
      tags = names(ww)[good]
      sppoly = sppoly[ sppoly$uid_internal %in% tags, ]
    }
    sppoly$uid_internal = NULL
    sppoly = as(sppoly, "sf")  # revert to sf
  }

  # --------------------

  sppoly = st_transform( sppoly, sp::CRS( areal_units_proj4string_planar_km ))
  sppoly[, "au_sa_km2"] = st_area(sppoly)

  # plot(st_geometry(sppoly))
  # plot(sppoly[,"au_sa_km2"])

  units( sa_threshold_km2 ) = units( sppoly$au_sa_km2 )
  toremove = which( c(sppoly$au_sa_km2) < c(sa_threshold_km2))
  if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?


  if (!exists("AUID", sppoly)) {
    sppoly[, "AUID"] = as.character( 1:nrow(sppoly) )
    row.names(sppoly) = sppoly$AUID
  }


  if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {
    # as a last pass, calculate surface areas of each subregion .. could be done earlier but it is safer done here due to deletions above
    message("Computing surface areas for each subarea ... can be slow if this is your first time")
    for (subarea in c("cfanorth", "cfasouth", "cfa23", "cfa24", "cfa4x" ) ) {
      print(subarea)
      csa = polygons_managementarea( species="snowcrab", area=subarea )
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


  # force saves as Spatial* data
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
  save(sppoly, file=fn, compress=TRUE)
  if (plotit) plot(sppoly)

  return( sppoly )

}
