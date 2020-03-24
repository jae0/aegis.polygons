areal_units = function( p=NULL, areal_units_source="lattice", areal_units_resolution_km=20, aegis_internal_resolution_km=1, areal_units_fn=NULL,
  spatial_domain="SSE", areal_units_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  timeperiod="default", plotit=FALSE, areal_units_overlay="none", sa_threshold_km2=0, areal_units_constraint="none", redo=FALSE, use_stmv_solution=FALSE  ) {

  if (0) {
    areal_units_source="lattice"
    areal_units_resolution_km=20
    aegis_internal_resolution_km=1
    areal_units_fn=NULL
    spatial_domain="SSE"
    areal_units_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km"
    timeperiod="default"
    plotit=FALSE
    areal_units_overlay="none"
    sa_threshold_km2=0
    areal_units_constraint="none"
    redo=FALSE
    use_stmv_solution=FALSE

    spatial_domain=p$spatial_domain
    areal_units_proj4string_planar_km=p$areal_units_proj4string_planar_km
    areal_units_source=p$areal_units_source
    areal_units_resolution_km=p$areal_units_resolution_km
    areal_units_overlay= ifelse(exists("areal_units_overlay",p), p$areal_units_overlay, areal_units_overlay)
    areal_units_constraint=areal_units_constraint
    areal_units_fn=p$areal_units_fn

  }

  if (is.null(areal_units_fn)) {
    areal_units_fn = paste(
      spatial_domain,
      paste0(areal_units_overlay, collapse="_"),
      areal_units_resolution_km,
      areal_units_source,
      timeperiod,
      sep="_"
    )
  }

  if ( !is.null(p) ) {
    # designed for operating with aegis data-based polygons adding neighbourhood structure as an attribute .. using only p and override everything
    if (exists("spatial_domain", p)) spatial_domain=p$spatial_domain
    if (exists("areal_units_proj4string_planar_km", p)) areal_units_proj4string_planar_km=p$areal_units_proj4string_planar_km
    if (exists("areal_units_source", p)) areal_units_source=p$areal_units_source
    if (exists("areal_units_resolution_km", p)) areal_units_resolution_km=p$areal_units_resolution_km
    if (exists("areal_units_overlay", p)) areal_units_overlay=p$areal_units_overlay
    if (exists("areal_units_constraint", p)) areal_units_constraint=areal_units_constraint
    if (exists("areal_units_fn", p)) areal_units_fn=p$areal_units_fn
  }


  fn = file.path( project.datadirectory("aegis", "polygons", "areal_units" ), paste(areal_units_fn, "rdata", sep="." ) )
  sppoly = NULL

  if (!redo) {
    if (file.exists(fn)) load(fn)
    if( !is.null(sppoly) ) return(sppoly)
  }


  if (areal_units_source == "lattice") {
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
      pn = spatial_parameters( spatial_domain=spatial_domain )  # geeneric defaults
      Z = bathymetry.db( p=pn, DS="aggregated_data" )
      names(Z)[which(names(Z)=="z.mean" )] = "z"
      Z = lonlat2planar(Z, pn$aegis_proj4string_planar_km)  # should not be required but to make sure
      Z = geo_subset( spatial_domain=spatial_domain, Z=Z )

      spdf0 = SpatialPoints(Z[, c("plon", "plat")], proj4string=sp::CRS(pn$aegis_proj4string_planar_km) )
      require(raster)
      raster_template = raster(extent(spdf0)) # +1 to increase the area
      res(raster_template) = areal_units_resolution_km  # in units of crs (which should be in  km)
      crs(raster_template) = projection(spdf0) # transfer the coordinate system to the raster

      sppoly = rasterize( Z[, c("plon", "plat")], raster_template, field=Z$z )
      sppoly = as(sppoly, "SpatialPolygonsDataFrame")
      sppoly$AUID = as.character(1:length(sppoly))  # row index
      row.names(sppoly) = sppoly$AUID

    }


    if ( grepl("groundfish_strata", areal_units_overlay) ) {

      timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod )

      ### next section is the same as "stratanal polygons" below .. copied here to avoid recursion
      ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
      # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
      # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
      gf = maritimes_groundfish_strata( timeperiod=timeperiod, returntype="polygons" )
      # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
      crs0 = proj4string( gf )
      gf$AUID = as.character(gf$AUID)
      row.names(gf) = gf$AUID
      gf$au_sa_km2 = gArea(spTransform(gf, sp::CRS(projection_proj4string("utm20")) ), byid=TRUE) # /10^3/10^3 # km^2  .. planar_crs_km should be in km units
      gf = sp::spChFIDs( gf, row.names(gf) )  #fix id's
      gf = gf[order(gf$AUID),]
      gf = spTransform(gf, sp::proj4string(sppoly) )
      o = over( sppoly, gf ) # match each datum to an area
      sppoly = sppoly[ which(!is.na(o$AUID)), ]

      shp = as( gf, "sf" )
      shp = st_simplify(shp)
      shp = st_buffer(shp, aegis_internal_resolution_km)
      shp = st_union(shp)
      shp = st_simplify(shp)

      oo = st_intersection( shp, as( sppoly, "sf") )
      qq = spTransform( as( oo, "Spatial" ), sp::proj4string(sppoly) )
      row.names(qq) = row.names(sppoly)
      sppoly = SpatialPolygonsDataFrame( qq, data=sppoly@data, match.ID=TRUE )

    }

    if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {

      pn = spatial_parameters( spatial_domain=spatial_domain )
      Z = lonlat2planar(Z, pn$aegis_proj4string_planar_km)
      spdf0 = SpatialPoints(Z[, c("plon", "plat")], proj4string=sp::CRS(pn$aegis_proj4string_planar_km) )

      require(raster)

      raster_template = raster(extent(spdf0)) # +1 to increase the area
      res(raster_template) = pn$pres/10  # in units of crs (which should be in 1 km .. so 100m)
      crs(raster_template) = projection(spdf0) # transfer the coordinate system to the raster

      ZZ = rasterize( Z[, c("plon", "plat")], raster_template, field=Z$z )
      ZZ = as(ZZ, "SpatialPolygonsDataFrame")
      ZZ = as(ZZ, "sf")
      ZZ = st_buffer(ZZ, aegis_internal_resolution_km)
      ZZ = st_simplify(ZZ)
      ZZ = st_union( ZZ)
      ZZ = st_simplify(ZZ)

      shp = polygons_managementarea( species="snowcrab", area="cfaall")
      shp = spTransform(shp, sp::proj4string(sppoly) )
      shp = as( shp, "sf" )
      shp = st_simplify(shp)
      shp = st_buffer(shp, aegis_internal_resolution_km)
      shp = st_union(shp)
      shp = st_simplify(shp)

      # domain = sf::st_join(shp, ZZ, join = st_intersects)
      # domain = sf::st_join(shp, as( sppoly, "sf"), join = st_intersects)

      domain = st_intersection( shp, ZZ )
      domain = st_union(domain)
      domain = st_simplify(domain)

      domain = spTransform(as(domain, "Spatial"), sp::proj4string(sppoly) )

      sppoly = raster::intersect( domain, sppoly  )

      sppoly = sppoly[order( sppoly$AUID ),]
      row.names(sppoly) = sppoly$AUID
      attr(sppoly, "region.id") =  sppoly$AUID
      sppoly = sp::spChFIDs( sppoly, sppoly$AUID  ) #fix id's

    }

    if (class( areal_units_constraint ) %in% c("data.frame", "matrix") ) {
      cst = SpatialPoints( coords=areal_units_constraint, CRS(projection_proj4string("lonlat_wgs84")) )
      cst = spTransform( cst, sp::proj4string(sppoly) )
      oo = over( sppoly, cst )
      sppoly = sppoly[ which(!is.na(oo) ), ]
      if (exists("areal_units_constraint_nmin", p)) {
        vv = over( cst, sppoly  )
        ww = tapply( rep(1, nrow(vv)), vv$AUID, sum, na.rm=T )
        good = which(ww > p$areal_units_constraint_nmin )
        tags = names(ww)[good]
        sppoly = sppoly[ sppoly$AUID %in% tags, ]
      }
    }

    sppoly$au_sa_km2 = gArea(sppoly, byid=TRUE)
    toremove = which(sppoly$au_sa_km2 < sa_threshold_km2)
    if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?

    sppoly = sppoly[order(sppoly$AUID),]
    row.names(sppoly) = sppoly$AUID
    attr(sppoly, "region.id") = sppoly$AUID
    sppoly = sp::spChFIDs( sppoly, sppoly$AUID ) #fix id's


    if ( grepl("snowcrab_managementareas", areal_units_overlay) ) {
      # as a last pass, calculate surface areas of each subregion .. could be done earlier but it is safer done here due to deletions above
      message("Computing surface areas for each subarea ... can be slow if this is your first time")
      csa_all = as(sppoly, "sf")
      for (subarea in c("cfanorth", "cfasouth", "cfa23", "cfa24", "cfa4x" ) ) {
        print(subarea)
        csa = polygons_managementarea( species="snowcrab", area=subarea )
        csa = spTransform(csa, sp::proj4string(sppoly) )
        csa = as(csa, "sf")
        ooo = st_intersection( csa, csa_all)
        ooo$surfacearea = st_area( ooo )
        vn = paste(subarea, "surfacearea", sep="_")
        sppoly[[ vn ]] = 0
        j = match( ooo$AUID, sppoly$AUID )
        if (length(j) > 0)  sppoly[[ vn ]][j] = ooo$surfacearea
      }
    }
  }

  # ------------------------------------------------

  if (areal_units_source == "stratanal_polygons") {
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    sppoly = maritimes_groundfish_strata( timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod ), returntype="polygons" )
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
    crs0 = proj4string( sppoly )
    sppoly$AUID = as.character(sppoly$AUID)
    row.names(sppoly) = sppoly$AUID
    sppoly$au_sa_km2 = gArea(spTransform(sppoly, sp::CRS(projection_proj4string("utm20")) ), byid=TRUE) # /10^3/10^3 # km^2  .. planar_crs_km should be in km units
    sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
    sppoly = sppoly[order(sppoly$AUID),]
    sppoly = spTransform(sppoly, sp::CRS(crs0) )  # for plot in UTM ccordinates .. default isobath and coastlines
  }


  if (areal_units_source %in% c("groundfish_polygons_inla_mesh",  "groundfish_polygons_tesselation") ) {

    ## method 1: Voronoi triangulation
    spset = survey.db( DS="set.base", p=p )
    spset = lonlat2planar(spset, areal_units_proj4string_planar_km)  # should not be required but to make sure
    spset = geo_subset( spatial_domain=spatial_domain, Z=spset )

    spset$AUID = spset$id
    coordinates(spset) = ~ lon+lat
    rownames(spset@data) = spset$AUID
    sp::proj4string(spset) = projection_proj4string("lonlat_wgs84")

    # crs_planar = sp::CRS(projection_proj4string("utm20"))
    crs_planar = sp::CRS( areal_units_proj4string_planar_km )
    spset = spTransform( spset, crs_planar )  # in km
    locs = coordinates( spset )
    locs = locs + runif( nrow(locs)*2, min=-1e-3, max=1e-3 ) # add  noise  to prevent a race condition

    spbuffer = 5
    hull_multiplier = 6
    bnd = non_convex_hull( locs, alpha=spbuffer*hull_multiplier  )
    bnd = list( Polygons(list( Polygon( as.matrix( bnd ) )), ID="boundary" ))
    bnd = SpatialPolygons( bnd, proj4string=sp::CRS(projection_proj4string("utm20")) )
    bnd = gBuffer( gUnaryUnion( gBuffer( bnd, width=spbuffer, byid=TRUE) ), width=spbuffer)
          # plot(bnd)

    if (!exists("areal_units_resolution_km", p)) areal_units_resolution_km = max( diff(range( locs[,1])), diff(range( locs[,2]) )) / 25  # in absence of range estimate take 1/10 of domain size

    if (areal_units_source == "groundfish_polygons_tesselation" ) {

      spmesh = aegis_mesh( SPDF=spset, resolution=areal_units_resolution_km, output_type="grid.count" )  # coarse grid representation
      spmesh = aegis_mesh( SPDF=spmesh, resolution=areal_units_resolution_km, output_type="polygons" )  # voroni tesslation and delaunay triagulation

    }

    if (areal_units_source == "groundfish_polygons_inla_mesh" ) {

      require(INLA)
      spmesh = inla.mesh.2d (
        loc=locs,
        max.edge = c( 0.5, 5 ) * areal_units_resolution_km,  #   # max size of a triange (in, out)
        offset = c( 0.1, 1 ) * areal_units_resolution_km , # how much to extend inside and outside of boundary,
        cutoff = c( 0.5, 5 ) * areal_units_resolution_km # min distance allowed between points #,
        # boundary =  inla.mesh.segment(st_coordinates( as(bnd, "sf") )[,c(1,2)])
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
        proj4string = crs_planar
        ),
        data = as.data.frame(spmesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
        match.ID = FALSE
      )
    }

    require(aegis.coastline)
    coast = (
        as( coastline.db( p=p, DS="eastcoast_gadm" ), "sf")
        %>% st_transform( crs_planar)
        %>% st_simplify()
        %>% st_buffer(0.1)
        %>% st_union()
    )

    sppoly = (
      as(bnd, "sf")
      %>% st_union()
      %>% st_buffer(0.1)
      %>% st_difference( coast)
      %>% st_cast( "MULTIPOLYGON" )
      %>% st_cast( "POLYGON" )
    )

    # must be done separately
    sppoly = (
      st_intersection( as( spmesh, "sf"), sppoly )
      %>% st_cast( "MULTIPOLYGON" )
      %>% st_cast( "POLYGON" )
    )

    sppoly[, "au_sa_km2"] = st_area(sppoly)
    # plot(st_geometry(sppoly))
    # plot(sppoly[,"au_sa_km2"])
    units( sa_threshold_km2 ) = units( sppoly$au_sa_km2 )
    toremove = which(sppoly$au_sa_km2 < sa_threshold_km2)
    if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?
    sppoly[, "AUID"] = as.character( 1:nrow(sppoly) )
    row.names(sppoly) = sppoly$AUID
    sppoly = as(sppoly, "Spatial")

  }



  # ------------------------------------------------

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
