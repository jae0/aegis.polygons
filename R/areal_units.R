areal_units = function( p=NULL, areal_units_source="lattice", areal_units_resolution_km=20, aegis_internal_resolution_km=1, areal_unit_type=NULL,
  spatial_domain="SSE", areal_units_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  timeperiod="default", plotit=FALSE, areal_units_overlay="none", sa_threshold_km2=0, areal_units_constraint="none", redo=FALSE, use_stmv_solution=FALSE  ) {

  if (0) {
    areal_units_source="lattice"
    areal_units_resolution_km=20
    aegis_internal_resolution_km=1
    areal_unit_type=NULL
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
    areal_unit_type=p$areal_unit_type

  }

  if (is.null(areal_unit_type)) {
    areal_unit_type = paste(
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
    if (exists("areal_unit_type", p)) areal_unit_type=p$areal_unit_type
  }


  fn = file.path( project.datadirectory("aegis", "polygons", "areal_units" ), paste(areal_unit_type, "rdata", sep="." ) )
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
      gf =  areal_units( areal_units_source="stratanal_polygons", areal_units_proj4string_planar_km=areal_units_proj4string_planar_km, timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod ) )
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

      sppoly = intersect( domain, sppoly )

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
