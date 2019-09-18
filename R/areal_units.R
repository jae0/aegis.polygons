areal_units = function( strata_type="lattice", resolution=20, resolution_aegis_internal=1,
  spatial.domain="SSE", proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km", proj4string_planar_km_aegis="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  timeperiod="default", plotit=FALSE, overlay="groundfish_strata", sa_threshold_km2=0, constraint="none", redo=FALSE  ) {

  fn = file.path( project.datadirectory("aegis", "polygons", "areal_units" ), paste( strata_type, overlay, spatial.domain, resolution, timeperiod, "rdata", sep="." ) )
  sppoly = NULL

  if (!redo) {
    if (file.exists(fn)) load(fn)
    if( !is.null(sppoly) ) return(sppoly)
  }

  if (strata_type == "lattice") {
    # res based on grids ... rather than arbitrary polygons
    # static features only so far
    # resolution = 20 # in units of crs (km)
    sppoly = aegis_db_spatial_object( spatial.domain=spatial.domain, proj4string=proj4string_planar_km, resolution=resolution, returntype="SpatialPolygonsDataFrame")
    sppoly$StrataID = as.character(sppoly$StrataID)


    if (overlay=="groundfish_strata") {
      gf =  areal_units( strata_type="stratanal_polygons", proj4string_planar_km=proj4string_planar_km, timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod ) )
      gf = spTransform(gf, sp::proj4string(sppoly) )
      o = over( sppoly, gf ) # match each datum to an area
      sppoly = sppoly[ which(!is.na(o$StrataID)), ]

      shp = as( gf, "sf" )
      shp = st_simplify(shp)
      shp = st_buffer(shp, resolution_aegis_internal)
      shp = st_union(shp)
      shp = st_simplify(shp)

      oo = st_intersection( shp, as( sppoly, "sf") )
      qq = spTransform( as( oo, "Spatial" ), sp::proj4string(sppoly) )
      row.names(qq) = row.names(sppoly)
      sppoly = SpatialPolygonsDataFrame( qq, data=sppoly@data, match.ID=TRUE )

    }

    if (overlay=="snowcrab") {
      cfaall = polygons_managementarea( species="snowcrab", area="cfaall")
      cfaall = spTransform(cfaall, sp::proj4string(sppoly) )
      o = over( sppoly, cfaall ) # match each datum to an area
      sppoly = sppoly[ which(!is.na(o)), ]

      shp = as( cfaall, "sf" )
      shp = st_simplify(shp)
      shp = st_buffer(shp, resolution_aegis_internal)
      shp = st_union(shp)
      shp = st_simplify(shp)

      oo = st_intersection( shp, as( sppoly, "sf") )
      qq = spTransform( as( oo, "Spatial" ), sp::proj4string(sppoly) )
      row.names(qq) = row.names(sppoly)
      sppoly = SpatialPolygonsDataFrame( qq, data=sppoly@data, match.ID=TRUE )

    }

    if (constraint != "none" ) {
      cst = SpatialPoints( coords=constraint, CRS("+proj=longlat +datum=WGS84") )
      cst = spTransform( cst, sp::proj4string(sppoly) )
      oo = over( sppoly, cst )
      sppoly = sppoly[ which(!is.na(oo) ), ]
    }

    sppoly$sa_strata_km2 = gArea(sppoly, byid=TRUE)
    toremove = which(sppoly$sa_strata_km2 < sa_threshold_km2)
    if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?

    sppoly$StrataID = factor( as.character(sppoly$StrataID) )
    sppoly$strata = as.numeric( sppoly$StrataID )
    sppoly = sppoly[order(sppoly$strata),]
    row.names(sppoly) = as.character(sppoly$StrataID)
    attr(sppoly, "region.id") = as.character( sppoly@data$StrataID )
    sppoly = sp::spChFIDs( sppoly, as.character(sppoly$StrataID) ) #fix id's


    if (overlay=="snowcrab") {
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
        j = match( ooo$StrataID, sppoly$StrataID )
        if (length(j) > 0)  sppoly[[ vn ]][j] = ooo$surfacearea
      }

    }


    if (plotit) plot(sppoly)
    save( sppoly, file=fn, compress=TRUE)
    return( sppoly )
  }

  # ------------------------------------------------

  if (strata_type == "stratanal_polygons") {
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
    # Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
    # This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)
    sppoly = maritimes_groundfish_strata( timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod ), returntype="polygons" )
    # prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
    crs0 = proj4string( sppoly )
    row.names(sppoly) = as.character(sppoly$StrataID)
    sppoly$sa_strata_km2 = gArea(spTransform(sppoly, sp::CRS("+proj=utm +ellps=WGS84 +zone=20 +units=km") ), byid=TRUE) # /10^3/10^3 # km^2  .. planar_crs_km should be in km units
    sppoly$StrataID = factor(sppoly$StrataID, levels=levels(sppoly$StrataID) ) # make sure in case not all strata are represented in set
    sppoly$strata = as.numeric( sppoly$StrataID )
    sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
    sppoly = sppoly[order(sppoly$strata),]
    sppoly = spTransform(sppoly, sp::CRS(crs0) )  # for plot in UTM ccordinates .. default isobath and coastlines

    if (plotit) plot(sppoly)
    save( sppoly, file=fn, compress=TRUE)
    return( sppoly )
  }


}
