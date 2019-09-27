areal_units = function( p=NULL, areal_units_strata_type="lattice", areal_units_resolution_km=20, aegis_internal_resolution_km=1, auid=NULL,
  spatial_domain="SSE", areal_units_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  timeperiod="default", plotit=FALSE, areal_units_overlay="none", sa_threshold_km2=0, areal_units_constraint="none", redo=FALSE  ) {

  if (is.null(auid)) auid = paste( spatial_domain, paste0(areal_units_overlay, collapse="_"), areal_units_resolution_km,     areal_units_strata_type,  areal_units_constraint, timeperiod, "rdata", sep="." )

  fn = file.path( project.datadirectory("aegis", "polygons", "areal_units" ), auid )
  sppoly = NULL

  if (!redo) {
    if (file.exists(fn)) load(fn)
    if( !is.null(sppoly) ) return(sppoly)
  }

  if ( !is.null(p) ) {
    # designed for operating with aegis data-based polygons adding neighbourhood structure as an attribute .. using only p
    sppoly = areal_units(
      spatial_domain=p$spatial_domain,
      areal_units_proj4string_planar_km=p$areal_units_proj4string_planar_km,
      areal_units_strata_type=p$areal_units_strata_type,
      areal_units_resolution_km=p$areal_units_resolution_km,
      areal_units_overlay= ifelse(!exists("areal_units_overlay", p), "none", p$areal_units_overlay),
      areal_units_constraint=ifelse(!exists("areal_units_constraint", p), "none", p$areal_units_constraint),
      redo=redo
    )

    W.nb = poly2nb(sppoly, row.names=sppoly$StrataID, queen=TRUE)  # slow .. ~1hr?
    W.remove = which(card(W.nb) == 0)

    if ( length(W.remove) > 0 ) {
      # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
      W.keep = which(card(W.nb) > 0)
      W.nb = nb_remove( W.nb, W.remove )
      sppoly = sppoly[W.keep,]
      row.names(sppoly) = as.character(sppoly$StrataID)
      sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
      sppoly$StrataID = factor( as.character(sppoly$StrataID) )
      sppoly$strata = as.numeric( sppoly$StrataID )
      sppoly = sppoly[order(sppoly$strata),]
    }

    attr(sppoly, "nb") = W.nb  # adding neighbourhood as an attribute to sppoly
    save(sppoly, file=fn, compress=TRUE)
    return( sppoly )
  }


  if (areal_units_strata_type == "lattice") {
    # res based on grids ... rather than arbitrary polygons
    # static features only so far

    sppoly = aegis_db_spatial_object( spatial_domain=spatial_domain, proj4string=areal_units_proj4string_planar_km, areal_units_resolution_km=areal_units_resolution_km, returntype="SpatialPolygonsDataFrame")
    sppoly$StrataID = as.character(sppoly$StrataID)


      if ( grepl("groundfish_strata", areal_units_overlay) ) {
        gf =  areal_units( areal_units_strata_type="stratanal_polygons", areal_units_proj4string_planar_km=areal_units_proj4string_planar_km, timeperiod=ifelse( timeperiod=="default", "pre2014", timeperiod ) )
        gf = spTransform(gf, sp::proj4string(sppoly) )
        o = over( sppoly, gf ) # match each datum to an area
        sppoly = sppoly[ which(!is.na(o$StrataID)), ]

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

        cfaall = polygons_managementarea( species="snowcrab", area="cfaall")
        cfaall = spTransform(cfaall, sp::proj4string(sppoly) )
        o = over( sppoly, cfaall ) # match each datum to an area
        sppoly = sppoly[ which(!is.na(o)), ]

        shp = as( cfaall, "sf" )
        shp = st_simplify(shp)
        shp = st_buffer(shp, aegis_internal_resolution_km)
        shp = st_union(shp)
        shp = st_simplify(shp)

        oo = st_intersection( shp, as( sppoly, "sf") )
        qq = spTransform( as( oo, "Spatial" ), sp::proj4string(sppoly) )
        row.names(qq) = row.names(sppoly)
        sppoly = SpatialPolygonsDataFrame( qq, data=sppoly@data, match.ID=TRUE )

      }


    if (areal_units_constraint != "none" ) {
      cst = SpatialPoints( coords=areal_units_constraint, CRS("+proj=longlat +datum=WGS84") )
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
        j = match( ooo$StrataID, sppoly$StrataID )
        if (length(j) > 0)  sppoly[[ vn ]][j] = ooo$surfacearea
      }

    }


    if (plotit) plot(sppoly)
    save( sppoly, file=fn, compress=TRUE)
    return( sppoly )
  }

  # ------------------------------------------------

  if (areal_units_strata_type == "stratanal_polygons") {
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
