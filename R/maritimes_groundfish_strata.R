
maritimes_groundfish_strata = function( W.nb=NULL, areal_units_timeperiod="pre2014", shapefilelocation=NULL, returntype="polygons") {
  # create a sf SpatialPolygons of the groundfish strata used in Martimes Region of Canada
  require(maptools)
  require(rgdal)
  require(sf)

  # if (is.null(shapefilelocation)) shapefilelocation=system.file("data", package="carstm", mustWork=TRUE)
  if (is.null(shapefilelocation)) shapefilelocation=project.datadirectory("aegis", "polygons", "Management_Areas", "Fisheries", "Groundfish")

  if (returntype=="polygons") {

    if (areal_units_timeperiod=="pre2014") {
      shapefilename = "GB_STRATA_VDC"
      crswgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      o1 = st_read( dsn=shapefilelocation, layer=shapefilename, crs=crswgs84)
      names(o1) = c( "AUID", "geometry")

      todrop = which(o1$AUID %in% c("5Z1", "5Z2") ) # drop as they are also in o2
      o1 = o1[-todrop,]

      shapefilename = "GF_SUMMER_STRATA_VDC"
      crswgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      o2 = st_read( dsn=shapefilelocation, layer=shapefilename, crs=crswgs84 )
      names(o2) = c( "AUID", "geometry" )

      groundfish_strata = rbind(o1, o2)

      # i = which( ! st_is_valid(groundfish_strata) )
      groundfish_strata = st_make_valid(groundfish_strata)
      st_geometry(groundfish_strata[40,]) = st_difference( st_geometry(groundfish_strata[40,]), st_geometry(groundfish_strata[38,]))
      groundfish_strata = st_make_valid(groundfish_strata)

      attr(groundfish_strata, "space.id") = as.character( groundfish_strata[,"AUID"] )

      if (0) {
        plot(st_geometry( groundfish_strata) )
        plot(st_geometry(groundfish_strata[38,]), col="blue", add=T)
        plot(st_geometry(groundfish_strata[40,]), col="red", add=T)
      }

      return(groundfish_strata)
    }

    if (areal_units_timeperiod=="post2014") {
      shapefilename = "MaritimesRegionEcosystemAssessmentStrata"
      crswgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      groundfish_strata = st_read( dsn=shapefilelocation, layer=shapefilename, crs=crswgs84 )
      names(groundfish_strata) = c( "AUID", "geometry")

      message( "Strata 70, 71, 72 in Misaine Banks are not contiguous .. this can be problematic")

      # i = which( ! st_is_valid(groundfish_strata) )
      groundfish_strata = st_make_valid(groundfish_strata)
      attr(groundfish_strata, "space.id") = as.character( groundfish_strata[,"AUID"] )

      return( groundfish_strata)
    }
  }

  if (returntype=="neighbourhoods") {

    if (is.null(W.nb)) {
      sppoly = maritimes_groundfish_strata( areal_units_timeperiod=areal_units_timeperiod, returntype="polygons" )
      W.nb = spdep::poly2nb(sppoly, row.names=sppoly[,"AUID"], queen=TRUE )
    }

    if (areal_units_timeperiod=="pre2014") {
      if (0){
        # polygons seem malformed .. manually modify connectivity
        isolated = which( card(W.nb) == 0)
        plot(sppoly)
        plot(sppoly[isolated,], add=T, col="red")
        dev.new()
        edit(W.nb, polys=as(sppoly, "Spatial"))
        card(W.nb) # last check if any more  isolated areas
      }

      # these are the results from the above ..
      # hand craft connectivity  .. these should be looked at more carefully a bit ..
      polyname = as.character( attr(W.nb, "region.id" ) )

      W.nb = nb_add( W.nb, which(polyname=="5Z1"), which(polyname=="478"))
      W.nb = nb_add( W.nb, which(polyname=="5Z1"), which(polyname=="482"))
      W.nb = nb_add( W.nb, which(polyname=="5Z1"), which(polyname=="483"))
      W.nb = nb_add( W.nb, which(polyname=="5Z1"), which(polyname=="5Z8"))
      W.nb = nb_add( W.nb, which(polyname=="5Z2"), which(polyname=="5Z3"))
      W.nb = nb_add( W.nb, which(polyname=="5Z2"), which(polyname=="483"))
      W.nb = nb_add( W.nb, which(polyname=="5Z3"), which(polyname=="5Z5"))
      W.nb = nb_add( W.nb, which(polyname=="5Z3"), which(polyname=="483"))
      W.nb = nb_add( W.nb, which(polyname=="5Z4"), which(polyname=="5Z1"))
      W.nb = nb_add( W.nb, which(polyname=="5Z4"), which(polyname=="5Z2"))
      W.nb = nb_add( W.nb, which(polyname=="5Z4"), which(polyname=="5Z5"))
      W.nb = nb_add( W.nb, which(polyname=="5Z4"), which(polyname=="5Z7"))
      W.nb = nb_add( W.nb, which(polyname=="474"), which(polyname=="476"))
      W.nb = nb_add( W.nb, which(polyname=="475"), which(polyname=="476"))
      W.nb = nb_add( W.nb, which(polyname=="476"), which(polyname=="481"))
      W.nb = nb_add( W.nb, which(polyname=="477"), which(polyname=="481"))
      W.nb = nb_add( W.nb, which(polyname=="478"), which(polyname=="482"))
      W.nb = nb_add( W.nb, which(polyname=="480"), which(polyname=="481"))
      W.nb = nb_add( W.nb, which(polyname=="480"), which(polyname=="482"))
      W.nb = nb_add( W.nb, which(polyname=="481"), which(polyname=="482"))
      W.nb = nb_add( W.nb, which(polyname=="481"), which(polyname=="484"))
      W.nb = nb_add( W.nb, which(polyname=="481"), which(polyname=="485"))
      return( W.nb )
    }

    if (areal_units_timeperiod=="post2014") {
      if (0){
        # polys look ok
        plot(sppoly)
        dev.new()
        edit(W.nb, polys=as(sppoly, "Spatial") )
        card(W.nb) # last check if any more  isolated areas
      }
      return( W.nb )
    }
  }



}
