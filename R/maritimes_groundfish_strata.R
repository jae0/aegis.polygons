
maritimes_groundfish_strata = function( W.nb=NULL, areal_units_timeperiod="pre2014", shapefilelocation=NULL, returntype="polygons") {
  # create a SpatialPolygons of the groundfish strata used in Martimes Region of Canada
  require(maptools)
  require(rgdal)

  # if (is.null(shapefilelocation)) shapefilelocation=system.file("data", package="carstm", mustWork=TRUE)
  if (is.null(shapefilelocation)) shapefilelocation=project.datadirectory("aegis", "polygons", "Management_Areas", "Fisheries", "Groundfish")

  if (returntype=="polygons") {

    if (areal_units_timeperiod=="pre2014") {
      shapefilename = "GB_STRATA_VDC"
      crswgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      o1 = rgdal::readOGR( dsn=shapefilelocation, layer=shapefilename, p4s=crswgs84 )
      names(o1) = "AUID"

      todrop = which(o1$AUID %in% c("5Z1", "5Z2") ) # drop as they are also in o2
      o1 = o1[-todrop,]

      shapefilename = "GF_SUMMER_STRATA_VDC"
      crswgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      o2 = rgdal::readOGR( dsn=shapefilelocation, layer=shapefilename, p4s=crswgs84 )
      names(o2) = "AUID"

      groundfish_strata = rbind(o1, o2)
      class(groundfish_strata) = class(o1)

      # i = which( ! st_is_valid(groundfish_strata) )
      # 38 overlaps with 40
      require(sf)
      groundfish_strata = as( groundfish_strata, "sf")
      groundfish_strata = st_make_valid(groundfish_strata)
      st_geometry(groundfish_strata[40,]) = st_difference( st_geometry(groundfish_strata[40,]), st_geometry(groundfish_strata[38,]))
      groundfish_strata = st_make_valid(groundfish_strata)
      groundfish_strata = as( groundfish_strata, "Spatial")

      attr(groundfish_strata, "region.id") = as.character( slot(groundfish_strata, "data")[,"AUID"] )

      # plot(groundfish_strata)
      return(groundfish_strata)
    }

    if (areal_units_timeperiod=="post2014") {
      shapefilename = "MaritimesRegionEcosystemAssessmentStrata"
      crswgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      groundfish_strata = rgdal::readOGR( dsn=shapefilelocation, layer=shapefilename, p4s=crswgs84 )
      names(groundfish_strata) = "AUID"

      require(sf)
      groundfish_strata = as( groundfish_strata, "sf")
      groundfish_strata = st_make_valid(groundfish_strata)
      # i = which( ! st_is_valid(groundfish_strata) )
      # groundfish_strata = st_make_valid(groundfish_strata)
      groundfish_strata = as( groundfish_strata, "Spatial")

      # plot(groundfish_strata)
      attr(groundfish_strata, "region.id") = as.character( slot(groundfish_strata, "data")[,"AUID"] )

      return( groundfish_strata)
    }
  }

  if (returntype=="neighbourhoods") {

    if (areal_units_timeperiod=="pre2014") {
      if (0){
        # polygons seem malformed .. manually modify connectivity
        isolated = which( card(W.nb) == 0)
        plot(sppoly)
        plot(sppoly[isolated,], add=T, col="red")
        dev.new()
        edit(W.nb, polys=sppoly)
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
      sppoly = maritimes_groundfish_strata( areal_units_timeperiod=areal_units_timeperiod, returntype="polygons" )
      W.nb = spdep::poly2nb(sppoly, row.names=rownames(slot(sppoly, "data")[,"AUID"] ) )
      if (0){
        # polys look ok
        plot(sppoly)
        dev.new()
        edit(W.nb, polys=sppoly)
        card(W.nb) # last check if any more  isolated areas
      }
      return( W.nb )
    }
  }



}
