
  maritimes_groundfish_strata_identify = function( Y, sppoly=NULL, xyvars=c("lon", "lat"),  plotdata=FALSE, planar_crs_km="+proj=utm +ellps=WGS84 +zone=20 +units=km"){
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
      # "pre2014" for older

      message( "NOTE: some groundfish strata seem to be either miscoded or set positions are incorrect in some cases ... ")
      message( "NOTE: fixing this will require some effort in looking at paper records! ")
      message( "NOTE: e.g.,: strata 455 453 451 484 on NED2012022.224, NED2016016.241, NED2015017.215, NED2017020.70 trip.set") #

          # Load file from groundfish.stratum where area is converted to trawlable units divide area by 0.011801 (41ft by 1.75 nm)
          # st = aegis.survey::groundfish.db(DS="gsstratum")
          # names(st) = paste("stratum", names(st), sep=".")
          # set = merge(all.sets, st, by.x="strat", by.y="stratum.strat", all.x=TRUE, all.y=FALSE, suffixes=c("", ".gsstratum") )


    pts <- SpatialPoints(Y[,xyvars], proj4string=CRS(proj4string( sppoly)) )

    # match each datum to an area
    o = over( pts, sppoly)
    Y$StrataID = as.character( o$StrataID )

    sppoly = spTransform(sppoly, sp::CRS(planar_crs_km) )  # "+proj=utm +ellps=WGS84 +zone=20 +units=km"
    sppoly$sa_strata_km2 = gArea(sppoly, byid=TRUE) # km^2  .. planar_crs_km should be in km units

    Y = merge(Y, slot(sppoly, "data"), all.x=TRUE, all.y=FALSE)

    st = aegis.survey::groundfish.db(DS="gsstratum")
    names(st) = paste("stratum", names(st), sep=".")
    st$stratum.dmin.fathoms = as.numeric(st$stratum.dmin) #units? 0.54680665 fathoms = 1 meter
    st$stratum.dmax.fathoms = as.numeric(st$stratum.dmax) #units?
    st$stratum.area.nauticalmiles = as.numeric(st$stratum.area)  # stratum.area is in sq nautical miles
    st$stratum.depth.fathoms = as.numeric(st$stratum.depth)  # units ?
    st$stratum.dmin = NULL
    st$stratum.dmax = NULL
    st$stratum.area = NULL
    st$stratum.depth = NULL


    Y = merge(Y, st, by.x="StrataID", by.y="stratum.strat", all.x=TRUE, all.y=FALSE, suffixes=c("", ".gsstratum") )

    if (plotdata) {
      nomatches = which(is.na(Y$StrataID))
      plot(sppoly)
      plot( pts, pch=20, col="green", add=TRUE)
      if (length(nomatches) >0) {
        pts2 <- SpatialPoints(Y[nomatches, xyvars], proj4string=CRS(proj4string( sppoly)) )
        plot( pts2, col="red", add=TRUE )
      }
    }
    return(Y)
  }
