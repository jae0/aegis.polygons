
  maritimes_groundfish_strata_identify = function( Y, sppoly=NULL, xyvars=c("lon", "lat"),  planar_crs_km="+proj=utm +ellps=WGS84 +zone=20 +units=km"){
    ## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
      # "pre2014" for older

      message( "NOTE: some groundfish strata seem to be either miscoded or set positions are incorrect in some cases ... ")
      message( "NOTE: fixing this will require some effort in looking at paper records! ")
      message( "NOTE: e.g.,: strata 455 453 451 484 on NED2012022.224, NED2016016.241, NED2015017.215, NED2017020.70 trip.set") #

          # Load file from groundfish.stratum where area is converted to trawlable units divide area by 0.011801 (41ft by 1.75 nm)
          # st = aegis.survey::groundfish_survey_db(DS="gsstratum")
          # names(st) = paste("stratum", names(st), sep=".")
          # set = merge(all.sets, st, by.x="strat", by.y="stratum.strat", all.x=TRUE, all.y=FALSE, suffixes=c("", ".gsstratum") )



    crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
    sppoly = st_transform(sppoly, crs=crs_lonlat )

    Y$AUID = st_points_in_polygons(
      pts = st_as_sf( Y[,xyvars], coords=c("lon","lat"), crs=crs_lonlat ),
      polys = sppoly[, "AUID"],
      varname="AUID"
    )

    sppoly = st_transform(sppoly, st_crs(planar_crs_km) )  # "+proj=utm +ellps=WGS84 +zone=20 +units=km"
    sppoly$au_sa_km2 = st_area(sppoly) # km^2  .. planar_crs_km should be in km units

    Y = merge(Y, st_drop_geometry(sppoly), all.x=TRUE, all.y=FALSE)

    st =  groundfish_survey_db(DS="gsstratum")
    names(st) = paste("stratum", names(st), sep=".")
    st$stratum.dmin.fathoms = as.numeric(st$stratum.dmin) #units? 0.54680665 fathoms = 1 meter
    st$stratum.dmax.fathoms = as.numeric(st$stratum.dmax) #units?
    st$stratum.area.nauticalmiles = as.numeric(st$stratum.area)  # stratum.area is in sq nautical miles
    st$stratum.depth.fathoms = as.numeric(st$stratum.depth)  # units ?
    st$stratum.dmin = NULL
    st$stratum.dmax = NULL
    st$stratum.area = NULL
    st$stratum.depth = NULL


    Y = merge(Y, st, by.x="AUID", by.y="stratum.strat", all.x=TRUE, all.y=FALSE, suffixes=c("", ".gsstratum") )

    return(Y)
  }
