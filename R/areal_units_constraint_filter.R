
areal_units_constraint_filter = function( sppoly, areal_units_constraint_nmin, areal_units_type, areal_units_proj4string_planar_km, constraintdata=NULL, sa_threshold_km2=0 ) {


  if ( !is.null(constraintdata)) {
    sppoly$npts  = 0
    sppoly$internal_id = 1:nrow(sppoly)
    row.names( sppoly) = sppoly$internal_id


    sppoly = st_transform( sppoly, st_crs( areal_units_proj4string_planar_km ))
    sppoly$au_sa_km2 = st_area(sppoly)
    attributes( sppoly$au_sa_km2 ) = NULL

    constraintdata = sf::st_as_sf( constraintdata, coords = c("lon","lat"), crs=st_crs(projection_proj4string("lonlat_wgs84")) )
    constraintdata = st_transform( constraintdata, st_crs(areal_units_proj4string_planar_km ))
    # constraintdata = st_join( constraintdata, sppoly, join=st_within )
    constraintdata$internal_id = st_points_in_polygons( constraintdata, sppoly, varname="internal_id" )
    ww = tapply( rep(1, nrow(constraintdata)), constraintdata$internal_id, sum, na.rm=T )
    sppoly$npts[ as.numeric(names(ww))] = ww
    constraintdata = NULL

    zeros = which( sppoly$npts == 0 )
    if ( length(zeros) > 0 ) sppoly = sppoly[-zeros,]
    todrop = which( sppoly$npts < areal_units_constraint_nmin )

    if (length(todrop) > 0 ) {

      if (areal_units_type == "lattice" ) {
        # laatice structure is required, simply drop where there is no data
        sppoly = sppoly[ -todrop , ]

      } else {
        # already done in tessilation method so not really needed but
        # try to join to adjacent au's
        # for other methods, this ensures a min no samples in each AU
        sppoly$dropflag = FALSE
        sppoly$nok = TRUE
        sppoly$nok[todrop] = FALSE
        W.nb = poly2nb(sppoly, row.names=sppoly$internal_id, queen=TRUE)  
        for (i in order(sppoly$au_sa_km2) ) {
          if ( sppoly$nok[i]) next()
          v = setdiff( 
            intersect( W.nb[[ i ]], which(sppoly$nok) ), ## AU neighbours that are OK and so can consider dropping
            which(sppoly$dropflag)                       ## AUs confirmed already to drop
          )
          if (length(v) > 0) {
            j = v[ which.min( sppoly$npts[v] )]
            g_ij = try( st_union( st_geometry(sppoly)[j] , st_geometry(sppoly)[i] ) )
            if ( !inherits(g_ij, "try-error" )) {
              st_geometry(sppoly)[j] = g_ij
              sppoly$npts[j] = sppoly$npts[j] + sppoly$npts[i]
              sppoly$nok[i] = FALSE
              sppoly$dropflag[i] = TRUE
              if ( sppoly$npts[j] >= areal_units_constraint_nmin) sppoly$nok[j] = TRUE
            }
          }
        }
        sppoly = sppoly[ - which( sppoly$dropflag ), ]
        sppoly$nok =NULL
      }
      
      sppoly$internal_id = NULL
      message( "Dropping due to areal_units_constraint_nmin, now there are : ", nrow(sppoly), " areal units." )
    }

  }


  # SA check
  sppoly = st_transform( sppoly, st_crs( areal_units_proj4string_planar_km ))
  sppoly$au_sa_km2 = st_area(sppoly)
  attributes( sppoly$au_sa_km2 ) = NULL
  toremove = which( sppoly$au_sa_km2 < sa_threshold_km2 )
  if ( length(toremove) > 0 ) sppoly = sppoly[-toremove,]  # problematic as it is so small and there is no data there?
  
  if ( nrow( sppoly) == 0 ) {
    message ( "No polygons meet the specified criteria (check sa_threshold_km2 ?)." )
    browser()
  }

  message( "Applying constraints leaves: ",  nrow(sppoly), " areal units." )

  return(sppoly)
}
