areal_units_neighbourhood_reset = function ( sppoly, snap=1 ) {
  # after modifying polygons, make sure neighbourhoods are consistent 
  # expected aegis.polygons attributes  

  require(spdep)
  
  W.nb = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE, snap=snap )  # slow .. ~1hr?
  W.remove = which(card(W.nb) == 0)
  if ( length(W.remove) > 0 ) {
    # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
    W.keep = which(card(W.nb) > 0)
    W.nb = nb_remove( W.nb, W.remove )
    sppoly = sppoly[W.keep,]
    row.names(sppoly) = sppoly$AUID
    sppoly = sppoly[order(sppoly$AUID),]
    sppoly = st_make_valid(sppoly)
  }

  attr(nb, "region.id") = sppoly$AUID
  attr(sppoly, "nb") = INLA::inla.read.graph( spdep::nb2mat( W.nb ))  # adding neighbourhood as an attribute to sppoly
  attr(sppoly, "W.nb") = W.nb  # adding neighbourhood as an attribute to sppoly

  return(sppoly)
}


