areal_units_neighbourhood_reset = function ( sppoly, snap=1 ) {
  # after modifying polygons, make sure neighbourhoods are consistent 
  # expected aegis.polygons attributes  

  require(spdep)
  
  NB_graph = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE, snap=snap )  # slow .. ~1hr?
  NB_remove = which(card(NB_graph) == 0)
  if ( length(NB_remove) > 0 ) {
    # remove isolated locations and recreate sppoly .. alternatively add links to NB_graph
    NB_keep = which(card(NB_graph) > 0)
    NB_graph = nb_remove( NB_graph, NB_remove )
    sppoly = sppoly[NB_keep,]
    row.names(sppoly) = sppoly$AUID
    sppoly = sppoly[order(sppoly$AUID),]
    sppoly = st_make_valid(sppoly)
  }

  attr(nb, "region.id") = sppoly$AUID
  attr(sppoly, "nb") = INLA::inla.read.graph( spdep::nb2mat( NB_graph ))  # adding neighbourhood as an attribute to sppoly
  attr(sppoly, "NB_graph") = NB_graph  # adding neighbourhood as an attribute to sppoly

  return(sppoly)
}


