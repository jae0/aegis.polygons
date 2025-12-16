
polygon_tabulate = function( plgn, pts ) {
    plgn$internal_id = 1:nrow( plgn)
    row.names( plgn) = plgn$internal_id
    cdd = setDT( st_drop_geometry( pts ) )
    cdd$internal_id = st_points_in_polygons( pts, plgn, varname="internal_id" )
    cdd = cdd[,.(npts=.N), .(internal_id) ]
    plgn$npts = cdd$npts[ match( as.character( plgn$internal_id),  as.character(cdd$internal_id) )]
    return( plgn$npts)
}
