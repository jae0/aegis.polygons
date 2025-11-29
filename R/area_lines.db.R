
area_lines.db = function( DS, returntype="list", project_to=st_crs( projection_proj4string("lonlat_wgs84") ) ) {

  if (DS=="cfa.regions") {
      
      nafo4Vn  = data.frame( rbind(
       
        
        c(-59.866046,47.717237),    
        c(-60.40917,47.035833)
      ))

      cfa.nens.23 = data.frame( rbind(
        c(-59.85, 46),
        c(-58.40, 46)
      ))
      cfa.23.24 = data.frame( rbind(
        c(-59.065, 43.5),
        c(-59.959007, 44.829624),
        c(-60.51667, 45.61667)
      ))
      cfa.4x.24 = data.frame( rbind(
        c( -63.333333, 42.61379),
        c( -63.333333,	44.332904),
        c( -63.50242, 44.502358)
      ))

      names( nafo4Vn ) =   c("lon", "lat")
      names( cfa.nens.23 ) =   c("lon", "lat")
      names( cfa.23.24 ) =  c("lon", "lat") 
      names( cfa.4x.24 ) = c("lon", "lat")

      out = list(
        nafo4Vn = nafo4Vn,
        cfa.nens.23=cfa.nens.23,
        cfa.23.24=cfa.23.24,
        cfa.4x.24=cfa.4x.24
      )

      if (returntype=="list") return(out)
      if (returntype=="sf") {
        require(sf)
        out = st_sfc( st_multilinestring( sapply( out, as.matrix) ) )
        st_cast(out, "MULTILINESTRING")
        st_crs(out) = st_crs( projection_proj4string("lonlat_wgs84") )
        out = st_transform(out, project_to) 
        return(out)
      }

      return(out)

  }

}
