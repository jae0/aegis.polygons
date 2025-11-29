
area_lines.db = function( DS, returntype="list", project_to=st_crs( projection_proj4string("lonlat_wgs84") ) ) {

  if (DS=="cfa.regions") {
      
      nafo4Vn  = data.frame( rbind(
    c(-60.25,45.7),    
        c(-60,45.666667),    
        c(-57.50092,45.666667),    
        c(-57.51692,45.681329),    
        c(-57.684715,45.826751),    
        c(-57.85251,45.972173),    
        c(-58.020304,46.117595),    
        c(-58.188099,46.263017),    
        c(-58.355894,46.408439),    
        c(-58.523688,46.553861),    
        c(-58.691483,46.699283),    
        c(-58.84952,46.836615),    
        c(-58.859278,46.844705),    
        c(-59.027072,46.990127),    
        c(-59.194867,47.135549),    
        c(-59.362662,47.280971),    
        c(-59.530457,47.426393),    
        c(-59.57857,47.4677),    
        c(-59.698251,47.571815),    
        c(-59.866046,47.717237),    
        c(-60,47.833333),    
        c(-60.40917,47.035833),    
        c(-60.68293,46.145712),    
        c(-61.11908,46.125),    
        c(-61.40417,45.641667),    
        c(-60.25,45.7)    
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
