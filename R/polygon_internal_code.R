
polygon_internal_code = function (area=NULL) {
  out = NULL

  if (length(area) > 0) {
    r = NULL
    for (a in area) {
      if (a=="cfanorth") r=c( "cfa20", "cfa21", "cfa22")
      if (a=="cfasouth") r=c( "cfa23", "cfa24")
      if (a=="cfaslope") r=c( "cfa24slope", "cfa23slope")
      if (a=="cfaall") r=c( "cfa20", "cfa21", "cfa22","cfa23", "cfa24","cfa4x")
      if (a=="cfa4x") r=c( "cfa4x")
      if (a=="cfanorth.not.glace.bay") r=c("cfa20", "cfa21", "cfa22inner")
      if (a=="cfa.23ab.24ab") r=c("cfa23a", "cfa23b", "cfa24a", "cfa24b")
      if (a=="nafo.2j3kl") r=c(r, "nafo.2j", "nafo.3k", "nafo.3l")
      if (a=="nafo.3no") r=c(r, "nafo.3n", "nafo.3o")
      if (a=="nafo.4rs") r=c(r, "nafo.4r", "nafo.4s")
      if (a=="nafo.4vw") r=c(r, "nafo.4v", "nafo.4w")
      if (a=="nafo.5zew") r=c(r, "nafo.5ze", "nafo.5zw")
      if (a=="4v") r="nafo.4v"
      if (a=="4w") r="nafo.4w"
      if (a=="4x") r="nafo.4x"
      if (a=="4vw") r=c("nafo.4v", "nafo.4w")
      if (a=="4vwx") r=c("nafo.4v", "nafo.4w", "nafo.4x")
      if (a=="5yz") r=c("nafo.5yz")
      if (a=="5y") r=c("nafo.5y")
      if (a=="5z") r=c("nafo.5z")
			if (a=="canada.east") {
				r=c("nafo.2j", "nafo.3k", "nafo.3l", "nafo.3n", "nafo.3o", 
						"nafo.4r", "nafo.4s", "nafo.4v", "nafo.4w", "nafo.4x",
						"nafo.5ze", "nafo.5zw"
			  )
			}
      if (is.null(r)) r = area
      out = c(out, r)
    }
  }

  return( unique(out) )
}



