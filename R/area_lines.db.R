
area_lines.db = function( DS ) {

  if (DS=="cfa.regions") {

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

      names( cfa.nens.23 ) = names( cfa.23.24 ) = names( cfa.4x.24 ) = c("lon", "lat")

      out = list(
        cfa.nens.23=cfa.nens.23,
        cfa.23.24=cfa.23.24,
        cfa.4x.24=cfa.4x.24
      )

      return(out)

  }

}
