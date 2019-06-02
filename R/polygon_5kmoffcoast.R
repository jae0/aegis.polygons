
polygon_5kmoffcoast = function( plotdata=FALSE ) {

  # Creates polygons 5km away from the coast
  #	that will be used to flag coastal and ocean data
  # before bottle data filtering (chlorophyll and nutrients).

  require(maps)
  require(mapdata)
  library(sp)
  #install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/stable")
  require(INLA)
  require(splancs)

  pd = project.datadirectory("aegis", "biochem")

  pddata = file.path( pd, "data")

  #construct a hull

  # get the coastline for the region of interest
  coastline = maps::map( "worldHires", regions=c("Canada", "USA"), xlim=c(-71, -48 ), ylim=c(37, 50), fill=FALSE, plot=FALSE )
  plot( coastline$x, coastline$y, pch="." )

  # Make bounding box interactively by placing points on the map.
  # Done once, no need to repeat. Results are bb, saved in fnbb.
  # This box will be used to select coastline part of interest.
  # bb = locator( type="o")  # bounding box manually determined

  fnbb = file.path( pddata, "boundingbox.coastline.ns.rdata" )
  # save( bb, file=fnbb)
  load( fnbb)

  # Convert to data frame
  aoi = data.frame( cbind( x=bb$x, y=bb$y))

  # Close the polygon
  aoi = rbind(aoi, aoi[1,] )

  # Select coastline within bounding box
  a = which( point.in.polygon( coastline$y, coastline$x, aoi$y, aoi$x) != 0 )
  cl = data.frame( cbind(x=coastline$x, y=coastline$y))
  cl = cl[a,]

  # Remove Canada-US border. Again, create bounding box interactively,
  # no need to repeat.The points of the box are in bborder and saved in fnborder.
  # bborder = locator( type="o")  # bounding box manually determined

  fnborder = file.path( pddata, "boundingbox.can_us_border.rdata" )
  # save( bborder, file=fnborder)
  load( fnborder)

  # Convert to data frame
  bborder = as.data.frame( bborder )

  # Close the polygon
  bborder = rbind(bborder, bborder[1,])

  # Find the points in the polygon i.e indices for US border
  a = which( point.in.polygon( cl$y, cl$x, bborder$y, bborder$x) != 0 )

  # Remove border from the coastline
  cl = cl[-a,]

  # Remove Cape Cod, the same way
  # bcapecod = locator( type="o")  # bounding box manually determined

  fnbcapecod = file.path( pddata, "boundingbox.bcapecod.rdata" )
  # save( bcapecod, file=fnbcapecod)
  load( fnbcapecod)

  bcapecod = as.data.frame( bcapecod )
  bcapecod = rbind(bcapecod, bcapecod[1,])

  a = which( point.in.polygon( cl$y, cl$x, bcapecod$y, bcapecod$x) != 0 )
  cl = cl[-a,]

  # Find end points of the coastline
  pt.top = which.min( cl$x )
  pt.bot = which.min( cl$y )

  # Plot the end points of the coastline
  points( y~x, cl[pt.top,], pch=20, col="red")
  points( y~x, cl[pt.bot,], pch=20, col="red")

  # Final polygons
  cl.fin = rbind( cl[1:pt.bot,], cl[pt.top,], cl[(pt.bot+1):nrow(cl) , ] )

  # convert to planar coordinates
  names(cl.fin) = c("lon", "lat")
  k = lonlat2planar( cl.fin, proj.type="utm20")
  k = k[, c("plon", "plat")] # k is coastline in planar

  ok = as.matrix(k)

  # Create a square with corners c1-c4
  c1 = c(-5,+5)
  c2 = c(5,  5)
  c3 = c(5, -5)
  c4 = c(-5, -5)

  # make a square around each point of the coastline k
  out = NULL
  out = rbind( out,  t(t(k)+c1) )
  out = rbind( out,  t(t(k)+c2) )
  out = rbind( out,  t(t(k)+c3) )
  out = rbind( out,  t(t(k)+c4) )

  # create non-convex hull using inla package

  coast5km = inla.nonconvex.hull( out, convex=5, resolution=500 )

  # output is a list. save
  fn5km= file.path( pddata, "coast5km.all.polygons.list.rdata" )
  save( coast5km, file=fn5km)
  load(fn5km)

  # extract x and y from the list and define flag for flagging polygons
  x=coast5km$loc[,1]
  y=coast5km$loc[,2]
  flag=NULL

  # flag the ploygons by looking at discontinuties
  d1=diff(x)
  d2=diff(y)

  # points of discontinuity (d1>27 to resolve Newfounland)
  dd=which(abs(d1)>27 | abs(d2)>30)

  # define polygons using "flag" field
  flag[1:dd[1]]=1 # first polygon

  for (i in 1:(length(dd)-1)) {
  flag[(dd[i]+1):dd[i+1]]=i+1
  }

  # last polygon
  flag[tail(dd,n=1):length(x)]=length(dd)+1


  # plot polygons in different colors
  colo=rainbow(length(dd)+1)
  plot(x,y) #plot all polygons in black

  for (i in 1:(length(dd)+1)) {
    p=which(flag==i)
    points(x[p],y[p],pch=".",cex=3,col=colo[i])
  }

  legend("topleft", legend = 1:length(dd), col=colo, pch=19)


  # convert the list into dataframe
  df5km=cbind.data.frame(x,y,flag)

  # Polygons of interest (5km off the coast at sea):
  # Nova Scotia and mainland: 1
  # Anticosti Island: 4
  # Magdalen Island: 8
  # Sable Island: 10
  # Newfoundland: 11

  co5k=df5km[which(df5km$flag %in% c(1,4,8,10,11)),]

  # convert from planar back to lon and lat
  tt=co5k
  tt$flag=NULL
  names(tt)=c("plon","plat")
  c5k=planar2lonlat(tt,proj.type="utm20")
  c5k$flag=co5k$flag

  # re-number flags for polygons:
  c5k$flag[which(c5k$flag==4)]=2
  c5k$flag[which(c5k$flag==8)]=3
  c5k$flag[which(c5k$flag==10)]=4
  c5k$flag[which(c5k$flag==11)]=5

  # Close each polygon
  c5kf=NULL
  for (i in 1:5) {
   pol=which(c5k$flag==i)
   c5kp = rbind(c5k[pol,], c5k[pol[1],] )
   c5kf=rbind(c5kf,c5kp)
  }



  # save polygons 5km from coastline in lat and lon
  fn2= file.path( pddata, "coast5km.polygons.4filtering.rdata" )
  fn3= file.path( pddata, "coast5km.polygons.4filtering.csv" )

  save(c5kf, file=fn2)
  # write.table( c5kf, file=fn3 )


  if (plotdata=="final") {

    # plot the coastline and final polygons
    plot(coastline$x, coastline$y, pch=".",xlab="Longitude", ylab="Latitude") #coastline

    uf=unique(c5k$flag) # find polygons
    colo=rainbow(length(uf)) #define colors

    #plot polygons
    for (i in 1:length(uf)) {
      p=which(c5k$flag==uf[i])
      points(c5k$lon[p],c5k$lat[p],pch=".",cex=3,col=colo[i])
    }

    # add legend
    leg=c("1: Mainland","2: Anticosti Island","3: Magdalen Island", "4: Sable Island", "5: Newfoundland")
    legend("bottomright", legend =leg, col=colo, pch=19, title="Polygons:")

  }



}
