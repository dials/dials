library(rgl)
source <- c(-0.002683765, 0.001318010, 0.666659962)
rot_ax <- c(1.000000000, 0.000000000, 0.000000000)
dstarmax <- 1.000000000

# draw the Ewald and limiting spheres
open3d()
spheres3d(source,radius=sqrt(sum(source*source)),color='#CCCCFF',alpha=0.3)
spheres3d(c(0,0,0),radius=dstarmax,color='red',alpha=0.1)

# draw the source vector and rotation axis
lines3d(rbind(c(0,0,0),source), col='red')
lines3d(rbind(c(0,0,0),rot_ax))
ap <- read.csv('ap_indices.dat', header=FALSE, col.names=c('h','k','l'), colClasses='numeric')
ra <- read.csv('ra_indices.dat', header=FALSE, col.names=c('h','k','l'), colClasses='numeric')
