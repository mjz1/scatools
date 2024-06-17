.get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


.normalize <- function(x, na.rm = TRUE) {
  return((x - min(x)) / (max(x) - min(x)))
}
