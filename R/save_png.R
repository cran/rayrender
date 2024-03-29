#'@title Save PNG
#'
#'@description Writes the hillshaded map to file.
#'
#'@param hillshade Array (or matrix) of hillshade to be written.
#'@param filename String with the filename. If `.png` is not at the end of the string, it will be appended automatically.
#'@param rotate Default 0. Rotates the output. Possible values: 0, 90, 180, 270.
#'@keywords internal
#'@examples
#'#internal
save_png = function(hillshade,filename ,rotate=0) {
  if(is.null(filename)) {
    stop("save_png requires a filename")
  }
  rotatef = function(x) t(apply(x, 2, rev))
  if(!(rotate %in% c(0,90,180,270))) {
    warning(paste0("Rotation value ",rotate," not in c(0,90,180,270). Ignoring"))
    number_of_rots = 0
  } else {
    number_of_rots = rotate/90
  }
  if(substring(filename, nchar(filename)-3,nchar(filename)) != ".png") {
    filename = paste0(filename,".png")
  }
  if(length(dim(hillshade)) == 2) {
    if(number_of_rots != 0) {
      for(j in 1:number_of_rots) {
        hillshade = rotatef(hillshade)
      }
    }
    final = array(t(hillshade[,ncol(hillshade):1]),dim=c(ncol(hillshade),nrow(hillshade),3))
    png::writePNG(final,filename)
  } else {
    if(number_of_rots != 0) {
      newarray = hillshade
      n_layers = dim(hillshade)[3]
      newarrayt = array(1,dim=c(ncol(hillshade),nrow(hillshade),n_layers))
      for(i in 1:number_of_rots) {
        for(j in 1:n_layers) {
          if(i == 2) {
            newarray[,,j] = rotatef(newarrayt[,,j])
          } else {
            newarrayt[,,j] = rotatef(newarray[,,j])
          }
        }
      }
      if(number_of_rots == 2) {
        hillshade = newarray
      } else {
        hillshade = newarrayt
      }
    }
    png::writePNG(hillshade,filename)
  } 
}
