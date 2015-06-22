setMethod("computeDistanceMatrix", "RasterLayer",
          function(object){
            # Extract coordinates from raster object
            coords = xyFromCell(object = object, cell = 1:ncell(object), spatial=FALSE)
            
            # Compute distance matrix
            dist <- apply(X = coords,
                          MARGIN = 1,
                          FUN = function(x){ spDistsN1(as.matrix(coords), x, longlat=TRUE) }
            )
            
            return(dist) 
          }
)
