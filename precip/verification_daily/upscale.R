## Funtion that computes sum/averages of a given number of columns
## Input (matrix,operation,scale)
## matrix <- a matrix of obs or ens
## operation -> [1,2] for [average,sum]

upscale <- function(matrix,operation,scale){
  if (operation == 1){ # For Temperature
    if (scale == 2){ # Weekly
      a0 <- rowMeans(matrix[,1:7],na.rm = TRUE) 
      a1 <- rowMeans(matrix[,8:14],na.rm = TRUE)
      a2 <- rowMeans(matrix[,15:21],na.rm = TRUE)
      a3 <- rowMeans(matrix[,22:28],na.rm = TRUE)
      a4 <- rowMeans(matrix[,29:35],na.rm = TRUE) 
      a5 <- rowMeans(matrix[,36:42],na.rm = TRUE)
      a6 <- rowMeans(matrix[,43:49],na.rm = TRUE) 
      a7 <- rowMeans(matrix[,50:56],na.rm = TRUE) 
      a8 <- rowMeans(matrix[,57:63],na.rm = TRUE) 
      a9 <- rowMeans(matrix[,64:70],na.rm = TRUE) 
      a10 <-rowMeans(matrix[,71:77],na.rm = TRUE) 
      a11 <- rowMeans(matrix[,78:90],na.rm = TRUE)
      b0 <- cbind(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)
      matrix_out <- b0
      return(matrix_out)
    }# End weekly
    if (scale == 3){ # For monthly
      a0 <- rowMeans(matrix[,1:30],na.rm = TRUE) 
      a1 <- rowMeans(matrix[,31:60],na.rm = TRUE)
      a2 <- rowMeans(matrix[,61:90],na.rm = TRUE)
      b0 <- cbind(a0,a1,a2)
      matrix_out <- b0
      return(matrix_out)
    }# End monthly
    if (scale == 4){ # For season
      a0 <- rowMeans(matrix,na.rm = TRUE) 
      matrix_out <- a0
      return(matrix_out)
    }# End season
    
  } else {  # For precipitation and et0
    if (scale == 2){ # Weekly
      a0 <- rowSums(matrix[,1:7],na.rm = TRUE) 
      a1 <- rowSums(matrix[,8:14],na.rm = TRUE)
      a2 <- rowSums(matrix[,15:21],na.rm = TRUE)
      a3 <- rowSums(matrix[,22:28],na.rm = TRUE)
      a4 <- rowSums(matrix[,29:35],na.rm = TRUE) 
      a5 <- rowSums(matrix[,36:42],na.rm = TRUE)
      a6 <- rowSums(matrix[,43:49],na.rm = TRUE) 
      a7 <- rowSums(matrix[,50:56],na.rm = TRUE) 
      a8 <- rowSums(matrix[,57:63],na.rm = TRUE) 
      a9 <- rowSums(matrix[,64:70],na.rm = TRUE) 
      a10 <-rowSums(matrix[,71:77],na.rm = TRUE) 
      a11 <- rowSums(matrix[,78:90],na.rm = TRUE)
      b0 <- cbind(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)
      matrix_out <- b0
      return(matrix_out)
    }# End weekly
    if (scale == 3){ # For monthly
      a0 <- rowSums(matrix[,1:30],na.rm = TRUE) 
      a1 <- rowSums(matrix[,31:60],na.rm = TRUE)
      a2 <- rowSums(matrix[,61:90],na.rm = TRUE)
      b0 <- cbind(a0,a1,a2)
      matrix_out <- b0
      return(matrix_out)
    }# End monthly
    if (scale == 4){ # For season
      a0 <- rowSums(matrix,na.rm = TRUE) 
      matrix_out <- a0
      return(matrix_out)
    }# End season
  } ## End of if condition
} # End of function