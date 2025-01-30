# Function to calculate Euclidean distance on a torus
euclidean_distance_torus <- function(x1, y1, x2, y2, width, height) {
  dx <- abs(x1 - x2)
  dy <- abs(y1 - y2)
  
  dx <- min(dx, width - dx)
  dy <- min(dy, height - dy)
  
  return(sqrt(dx^2 + dy^2))
}

# Function to calculate distances from a certain index
calculate_distances_torus <- function(matrix, start_row, start_col) {
  num_rows <- nrow(matrix)
  num_cols <- ncol(matrix)
  
  distances <- matrix(0, nrow = num_rows, ncol = num_cols)
  
  for (i in 1:num_rows) {
    for (j in 1:num_cols) {
      distances[i, j] <- euclidean_distance_torus(start_row, start_col, i, j, num_rows, num_cols)
    }
  }
  
  return(distances)
}