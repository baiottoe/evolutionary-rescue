# Function to generate binary sequences
generate_binary_sequence <- function(length, ones_count) {
  # Create a sequence with 'ones_count' ones and 'length - ones_count' zeros
  sequence <- c(rep(1, ones_count), rep(0, length - ones_count))
  sequence <- sample(sequence)  # Shuffling in R
  return(paste(sequence, collapse = ''))
}