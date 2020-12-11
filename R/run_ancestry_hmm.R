#' run Ancestry_HMM
#' @param sample_matrix matrix with sample names and ploidy
#' @param genetic_data matrix with the genetic data (see manual)
#' @param cmd_line_options cmd line entries in a vector, where each element
#' separated by a space is in a different position, e.g. c("-i", "input.txt")
#' @param viterbi should viterbi output be generated?
#' @param pulse1 if no cmd_line_options are given, a pulse can be added manually
#' with 3 properties: c(ancestor, time, proportion)
#' @param pulse2 if no cmd_line_options are given, a pulse can be added manually
#' with 3 properties: c(ancestor, time, proportion)
#' @param use_genetic_data is the genetic data in reads, or alleles?
#' @return list with posterior probabilities, age estimate, and likelihood
#' @export
run_ancestry_hmm <- function(sample_matrix,
                             genetic_data,
                             cmd_line_options,
                             viterbi = FALSE,
                             pulse1 = numeric(0),
                             pulse2 = numeric(0),
                             use_genetic_data = TRUE) {

  # add input data checks.

  results <- run_ancestry_hmm_cpp(sample_matrix,
                                  cmd_line_options,
                                  genetic_data,
                                  viterbi,
                                  pulse1,
                                  pulse2,
                                  use_genetic_data)

  return(results)
}