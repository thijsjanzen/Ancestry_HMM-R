context("compare methods")

test_that("compare methods", {

  sim_data <- junctions::sim_phased_unphased(pop_size = 1000,
                                             freq_ancestor_1 = 0.5,
                                             total_runtime = 50,
                                             time_points = 10,
                                             markers = 1000,
                                             num_indiv_sampled = 10,
                                             seed = 42)

  anchmm_data <- junctions_to_anchmm(sim_data)
  which(anchmm_data[, 7] < 0)


  sample_matrix <- matrix(nrow = length(unique(sim_data$individual)), ncol = 2)
  for (i in 1:length(sample_matrix[, 1])) {
    sample_matrix[i, ] <- c(i, 2)
  }

  cmd_line_options <- c("-i", "input_file.txt",
                        "-s", "sample.txt",
                        "-a", "2", "0.5", "0.5",
                        "-p", "0", "100000", "0.5",
                        "-p", "1", "-10", "0.5",
                        "-g")


  result <- ancestryhmmR::run_ancestry_hmm(sample_matrix = sample_matrix,
                                           genetic_data = anchmm_data,
                                           cmd_line_options = cmd_line_options,
                                           viterbi = FALSE,
                                           use_genetic_data = TRUE)


  result2 <- ancestryhmmR::run_ancestry_hmm(sample_matrix = sample_matrix,
                                            genetic_data = anchmm_data,
                                            cmd_line_options = cmd_line_options,
                                            viterbi = FALSE,
                                            pulse1 = c(0, 100000, 0.5),
                                            pulse2 = c(1, -10, 0.5),
                                            use_genetic_data = TRUE)

  testthat::expect_true(all.equal(result, result2))
})
