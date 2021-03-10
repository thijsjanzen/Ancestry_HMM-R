get_emp_data_entry <- function(rel_pos, rel_emp_pos) {
  a <- max(which(rel_emp_pos < rel_pos))
  return(a)
}

#' converts output from simulations of the junctions package to input for
#' Ancestry_HMM
#' @param sim_data  output from junctions::sim_phased_unphased
#' @param emp_data empirical data.
#' @return list with posterior probabilities, age estimate, and likelihood
#' @export
junctions_to_anchmm_using_data <- function(sim_data,
                                           emp_data) {

  # we need per line:
  # chrom, pos, A_1, a_1, A_2, a_2, recomdiff, A_indiv1, a_indiv1, A_indiv2 etc
  num_indiv <- length(unique(sim_data$individual))

  rel_emp_pos <- emp_data$V2
  min_pos <- min(rel_emp_pos)
  max_pos <- max(rel_emp_pos)
  rel_emp_pos <- (rel_emp_pos - min_pos) /
                    (max_pos - min_pos)


  a1 <- subset(sim_data, sim_data$individual == unique(sim_data$individual)[1])

  output_matrix <- matrix(NA, nrow = length(a1$time), ncol = 7 + 2 * num_indiv)

  for (i in unique(sim_data$individual)) {
    focal_data <- subset(sim_data, sim_data$individual == i)
    testit::assert(length(focal_data$time) == nrow(output_matrix))
    size_chrom <- 10000 / min(focal_data$location)

    for (j in 1:length(focal_data$anc_chrom_1)) {

      location <- focal_data$location[j]
      testit::assert(location <= 1)
      emp_data_index <-  get_emp_data_entry(location, rel_emp_pos)
      emp_data_entry <-  emp_data[emp_data_index, ]

      if (i == unique(sim_data$individual[1])) {
        output_matrix[j, 1] <- 1 # chrom
        output_matrix[j, 2] <- floor(focal_data$location[j] * size_chrom)
        output_matrix[j, 3] <- emp_data_entry[[3]]
        output_matrix[j, 4] <- emp_data_entry[[4]] # a_1
        output_matrix[j, 5] <- emp_data_entry[[5]] # A_2
        output_matrix[j, 6] <- emp_data_entry[[6]]  # a_2

        recomdiff <- -1
        if (j > 1) {
          recomdiff <- focal_data$location[j] - focal_data$location[j - 1]
        }

        output_matrix[j, 7] <- recomdiff  # recompos
      }

      # now we have to count alleles
      genotype <- as.numeric(c(focal_data$anc_chrom_1[j],
                               focal_data$anc_chrom_2[j]))

      alleles <- sample_from_emp(emp_data_entry, genotype)

      location_1 <- 8 + 1 * (i * 2)
      location_2 <- location_1 + 1
      output_matrix[j, c(location_1, location_2)] <- alleles
    }
  }
  return(output_matrix)
}

#' @keywords internal
sample_from_emp <- function(emp_data_entry, genotype) {
  allele_dist <- list(c(emp_data_entry[[3]], emp_data_entry[[4]]),
                      c(emp_data_entry[[5]], emp_data_entry[[6]]))
  alleles <- c(0, 0) # options:  0 2  //  1 1 // 2 0
  for (i in 1:2) {
    prob_alleles <- allele_dist[[ 1 + genotype[i]]]
    chosen_allele <- sample(x = 1:2, size = 1, prob = prob_alleles)
    alleles[chosen_allele] <- alleles[chosen_allele] + 1
  }
  return(alleles)
}

