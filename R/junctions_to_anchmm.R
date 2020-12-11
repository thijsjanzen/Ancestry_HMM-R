
#' converts output from simulations of the junctions package to input for
#' Ancestry_HMM

junctions_to_anchmm <- function(sim_data) {

  # we need per line:
  # chrom, pos, A_1, a_1, A_2, a_2, recomdiff, A_indiv1, a_indiv1, A_indiv2 etc
  num_indiv <- length(unique(sim_data$individual))

  a1 <- subset(sim_data, sim_data$individual == unique(sim_data$individual)[1])

  output_matrix <- matrix(NA, nrow = length(a1$time), ncol = 7 + 2 * num_indiv)

  for(i in unique(sim_data$individual)) {
    focal_data <- subset(sim_data, sim_data$individual == i)
    testit::assert(length(focal_data$time) == nrow(output_matrix))

    for(j in 1:length(focal_data$time)) {
      pos <- focal_data$location[j] * 1e12
      if (i == unique(sim_data$individual[1])) {
        output_matrix[j, 1] <- 1 # chrom
        output_matrix[j, 2] <- floor(focal_data$location[j] * 1e7)
        output_matrix[j, 3] <- 100  # A_1
        output_matrix[j, 4] <- 0 # a_1
        output_matrix[j, 5] <- 0 # A_2
        output_matrix[j, 6] <- 100  # a_2

        recomdiff <- -1
        if (j > 1) {
          recomdiff <- focal_data$location[j] - focal_data$location[j-1]
        }

        output_matrix[j, 7] <- recomdiff  # recompos
      }

      # now we have to count alleles
      genotype <- as.numeric(c(focal_data$anc_chrom_1[j],
                               focal_data$anc_chrom_2[j]))
      alleles <- c(sum(genotype == 0), sum(genotype == 1))
      location_1 <- 8 + 1 * (i * 2)
      location_2 <- location_1 + 1
      output_matrix[j, c(location_1, location_2)] <- alleles
    }
  }
  return(output_matrix)
}
