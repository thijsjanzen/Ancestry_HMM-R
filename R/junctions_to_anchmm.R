#' converts output from simulations of the junctions package to input for
#' Ancestry_HMM
#' @param sim_data  output from junctions::sim_phased_unphased
#' @return list with posterior probabilities, age estimate, and likelihood
#' @rawNamespace useDynLib(ancestryhmmR)
#' @rawNamespace import(Rcpp)
#' @export
junctions_to_anchmm <- function(sim_data) {

  # we need per line:
  # chrom, pos, A_1, a_1, A_2, a_2, recomdiff, A_indiv1, a_indiv1, A_indiv2 etc
  num_indiv <- length(unique(sim_data$individual))

  a1 <- subset(sim_data, sim_data$individual == unique(sim_data$individual)[1])

  output_matrix <- matrix(NA, nrow = length(a1$time), ncol = 7 + 2 * num_indiv)

  for (i in unique(sim_data$individual)) {
    focal_data <- subset(sim_data, sim_data$individual == i)
    testit::assert(length(focal_data$time) == nrow(output_matrix))

    for (j in 1:length(focal_data$time)) {
      if (i == unique(sim_data$individual[1])) {
        output_matrix[j, 1] <- 1 # chrom
        output_matrix[j, 2] <- floor(focal_data$location[j] * 1e7)
        output_matrix[j, 3] <- 100  # A_1
        output_matrix[j, 4] <- 0 # a_1
        output_matrix[j, 5] <- 0 # A_2
        output_matrix[j, 6] <- 100  # a_2

        recomdiff <- -1
        if (j > 1) {
          recomdiff <- focal_data$location[j] - focal_data$location[j - 1]
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

get_alleles <- function(d_value) {

  px <- runif(n = 1, min = d_value, max = 1)
  py <- px - d_value

  # now that we have px and py, we only have to create the associated
  # alleles
  alleles_x <- round(c(px * 100, (1 - px) * 100))
  alleles_y <- round(c(py * 100, (1 - py) * 100))

  return(c(alleles_x, alleles_y))
}

#' converts output from simulations of the junctions package to input for
#' Ancestry_HMM
#' @param sim_data  output from junctions::sim_phased_unphased
#' @param d_value average discriminatory value of the markers used, where the
#' discriminatory value is given by: p_A / (p_A + p_a), e.g. the d value
#' indicates the probability of getting the ancestry right.
#' @return list with posterior probabilities, age estimate, and likelihood
junctions_to_anchmm_imperfect <- function(sim_data,
                                          d_value) {

  # we need per line:
  # chrom, pos, A_1, a_1, A_2, a_2, recomdiff, A_indiv1, a_indiv1, A_indiv2 etc
  num_indiv <- length(unique(sim_data$individual))

  a1 <- subset(sim_data, sim_data$individual == unique(sim_data$individual)[1])

  output_matrix <- matrix(NA, nrow = length(a1$time), ncol = 7 + 2 * num_indiv)

  for (i in unique(sim_data$individual)) {
    focal_data <- subset(sim_data, sim_data$individual == i)
    testit::assert(length(focal_data$time) == nrow(output_matrix))

    for (j in 1:length(focal_data$time)) {
      if (i == unique(sim_data$individual[1])) {
        output_matrix[j, 1] <- 1 # chrom
        output_matrix[j, 2] <- floor(focal_data$location[j] * 1e7)

        alleles <- get_alleles(d_value)

        output_matrix[j, 3] <- alleles[1]  # A_1
        output_matrix[j, 4] <- alleles[2] # a_1
        output_matrix[j, 5] <- alleles[3] # A_2
        output_matrix[j, 6] <- alleles[4]  # a_2

        recomdiff <- -1
        if (j > 1) {
          recomdiff <- focal_data$location[j] - focal_data$location[j - 1]
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
