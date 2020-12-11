#ifndef __READ_INPUT_R_H
#define __READ_INPUT_R_H

#include <sstream>


void read_file_r( cmd_line &options,
                std::vector<markov_chain> &markov_chain_information,
                std::map<int, std::vector< std::vector<int> > > &state_list,
                std::vector<int> &position,
                std::vector<double> &recombination_rate,
                std::vector< std::string > &chromosomes,
                const Rcpp::NumericMatrix& ancestry_matrix) {

  /// vector to hold index of inbred path if we have variable ploidy
  std::vector<int> path_index( markov_chain_information.size(), 0 ) ;

  //// since the first site transition matrix does not matter, we can print anything
  double extra_recombination = 1 ;
  std::string last_chrom = "" ;

  for(int i = 0; i < ancestry_matrix.nrow(); ++i) {

    Rcpp::NumericVector focal_line = ancestry_matrix(i, Rcpp::_);

    input_line new_line;
    int index = 0;
    //  Rcpp::Rcout << focal_line[index] << " " << "\n"; force_output();

    std::stringstream ss;
    ss << focal_line[index];

    new_line.chrom = ss.str();

    new_line.pos   = focal_line[++index];         // in >> new_line.chrom >> new_line.pos ;
    // Rcpp::Rcout << new_line.chrom << " " << new_line.pos << " ";
    /// if two adjacent sites have the same positions, skip second
    if ( ( position.size() > 0 && new_line.pos == position.back() ) ) {
      //getline( in, new_line.chrom ) ;
      i++;
      continue ;
    }

    // read reference panel genotype counts
    new_line.reference_counts.resize( options.ancestry_proportion.size() ) ;


    for ( int p = 0 ; p < options.ancestry_proportion.size() ; p ++ ) {
      double count1 = focal_line[++index];
      double count2 = focal_line[++index];
      new_line.reference_counts[p].push_back(count1) ;
      new_line.reference_counts[p].push_back(count2) ;
      new_line.reference_counts[p].push_back(count1+count2) ;
      //  Rcpp::Rcout << count1 << " " << count2 << " ";
    }

    /// read recombination rate
    //in >> new_line.recombination_rate ;
    new_line.recombination_rate = focal_line[++index];
    // Rcpp::Rcout << new_line.recombination_rate << " ";

    /// if line specific error rates are provided
    new_line.error_1 = options.error_rate ;
    new_line.error_2 = options.error_rate ;

    if ( options.error_rates == true ) {
      new_line.error_1 = focal_line[++index];
      new_line.error_2 = focal_line[++index];
      // in >> new_line.error_1 >> new_line.error_2 ;
    }

    // read sample panel read counts
    new_line.sample_counts.resize( markov_chain_information.size() ) ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
      double count1 = focal_line[++index];
      double count2 = focal_line[++index];
      // Rcpp::Rcout << count1 << " " << count2 << " ";
      //in >> count1 >> count2 ;

      /// subsample reads to a maximum depth so we can compute multinomial probs without overflow errors
      if ( count1 + count2 > 170 ) {
        subsample_reads( count1, count2 ) ;
      }

      /// now store counts and total for sample
      new_line.sample_counts[m].push_back(count1) ;
      new_line.sample_counts[m].push_back(count2) ;
      new_line.sample_counts[m].push_back(count1+count2) ;
    }

    if ( new_line.chrom != last_chrom ) {
      recombination_rate.push_back( 0.5 ) ;
      last_chrom = new_line.chrom ;
      extra_recombination = 0 ;
    }

    /// ignore lines where recombination may not be sufficient to make sites independent
    /// this might be useful in place of LD pruning
    else {

      extra_recombination += new_line.recombination_rate ;
      if ( extra_recombination < options.minimum_distance ) {
        continue ;
      }
      new_line.recombination_rate = extra_recombination ;
      extra_recombination = 0 ;

      recombination_rate.push_back( new_line.recombination_rate/ ( new_line.pos - position.back() ) ) ;
    }

    /// record position
    position.push_back( new_line.pos ) ;
    chromosomes.push_back( new_line.chrom ) ;

    /// check all path indexes and update as needed
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {

      if ( markov_chain_information[m].path_file != "null" ) {

        /// record previous ploidy
        int previous_ploidy = markov_chain_information[m].sample_ploidy_path[path_index[m]].ploidy ;

        /// check to make sure we're on the right ploidy tract
        while ( new_line.chrom != markov_chain_information[m].sample_ploidy_path[path_index[m]].chrom ) {
          path_index[m] ++ ;
        }
        while ( new_line.pos > markov_chain_information[m].sample_ploidy_path[path_index[m]].stop ) {
          path_index[m] ++ ;
        }

        /// record switches
        if ( previous_ploidy != markov_chain_information[m].sample_ploidy_path[path_index[m]].ploidy ) {

          markov_chain_information[m].ploidy_switch_position.push_back( position.size() - 1 ) ;
          markov_chain_information[m].ploidy_switch.push_back( markov_chain_information[m].sample_ploidy_path[path_index[m]].ploidy ) ;
        }
      }
    }

    ///
    if ( options.genotype == false ) {
      for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        vec emissions ;
        create_emissions_matrix( markov_chain_information[m].sample_ploidy_path[path_index[m]].ploidy, new_line, options.ancestral_fixed, state_list[markov_chain_information.at(m).sample_ploidy_path[path_index[m]].ploidy], m, options.ancestry_pulses, emissions ) ;
        markov_chain_information[m].emission_probabilities.push_back( emissions ) ;
      }
    }

    /// create emissions matrix with genotypes
    else {
      for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        vec emissions ;
        create_emissions_matrix_genotype( markov_chain_information[m].sample_ploidy_path[path_index[m]].ploidy, new_line, options.ancestral_fixed, state_list[markov_chain_information.at(m).sample_ploidy_path[path_index[m]].ploidy], m, options.ancestry_pulses, emissions ) ;
        markov_chain_information[m].emission_probabilities.push_back( emissions ) ;
      }
    }
    // Rcpp::Rcout << "\n";
  }

  /// to avoid lookahead errors
  for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
    markov_chain_information[m].ploidy_switch_position.push_back( position.size() ) ;
  }

}

#endif