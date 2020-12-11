#ifndef __READ_SAMPLES_R_H
#define __READ_SAMPLES_R_H

void read_samples_r(std::vector<markov_chain> &markov_chain_information,
                  const Rcpp::NumericMatrix& m,
                  bool viterbi) {

  for(int i = 0; i < m.nrow(); ++i) {

    Rcpp::NumericVector local_sample = m(i, Rcpp::_);

    markov_chain new_sample;

    new_sample.output_file = local_sample[0];

    if ( new_sample.output_file == "" ) {
      continue ;
    }

    if ( viterbi == false ) {
      new_sample.output_file.append( ".posterior" ) ;
    }
    else {
      new_sample.output_file.append( ".viterbi" ) ;
    }

    new_sample.number_chromosomes = local_sample[1];

    new_sample.path_file = "null" ;

    /// store new samples
    markov_chain_information.push_back( new_sample ) ;
  }

  /// read or generate ploidy paths
  for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
    if ( markov_chain_information[m].path_file != "null" ) {
      read_ploidy_file( markov_chain_information[m].path_file, markov_chain_information[m].sample_ploidy_path ) ;
      markov_chain_information[m].ploidy_switch_position.push_back( 0 ) ;
      markov_chain_information[m].ploidy_switch.push_back( markov_chain_information[m].sample_ploidy_path[0].ploidy ) ;
    }
    else {
      markov_chain_information[m].ploidy_switch_position.push_back( 0 ) ;
      markov_chain_information[m].ploidy_switch.push_back( markov_chain_information[m].number_chromosomes ) ;
      ploidy_entry new_entry ;
      new_entry.ploidy = markov_chain_information[m].number_chromosomes ;
      markov_chain_information[m].sample_ploidy_path.push_back( new_entry ) ;
    }
  }

  for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
    markov_chain_information[m].end_prob = 1 ;
    markov_chain_information[m].start_prob = 1 ;
  }

  return ;
}

#endif