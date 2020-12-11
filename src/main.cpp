/*

 copyright: Russ Corbett-Detig
 rucorbet@ucsc.edu

 This is software distributed under the gnu public license version 3.

 */

/// headers
#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <unistd.h>
using namespace std ;

/// linear algebra library is armadillo
#define ARMA_NO_DEBUG
#include "RcppArmadillo.h"
using namespace arma ;
using namespace Rcpp;

#include <chrono>
#include <thread>

void force_output() {
//  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
 // R_FlushConsole();
//  R_ProcessEvents();
//  R_CheckUserInterrupt();
}

// forward declaration
Rcpp::NumericMatrix convert_double_vector(const std::vector< std::vector< double >>& invec);

/// our header files in /src directory
#include "print_usage.h"
#include "factorial.h"
#include "nchoosek.h"
#include "subsample.h"
#include "multichoose.h"
#include "multipermute.h"
#include "normalize.h"
#include "ancestry_pulse.h"
#include "ploidy_path.h"
#include "markov_chain.h"
#include "read_samples.h"
#include "read_samples_r.h"
#include "pulses_to_ancestry.h"
#include "compute_forward.h"
#include "compute_backward.h"
#include "forward_backward.h"
#include "forward_backward_r.h"
#include "viterbi.h"
#include "transition_information.h"
#include "exponentiate_matrix.h"
#include "cmd_line.h"
#include "create_transition_rates.h"
#include "read_cmd_line.h"
#include "evaluate_vertex.h"
#include "check_vertex.h"
#include "sort_vertices.h"
#include "create_pulses.h"
#include "create_states.h"
#include "input_line.h"
#include "distribute_alleles.h"
#include "binomial.h"
#include "read_emissions.h"
#include "genotype_emissions.h"
#include "read_input.h"
#include "read_input_r.h"
#include "nelder_mead.h"
#include "golden_search.h"
#include "bootstrap.h"



pulse create_pulse(const Rcpp::NumericVector& p);


// [[Rcpp::export]]
Rcpp::List run_ancestry_hmm_cpp(const Rcpp::NumericMatrix& sample_matrix,
                            Rcpp::StringVector& cmd_line_options,
                            const Rcpp::NumericMatrix& genetic_data,
                            bool viterbi,
                            const Rcpp::NumericVector& pulse1,
                            const Rcpp::NumericVector& pulse2,
                            bool use_genome_data) {

  const int argc = cmd_line_options.size();
  char* argv[argc];

  Rcout << argc << "\n"; force_output();

  if (argc > 1) {
    for(int i = 0; i < cmd_line_options.size(); ++i) {
      std::string input = Rcpp::as<std::string>(cmd_line_options[i]);
      const char* s = input.c_str();
      std::size_t sz = std::strlen(s);
      char *p = new char[sz];
      std::strcpy(p, s);
      argv[i] = p;
    }
  }

  clock_t t = clock();
  clock_t total = clock() ;

  // read cmd line
  cmd_line options ;
  Rcout << "reading command line\n" ; t = clock(); force_output();

  options.read_cmd_line( argc, &argv[0]);

  if (argc <= 1) {
    Rcout << "using manual entries\n"; force_output();
    options.ancestry_pulses.clear() ;
    // set pulses from numeric Vectors
    pulse new_ancestry_pulse  = create_pulse(pulse1);

    options.ancestry_pulses.push_back( new_ancestry_pulse ) ;
    options.ancestry_pulses.back().entry_order = options.ancestry_pulses.size() - 1 ;

    pulse new_ancestry_pulse2 = create_pulse(pulse2);
    options.ancestry_pulses.push_back( new_ancestry_pulse2 ) ;
    options.ancestry_pulses.back().entry_order = options.ancestry_pulses.size() - 1 ;

    options.genotype = use_genome_data;
  }

  /// chain objects for each sample
  vector<markov_chain> markov_chain_information ;

  /// get sample ids and ploidy from input file
  Rcout << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading sample ids and ploidy" ; t = clock();
  read_samples_r(markov_chain_information,
                 sample_matrix,
                 viterbi) ;

  /// create states matrix
  Rcout << "\t\t\t" << (double) (clock() - t) << " ms\n" << "creating states matrix" ; t = clock();
  /// store all possible state space arranged by ploidy and then vector of state counts

  map<int,vector<vector<int> > > state_list ;
  /// now create initial state list
  for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
    for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
      create_initial_states( markov_chain_information.at(m).sample_ploidy_path[p].ploidy,
                             options.ancestry_pulses,   // this one!
                             state_list ) ;
    }
  }

  /// read in panels and update matrices
  Rcout << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading data and creating emissions matrices\t" ; t = clock() ;
  /// store recombination rates and positions
  vector<int> position ;
  vector<double> recombination_rate ;
  vector<string> chromosomes ;

  read_file_r( options,
             markov_chain_information,
             state_list,
             position,
             recombination_rate,
             chromosomes,
             genetic_data);

  /// create basic transition information
  Rcout << (double) (clock() - t) << " ms" << endl << "computing transition routes\t\t\t" ; t = clock() ;
  /// 3d map to look up by ploidy, start state, end state, and then relevant transition information
  map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information ;
  for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
    for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
      create_transition_information( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, transition_matrix_information, state_list[markov_chain_information.at(m).sample_ploidy_path[p].ploidy] ) ;
    }
  }

  /// create admixture model(s)
  Rcout << (double) (clock() - t) << " ms" << endl << "creating initial admixture model(s)\t\t" ; t = clock();
  vector<vector<pulse> > vertices ;
  int nparams = create_pulses( vertices, options ) ;

  /// set number of restarts if unspecified default is factorial * 2
  Rcout << (double) (clock() - t) << " ms" << endl << "estimating " << nparams << " parameters\n" ;
  if ( options.n_restarts < 0 ) {
    options.n_restarts = factorial_vec[nparams] * 2 ;
  }

  /// vector of models to be evaluated and optimized
  vector<pulse> optimum ;

  /// if there are params to estimate, do amoeba search
  if ( nparams > 1 ) {
    Rcout << "starting nelder-mead search\t\t" << endl ;
    optimum = nelder_mead_search( vertices, options, markov_chain_information, transition_matrix_information, recombination_rate, position, state_list ) ;
    Rcout << "\n\t\t\t\t\tSEARCH TIME: " << (double) (clock() - t) << " ms" << endl << endl << "optimal model found:\n\n" ;
  }

  /// or do golden section line search for single parameter optimization
  else if ( nparams == 1 ) {
    Rcout << "starting golden section search\t\t" << endl ;
    optimum = golden_search( options, markov_chain_information, transition_matrix_information, recombination_rate, position, state_list ) ;
    Rcout << "\n\t\t\t\t\tSEARCH TIME: " << (double) (clock() - t) << " ms" << endl << endl << "optimal model found:\n\n" ;
  }

  /// otherwise just evaluate the supplied model
  else {
    optimum = options.ancestry_pulses ;
    Rcout << endl << endl << "evaluating supplied model:\n\n" ;
  }

  /// print model
  Rcout << "\ttype\ttime\tproportion\n" ;
  Rcout << "optimum: \n" ;
  Rcout << "\ttype\ttime\tproportion\n" ;

  vector<double> a = options.ancestry_proportion ;
  for ( int p = 0 ; p < optimum.size() ; p ++ ) {
    optimum[p].proportion = a[optimum[p].type] * optimum[p].fraction_of_remainder ;
    a[optimum[p].type] -= optimum[p].proportion ;
    Rcout << "\t" << optimum[p].type << "\t" << optimum[p].time << "\t" << optimum[p].proportion << endl ;
    //    cout << "\t" << optimum[p].type << "\t" << optimum[p].time << "\t" << optimum[p].proportion << endl ;
  }

  /// bootstrap models as necessary
  if ( options.n_bootstraps > 0 ) {
    Rcout << "computing " << options.n_bootstraps << " bootstrap models" << endl ;

    vector<vector<pulse> > bootstrap = bootstraps( vertices, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_list, chromosomes ) ;

    //// print out bootstrapped admixture models
    for ( int b = 0 ; b < options.n_bootstraps ; b ++ ) {

      Rcout << "bootstrap: " << b << endl ;
      Rcout << "\ttype\ttime\tproportion\n" ;

      // cerr << endl << "bootstrap: " << b << endl ;
      // cerr << "\ttype\ttime\tproportion\n" ;

      vector<double> a = options.ancestry_proportion ;
      for ( int p = 0 ; p < optimum.size() ; p ++ ) {
        bootstrap[b][p].proportion = a[optimum[p].type] * bootstrap[b][p].fraction_of_remainder ;
        a[bootstrap[b][p].type] -= bootstrap[b][p].proportion ;
        // cerr << "\t" << bootstrap[b][p].type << "\t" << bootstrap[b][p].time << "\t" << bootstrap[b][p].proportion << endl ;
        Rcout << "\t" << bootstrap[b][p].type << "\t" << bootstrap[b][p].time << "\t" << bootstrap[b][p].proportion << endl ;
      }
    }
  }

  /// create transition rates for the optimal or supplied set of pulses
  Rcout << endl << "creating per morgan transition rates\t\t" ; t = clock();
  mat transition_rates = create_transition_rates( optimum, options.ne, options.ancestry_proportion ) ;

  /// create transition information
  Rcout << (double) (clock() - t) << " ms" << endl << "creating transition matrices\t\t\t" ; t = clock();
  map<int,vector<mat> > transition_matrix ;
  for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
    create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information.at(m).number_chromosomes], recombination_rate, position, markov_chain_information.at(m).number_chromosomes, transition_rates ) ;
    for ( int p = 0 ; p < markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
      create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information[m].ploidy_switch[p]], recombination_rate, position, markov_chain_information[m].ploidy_switch[p], transition_rates ) ;
    }
  }

  //// create interploidy transition matrix
  vector<mat> interploidy_transitions = create_interploidy_transitions ( state_list, optimum, options.ancestry_proportion ) ;

  NumericMatrix probabilities;
  double lnl = 0 ;
  /// output viterbi path for optimized model
  if ( options.viterbi == true ) {
    Rcout << (double) (clock() - t) << " ms" << endl << "viterbi posterior decoding and printing\t" ; t = clock() ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
      markov_chain_information[m].viterbi( position, recombination_rate, state_list, chromosomes, transition_matrix, interploidy_transitions, options.output_pulses, optimum ) ;
    }
  }

  /// output forward-backward full probability distribution by default
  else {
    Rcpp::Rcout << (double) (clock() - t) << " ms" << endl << "computing forward probabilities\t" ; t = clock() ;

    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
      lnl += markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions ) ;
    }

    Rcpp::Rcout << "lnl: " << lnl << "\t\t" << (double) (clock() - t) << " ms" << endl << "computing backward probabilities\t\t\t" ; t = clock() ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
      markov_chain_information[m].compute_backward_probabilities( transition_matrix, interploidy_transitions ) ;
    }

    Rcpp::Rcout << (double) (clock() - t) << " ms" << endl << "forward-backward posterior decoding and printing\t\t\t" ; t = clock() ;
    std::vector< std::vector< double > > output;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
      markov_chain_information[m].combine_prob( position,
                                                state_list,
                                                chromosomes,
                                                options.output_pulses,
                                                optimum,
                                                output,
                                                m) ;
    }
    Rcpp::Rcout << "converting to NumericMatrix\n"; force_output();
    probabilities = convert_double_vector(output);
  }

  Rcout << (double) (clock() - t) << " ms" << endl ;
  Rcout << "total run time:\t\t\t" << (double) (clock() - total) << " ms" << endl ;

  NumericVector time_estimates;
  for(int i = 0; i < optimum.size(); ++i) {
    time_estimates.push_back(optimum[i].time);
  }

  return Rcpp::List::create(Named("posterior") = probabilities,
                            Named("age_estimate") = time_estimates,
                            Named("age_ll") = lnl);
}

Rcpp::NumericMatrix convert_double_vector(const std::vector< std::vector< double >>& invec) {
  int nrow = invec.size();
  int ncol = invec[0].size();
  Rcpp::NumericMatrix output(nrow, ncol);

  for(int i = 0; i < invec.size(); ++i) {
    for(int j = 0; j < invec[i].size(); ++j) {
      output(i, j) = invec[i][j];
    }
  }
  return output;
}

pulse create_pulse(const Rcpp::NumericVector& p) {
  pulse new_ancestry_pulse;
  new_ancestry_pulse.type = p[1];
  new_ancestry_pulse.time = p[2];
  new_ancestry_pulse.proportion = p[3];

  // if time is set, we are not estimating it
  ///// set time with a negative number to provide the starting guess for this parameter
  if ( new_ancestry_pulse.time > 0 ) {
    new_ancestry_pulse.time_fixed = true ;
  }
  else {
    new_ancestry_pulse.time = new_ancestry_pulse.time * -1 ;
    new_ancestry_pulse.time_fixed = false ;
  }

  // if proporion is set, we are not estimating it
  ////// set proporiton with a negative number to provide the starting guess for this parameter
  if ( new_ancestry_pulse.proportion > 0 ) {
    new_ancestry_pulse.proportion_fixed = true ;
  }
  else {
    new_ancestry_pulse.proportion_fixed = false ;
    new_ancestry_pulse.proportion = -1 * new_ancestry_pulse.proportion ;
  }
  return new_ancestry_pulse;
}