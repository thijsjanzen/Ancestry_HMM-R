#ifndef __PRINT_USAGE_H
#define __PRINT_USAGE_H

void print_usage() {

    Rcpp::Rcout << endl << endl << "ancestry_hmm usage:" << endl << endl ;
    Rcpp::Rcout << "\trequired:" << endl ;
    Rcpp::Rcout << "\t\t-i [string]\t\tinput file name" << endl ;
    Rcpp::Rcout << "\t\t-s [string]\t\tsample id and ploidy file" << endl ;
    Rcpp::Rcout << "\t\t-a [int] [float] [float] ..." << endl ;
    Rcpp::Rcout << "\t\t\tnumber of ancestral populations and ancestry proportion attributable to each" << endl ;
    Rcpp::Rcout << "\t\t-p [int] [int] [float]" << endl ;
    Rcpp::Rcout << "\t\t\tancestry pulse with format, ancestral population, time," << endl ;
    Rcpp::Rcout << "\t\t\tand proportion of final ancestry from this pulse" << endl ;
    Rcpp::Rcout << "\t\t\tnegative time or proportions indicate that parameters are to be estimated" << endl << endl ;

    Rcpp::Rcout << "\toptional:" << endl ;
    Rcpp::Rcout << "\t\t--help\t\t\tprint this help statement" << endl ;
    Rcpp::Rcout << "\t\t--ne [int]\t\teffective population size of the admixed population" << endl ;
    Rcpp::Rcout << "\t\t-g\t\t\tsamples are specified with genotypes rather than read counts" << endl ;
    Rcpp::Rcout << "\t\t--precision [int]\tmodify float and double precision to int" << endl ;
    Rcpp::Rcout << "\t\t-v\t\t\tviterbi decoding" << endl ;
    Rcpp::Rcout << "\t\t-b [int] [int]\t\tnumber of bootstraps and bootstrap block size in number of SNPs" << endl ;
    Rcpp::Rcout << "\t\t--tmax [int]\t\tmaximum time of an admixture pulse" << endl ;
    Rcpp::Rcout << "\t\t--tmin [int]\t\tminimum time of an admixture pulse" << endl ;
    Rcpp::Rcout << "\t\t--tolerance [float]\tdistance in lnL units to just convergence" << endl ;
    Rcpp::Rcout << "\t\t-e [float]\t\terror rates" << endl ;
    Rcpp::Rcout << "\t\t-E\t\t\tsite specific error rates are included" << endl ;
    Rcpp::Rcout << "\t\t--fix\t\t\tancestral allele frequencies are certain" << endl << endl ;

    Rcpp::Rcout << "\toptional and relevant only for multiple pulse models:" << endl ;
    Rcpp::Rcout << "\t\t--output-ancestry\toutput ancestry posteriors rather than pulses" << endl ;
    Rcpp::Rcout << "\t\t-r [int]\t\tnumber of random restarts during nelder-mead optimization" << endl ;
    Rcpp::Rcout << "\t\t--pmax [int]\t\tmaximum proportion ancestry in an admixture pulse" << endl ;
    Rcpp::Rcout << "\t\t--pmin [int]\t\tminimum proportion ancestry in an admixture pulse" << endl ;
}

#endif

