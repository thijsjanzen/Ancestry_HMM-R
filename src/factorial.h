#ifndef __FACTORIAL_H
#define __FACTORIAL_H

#include <vector>

std::vector<double> create_factorial () {

  std::vector<double> facto(1755) ;
	facto[0] = 1 ;
	double fact = 1 ;
	for ( double base = 1 ; base < 1755 ; base ++ ) {
		fact *= base ;
		facto[base] = fact ;
	}
	return facto ;
}

const std::vector<double> factorial_vec = create_factorial() ;

#endif
