#ifndef __SUBSAMPLE_H
#define __SUBSAMPLE_H

void subsample_reads ( double &c1, double &c2 ) {
    std::random_device rd;
    std::mt19937 rndgen(rd());
    std::uniform_real_distribution<double> rand_dist(0, 1);

    while ( c1 + c2 > 170 ) {
        double r = rand_dist(rndgen) ;
        if ( r < c1/(c1+c2) ) {
            c1 -- ;
        }
        else {
            c2 -- ;
        }
    }
}

#endif
