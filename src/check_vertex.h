#ifndef __CHECK_VERTEX_H
#define __CHECK_VERTEX_H

/// check vertex to make sure it's within the parameter bounds
/// then use random box method to select points within the admissable parameter space
/// hope to prevent collapsing the edge by using random box, but boundary cases should be checked carefully

/// rng stuff
double range_term = 0.05 ;

/// perform check
void check_vertex( vector<pulse> &vertex, cmd_line &options ) {

    std::random_device rd;
    std::mt19937 rndgen(rd());
    std::uniform_real_distribution<double> rand_dist(0, 1);


    for ( int p = 0 ; p < vertex.size() ; p ++ ) {
        if ( vertex[p].time_fixed == false ) {
            if ( vertex[p].time > options.t_max ) {
                vertex[p].time = options.t_max - (rand_dist(rndgen))*range_term * ( options.t_max - options.t_min ) ;
            }
            else if ( vertex[p].time < options.t_min ) {
                vertex[p].time = options.t_min + (rand_dist(rndgen))*range_term * ( options.t_max - options.t_min ) ;
            }
        }
        if ( vertex[p].proportion_fixed == false ) {
            if ( vertex[p].fraction_of_remainder > options.p_max ) {
                vertex[p].fraction_of_remainder = options.p_max - (rand_dist(rndgen))*range_term * ( options.p_max - options.p_min ) ;
            }
            else if ( vertex[p].fraction_of_remainder < options.p_min ) {
                vertex[p].fraction_of_remainder = options.p_min + (rand_dist(rndgen))*range_term * ( options.p_max - options.p_min ) ;
            }
        }
    }
}

#endif