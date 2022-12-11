#include "initialization.hpp"

// The calculation of neumann F value
void dFdn_calc(element& elm){
    double _guess_value = 0.0e0;            // Still can be changed
    for (int i = 0; i < elm.num; i++){
        elm.dFdn[i] = _guess_value;
    }
}

// The calculation of dirichlet F value
void F_calc(element& elm){
    //
}

// The calculation of neuman phi value
void dpdn_calc(element& elm){
    double _guess_value = 0.0e0;            // Still can be changed
    for (int i = 0; i < elm.num; i++){
        elm.dFdn[i] = _guess_value;
    }
}

// The calculation of dirichelt phi value
void p_calc(element& elm){
    //
}