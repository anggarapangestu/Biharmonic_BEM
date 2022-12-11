#include "initialization.hpp"

// ================================================================================
// ================================================================================
// The calculation of F boundary value
void initialization::F_val_calc(element& elm, int ID){
    double _guess_value = 0.0e0;            // Still can be changed
    for (int i = 0; i < elm.num; i++){
        elm.dFdn[i] = _guess_value;
    }

    // Other case if dirichlet obtained
    /* Code lies here ... */
}

// ================================================================================
// ================================================================================
// The calculation of phi boundary value
void initialization::p_val_calc(element& elm, int ID){
    // Calculate the dpdn
    double dpdx_base = 0.0e0;            // Still can be changed
    double dpdy_base = 0.0e0;            // Still can be changed
    
    // Marking whether a base or inner geometry
    int dir;
    bool base;
    if (ID >= 0){
        base = false;
        dir = -1;
    }else{
        base = true;
        dir = 1;
    }

    // Initialize the internal variable
    std::vector<double> dpdx_node(elm.num + 1,0.0e0);
    std::vector<double> dpdy_node(elm.num + 1,0.0e0);
    dpdx_node[0] = {dpdx_base};
    dpdy_node[0] = {dpdy_base};
    double _dpdx;
    double _dpdy;
    
    // Calculate the partial differential of phi at each node
    for (int i = 0; i < elm.num; i++){
        dpdx_node[i+1] = dpdx_node[i] - (dir * elm.L[i] * elm.Ty[i]);
        dpdy_node[i+1] = dpdy_node[i] + (dir * elm.L[i] * elm.Tx[i]);
    }
    
    // Check for the value of final and initial
    std::cout << "[DEBUGGING] Comparison\n";
    std::cout << "The initial value: \n";
    std::cout << " - dpdx init value: " << dpdx_node[0] << std::endl;
    std::cout << " - dpdy init value: " << dpdy_node[0] << std::endl;
    std::cout << "\nThe final value: \n";
    std::cout << " - dpdx final value: " << dpdx_node[dpdx_node.size() - 1] << std::endl;
    std::cout << " - dpdy final value: " << dpdy_node[dpdy_node.size() - 1] << std::endl;

    // Assign the value of dpdn
    for (int i = 0; i < elm.num; i++){
        _dpdx = 0.5 * (dpdx_node[i+1] + dpdx_node[i]);
        _dpdy = 0.5 * (dpdy_node[i+1] + dpdy_node[i]);
        elm.dpdn[i] = _dpdx * elm.xn[i] + _dpdy * elm.yn[i];
    }

    // Calculate the phi (If the algorithm is known)
    /* Code lies here ... */
}