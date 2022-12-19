#include "initialization.hpp"

// ================================================================================
// ================================================================================
// The calculation of F boundary value
void initialization::F_val_calc(element& elm, int ID){
    // Define the calculation type
    int F_calc_type = 1;
    
    if (F_calc_type == 1){
        // Constant neumann value dFdn = 0
        double _guess_value = 1.0e-12;            // Still can be changed
        for (int i = 0; i < elm.num; i++){
            elm.dFdn[i] = _guess_value;
        }
    }
    else if (F_calc_type == 2){
        // Method by Dewangga [2021]
        for (int i = 0; i < elm.num; i++){
            elm.dFdn[i] = elm.Tx[i] * elm.xn[i] + elm.Ty[i] * elm.yn[i];
        }
    }
    else if(F_calc_type == 3){
        // The case of dirichlet F = Tx + Ty
        double _guess_value = Par::trac_t_y + Par::trac_r_x;            // Still can be changed
        for (int i = 0; i < elm.num; i++){
            elm.F[i] = _guess_value;
            elm.F_type[i] = false;
        }
    }
    else if(F_calc_type == 4){
        // The case of dirichlet F = Tx + Ty
        double _guess_value = Par::trac_t_y + Par::trac_r_x;            // Still can be changed
        for (int i = 0; i < elm.num; i++){
            elm.F[i] = 0.0e0;
            elm.dFdn[i] = 0.0e0;
        }
    }
    
    // Other Case 
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
        // std::cout << "The node x diff value :" << dpdx_node[i+1] << std::endl;
        // std::cout << "The node y diff value :" << dpdy_node[i+1] << std::endl;
    }

    
    // Update for symmetric extremes (work as well, but still need to change the value of F and phi)
    double max_x = 0, max_y = 0, min_x = 0, min_y = 0;
    for (int i = 0; i < dpdx_node.size(); i++){
        max_x = dpdx_node[i] > max_x ? dpdx_node[i] : max_x;
        max_y = dpdy_node[i] > max_y ? dpdy_node[i] : max_y;
        min_x = dpdx_node[i] < min_x ? dpdx_node[i] : min_x;
        min_y = dpdy_node[i] < min_y ? dpdy_node[i] : min_y;
    }

    for (int i = 0; i < dpdx_node.size(); i++){
        dpdx_node[i] += (max_x + min_x)/(-2.0);
        dpdy_node[i] += (max_y + min_y)/(-2.0);
    }

    // Assign the value of dpdn
    for (int i = 0; i < elm.num; i++){
        _dpdx = 0.5 * (dpdx_node[i+1] + dpdx_node[i]);
        _dpdy = 0.5 * (dpdy_node[i+1] + dpdy_node[i]);
        elm.dpdn[i] = _dpdx * elm.xn[i] + _dpdy * elm.yn[i];
        // std::cout << "The node x value :" << _dpdx << std::endl;
        // std::cout << "The normal x :" << elm.xn[i] << std::endl;
        // std::cout << "The node y value :" << _dpdy << std::endl;
        // std::cout << "The normal y :" << elm.yn[i] << std::endl;
        // if (elm.dpdn[i] == 0){
        //     elm.p[i] = 0;
        //     elm.p_type[i] = false;
        // }
    }

    
    // Calculate the phi (If the algorithm is known)
    /* Code lies here ... */

    bool given_P = false;
    if (given_P == true){
        for (int i = 0; i < elm.num; i++){
            elm.p[i] = Par::trac_t_y/2.0 * elm.xm[i] * elm.xm[i] + Par::trac_r_x/2.0 * elm.ym[i] * elm.ym[i];
        }
    }
}


// ================================================================================
// ================================================================================
// The calculation of F boundary value
void initialization::T_val_calc(element& elm, int ID){
    // Marking whether a base or inner geometry
    bool base;
    if (ID >= 0){
        base = false;
    }else{
        base = true;
    }

    // Define the calculation type
    int T_calc_type = 2;
    
    if (T_calc_type == 1){
        // Constant neumann value dTdn = 0
        double _guess_value = 0.0e0;            // Still can be changed
        for (int i = 0; i < elm.num; i++){
            elm.dTdn[i] = _guess_value;
        }
    }
    else if (T_calc_type == 2){
        // Constant dirichlet value T = prescribed
        if (base == true){
            for (int i = 0; i < elm.num; i++){
                elm.T[i] = Par::Temp;
                elm.T_type[i] = false;
            }
        }else{
            for (int i = 0; i < elm.num; i++){
                elm.T[i] = Par::In_temp[ID];
                elm.T_type[i] = false;
            }
        }
    }
    else if (T_calc_type == 3){
        // Find the position at each boundary
        double max_x = 0, max_y = 0, min_x = 0, min_y = 0;
        for (int i = 0; i < elm.num; i++){
            max_x = elm.xm[i] > max_x ? elm.xm[i] : max_x;
            max_y = elm.ym[i] > max_y ? elm.ym[i] : max_y;
            min_x = elm.xm[i] < min_x ? elm.xm[i] : min_x;
            min_y = elm.ym[i] < min_y ? elm.ym[i] : min_y;
        }

        // Constant dirichlet value T = prescribed
        for (int i = 0; i < elm.num; i++){
            if (elm.xm[i] == max_x){
                // Right boundary
                elm.T[i] = Par::temp_r;
                elm.T_type[i] = false; // (dirichet)
            }
            else if (elm.xm[i] == min_x){
                // Left boundary
                elm.T[i] = Par::temp_l;
                elm.T_type[i] = false; // (dirichet)
            }
            else if (elm.ym[i] == min_y){
                // Bottom boundary
                elm.T[i] = Par::temp_b;
                elm.T_type[i] = true; // (neumann)
            }
            else if (elm.ym[i] == max_y){
                // Top boundary
                elm.T[i] = Par::temp_t;
                elm.T_type[i] = true; // (neumann)
            }
            
        }
    }

    // Other Case 
    /* Code lies here ... */
}
