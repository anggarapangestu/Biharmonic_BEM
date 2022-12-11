#include "physical_prop.hpp"

void propertyCalc::calculate_property(intElement& intNode){
    // Calculate the stresses
    this->calculate_stress(intNode);
    if (Par::opt_sim_type == 1){    // Plane strain
        // Calculate the stress at z direction
        intNode.s_zz.resize(intNode.num,0.0e0);
        for (int i = 0; i < intNode.num; i++){
            intNode.s_zz[i] = Par::nu * (intNode.s_xx[i] + intNode.s_yy[i]);
        }
    }
    else{   // Plane stress
        // No stresses at z direction
    }

    // Calculate the strains
    this->calculate_strain(intNode);

    // Calculate the displacements
    this->calculate_disp(intNode);
}

// ======================================================================
// ======================================================================
// Calculation of stress by using LSMPS
void propertyCalc::calculate_stress(intElement& intNode){
    // The main code
}

void propertyCalc::calculate_strain(intElement& intNode){
    // Resize the strain properties
    intNode.e_xx.resize(intNode.num, 0.0e0);
    intNode.e_yy.resize(intNode.num, 0.0e0);
    intNode.e_xy.resize(intNode.num, 0.0e0);
    
    // Strain are calculated base on the simulation type

    // Plane strain
    if (Par::opt_sim_type == 1){
        double A = Par::lambda + 2 * Par::mu;
        double B = Par::lambda;
        double D = 1 / (A*A - B*B);

        for (int i = 0; i < intNode.num; i++){
            intNode.e_xx[i] = D * (A * intNode.s_xx[i] - B * intNode.s_yy[i]);
            intNode.e_yy[i] = D * (A * intNode.s_yy[i] - B * intNode.s_xx[i]);
            intNode.e_xy[i] = intNode.t_xy[i] / (2.0 * Par::mu);
        }
    }
    
    // Plane stress
    else if (Par::opt_sim_type == 2){
        // Strain at z direction
        intNode.e_zz.resize(intNode.num, 0.0e0);

        // Strain calculation at z direction
        for (int i = 0; i < intNode.num; i++){
            intNode.e_xx[i] = (1/Par::E) * (intNode.s_xx[i] - Par::nu * intNode.s_yy[i]);
            intNode.e_yy[i] = (1/Par::E) * (intNode.s_yy[i] - Par::nu * intNode.s_xx[i]);
            intNode.e_zz[i] = (Par::nu/(Par::nu - 1)) * (intNode.e_xx[i] + intNode.e_yy[i]);
            intNode.e_xy[i] = ((1 + Par::nu)/Par::E) * intNode.t_xy[i];
        }
    }
}

void propertyCalc::calculate_disp(intElement& intNode)
{
    // In development
}
