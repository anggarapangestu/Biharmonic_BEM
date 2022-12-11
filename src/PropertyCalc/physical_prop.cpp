#include "physical_prop.hpp"

void propertyCalc::calculate_property(intElement& intNode){
    // Calculate the stresses
    // Initialization stress calculation log
    printf("\nStress calculation ...\n");
    clock_t _time = clock();
    this->calculate_stress(intNode);
    
    // Out plane stress calculation (z-direction)
    if (Par::opt_sim_type == 1){        // Plane strain
        // Calculate the stress at z direction
        intNode.s_zz.resize(intNode.num,0.0e0);
        for (int i = 0; i < intNode.num; i++){
            intNode.s_zz[i] = Par::nu * (intNode.s_xx[i] + intNode.s_yy[i]);
        }
    }
    else if (Par::opt_sim_type == 2){   // Plane stress
        // No stresses at z direction
    }

    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Stress calculation comp. time      [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
    
    // Calculate the strains
    // Initialization strain calculation log
    printf("\nStrain calculation ...\n");
    _time = clock();
    this->calculate_strain(intNode);
    
    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Strain calculation comp. time      [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);

    // Calculate the displacements
    // Initialization strain calculation log
    printf("\nDisplacement calculation ...\n");
    _time = clock() - _time;
    this->calculate_disp(intNode);
    
    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Displacement calculation\n");
    printf("    comp. time                         [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}

// ======================================================================
// ======================================================================
// Calculation of stress by using LSMPS
void propertyCalc::calculate_stress(intElement& intNode){
    // Resize the stress properties
    intNode.s_xx.resize(intNode.num, 0.0e0);
    intNode.s_yy.resize(intNode.num, 0.0e0);
    intNode.t_xy.resize(intNode.num, 0.0e0);

    // Evaluate the neigbor
	std::vector<std::vector<int>> ngh_ID = this->neighEval.link_list(intNode.num, intNode.x, intNode.y, Par::R_s);
	
	// Performing LSMPS calculation
	lsmpsa_phi.set_LSMPS(intNode.x, intNode.y, intNode.s, intNode.phi, ngh_ID);
	std::vector<double> d2pd2x  = lsmpsa_phi.get_d2d2x();
	std::vector<double> d2pdxdy = lsmpsa_phi.get_d2dxdy();
    std::vector<double> d2pd2y  = lsmpsa_phi.get_d2d2y();

    // Assign the stresses value
    for (int i = 0; i < intNode.num; i++){
        intNode.s_xx[i] = d2pd2y[i];
        intNode.s_yy[i] = d2pd2x[i];
        intNode.t_xy[i] = - d2pdxdy[i];
    }
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
    // In development ...
}
