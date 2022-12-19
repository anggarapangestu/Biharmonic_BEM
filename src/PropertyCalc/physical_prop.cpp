#include "physical_prop.hpp"
#include <fstream>

// =====================================================================================
// =====================================================================================
// The property calculation manager
void propertyCalc::calculate_property(intElement& intNode){
    // ====================== STRESS CALCULATION ======================
    // ****************************************************************
    // Initialization stress calculation log
    printf("\nStress calculation ...\n");
    clock_t _time = clock();

    // Stress calculation
    this->calculate_stress(intNode);
    
    // Out plane stress calculation (z-direction)
    // Plane strain
    if (Par::opt_biharmonic_type == 1){
        // Calculate the stress at z direction
        intNode.s_zz.resize(intNode.num,0.0e0);
        for (int i = 0; i < intNode.num; i++){
            intNode.s_zz[i] = Par::nu * (intNode.s_xx[i] + intNode.s_yy[i]);
        }
    }
    // Plane stress
    else if (Par::opt_biharmonic_type == 2){
        // No stresses at z direction
    }

    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Stress calculation comp. time      [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
    

    // ====================== STRAIN CALCULATION ======================
    // ****************************************************************
    // Initialization strain calculation log
    printf("\nStrain calculation ...\n");
    _time = clock();

    // Strain calculation
    this->calculate_strain(intNode);
    
    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Strain calculation comp. time      [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);


    // =================== DISPLACEMENT CALCULATION ===================
    // ****************************************************************
    // Initialization strain calculation log
    printf("\nDisplacement calculation ...\n");
    _time = clock();
    
    // Displacement calculation
    this->calculate_disp(intNode);
    
    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Displacement calculation\n");
    printf("    comp. time                         [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);

    // ============ COORDINATE TRANSFORMATION CALCULATION =============
    // ****************************************************************
    // Calculate the properties transformation
    if (Par::flag_cylinder == true){
        // Initialization of transformation calculation log
        printf("\nTransformation calculation ...\n");
        _time = clock();
        
        // Coordinate transformation calculation
        this->coor_transform(intNode);
        
        // Displaying the computational time
        _time = clock() - _time;
        printf("<-> Transformation calculation\n");
        printf("    comp. time                         [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
    }
}

// For TESTING !!!
// The analytical soltuion of phi biaxial
void propertyCalc::phi_analytic_biaxial(intElement& intNode){
    /* The biaxial analytical solution for phi
        phi(x,y) = Ax^2 + Bxy + Cy^2
        With:
            A = s_yy/2
            B = -s_xy
            C = s_xx/2
    */
    
    // Internal variables definition
    double A, B, C;
    A = Par::trac_t_y / 2.0;
    B = 0.0;
    C = Par::trac_r_x / 2.0;
    
    // The phi analytical solution calculation and assignment
    intNode.phi_an.resize(intNode.num,0.0e0);
    for (size_t i = 0; i < intNode.num; i++){
        intNode.phi_an[i] =   A * std::pow(intNode.x[i],2) 
                            + B * intNode.x[i]*intNode.y[i]
                            + C * std::pow(intNode.y[i],2);
    }

}

// ======================================================================
// ======================================================================
// Calculation of stress from airy stress by using LSMPS
void propertyCalc::calculate_stress(intElement& intNode){
    // Resize the stress properties
    intNode.s_xx.resize(intNode.num, 0.0e0);
    intNode.s_yy.resize(intNode.num, 0.0e0);
    intNode.t_xy.resize(intNode.num, 0.0e0);

    // Evaluate the neigbor
	std::vector<std::vector<int>> ngh_ID = this->neighEval.link_list(intNode.num, intNode.x, intNode.y, Par::R_s);

    // WRITE the neighbor data
    if(Par::flag_save_Neigh == true){
        std::ofstream save;
        save.open("output/neighbor_dat.csv");

        // Write the data header
        save << "" << "x" 
            << "," << "y"
            << "," << "ngh"
            << "\n";

        // Write each node data
        for (int i = 0; i < intNode.num; i++){
            save << "" << intNode.x[i]
                << "," << intNode.y[i]
                << "," << ngh_ID[i][0];
            for (int j = 1; j < ngh_ID[i].size(); j++){
            save << ";" << ngh_ID[i][j];
            }
            save<< "\n";
        }
    }
	
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

// Calculation of srain using plane strain/stress relation
void propertyCalc::calculate_strain(intElement& intNode){
    // Resize the strain properties
    intNode.e_xx.resize(intNode.num, 0.0e0);
    intNode.e_yy.resize(intNode.num, 0.0e0);
    intNode.e_xy.resize(intNode.num, 0.0e0);
    
    // Strain are calculated base on the simulation type
    // *************************************************
    // PLANE STRAIN
    if (Par::opt_biharmonic_type== 1){
        double A = Par::lambda + 2 * Par::mu;
        double B = Par::lambda;
        double D = 1 / (A*A - B*B);

        for (int i = 0; i < intNode.num; i++){
            intNode.e_xx[i] = D * (A * intNode.s_xx[i] - B * intNode.s_yy[i]);
            intNode.e_yy[i] = D * (A * intNode.s_yy[i] - B * intNode.s_xx[i]);
            intNode.e_xy[i] = intNode.t_xy[i] / (2.0 * Par::mu);
        }
    }
    
    // PLANE STRESS
    else if (Par::opt_biharmonic_type == 2){
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

// Calculation of displacement
void propertyCalc::calculate_disp(intElement& intNode)
{
    // Still in development ...
    // See you soon :)
}

// =====================================================================================
// =====================================================================================
// Transformation stress and strain tensor into cylindrical coordinate
void propertyCalc::coor_transform(intElement& intNode){
    // Resize each vector at new coordinate
    intNode.s_rr.resize(intNode.num,0.0e0);
    intNode.s_tt.resize(intNode.num,0.0e0);
    intNode.t_rt.resize(intNode.num,0.0e0);
    intNode.e_rr.resize(intNode.num,0.0e0);
    intNode.e_tt.resize(intNode.num,0.0e0);
    intNode.e_rt.resize(intNode.num,0.0e0);
    
    // Variable for internal calculation
    double s11,s12,s22,e11,e12,e22;
    double theta, _sin2t, _cos2t;
    for (int i = 0; i < intNode.num; i++){
        // Assign the value at each evaluation
        theta = std::atan2(intNode.y[i],intNode.x[i]);
        _sin2t = std::sin(2 * theta);
        _cos2t = std::cos(2 * theta);
        s11 = intNode.s_xx[i];
        s12 = intNode.t_xy[i];
        s22 = intNode.s_yy[i];
        e11 = intNode.e_xx[i];
        e12 = intNode.e_xy[i];
        e22 = intNode.e_yy[i];

        // Calculation of coordinate transformation
        intNode.s_rr[i] = (s11 + s22)/2.0 + ((s11 - s22)/2.0)*_cos2t + s12*_sin2t;
        intNode.t_rt[i] = - ((s11 - s22)/2.0)*_sin2t + s12*_cos2t;
        intNode.s_tt[i] = (s11 + s22)/2.0 - ((s11 - s22)/2.0)*_cos2t - s12*_sin2t;
        intNode.e_rr[i] = (e11 + e22)/2.0 + ((e11 - e22)/2.0)*_cos2t + e12*_sin2t;;
        intNode.e_rt[i] = - ((e11 - e22)/2.0)*_sin2t + e12*_cos2t;;
        intNode.e_tt[i] = (e11 + e22)/2.0 - ((e11 - e22)/2.0)*_cos2t - e12*_sin2t;;
    }
}

