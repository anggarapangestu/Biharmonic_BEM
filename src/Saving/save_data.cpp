#include "save_data.hpp"

// Method to display and write the simulation log
void dataSaving::simulation_log(){
    // Simulation Log Initial Parameter Summary Data
    std::cout << std::endl;
    std::cout << "#=================================================#\n";
    std::cout << "+---------------- SIMULATION LOG -----------------+\n";
    std::cout << "#=================================================#\n";
    std::cout << "<!> Computation Message Log\n";
	
	// Initial flow parameters summary data
	printf("\n+-------------- Material Properties --------------+\n");
	printf("Elasticity (E)                      : %8.2f GPA\n", Par::E/1.0e9);
	printf("Poisson Ratio (nu)                  : %8.2f [-]\n", Par::nu);
	printf("Lame's constant mu                  : %8.2f GPA\n", Par::mu/1.0e9);
	printf("Lame's constant lambda              : %8.2f GPA\n", Par::lambda/1.0e9);
	printf("+-------------------------------------------------+\n");

    // Simulation parameter summary data
	printf("\n+---------- Simulation Parameters Data -----------+\n");
    if (Par::opt_sim_type == 1){
        printf("Simulation type                     :   Biharmonic \n");
        if (Par::opt_biharmonic_type == 1){
            printf("Biharmonic simulation type          : Plane Strain \n");
        }else if (Par::opt_biharmonic_type == 2){
            printf("Biharmonic simulation type          : Plane Stress \n");
        }
    }else if (Par::opt_sim_type == 2){
        printf("Simulation type                     :  Temperature \n");
    }
    printf("Body option                         :       type %d \n", Par::G_type);
	printf("Internal node distribution option   :       type %d \n", Par::opt_int_init);
    printf("BEM calculation type                :       type %d \n", Par::opt_BEM);
	if (Par::flag_cylinder == true){
        printf("Cylindrical coordinate calculation  :    turned on\n");
    }
	printf("+-------------------------------------------------+\n");

	// Additional calculation data
	printf("\n+---------------- Additional Data ----------------+\n");
	printf("Geometry x length                   :   %8.4f m\n", Par::dom_Lx);
    printf("Geometry y length                   :   %8.4f m\n", Par::dom_Ly);
    printf("Base internal node spacing          :   %8.4f m\n", Par::spc);
    printf("Base panel length                   :   %8.4f m\n", Par::len);
    printf("Number of boundary surface          :   %8d \n", Par::N_Gin + 1);
    printf("Support domain radius factor        :   %8.2f \n", Par::R_s);
	printf("+-------------------------------------------------+\n");

    // Cancel the saving procedure if flag is closed
	if (Par::flag_save_log == false){return;}
	
	// open the write file
	this->save.open("output/simulation_log.dat");
	this->save << "#=================================================#\n"
               << "+---------------- SIMULATION LOG -----------------+\n"
               << "#=================================================#\n\n";
	
	this->save << "+-------------- Material Properties --------------+\n";
	this->save << std::fixed << std::setprecision(2)
               << "Elasticity (E)                      : "; save.width(8); save << std::right << Par::E/1.0e9  << " GPa\n"
               << "Poisson Ratio (nu)                  : "; save.width(8); save << std::right << Par::nu << " [-]\n"
               << "Lame's constant mu                  : "; save.width(8); save << std::right << Par::mu/1.0e9 << " GPa\n"
               << "Lame's constant lambda              : "; save.width(8); save << std::right << Par::lambda/1.0e9 << " GPa\n"
               << "+-------------------------------------------------+\n\n"
               ;

	this->save << "+----------- Simulation Setting Option -----------+\n";
    if (Par::opt_sim_type == 1){
        this->save << "Simulation type                     :   " << "Biharmonic\n";
        if (Par::opt_biharmonic_type == 1){
            this->save << "Biharmonic simulation type          : " << "Plane Strain\n";
        }else if (Par::opt_biharmonic_type == 2){
            this->save << "Biharmonic simulation type          : " << "Plane Stress\n";
        }
    }else if (Par::opt_sim_type == 2){
        this->save << "Simulation type                     :  " << "Temperature\n";
    }
    this->save << "Body option                         :       " << "type " << Par::G_type << "\n"
               << "Internal node distribution option   :       " << "type " << Par::opt_int_init << "\n"
               << "BEM calculation option              :       " << "type " << Par::opt_BEM << "\n"
               ;
    if (Par::flag_cylinder == true){
        this->save << "Cylindrical coordinate calculation  :    " << "turned on\n";
    }
    this->save << "+-------------------------------------------------+\n\n";

	this->save << std::fixed << std::setprecision(4)
               << "+---------------- Additional Data ----------------+\n"
               << "Domain x length                     :   "; save.width(8); save << std::right << Par::dom_Lx << " m\n"
               << "Domain x length                     :   "; save.width(8); save << std::right << Par::dom_Ly << " m\n"
               << "Base internal node spacing          :   "; save.width(8); save << std::right << Par::spc << " m\n"
               << "Base panel length                   :   "; save.width(8); save << std::right << Par::len << " m\n"
               << "Number of boundary surface          :   "; save.width(8); save << std::right << Par::N_Gin + 1 << "\n"
               << "Support domain radius factor        :   "; save.width(8); save << std::right << Par::R_s << "\n"
               << "+-------------------------------------------------+\n\n"
               ;
	
	this->save << std::setprecision(-1);
	
	// close the write file
	this->save.close();
}

// =====================================================================================
// =====================================================================================
// Method to write the internal data properties
void dataSaving::write_internal_data(const intElement& intElm){
	// Cancel the saving procedure if flag is closed
    if (Par::flag_save_Int_Node == false){return;}

    // saving starting log
    printf("\nSaving the internal node data ...\n");
    
    // write file name
	std::string name;
		
	// open the write file
	name.append("output/internal_node_data.csv");
	this->save.open(name.c_str());

    // write the data header
	this->save << "" << "x" 
               << "," << "y"
               << "," << "R"
               << "," << "s"
               << "," << "phi"
               << "," << "phi_an"
               << "," << "sigma_xx"
               << "," << "sigma_yy" 
               << "," << "tau_xy" 
               << "," << "epsilon_xx" 
               << "," << "epsilon_yy"
               << "," << "epsilon_xy";
            //    << "," << "u"
            //    << "," << "v";
    if (Par::flag_cylinder == true){
    this->save << "," << "sigma_rr"
               << "," << "sigma_tt" 
               << "," << "tau_rt" 
               << "," << "epsilon_rr" 
               << "," << "epsilon_tt"
               << "," << "epsilon_rt";
    }
    this->save << "\n";

	// write each node data
    for (int i = 0; i < intElm.num; i++){
        this->save << "" << intElm.x[i]
                   << "," << intElm.y[i]
                   << "," << intElm.R[i]
                   << "," << intElm.s[i]
                   << "," << intElm.phi[i]
                   << "," << intElm.phi_an[i]
                   << "," << intElm.s_xx[i]
                   << "," << intElm.s_yy[i]
                   << "," << intElm.t_xy[i]
                   << "," << intElm.e_xx[i]
                   << "," << intElm.e_yy[i]
                   << "," << intElm.e_xy[i];
                //    << "," << intElm.u[i]
                //    << "," << intElm.v[i];
    if (Par::flag_cylinder == true){
    this->save << "," << intElm.s_rr[i]
               << "," << intElm.s_tt[i]
               << "," << intElm.t_rt[i]
               << "," << intElm.e_rr[i]
               << "," << intElm.e_tt[i]
               << "," << intElm.e_rt[i];
    }
    this->save << "\n";
	}
    
    // close the write file
	this->save.close();

    printf("<+> Done saving internal node data\n");
}

// Method to write the internal data properties of temperature
void dataSaving::write_internal_data_temp(const intElement& intElm){
	// Cancel the saving procedure if flag is closed
    if (Par::flag_save_Int_Node == false){return;}

    // saving starting log
    printf("\nSaving the internal node data ...\n");
    
    // write file name
	std::string name;
		
	// open the write file
	name.append("output/internal_node_data_temp.csv");
	this->save.open(name.c_str());

    // write the data header
	this->save << "" << "x" 
               << "," << "y"
               << "," << "R"
               << "," << "s"
               << "," << "T"
               << "\n";

	// write each node data
    for (int i = 0; i < intElm.num; i++){
        this->save << "" << intElm.x[i]
                   << "," << intElm.y[i]
                   << "," << intElm.R[i]
                   << "," << intElm.s[i]
                   << "," << intElm.T[i]
                   << "\n";
	}
    
    // close the write file
	this->save.close();

    printf("<+> Done saving internal node data\n");
}

// Method to write the boundary element data
void dataSaving::write_BEM_data(const element& elm, const std::vector<element>& in_elm){
    // Cancel the saving procedure if flag is closed
    if (Par::flag_save_BEM == false){return;}
    
    // saving starting log
    printf("\nSaving the boundary element data ...\n");
    
    // write file name
	std::string name;
		
	// open the write file
	name.append("output/boundary_panel_data.csv");
	this->save.open(name.c_str());

    // write the data header
	this->save << "" << "xm" 
               << "," << "ym"
               << "," << "length"
               << "," << "xnormal"
               << "," << "ynormal"
               << "," << "Tx"
               << "," << "Ty" 
               << "," << "F"
               << "," << "phi" 
               << "," << "dF_dn" 
               << "," << "dphi_dn" 
               << "\n";

	// write each element data

    // Base geometry
    for (int i = 0; i < elm.num; i++){
        this->save << "" << elm.xm[i]
                   << "," << elm.ym[i]
                   << "," << elm.L[i]
                   << "," << elm.xn[i]
                   << "," << elm.yn[i]
                   << "," << elm.Tx[i]
                   << "," << elm.Ty[i]
                   << "," << elm.F[i]
                   << "," << elm.p[i]
                   << "," << elm.dFdn[i]
                   << "," << elm.dpdn[i]
                   << "\n";
	}

    // Inner geometry
    for(int ID = 0; ID < Par::N_Gin; ID++){
        for (int i = 0; i < in_elm[ID].num; i++){
            this->save << "" << in_elm[ID].xm[i]
                    << "," << in_elm[ID].ym[i]
                    << "," << in_elm[ID].L[i]
                    << "," << in_elm[ID].xn[i]
                    << "," << in_elm[ID].yn[i]
                    << "," << in_elm[ID].Tx[i]
                    << "," << in_elm[ID].Ty[i]
                    << "," << in_elm[ID].F[i]
                    << "," << in_elm[ID].p[i]
                    << "," << in_elm[ID].dFdn[i]
                    << "," << in_elm[ID].dpdn[i]
                    << "\n";
        }
    }
    
    // close the write file
	this->save.close();

    printf("<+> Done saving boundary element data\n");
}

// Method to write the boundary element data
void dataSaving::write_BEM_data_temp(const element& elm, const std::vector<element>& in_elm){
    // Cancel the saving procedure if flag is closed
    if (Par::flag_save_BEM == false){return;}
    
    // saving starting log
    printf("\nSaving the boundary element data ...\n");
    
    // write file name
	std::string name;
		
	// open the write file
	name.append("output/boundary_panel_data_temp.csv");
	this->save.open(name.c_str());

    // write the data header
	this->save << "" << "xm" 
               << "," << "ym"
               << "," << "length"
               << "," << "xnormal"
               << "," << "ynormal"
               << "," << "T"
               << "," << "dT_dn" 
               << "\n";

	// write each node data

    // Base geometry
    for (int i = 0; i < elm.num; i++){
        this->save << "" << elm.xm[i]
                   << "," << elm.ym[i]
                   << "," << elm.L[i]
                   << "," << elm.xn[i]
                   << "," << elm.yn[i]
                   << "," << elm.T[i]
                   << "," << elm.dTdn[i]
                   << "\n";
	}

    // Inner geometry
    for(int ID = 0; ID < Par::N_Gin; ID++){
        for (int i = 0; i < in_elm[ID].num; i++){
            this->save << "" << in_elm[ID].xm[i]
                    << "," << in_elm[ID].ym[i]
                    << "," << in_elm[ID].L[i]
                    << "," << in_elm[ID].xn[i]
                    << "," << in_elm[ID].yn[i]
                    << "," << in_elm[ID].T[i]
                    << "," << in_elm[ID].dTdn[i]
                    << "\n";
        }
    }
    
    // close the write file
	this->save.close();

    printf("<+> Done saving boundary element data\n");
}

// =====================================================================================
// =====================================================================================
// Write matrix data
void dataSaving::write_Matrix(const Eigen::MatrixXd& MAT, std::string name){
    // saving starting log
    printf("<+> Saving matrix %s ...\n", name.c_str());
    
	// open the write file
    std::string _name = "output/matrix_";
    _name.append(name.c_str());
	_name.append(".csv");
	this->save.open(_name.c_str());

    // write the matrix data
    for (size_t i = 0; i < MAT.rows(); i++){
        this->save << MAT(i, 0);
        for (size_t j = 1; j < MAT.cols(); j++){
            this->save << "," << MAT(i, j);
        }
        this->save << "\n";
    }
    this->save << "\n";
    
    // close the write file
	this->save.close();
}

// Write vector data
void dataSaving::write_Matrix(const Eigen::VectorXd& VEC, std::string name){
    // saving starting log
    printf("<+> Saving vector %s ...\n", name.c_str());
    
	// open the write file
    std::string _name = "output/vector_";
    _name.append(name.c_str());
	_name.append(".csv");
	this->save.open(_name.c_str());

    // write the matrix data
    for (size_t i = 0; i < VEC.rows(); i++){
        this->save << VEC(i, 0);
        this->save << "\n";
    }
    
    // close the write file
	this->save.close();
}

// =====================================================================================
// =====================================================================================
// Saving method for testing
void dataSaving::save_Test(const intElement & elm){
    // saving starting log
    printf("<+> Saving testing data ... \n");
    
	// open the write file
	this->save.open("output/test.csv");
    this->save << "x,y,theta\n";

    // write the matrix data
    for (size_t i = 0; i < elm.x.size(); i++){
        this->save << elm.x[i] << "," << elm.y[i] << "," << elm.phi[i] << "\n";
    }
    
    // close the write file
	this->save.close();
}