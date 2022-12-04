#include "save_data.hpp"

// Method to display and write the simulation log
void dataSaving::simulation_log(){
    // Simulation Log Initial Parameter Summary Data
    std::cout << "#=================================================#\n";
    std::cout << "+---------------- SIMULATION LOG -----------------+\n";
    std::cout << "#=================================================#\n";
    std::cout << "<!> Computation Message Log\n";
	
	// Simulation parameter summary data
	printf("\n+---------- Simulation Parameters Data -----------+\n");
	printf("Body option                           :   type %d \n", Par::opt_geom);
	printf("Initialization option                 :   type %d \n", Par::opt_int_init);
	// printf("Neighbor search option                :   type %d \n", Par::opt_neighbor);
	// printf("Penalization option                   :   type %d \n", Par::opt_pen);
	// printf("Maximum resolution level              :        %d \n", Par::max_level);
	// printf("Core size                             : %8.4f m\n", Par::sigma);
	// printf("Time step                             : %8.4f s\n", Par::dt);
	// printf("Total simulation time                 : %8.2f s\n", Par::simulation_time);
	printf("+-------------------------------------------------+\n");
	
	// Initial flow parameters summary data
	printf("\n+------------- Flow Parameters Data --------------+\n");
	// printf("Reynolds number (RE)           : %10.2f [-]\n", Par::RE);
	// printf("Freestream velocity (U)        : %10.2f m/s\n", Par::U_inf);
	// printf("Fluid density (rho)            : %10.2f kg/m^3\n", Par::RHO);
	// printf("Fluid viscosity (nu)           : %10f m^2/s\n", Par::NU);
	printf("+-------------------------------------------------+\n");

	// Additional calculation data
	printf("\n+---------------- Additional Data ----------------+\n");
	// printf("Courant number (C)                  : %12f \n", Par::Courant);
	// printf("Diffusion number (Phi)              : %12f \n", Par::Diffusion);
	// printf("Turbulent scaling                   : %12f \n", std::ceil(1.0e0 / Par::Courant));
	// printf("Maximum time step (Phi criteria)    : %12f \n", 0.25*Par::sigma*Par::sigma/Par::NU);
	printf("+-------------------------------------------------+\n");

    // Cancel the saving procedure if flag is closed
	if (Par::flag_save_log == false){return;}
	
	// open the write file
	this->save.open("output/Log.dat");
	this->save << "#=================================================#\n"
               << "+---------------- SIMULATION LOG -----------------+\n"
               << "#=================================================#\n\n";
	
	this->save << "+----------- Simulation Flow Parameter -----------+\n";
	this->save << std::fixed << std::setprecision(2)
               //  << "Reynolds number (RE)             : "; save.width(8); save << std::right << Par::RE    << " [-]" << "\n"
               //  << "Freestream velocity (U)          : "; save.width(8); save << std::right << Par::u_inf << " m/s"<< "\n"
               << "+-------------------------------------------------+\n\n"
               ;

	this->save << "+----------- Simulation Setting Option -----------+\n"
               << "Body option                             : " << "type " << Par::opt_geom << "\n"
               << "Initialization option                   : " << "type " << Par::opt_int_init << "\n"
               //  << "Neighbor search option                  : " << "type " << Par::opt_neighbor << "\n"
               //  << "Penalization option                     : " << "type " << Par::opt_pen << "\n"
               ;
	this->save << std::fixed << std::setprecision(4)
               //  << "Number of save data                     : "; save.width(6); save << std::right << Par::nt_data << "\n"
               //  << "Saving step interval                    : "; save.width(6); save << std::right << Par::step_inv << "\n"
               //  << "Core size                               : "; save.width(6); save << std::right << Par::sigma << " m\n"
               //  << "Time step                               : "; save.width(6); save << std::right << Par::dt << " s\n"
               ;
	this->save << std::fixed << std::setprecision(2)
               //  << "Total simulation time                   : "; save.width(6); save << std::right << Par::simulation_time << " s\n"
               << "+-------------------------------------------------+\n\n"
               ;

	this->save << "+----------- Simulation Parameter Data -----------+\n"
               //  << "Domain x length                       : "; save.width(8); save << std::right << Par::lxdom << " m\n"
               //  << "Domain y length                       : "; save.width(8); save << std::right << Par::lydom << " m\n"
               ;
	this->save << std::fixed << std::setprecision(4)
               //  << "Courant number (C)                    : "; save.width(8); save << std::right << Par::Courant << "\n"
               << "+-------------------------------------------------+\n\n"
               ;
	
	this->save << std::setprecision(-1);
	
	// close the write file
	this->save.close();
}

// Method to write the internal data properties
void dataSaving::write_internal_data(const intElement& intElm){
	// saving starting log
    printf("\nSaving the particle state ...\n");
    
    // write file name
	std::string name;
		
	// open the write file
	name.append("output/internal_node_data.csv");
	name.append(".csv");
	this->save.open(name.c_str());

    // write the data header
	this->save << "" << "xp" 
               << "," << "yp"
               << "," << "phi"
               << "," << "sigma_xx"
               << "," << "sigma_yy" 
               << "," << "tau_xy" 
               << "," << "epsilon_xx" 
               << "," << "epsilon_yy"
               << "," << "epsilon_xy"
               << "\n";

	// write each node data
    for (int i = 0; i < intElm.num; i++){
        this->save << "" << intElm.x[i]
                   << "," << intElm.y[i]
                   << "," << intElm.phi[i]
                   << "," << intElm.s_xx[i]
                   << "," << intElm.s_yy[i]
                   << "," << intElm.t_xy[i]
                   << "," << intElm.e_xx[i]
                   << "," << intElm.e_yy[i]
                   << "," << intElm.e_xy[i]
                   << "\n";
	}
    
    // close the write file
	this->save.close();
}

// Method to write the boundary element data
void dataSaving::write_BEM_data(const element& elm){
    // A command
}