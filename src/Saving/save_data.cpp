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
	printf("Reynolds number (RE)           : %10.2f [-]\n", Par::RE);
	printf("Freestream velocity (U)        : %10.2f m/s\n", Par::U_inf);
	printf("Fluid density (rho)            : %10.2f kg/m^3\n", Par::RHO);
	printf("Fluid viscosity (nu)           : %10f m^2/s\n", Par::NU);
	printf("+-------------------------------------------------+\n");

	// Additional calculation data
	printf("\n+---------------- Additional Data ----------------+\n");
	printf("Courant number (C)                  : %12f \n", Par::Courant);
	printf("Diffusion number (Phi)              : %12f \n", Par::Diffusion);
	printf("Turbulent scaling                   : %12f \n", std::ceil(1.0e0 / Par::Courant));
	printf("Maximum time step (Phi criteria)    : %12f \n", 0.25*Par::sigma*Par::sigma/Par::NU);
	printf("+-------------------------------------------------+\n");

	// Cancel the saving procedure if flag is closed
	if (Par::flag_save_log == false){return;}
	
	// Log saving before simulation
	std::ofstream data;
	data.open("output/Parameter.dat");
	data << "#=================================================#\n"
	     << "+---------------- SIMULATION LOG -----------------+\n"
		 << "#=================================================#\n\n";
	
	data << "+----------- Simulation Flow Parameter -----------+\n";
	data << std::fixed << std::setprecision(2)
	    //  << "Reynolds number (RE)             : "; data.width(8); data << std::right << Par::RE    << " [-]" << "\n"
		//  << "Freestream velocity (U)          : "; data.width(8); data << std::right << Par::u_inf << " m/s"<< "\n"
 		//  << "Fluid density (rho)              : "; data.width(8); data << std::right << Par::RHO   << " kg/m^3" << "\n"
		//  << "Fluid viscosity (nu)             : "; data.width(8); data << std::right << Par::NU    << " m^2/s" << "\n"
		 << "+-------------------------------------------------+\n\n"
         ;

	data << "+----------- Simulation Setting Option -----------+\n"
	     << "Body option                             : " << "type " << Par::opt_geom << "\n"
		 << "Initialization option                   : " << "type " << Par::opt_int_init << "\n"
		//  << "Neighbor search option                  : " << "type " << Par::opt_neighbor << "\n"
		//  << "Penalization option                     : " << "type " << Par::opt_pen << "\n"
         ;
	data << std::fixed << std::setprecision(4)
	    //  << "Number of save data                     : "; data.width(6); data << std::right << Par::nt_data << "\n"
		//  << "Saving step interval                    : "; data.width(6); data << std::right << Par::step_inv << "\n"
		//  << "Support radius factor                   : "; data.width(6); data << std::right << Par::r_sup << "\n"
		//  << "Buffer radius factor                    : "; data.width(6); data << std::right << Par::r_buff << "\n"
		//  << "Maximum resolution level                : "; data.width(6); data << std::right << Par::max_level << "\n"
		//  << "Core size                               : "; data.width(6); data << std::right << Par::sigma << " m\n"
		//  << "Time step                               : "; data.width(6); data << std::right << Par::dt << " s\n"
         ;
	data << std::fixed << std::setprecision(2)
		//  << "Total simulation time                   : "; data.width(6); data << std::right << Par::simulation_time << " s\n"
		 << "+-------------------------------------------------+\n\n"
         ;

	data << "+----------- Simulation Parameter Data -----------+\n"
	    //  << "Domain x length                       : "; data.width(8); data << std::right << Par::lxdom << " m\n"
		//  << "Domain y length                       : "; data.width(8); data << std::right << Par::lydom << " m\n"
		//  << "Origin gap length                     : "; data.width(8); data << std::right << Par::xdom  << " m\n"
		//  << "Reference length (D)                  : "; data.width(8); data << std::right << Par::Df    << " m\n"
		//  << "Plate_thick                           : "; data.width(8); data << std::right << Par::Df*Par::H_star   << " m\n"
         ;
	data << std::fixed << std::setprecision(4)
		//  << "Courant number (C)                    : "; data.width(8); data << std::right << Par::Courant << "\n"
		//  << "Diffusion number (Phi)                : "; data.width(8); data << std::right << Par::Diffusion << "\n"
		//  << "Turbulent scaling                     : "; data.width(8); data << std::right << std::ceil(1.0e0 / Par::Courant) << "\n"
		//  << "Maximum time step (Phi criteria)      : "; data.width(8); data << std::right << 0.25*Par::sigma*Par::sigma/Par::NU << " s\n"
		 << "+-------------------------------------------------+\n\n"
         ;
	
	data << std::setprecision(-1);
	
	// End of writting simulation setting
	data.close();
}

// Method to write the internal data properties
void dataSaving::write_internal_data(const internalElement& intElm){
    // -- accessing struct data, variables that are not commented are an input only
	int &np = p.num;
	double limdom = 1;
	
	// Store neighbor ID flag
	bool _ngh = false;

	// internal variables
	std::string name;
	std::ofstream ofs;

	printf("\nSaving the particle state ...\n");
		
	// Save common data (information of free particles):
	name.append("output/particle_state_");
	name.append(s);
	name.append(".csv");
	ofs.open(name.c_str());
	ofs << "" << "xp" 
		<< "," << "yp" 
		<< "," << "gpz"
		<< "," << "vor" 
		<< "," << "up" 
		<< "," << "vp" 
		<< "," << "sp" 
		<< "," << "active" 
		<< "," << "chi" 
		// << "," << "R" 
		// << "," << "Basis_CELL_ID" 
		// << "," << "CELL_ID"
		;
	if (_ngh){
		ofs << "," << "ngh_ID";
	}
	ofs << "\n";

	for (int i = 0; i < np; i++)
	{
		if (type == 1)	// Only save the interest domain
		{
			if (p.x[i] > -limdom && p.x[i] < Par::Df + limdom && p.y[i] > -limdom && p.y[i] < limdom)
			{
				ofs << "" << p.x[i]
					<< "," << p.y[i]
					<< "," << p.gz[i]
					<< "," << p.vorticity[i]
					<< "," << p.u[i] - Pars::ubody 
					<< "," << p.v[i] - Pars::vbody
					<< "," << p.s[i]
					<< "," << p.isActive[i]
					<< "," << p.chi[i]
					// << "," << p.R[i]
					// << "," << p.basis_label[i]
		            // << "," << p.cell_label[i]
					;
				if (_ngh){
					ofs << "," ;
					for (int j = 0; j < p.neighbor[i].size(); j++)
					{
						if(j == 0)
							ofs << p.neighbor[i][j];
						else
							ofs << "." << p.neighbor[i][j];
					}
				}
				ofs << "\n";
			}
		}
		else if (type == 2)	// Save whole domain
		{
			ofs << "" << p.x[i]
				<< "," << p.y[i]
				<< "," << p.gz[i]
				<< "," << p.vorticity[i]
				<< "," << p.u[i] - Pars::ubody 
				<< "," << p.v[i] - Pars::vbody
				<< "," << p.s[i]
				<< "," << p.isActive[i]
				<< "," << p.chi[i]
				// << "," << p.R[i]
				// << "," << p.basis_label[i]
		    	// << "," << p.cell_label[i]
				;
			if (_ngh){
				ofs << "," ;
				for (int j = 0; j < p.neighbor[i].size(); j++)
				{
					if(j == 0)
						ofs << p.neighbor[i][j];
					else
						ofs << "." << p.neighbor[i][j];
				}
			}
			ofs << "\n";
		}
	}
	ofs.close();
}

// Method to write the boundary element data
void dataSaving::write_BEM_data(const element& elm){
    // A command
}