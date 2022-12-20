#ifndef PHYSICAL_PROPERTIES_CALCULATION
#define PHYSICAL_PROPERTIES_CALCULATION

#ifndef INCLUDE_VARIABLE
#include "../../variable.hpp"
#endif

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

#ifndef INCLUDED_NEIGHBOR
#include "../Neighbor/neighbor.hpp"
#endif

class propertyCalc
{
private:
    // The private method of property calculation
    // ******************************************
    // The internal object variable
    neighbor neighEval; // Method to calculate neighbor
    LSMPSa lsmpsa_phi;  // Method to calculate LSMPS A
    
    // The property calculation
    void calculate_stress(intElement& intNode);
    void calculate_strain(intElement& intNode);
    void calculate_disp(intElement& intNode);
    
    // The property calculation by transformation to cylindrical
    void coor_transform(intElement& intNode);
public:
    // The public method of property calculation
    // *****************************************
    // The property calculation manager
    void calculate_property(intElement& intNode);

    // For TESTING !!!
    // The analytical soltuion of phi biaxial
    void phi_analytic_biaxial(intElement& intNode);

    // The analytical solution of polynomial laplace and biharmonic
    void lap_analytic(intElement& intNode);
    void bhm_analytic(intElement& intNode);
};

#endif