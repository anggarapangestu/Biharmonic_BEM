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
    neighbor neighEval;
    LSMPSa lsmpsa_phi;
    void calculate_stress(intElement& intNode);
    void calculate_strain(intElement& intNode);
    void calculate_disp(intElement& intNode);
public:
    void calculate_property(intElement& intNode);
};

#endif