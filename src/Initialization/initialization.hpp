#ifndef INITIALIZATION
#define INITIALIZATION

#ifndef INCLUDE_VARIABLE
#include "../../variable.hpp"
#endif

class initialization
{
private:
    // The class of private method belong here
    // Internal node generation
    void internal_regular(intElement& intElm, const element& elm, const std::vector<element>& in_elm);
    void internal_finer_near_panel(intElement& intElm, const element& elm, const std::vector<element>& in_elm);

    // Boundary element generation
    void element_rectangular(element& elm, int ID);
    void element_circular(element& elm, int ID);

    // Panel boundary value calculation
    void F_val_calc(element& elm, int ID);
    void p_val_calc(element& elm, int ID);
    void T_val_calc(element& elm, int ID);
    
    // Given function boundary condition
    void F_bihar_val_calc(element& elm, int ID);
    void p_bihar_val_calc(element& elm, int ID);
    void p_lap_val_calc(element& elm, int ID);

public:
    // The class of public method belong here
    void generate_internal_node(intElement& intElm, const element& elm, const std::vector<element>& in_elm);
    void generate_boundary_element(element& elm, std::vector<element>& in_elm);
    void calculate_boundary_condition(element& elm, std::vector<element>& in_elm);
};

#endif
