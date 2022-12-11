#include "initialization.hpp"

// Generation of internal node
void initialization::generate_internal_node(intElement& intElm, const element& elm, const std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nGenerate the internal node ...\n");

    // Procedure: \
       1. Generate the regular node distribution \
       2. Eleminate the node outside the domain\
       3. Save the position data into intElm\
       4. Neighbor searching

    // Generate the internal particle data
    if (Par::opt_int_init == 1){
        printf("<+> Regular internal node\n");
        this->internal_regular(intElm, elm, in_elm);
    }else if (Par::opt_int_init == 2){
        printf("<+> Finer near panel internal node\n");
        this->internal_finer_near_panel(intElm, elm, in_elm);
    }
    
    // Calculate the neighbor
        // CODE for NEIGHBOR EVAL ...
}

// ===================================================
// Generation of boundary panel element
void initialization::generate_boundary_element(element& elm, std::vector<element>& in_elm){
    // Procedure: \
       1. Generate the base geometry panel element\
       2. Generate the inner domain geometry panel element

    // initialization generate boundary panel starting log
    printf("\nGenerate the base geometry boundary panel ...\n");
    
    // Generate the base geometry panel element
    if (Par::G_type == 1){
        printf("<+> Rectangular geometry type\n");
        this->element_rectangular(elm, -1);
    }else if (Par::G_type == 2){
        printf("<+> Circular geometry type\n");
        this->element_circular(elm, -1);
    }

    // initialization generate boundary panel starting log
    printf("\nGenerate the inner geometry boundary panel ...\n");

    // Generate the base geometry panel element
    for (int ID = 0; ID < Par::N_Gin; ID++){
        // Create the panel for each inner geometry
        element innerPanel;
        
        // Create the panel
        if (Par::Gin_type[ID] == 1){
            printf("<+> Rectangular type for inner geometry %d \n", ID+1);
            this->element_rectangular(innerPanel, ID);
        }else if (Par::Gin_type[ID] == 2){
            printf("<+> Circular type for inner geometry %d \n", ID+1);
            this->element_circular(innerPanel, ID);
        }

        // Insert the panel into the list
        in_elm.emplace_back(innerPanel);
    }
}

// ===================================================
// Calculate the boundary condition at each panel
void initialization::calculate_boundary_condition(element& elm, std::vector<element>& in_elm){
    // Calculate the boundary condition for \
       1. Base geometry element and \
       2. Inner domain geometry element

    // initialization generate boundary panel starting log
    printf("\nCalculating boundary condition for base panel ...\n");
    
    // Generate the base geometry panel element
    printf("<+> Calculate dFdn value\n");
    this->element_rectangular(elm, -1);
    printf("<+> Calculate dpdn value\n");
    this->element_rectangular(elm, -1);

    // initialization generate boundary panel starting log
    printf("\nCalculating boundary condition for inner panel ...\n");

    // Generate the base geometry panel element
    for (int ID = 0; ID < Par::N_Gin; ID++){
        // Generate the base geometry panel element
        printf("<+> Calculate dFdn for inner geometry %d \n", ID+1);
        this->element_rectangular(elm, -1);
        printf("<+> Calculate dpdn for inner geometry %d \n", ID+1);
        this->element_rectangular(elm, -1);
    }
}