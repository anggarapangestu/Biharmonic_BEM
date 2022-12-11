#include "initialization.hpp"

// Generation of internal node
void initialization::generate_internal_node(intElement& intElm, const element& elm, const std::vector<element>& in_elm){
    // initialization generate internal node starting log
    printf("\nGenerate the internal node ...\n");
    clock_t _time = clock();

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
    printf("<+> Number of element                 : %8d\n", intElm.num);
    
    // Displaying the computational time
    _time = clock() - _time;
	printf("<-> Internal node generation\n");
    printf("    comp. time                         [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}

// ===================================================
// Generation of boundary panel element
void initialization::generate_boundary_element(element& elm, std::vector<element>& in_elm){
    // Procedure: \
       1. Generate the base geometry panel element\
       2. Generate the inner domain geometry panel element

    // initialization generate boundary panel starting log
    printf("\nGenerate the base geometry boundary panel ...\n");
    clock_t _time = clock();
    
    // Generate the base geometry panel element
    if (Par::G_type == 1){
        printf("<+> Rectangular geometry type\n");
        this->element_rectangular(elm, -1);
    }else if (Par::G_type == 2){
        printf("<+> Circular geometry type\n");
        this->element_circular(elm, -1);
    }
    printf("<+> Number of element                 : %8d\n", elm.num);


    // Generate the base geometry panel element
    for (int ID = 0; ID < Par::N_Gin; ID++){
        // initialization generate boundary panel starting log
        printf("\nGenerate boundary panel of inner geometry %d ...\n", ID+1);

        // Create the panel for each inner geometry
        element innerPanel;
        
        // Create the panel
        if (Par::Gin_type[ID] == 1){
            printf("<+> Rectangular geometry type\n");
            this->element_rectangular(innerPanel, ID);
        }else if (Par::Gin_type[ID] == 2){
            printf("<+> Circular geometry type\n");
            this->element_circular(innerPanel, ID);
        }
        printf("<+> Number of element                 : %8d\n", innerPanel.num);

        // Insert the panel into the list
        in_elm.emplace_back(innerPanel);
    }

    // Displaying the computational time
    _time = clock() - _time;
	printf("\n<-> All boundary element generation \n");
    printf("    comp. time                         [%8.4f s]\n", (double)_time/CLOCKS_PER_SEC);
}

// ===================================================
// Calculate the boundary condition at each panel
void initialization::calculate_boundary_condition(element& elm, std::vector<element>& in_elm){
    // Calculate the boundary condition for \
       1. Base geometry element and \
       2. Inner domain geometry element

    // ======== BASE GEOMETRY ========
    // initialization generate boundary panel starting log
    printf("\nCalculating boundary condition for base panel ...\n");
    
    // Resize the boundary value array
    elm.F.resize(elm.num, 0.0e0);
    elm.dFdn.resize(elm.num, 0.0e0);
    elm.p.resize(elm.num, 0.0e0);
    elm.dpdn.resize(elm.num, 0.0e0);
    elm.F_type.resize(elm.num, true);   // The basic known value (dFdn)
    elm.p_type.resize(elm.num, true);   // The basic known value (dpdn)

    // Generate the base geometry panel element
    printf("<+> Calculate dFdn value\n");
    this->F_val_calc(elm, -1);
    printf("<+> Calculate dpdn value\n");
    this->p_val_calc(elm, -1);

    // ======== INNER GEOMETRY ========
    // initialization generate boundary panel starting log
    printf("\nCalculating boundary condition for inner panel ...\n");

    // Generate the base geometry panel element
    for (int ID = 0; ID < Par::N_Gin; ID++){
        // Resize the boundary value array
        in_elm[ID].F.resize(in_elm[ID].num, 0.0e0);
        in_elm[ID].dFdn.resize(in_elm[ID].num, 0.0e0);
        in_elm[ID].p.resize(in_elm[ID].num, 0.0e0);
        in_elm[ID].dpdn.resize(in_elm[ID].num, 0.0e0);
        in_elm[ID].F_type.resize(in_elm[ID].num, true);   // The basic known value (dFdn)
        in_elm[ID].p_type.resize(in_elm[ID].num, true);   // The basic known value (dpdn)

        // Generate the base geometry panel element
        printf("<+> Calculate dFdn for inner geometry %d \n", ID+1);
        this->F_val_calc(in_elm[ID], ID);
        printf("<+> Calculate dpdn for inner geometry %d \n", ID+1);
        this->p_val_calc(in_elm[ID], ID);
    }
}