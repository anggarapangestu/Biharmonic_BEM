#	File name		: Makefile
#	Date			: December 2022
#	Version			: 1.0.0
#	Author			: Angga

# Useful Note:
# * Shorten the directory in terminal WSL: PROMPT_DIRTRIM=1
# * There are 4 method: (1)Compile, (2)Run, (3)Clean Files, (4)Delete Output Data

# Procedure:
# * Create an object file: g++ -c *.cpp -o *.o
# * Link the object file: g++ *1.o *2.o -o program
# * Target will run the prerequisite word by word

# Defining variables
DEST		= output/
COMPILER	= g++
PROGRAM		= program

# List of cpp source file path
SRCS		= main.cpp\
			setting.cpp\
			src/Initialization/initialization.cpp \
			src/Initialization/init_element.cpp \
			src/Initialization/init_internal.cpp \
			src/Initialization/init_assignBC.cpp \
			src/Saving/save_data.cpp \
			src/BEM/BEM.cpp \
			src/BEM/BEM_utils.cpp \
			# src/LSMPS/LSMPSa.cpp \
			# src/LSMPS/LSMPSb.cpp \
			# src/geometry/generateBody.cpp	\
			# src/geometry/2D_objects_generator.cpp\
			# src/geometry/3D_objects_generator.cpp\
			# src/geometry/u_var.cpp \
			# src/DC_operator/dc_gradient.cpp \

# List of target object file path
OBJS		= $(SRCS:%.cpp=%.o)            # Change from SRCS the ".cpp" into ".o"

# Target to create the object file from cpp source file
%.o: %.cpp
			@echo "\nBuild $* object ..."
			$(COMPILER) -c $*.cpp -o $*.o
			@echo "... $* object is built"

# Target to compile and link all program
$(PROGRAM):	OBPR $(OBJS)
			@echo "\n[DONE] Object file is completely built"
			@echo "\nStart linking the object ..."
			$(COMPILER) $(OBJS) -o $(PROGRAM)
			@echo "\n[DONE] Program compiled and linked successfully"

# Target to compile and link all program
compile:	$(PROGRAM)

# Target for displaying initial object building prompt
OBPR:
			@echo "Start building the object file ..."

# Target to run the program
run:
			@./$(PROGRAM)

# Target to delete the object and program files
clean:
			@while [ -z "$$CONTINUE" ]; do \
			read -r -p "<?> Do you really want to clean program files [Y/N] ? " CONTINUE; \
			done ; \
			if [ $$CONTINUE != "y" ] && [ $$CONTINUE != "Y" ]; then \
			echo "Exiting." ; exit 1; \
			fi
			@rm -f $(OBJS)
			@rm -f $(PROGRAM)
			@echo "+------------- FILE CLEANED -------------+"

# Target to delete the output data file
delete:
			@while [ -z "$$CONTINUE" ]; do \
			read -r -p "<?> Do you really want to delete data files [Y/N] ? " CONTINUE; \
			done ; \
			if [ $$CONTINUE != "y" ] && [ $$CONTINUE != "Y" ]; then \
			echo "Exiting." ; exit 1; \
			fi
			@rm -f output/*.dat output/*.jpg output/*.png output/*.csv
			@echo "+------------- DATA DELETED -------------+"


# github link 	:= https://github.com/anggarapangestu/NEW_VPM.git
# token 		:= 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