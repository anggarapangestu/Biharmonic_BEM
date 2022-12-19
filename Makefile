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
PROGRAM_LOG	= log

# ==========================================
# ============ PROGRAM COMPILE =============
# ==========================================
# List of cpp source file path
SRCS		= main.cpp\
			setting.cpp\
			src/Initialization/initialization.cpp\
			src/Initialization/init_element.cpp\
			src/Initialization/init_internal.cpp\
			src/Initialization/init_assignBC.cpp\
			src/Saving/save_data.cpp\
			src/BEM/BEM.cpp\
			src/BEM/BEM_utils.cpp\
			src/LSMPS/LSMPSa.cpp\
			src/Neighbor/link_list.cpp\
			src/Neighbor/link_list_utils.cpp\
			src/PropertyCalc/physical_prop.cpp

# List of target object file path
OBJS		= $(SRCS:%.cpp=%.o)           # Change from SRCS the ".cpp" into ".o"

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

# Target to rebuild the program
rebuild:	clean $(PROGRAM)

# Target to run the program
run:
			@./$(PROGRAM)

# ==========================================
# ============= SIMULATION LOG =============
# ==========================================
# List of cpp source file path for log
LOG			= log.cpp\
			setting.cpp\
			src/Saving/save_data.cpp

# Object file path for log
OBJS_LOG	= $(LOG:%.cpp=%.o)            # Change from SRCS the ".cpp" into ".o"

# Target to compile and link all log
$(PROGRAM_LOG):	OBPR_LOG $(OBJS_LOG)
			@echo "\n[DONE] Object file is completely built"
			@echo "\nStart linking the object ..."
			$(COMPILER) $(OBJS_LOG) -o $(PROGRAM_LOG)
			@echo "\n[DONE] Program compiled and linked successfully"

# Target for displaying initial object building prompt
OBPR_LOG:
			@echo "Start building the log object file ..."

# Target to run the log message
run_log: $(PROGRAM_LOG)
			@./$(PROGRAM_LOG)

# ==========================================
# ================== DATA ==================
# ==========================================
# Target to delete the log object program files
clean_log:
			@rm -f $(OBJS_LOG)
			@rm -f $(PROGRAM_LOG)
			@echo "+-------------- LOG CLEANED -------------+"

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


# github link 	:= https://github.com/anggarapangestu/Biharmonic_BEM.git
# token 		:= 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