# Copyright (C) 2005, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
# If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
# https://www.eclipse.org/legal/epl-2.0/

# Authors:  N.N
#           Daniel Doehring                  RWTH Aachen University 2023-05-12

DCO_PATH=$(HOME)/Software/dco

# This should be the name of your executable
EXE = Roots_Real Roots_RealImag

# Additional libraries. Here: dco
ADDLIBS = -L $(DCO_PATH)/lib -l dcoc

# Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I $(DCO_PATH)/include -I include/ \
							-DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING 

##########################################################################

# C++ Compiler command
CXX = g++ -std=c++17

# C++ Compiler options

CXXFLAGSRUN = -Ofast

CXXFLAGSDEBUG = -O0 -g -Wall -Wextra -Wpedantic
# Change only this
CXXFLAGS = ${CXXFLAGSRUN}

# additional C++ Compiler options for linking
#CXXLINKFLAGS =  -Wl,--rpath -Wl,/usr/local/lib
CXXLINKFLAGS =

# Include directories
INCL = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --libs ipopt`

SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

all: Roots_Real Roots_RealImag

.SUFFIXES: .cpp .o

Roots_Real: $(OBJ_DIR)/Roots_Real.o $(OBJ_DIR)/Main_Roots_Real.o | $(BIN_DIR)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@.exe $(OBJ_DIR)/Roots_Real.o $(OBJ_DIR)/Main_Roots_Real.o $(ADDLIBS) $(LIBS)

Roots_RealImag: $(OBJ_DIR)/Roots_RealImag.o $(OBJ_DIR)/Main_Roots_RealImag.o | $(BIN_DIR)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@.exe $(OBJ_DIR)/Roots_RealImag.o $(OBJ_DIR)/Main_Roots_RealImag.o $(ADDLIBS) $(LIBS)

# Object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCL) -c $< -o $@ 

$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

.PHONY: clean cleanout

clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR) *.out 

cleanout:
	rm *.txt
