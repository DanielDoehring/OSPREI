// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Carl Laird, Andreas Waechter     IBM                    2005-08-16
//           Daniel Doehring                  RWTH Aachen University 2022-11-20


#include "IpIpoptApplication.hpp"
#include "Roots_Real.hpp"

#include <iostream>

using namespace Ipopt;

int main(int argc, char** argv) {
   assert(argc >= 6);

   const int NumStages    = std::stoi(argv[1]);
   const int ConsOrder    = std::stoi(argv[2]);
   const int NumStagesRef = std::stoi(argv[3]);
   const Number dtRef     = std::stod(argv[4]);

   // Currently only order 1 and 2 implemented
   assert(ConsOrder == 1 || ConsOrder == 2);

   // Create a new instance of your nlp (use Ipopt::SmartPtr)
   SmartPtr<TNLP> mynlp;
   if(argc == 7)
      // Case for which hull is used
      mynlp = new Roots_Real(NumStages, ConsOrder, NumStagesRef, dtRef, std::string(argv[5]), std::string(argv[6]));
   else
      mynlp = new Roots_Real(NumStages, ConsOrder, NumStagesRef, dtRef, std::string(argv[5]));

   // Create a new instance of IpoptApplication (use Ipopt::SmartPtr)
   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

   app->Options()->SetStringValue("output_file", "Roots_Real.out");

   // The following overwrites the default name (ipopt.opt) of the options file
   app->Options()->SetStringValue("option_file_name", "Roots_Real.opt");

   if(ConsOrder == 1)
      // There are no equality constraints => constant Eq.-Constr. Jacobian (not sure of option does something)
      app->Options()->SetStringValue("jac_c_constant", "yes");

   // Initialize the IpoptApplication and process the options
   ApplicationReturnStatus status;
   status = app->Initialize();

   if( status != Solve_Succeeded ) {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return (int) status;
   }

   // Ask Ipopt to solve the problem
   status = app->OptimizeTNLP(mynlp);

   std::cout << std::endl << std::endl << "*** The problem terminated!" << std::endl;

   return (int) status;
}
