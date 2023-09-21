// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Carl Laird, Andreas Waechter     IBM                    2005-08-16
//           Daniel Doehring                  RWTH Aachen University 2022-11-20

#include "IpIpoptCalculatedQuantities.hpp" // To access current original, unscaled violations

#include "Roots_Real.hpp"

#include <cassert>

#include <iostream>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "IO_Funcs.hpp"
#include "RootDistribution.hpp"
#include "StabConstraints_Real.hpp"
#include "OrderConstraints_Real.hpp"

using namespace Ipopt;


Roots_Real::Roots_Real(
   const int NumStages_,
   const int ConsOrder_,
   const int NumStagesRef_, 
   const Number dtRef_,
   const std::string EigValFileName
) : NumStages(NumStages_), Degree(NumStages_), ConsOrder(ConsOrder_),
    NumStagesRef(NumStagesRef_), dtRef(dtRef_), dtExp((dtRef / NumStagesRef_) *  NumStages_)
{
  std::cout << std::endl << "Optimize a " << Degree << " degree stability polynomial" << std::endl
            << std::endl << "Optimize roots of the " << Degree - 1 << " \"lower\" degree polynomial" 
            << std::endl << std::endl;

  std::cout.precision(std::numeric_limits<Number>::max_digits10);
  std::cout << "The expected maximal stable timestep is: " << dtExp  << std::endl << std::endl;

  NumEigVals = -1;
  read_EigVals(EigValFileName, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
  // TODO: This assumes that the eigenvalues are sorted ascending with real part -> Implement sorting!

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEigValsScaled[i] *= dtExp;
    ImagEigValsScaled[i] *= dtExp;
  }

  ImagDiff_over_RealDiff.resize(NumEigVals - 1);
  for(size_t i = 0; i < NumEigVals-1; i++) {
      ImagDiff_over_RealDiff[i] = (ImagEigValsScaled[i+1] - ImagEigValsScaled[i]) / 
                                  (RealEigValsScaled[i+1] - RealEigValsScaled[i]);
  }

  RealUB  = std::min(RealEigValsScaled[NumEigVals-1], -1e-9); // Division by zero guard
  RealMin = *min_element(std::begin(RealEigValsScaled), std::end(RealEigValsScaled));
  // CARE: Assume that smallest real is @0, next @1
  RealMargin = std::abs(RealMin - std::real(RealEigValsScaled[1]));

  NumConstr = NumEigVals + ConsOrder - 1; // -1 Since one consistency order comes for free

  OddDegree = Degree % 2;
  NumRoots  = NumStages / 2; // Note: Integer division is here desired

  // Optimize timestep now as well
  NumUnknowns = NumRoots + 1;

   std::string PE_HalfStagesFileName = "./RealImag_Optimized_" + std::to_string(NumStages/2) + ".txt";
   if(std::filesystem::exists(PE_HalfStagesFileName)) {
      std::cout << "Use pseudo-extrema of " << std::to_string(NumStages/2) 
                << " stage RKM for initialization, stored in file" 
                << PE_HalfStagesFileName << std::endl << std::endl;

      std::vector<Number> Real_PE_HalfStagesScaled;
      read_PE(PE_HalfStagesFileName, NumStages/4, Real_PE_HalfStagesScaled);
      // Scale Pseudo-Extrema
      for(size_t i = 0; i < NumStages/4; i++)
         Real_PE_HalfStagesScaled[i] *= 2.0;

      x0 = InitialRootDistr(NumEigVals, NumUnknowns, NumStages, RealEigValsScaled, ImagEigValsScaled, 
                            NumStages/4, Real_PE_HalfStagesScaled);
   }
   else 
     x0 = InitialRootDistr(NumEigVals, NumUnknowns, NumStages, RealMin, RealEigValsScaled, ImagEigValsScaled);

  // Expect some slightly less optimal timestep
  x0[NumRoots] = 0.95 * dtExp; 
  
  std::cout << "Initial values are: " << std::endl;
  for(size_t i = 0; i < NumRoots; i++)
    std::cout << x0[i] << std::endl;
  std::cout << x0[NumRoots] << std::endl << std::endl;

  UseHull = false;
  std::cout << "Use hull for interpolation? " << UseHull << std::endl << std::endl;

  xMaxdt = new Number[NumUnknowns];
  Maxdt = 0.;
  InfPr = 42e6;
 }

 Roots_Real::Roots_Real(
   const int NumStages_,
   const int ConsOrder_,
   const int NumStagesRef_,
   const Number dtRef_,
   const std::string EigValFileName,
   const std::string HullPointPath
) : NumStages(NumStages_), Degree(NumStages_), ConsOrder(ConsOrder_),
    NumStagesRef(NumStagesRef_), dtRef(dtRef_), dtExp((dtRef / NumStagesRef_) *  NumStages_)
{
    std::cout << std::endl << "Optimize a " << Degree << " degree stability polynomial" << std::endl
              << std::endl << "Optimize roots of the " << Degree - 1 << " \"lower\" degree polynomial" 
              << std::endl << std::endl;

  std::cout.precision(std::numeric_limits<Number>::max_digits10);
  std::cout << "The expected maximal stable timestep is: " << dtExp  << std::endl << std::endl;

  NumEigVals = -1;
  read_EigVals(EigValFileName, NumEigVals, RealEigValsScaled, ImagEigValsScaled);

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEigValsScaled[i] *= dtExp;
    ImagEigValsScaled[i] *= dtExp;
  }

  RealUB  = std::min(RealEigValsScaled[NumEigVals-1], -1e-9); // Division by zero guard
  RealMin = *min_element(std::begin(RealEigValsScaled), std::end(RealEigValsScaled));
  // CARE: Assume that smallest real is @0, next @1
  RealMargin = std::abs(RealMin - std::real(RealEigValsScaled[1]));

  NumConstr = NumEigVals + ConsOrder - 1; // -1 Since one consistency order comes for free

  OddDegree = Degree % 2;
  NumRoots  = NumStages / 2; // Note: Integer division is here desired

  // Optimize timestep now as well
  NumUnknowns = NumRoots + 1;

  read_Hull(HullPointPath + "_real.txt", HullRealScaled);
  read_Hull(HullPointPath + "_imag.txt", HullImagScaled);
  assert(HullRealScaled.size() == HullImagScaled.size());

  for(size_t i = 0; i < HullRealScaled.size(); i++) {
      HullRealScaled[i] *= dtExp;
      HullImagScaled[i] *= dtExp;
  }

  ImagDiff_over_RealDiff.resize(HullRealScaled.size() - 1);
  for(size_t i = 0; i < HullRealScaled.size()-1; i++) {
      ImagDiff_over_RealDiff[i] = (HullImagScaled[i+1] - HullImagScaled[i]) / 
                                  (HullRealScaled[i+1] - HullRealScaled[i]);
  }

  UseHull = true;
  std::cout << "Use hull for interpolation? " << UseHull << std::endl << std::endl;

   std::string PE_HalfStagesFileName = "./RealImag_Optimized_" + std::to_string(NumStages/2) + ".txt";
   if(std::filesystem::exists(PE_HalfStagesFileName)) {
      std::cout << "Use pseudo-extrema of " << std::to_string(NumStages/2) 
               << " stage RKM for initialization, stored in file" 
               << PE_HalfStagesFileName << std::endl << std::endl;

      std::vector<Number> Real_PE_HalfStagesScaled;
      read_PE(PE_HalfStagesFileName, NumStages/4, Real_PE_HalfStagesScaled);
      // Scale Pseudo-Extrema
      for(size_t i = 0; i < NumStages/4; i++)
         Real_PE_HalfStagesScaled[i] *= 2.0;

      x0 = InitialRootDistr(HullImagScaled.size(), NumUnknowns, NumStages, HullRealScaled, HullImagScaled, 
                              NumStages/4, Real_PE_HalfStagesScaled);
   }
   else
      x0 = InitialRootDistr(HullImagScaled.size(), NumUnknowns, NumStages, RealMin, HullRealScaled, HullImagScaled);
   

   // Expect some slightly less optimal timestep
   x0[NumRoots] = 0.95 * dtExp; 
  
  std::cout << "Initial values are: " << std::endl;
  for(size_t i = 0; i < NumRoots; i++)
    std::cout << x0[i] << std::endl;
  std::cout << x0[NumRoots] << std::endl << std::endl;

  xMaxdt = new Number[NumUnknowns];
  Maxdt = 0.;
  InfPr = 42e6;
}

// destructor
Roots_Real::~Roots_Real()
{
   delete[] xMaxdt;
}

// [TNLP_get_nlp_info]
// returns the size of the problem
bool Roots_Real::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g, // NNZ: Number NonZero
   Index&          nnz_h_lag, // LAGrangian
   IndexStyleEnum& index_style
)
{
   // n = Number variables
   n = NumUnknowns;

   // m = Number constraints
   m = NumConstr;
   std::cout << "Total number of constraints is: " << m << std::endl;

   // Jacobian will in general be dense, dimensions: m * n
   nnz_jac_g = m * NumUnknowns;

   // The Hessian is also dense and has NumRoots*NumRoots total nonzeros, but we
   // only need the lower left corner (since it is symmetric)
   nnz_h_lag = (NumUnknowns * (NumUnknowns + 1))/2;

   // use the C style indexing (0-based)
   index_style = TNLP::C_STYLE;

   return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool Roots_Real::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{
   // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
   // If desired, we could assert to make sure they are what we think they are.
   // Omitted for performance
   /*
   assert(n == NumUnknowns);
   assert(m == NumEigVals+ConsOrder-1);
   */

   for( Index i = 0; i < NumRoots; i++ ) {
      // Tried stricter bounds (around initial values, for instance). In practice not (yet) necessary.
      /*
      x_l[i] = std::min(x0[i] - RealPartMargin, RealMin);
      x_u[i] = std::max(x0[i] + RealPartMargin, RealEigValsScaled[NumEigVals-1]);
      */
      //x_l[i] = RealMin - RealMargin;
      x_l[i] = RealMin;

      // 'RealUB': Division-by-zero guard
      x_u[i] = RealUB;
   }
   // Timestep
   x_l[NumRoots] = 1e-6; // Safety guard (no division by zero)
   //x_l[NumRoots] = 0.8 * dtExp; // Safety guard (no division by zero)

   x_u[NumRoots] = 1.1 * dtExp; // Allow for slightly larger timestep

   for( Index i = 0; i < NumEigVals; i++ ) {
     // Not clear what is easier for the optimizer to handle
     //g_l[i] = 0.; // Automatically ensured by abs
     g_l[i] = -2e19; // Two times ipopt default value for treating number as inf

     g_u[i] = 1.; // Actual stability bound
   }

   // Order constraints: Equality constraints
   for ( Index i = NumEigVals; i < NumEigVals + ConsOrder - 1; i++ ) {
      g_l[i] = 0.;
      g_u[i] = 0.;
   }

   return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool Roots_Real::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
)
{
   // Omitted for performance
   /*
   assert(n == NumUnknowns);
   assert(m == NumEigVals+ConsOrder-1);
   */

   // Here, we assume we only have starting values for x, if you code
   // your own NLP, you can provide starting values for the dual variables
   // if you wish
   assert(init_x == true);
   // initialize to the given starting point
   for( Index i = 0; i < NumUnknowns; i++ ) {
      x[i] = x0[i];
   }

   if(init_lambda)
      for( Index i = 0; i < m; i++ ) {
         lambda[i] = 0.;
      }

   return true;
}
// [TNLP_get_starting_point]

// [TNLP_eval_f]
// returns the value of the objective function
bool Roots_Real::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
   //assert(n == NumUnknowns); // Removed for performance

   obj_value = -x[NumRoots]; // This is the (negated timestep)

   return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool Roots_Real::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
   //assert(n == NumUnknowns); // Removed for performance

   for(size_t i = 0; i < NumRoots; i++)
     grad_f[i] = 0.;

   grad_f[NumRoots] = -1.;

   return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool Roots_Real::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
   // Removed for performance
   /*
   assert(n == NumUnknowns);
   assert(m == NumEigVals+ConsOrder-1);
   */

   if(OddDegree) {
      if(UseHull)
         StabConstr_Real(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                         HullRealScaled, HullImagScaled, dtExp);
      else
         StabConstr_Real(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, dtExp);

      if(ConsOrder >= 2) {
         if(UseHull)
            SecOrder(x, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
         else
            SecOrder(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);

         if(ConsOrder >= 3) {
            if(UseHull)
               ThirdOrder(x, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
            else
               ThirdOrder(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
         }
      }
   }
   else {
      i_min = 0;
      for(size_t i = 1; i < NumRoots; i++) {
         if(x[i] < x[i_min])
            i_min = i;
      }

      if(UseHull)
         StabConstr_Real(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                         HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
      else
         StabConstr_Real(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff, 
                         dtExp, i_min);

      if(ConsOrder >= 2) {
         if(UseHull)
            SecOrder(x, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, i_min);
         else
            SecOrder(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff, i_min);

         if(ConsOrder >= 3) {
            if(UseHull)
               ThirdOrder(x, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, i_min);
            else
               ThirdOrder(x, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff, i_min);
         }
      }
   }

   return true;
}
// [TNLP_eval_g]

// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool Roots_Real::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
   // Removed for performance
   /*
   assert(n == NumUnknowns);
   assert(m == NumEigVals+ConsOrder-1);
   */

   if( values == NULL ) {
      // return the structure of the Jacobian
      // Dimensions: m x n

      // Jacobian is dense
      size_t ind = 0;
      for(size_t i = 0; i < m; i++)
        for(size_t j = 0; j < n; j++) {
          iRow[ind] = i;
          jCol[ind] = j;

          ind++;
        }
   }
   else {
      using DCO_M  = typename dco::ga1s<Number>;  // Turn on adjoint mode
      using DCO_T  = typename DCO_M::type;   // Declare adjoint type
      using DCO_TT = typename DCO_M::tape_t; // Specify tape type

      // Removed for performance
      //assert(NumEigVals == RealEigValsScaled.size());

      std::vector<DCO_T> x_dco(NumUnknowns);
      for(size_t i = 0; i < NumUnknowns; i++)
         dco::value(x_dco[i]) = x[i];

      DCO_M::global_tape = DCO_TT::create(); // Prepare / "touch" tape
      DCO_M::global_tape->register_variable(x_dco);

      std::vector<DCO_T> g(m);
      
      if(OddDegree) {
         if(UseHull)
            StabConstr_Real(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                            HullRealScaled, HullImagScaled, dtExp);
         else
            StabConstr_Real(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                            dtExp);
      }
      else {
         i_min = 0;
         for(size_t i = 1; i < NumRoots; i++) {
            if(x[i] < x[i_min])
               i_min = i;
         }

         if(UseHull)
            StabConstr_Real(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                            HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
         else
            StabConstr_Real(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                            ImagDiff_over_RealDiff, dtExp, i_min);
      }

      DCO_M::global_tape->register_output_variable(g); // Record active output

      size_t ind = 0;
      for(size_t i = 0; i < NumEigVals; i++) {
         dco::derivative(g)[i] = 1.; // Seed component

         DCO_M::global_tape->interpret_adjoint(); // Interpret (stored) tape

         // Harvest
         for (size_t j = 0; j < NumUnknowns; j++) {
            values[ind] = dco::derivative(x_dco[j]);
            ind++;
         }

         // Reset adjoints
         //dco::derivative(g)[i] = 0.; // Unseed component
         DCO_M::global_tape->zero_adjoints();
      }

      if(ConsOrder >= 2) {
         if(OddDegree) {
            if(UseHull)
               SecOrder(x_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
            else
               SecOrder(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
         }
         else {
            if(UseHull)
               SecOrder(x_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                        ImagDiff_over_RealDiff, i_min);
            else
               SecOrder(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                        ImagDiff_over_RealDiff, i_min);
         }
         
         dco::derivative(g)[NumEigVals] = 1.; // Seed component

         DCO_M::global_tape->interpret_adjoint(); // Interpret (stored) tape

         // Harvest
         for (size_t j = 0; j < NumUnknowns; j++) {
            values[ind] = dco::derivative(x_dco[j]);
            ind++;
         }

         if(ConsOrder >= 3) {
            DCO_M::global_tape->zero_adjoints();

            if(OddDegree) {
               if(UseHull)
                  ThirdOrder(x_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
               else
                  ThirdOrder(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
            }
            else {
               if(UseHull)
                  ThirdOrder(x_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                           ImagDiff_over_RealDiff, i_min);
               else
                  ThirdOrder(x_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                           ImagDiff_over_RealDiff, i_min);
            }
            
            dco::derivative(g)[NumEigVals + 1] = 1.; // Seed component

            DCO_M::global_tape->interpret_adjoint(); // Interpret (stored) tape

            // Harvest
            for (size_t j = 0; j < NumUnknowns; j++) {
               values[ind] = dco::derivative(x_dco[j]);
               ind++;
            }
         }
      }

      DCO_TT::remove(DCO_M::global_tape); // Deallocate tape
   }

   return true;
}
// [TNLP_eval_jac_g]


// [TNLP_eval_h]
//return the structure or values of the Hessian
bool Roots_Real::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess, // I guess this is Number ELEments HESSian
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{
   // Removed for performance
   /*
   assert(n == NumUnknowns);
   assert(m == NumEigVals+ConsOrder-1);
   */

   if( values == NULL ) {
      // return the structure. 
      // This is a symmetric matrix, fill the lower left triangle only.

      // the hessian for this problem is actually dense
      Index idx = 0;
      for( Index row = 0; row < n; row++ ) {
         for( Index col = 0; col <= row; col++ ) {
            iRow[idx] = row;
            jCol[idx] = col;
            idx++;
         }
      }

      //assert(idx == nele_hess); // Removed for performance
   }
   else {
      // return the values. 
      // This is a symmetric matrix, fill the lower left triangle only

      // fill the objective portion
      // Not existent for this problem (feasibility problem)

      // Constraints

      /// First Eigenvalue ///
      using DCO_BM  = typename dco::ga1s<Number>; // adjoint (base) mode
      using DCO_BT  = typename DCO_BM::type; // adjoint (base) type
      using DCO_BTT = typename DCO_BM::tape_t; // base tape type
      using DCO_M   = typename dco::ga1s<DCO_BT>; // adjoint of base mode
      using DCO_T   = typename DCO_M::type; // adjoint of base type
      using DCO_TT  = typename DCO_M::tape_t; // tape type

      std::vector<DCO_T> x_dco(NumUnknowns);

      DCO_BM::global_tape = DCO_BTT::create(); // "touch" base tape
      DCO_M::global_tape  = DCO_TT::create(); // "touch" tape

      DCO_T g; // Scalar output

      for (size_t i = 0; i < NumUnknowns; i++) {
         DCO_BM::global_tape->register_variable(dco::value(x_dco[i]) ); // record active input
         DCO_BM::global_tape->register_variable(dco::derivative(x_dco[i]) ); // record active input

         dco::value(dco::value(x_dco[i]) ) = x[i];

         DCO_M::global_tape->register_variable(x_dco[i]);
      }
      
      if(OddDegree) {
         if(UseHull)
            g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, 
                                  HullRealScaled, HullImagScaled, dtExp);
         else
            g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, dtExp);
      }
      else {
         i_min = 0;
         for(size_t i = 1; i < NumRoots; i++) {
            if(x[i] < x[i_min])
               i_min = i;
         }
         if(UseHull)
            g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, 
                                  HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
         else
            g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, 
                                  ImagDiff_over_RealDiff, dtExp, i_min);
      }

      dco::value(dco::derivative(g) ) = 1.; // Seed
      DCO_M::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint

      int ind = 0;
      for(size_t i = 0; i < NumUnknowns; i++) {
         dco::derivative(dco::derivative(x_dco[i]) ) = 1.; // Seed

         DCO_BM::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint
         for(size_t j = 0; j <= i; j++) { // Fill lower left triangle only
            values[ind] = lambda[0] * dco::derivative(dco::value(x_dco[j]) );
            ind++;
         }

         //dco::derivative(dco::derivative(x_dco[i]) ) = 0.; // Unseed

         DCO_BM::global_tape->zero_adjoints();
      }

      // Clear tape only, do not remove
      DCO_M::global_tape->reset();
      DCO_BM::global_tape->reset();

      /// Remaining Eigenvalues ///
      for(size_t i = 1; i < NumEigVals; i++) {
         for (size_t j = 0; j < NumUnknowns; j++) {
            DCO_BM::global_tape->register_variable(dco::value(x_dco[j]) ); // record active input
            DCO_BM::global_tape->register_variable(dco::derivative(x_dco[j]) ); // record active input

            dco::value(dco::value(x_dco[j]) ) = x[j];

            DCO_M::global_tape->register_variable(x_dco[j]);
         }

         if(OddDegree) {
            if(UseHull)
               g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, 
                                     HullRealScaled, HullImagScaled, dtExp);
            else
               g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, dtExp);
         }
         else {
            if(UseHull)
               g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, 
                                     HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
            else
               g = StabConstr_Real_i(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, 
                                    ImagDiff_over_RealDiff, dtExp, i_min);
         }           

         dco::value(dco::derivative(g) ) = 1.; // Seed
         DCO_M::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint

         ind = 0;
         for(size_t k = 0; k < NumUnknowns; k++) {
            dco::derivative(dco::derivative(x_dco[k]) ) = 1.; // Seed

            DCO_BM::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint
            for(size_t j = 0; j <= k; j++) {
               values[ind] += lambda[i] * dco::derivative(dco::value(x_dco[j]) );
               ind++;
            }

            //dco::derivative(dco::derivative(x_dco[k]) ) = 0.; // Unseed
            DCO_BM::global_tape->zero_adjoints();
         }

         // Clear tape only, do not remove
         DCO_M::global_tape->reset();
         DCO_BM::global_tape->reset();
      }

      /// Order constraints ///

      if(ConsOrder >= 2) {
         for (size_t j = 0; j < NumUnknowns; j++) {
            DCO_BM::global_tape->register_variable(dco::value(x_dco[j]) ); // record active input
            DCO_BM::global_tape->register_variable(dco::derivative(x_dco[j]) ); // record active input

            dco::value(dco::value(x_dco[j]) ) = x[j];

            DCO_M::global_tape->register_variable(x_dco[j]);
         }

         if(OddDegree) {
            if(UseHull)
               g = SecOrder(x_dco, NumRoots, HullRealScaled, HullImagScaled);
            else
               g = SecOrder(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled);
         }
         else {
            if(UseHull)
               g = SecOrder(x_dco, NumRoots, HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, i_min);
            else
               g = SecOrder(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff, i_min);
         }

         dco::value(dco::derivative(g) ) = 1.; // Seed
         DCO_M::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint

         ind = 0;
         for(size_t k = 0; k < NumUnknowns; k++) {
            dco::derivative(dco::derivative(x_dco[k]) ) = 1.; // Seed

            DCO_BM::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint
            for(size_t j = 0; j <= k; j++) {
               values[ind] += lambda[NumEigVals] * dco::derivative(dco::value(x_dco[j]) );
               ind++;
            }

            DCO_BM::global_tape->zero_adjoints(); // Unseed
         }

         if(ConsOrder >= 3) {
            for (size_t j = 0; j < NumUnknowns; j++) {
               DCO_BM::global_tape->register_variable(dco::value(x_dco[j]) ); // record active input
               DCO_BM::global_tape->register_variable(dco::derivative(x_dco[j]) ); // record active input

               dco::value(dco::value(x_dco[j]) ) = x[j];

               DCO_M::global_tape->register_variable(x_dco[j]);
            }

            if(OddDegree) {
               if(UseHull)
                  g = ThirdOrder(x_dco, NumRoots, HullRealScaled, HullImagScaled);
               else
                  g = ThirdOrder(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled);
            }
            else {
               if(UseHull)
                  g = ThirdOrder(x_dco, NumRoots, HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, i_min);
               else
                  g = ThirdOrder(x_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff, i_min);
            }

            dco::value(dco::derivative(g) ) = 1.; // Seed
            DCO_M::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint

            ind = 0;
            for(size_t k = 0; k < NumUnknowns; k++) {
               dco::derivative(dco::derivative(x_dco[k]) ) = 1.; // Seed

               DCO_BM::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint
               for(size_t j = 0; j <= k; j++) {
                  values[ind] += lambda[NumEigVals + 1] * dco::derivative(dco::value(x_dco[j]) );
                  ind++;
               }

               DCO_BM::global_tape->zero_adjoints(); // Unseed
            }
         }
      }

      DCO_TT::remove(DCO_M::global_tape); // deallocate tape
      DCO_BTT::remove(DCO_BM::global_tape); // deallocate base tape
   }

   return true;
}
// [TNLP_eval_h]

// [TNLP_intermediate_callback]
bool Roots_Real::intermediate_callback(
   AlgorithmMode              mode,
   Index                      iter,
   Number                     obj_value,
   Number                     inf_pr,
   Number                     inf_du,
   Number                     mu,
   Number                     d_norm,
   Number                     regularization_size,
   Number                     alpha_du,
   Number                     alpha_pr,
   Index                      ls_trials,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{  
   //if(ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX) < std::numeric_limits<Number>::epsilon()) {  
   //if(ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX) < 1e-14) {  

   if(mode == RegularMode) { // Otherwise 'ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX)' returns erronous 0!
      // Check if current & best yet solution are of comparable feasibility
      if(fabs(ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX) - InfPr) < 1e-12) {
         // Then, check whether one has a larger timestep
         if(obj_value < Maxdt) {
            Maxdt = obj_value;
            InfPr = ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX);
            get_curr_iterate(ip_data, ip_cq, false, NumUnknowns, xMaxdt, NULL, NULL, NumConstr, NULL, NULL);
         }
      }
      else { // Case where current is significantly better then current best
         if(ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX) < InfPr) {
            InfPr = ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX);
            get_curr_iterate(ip_data, ip_cq, false, NumUnknowns, xMaxdt, NULL, NULL, NumConstr, NULL, NULL);
            // In that case: Always update timestep
            Maxdt = obj_value;
         }
      }
   }

   return true;
}
// [TNLP_intermediate_callback]

// [TNLP_finalize_solution]
void Roots_Real::finalize_solution(
   SolverReturn               status,
   Index                      n,
   const Number*              x,
   const Number*              z_L,
   const Number*              z_U,
   Index                      m,
   const Number*              g,
   const Number*              lambda,
   Number                     obj_value,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{
   std::cout.precision(std::numeric_limits<Number>::max_digits10);
   std::cout << std::endl << std::endl << std::endl << "### RESULTS ###" << std::endl;
   std::cout << std::endl << "Minimum primal infeasibility is: " 
             << InfPr << " with timestep: " << -Maxdt << std::endl;

   // For this example, we write the solution to the console
   std::cout << std::endl << "Solution:" << std::endl;
   for( Index i = 0; i < n-1; i++ ) {
      std::cout << "x[" << i << "] = " << xMaxdt[i] << std::endl;
   }

   std::cout << std::endl << std::endl << "Interpolated imaginary values:" << std::endl;
   for( Index i = 0; i < n-1; i++ ) {
      if(UseHull)
         std::cout << "y[" << i << "] = " << Lin_IntPol(xMaxdt[i], HullRealScaled, HullImagScaled) << std::endl;
      else
         std::cout << "y[" << i << "] = " << Lin_IntPol(xMaxdt[i], RealEigValsScaled, ImagEigValsScaled) 
                   << std::endl;
   }

   std::cout << std::endl << std::endl << "Optimal timestep: " << xMaxdt[NumRoots]
             << std::endl << "This corresponds to " << xMaxdt[NumRoots] / dtExp << " Efficiency" 
             << std::endl << "and a reference timestep of: " << dtExp 
             << std::endl << std::endl;

   /*
   std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
   for( Index i = 0; i < n; i++ )
   {
      std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
   }
   for( Index i = 0; i < n; i++ )
   {
      std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
   }

   std::cout << std::endl << "Lagrange multipliers:" << std::endl;
   for( Index i = 0; i < m; i++ )
   {
      std::cout << "lambda(" << i << ") = " << lambda[i] << std::endl;
   }
   */

   Number Constr[NumConstr];
   if(OddDegree)
      if(UseHull)
         StabConstr_Real(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                         HullRealScaled, HullImagScaled, dtExp);
      else 
         StabConstr_Real(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, dtExp);
   else {
      i_min = 0;
      for(size_t i = 1; i < NumRoots; i++) {
         if(xMaxdt[i] < xMaxdt[i_min])
            i_min = i;
      }

      if(UseHull)
         StabConstr_Real(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                         HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
      else 
         StabConstr_Real(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                         ImagDiff_over_RealDiff, dtExp, i_min);
   }

   std::cout << std::endl << "Final value of violating eigenvalue constraints:" << std::endl;
   for( Index i = 0; i < NumEigVals; i++ ) {
      if(Constr[i] > 1.)
         std::cout << "g(" << i << ") = " << Constr[i] << std::endl;
   }

   if(ConsOrder >= 2) {   
      if(OddDegree) {
         if(UseHull)
            SecOrder(xMaxdt, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
         else
            SecOrder(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
         }
      else {
         if(UseHull)
            SecOrder(xMaxdt, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                     ImagDiff_over_RealDiff, i_min);
         else
            SecOrder(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                     ImagDiff_over_RealDiff, i_min);
      }
      std::cout << std::endl << "Final value of 2nd order constraint: " << Constr[NumEigVals] << std::endl;

      if(ConsOrder >= 3) {   
         if(OddDegree) {
            if(UseHull)
               ThirdOrder(xMaxdt, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
            else
               ThirdOrder(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
            }
         else {
            if(UseHull)
               ThirdOrder(xMaxdt, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                          ImagDiff_over_RealDiff, i_min);
            else
               ThirdOrder(xMaxdt, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                          ImagDiff_over_RealDiff, i_min);
         }
         std::cout << "Final value of 3rd order constraint: " << Constr[NumEigVals+1] << std::endl;
      }
   }

   /// Write to file section ///

   std::ofstream RealOptFile("./Real_Optimized_" + std::to_string(NumStages) + ".txt");
   for(size_t i = 0; i < n; i++) {
      std::stringstream StringStr; // On purpose within loop (automatic reset)
      StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
      StringStr << xMaxdt[i];
      RealOptFile << StringStr.str();
      if(i != n-1)
         RealOptFile << "\n";
   }
   RealOptFile.close();

   // TODO: Write out with special treatment for min_real
   const std::string PE_FileName = "./PE_" + std::to_string(NumStages);
   std::ofstream PseudoExtremaFile(PE_FileName + ".txt");

   std::stringstream StringStr;
   StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
   StringStr << xMaxdt[i_min];
   PseudoExtremaFile << StringStr.str();
   PseudoExtremaFile << "\n";
   for(size_t i = 0; i < i_min; i++) {
      std::stringstream StringStr; // On purpose within loop (automatic reset)
      StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
      StringStr << xMaxdt[i] << "+";
      if(UseHull)
         StringStr << Lin_IntPol(xMaxdt[i], HullRealScaled, HullImagScaled) << "i";
      else
         StringStr << Lin_IntPol(xMaxdt[i], RealEigValsScaled, ImagEigValsScaled) << "i";
      PseudoExtremaFile << StringStr.str();
      if(i != n-2)
         PseudoExtremaFile << "\n";
   }
   for(size_t i = i_min + 1; i < n-1; i++) {
      std::stringstream StringStr; // On purpose within loop (automatic reset)
      StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
      StringStr << xMaxdt[i] << "+";
      if(UseHull)
         StringStr << Lin_IntPol(xMaxdt[i], HullRealScaled, HullImagScaled) << "i";
      else
         StringStr << Lin_IntPol(xMaxdt[i], RealEigValsScaled, ImagEigValsScaled) << "i";
      PseudoExtremaFile << StringStr.str();
      if(i != n-2)
         PseudoExtremaFile << "\n";
   }
   PseudoExtremaFile.close();
}
// [TNLP_finalize_solution]
