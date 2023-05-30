// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Carl Laird, Andreas Waechter     IBM                    2005-08-16
//           Daniel Doehring                  RWTH Aachen University 2022-11-20

//#include "IpIpoptCalculatedQuantities.hpp" // To access current original, unscaled violations

#include "Roots_RealImag.hpp"

#include <cassert>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "IO_Funcs.hpp"
#include "StabConstraints_RealImag.hpp"
#include "OrderConstraints_RealImag.hpp"

// Post-Processing
#include "RKCoeffs.hpp"

using namespace Ipopt;


// constructor
Roots_RealImag::Roots_RealImag(
   const int NumStages_,
   const int ConsOrder_,
   const int NumStagesRef_, 
   const Number dtRef_,
   const std::string EigValFileName
) : NumStages(NumStages_), Degree(NumStages_), ConsOrder(ConsOrder_), NumStagesRef(NumStagesRef_), 
    dtRef(dtRef_), dtExp((dtRef / NumStagesRef_) *  NumStages_)
{
    std::cout << std::endl << "Optimize a " << Degree << " degree stability polynomial" << std::endl
              << std::endl << "Optimize roots of the " << Degree - 1 << " \"lower\" degree polynomial" 
              << std::endl << std::endl;

  std::cout << "The expected maximal stable timestep is: " << dtExp  << std::endl << std::endl;

  NumEigVals = -1;
  read_EigVals(EigValFileName, NumEigVals, RealEigValsScaled, ImagEigValsScaled);

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEigValsScaled[i] *= dtExp;
    ImagEigValsScaled[i] *= dtExp;
  }

  ImagDiff_over_RealDiff.resize(NumEigVals - 1);
  for(size_t i = 0; i < NumEigVals-1; i++) {
      ImagDiff_over_RealDiff[i] = (ImagEigValsScaled[i+1] - ImagEigValsScaled[i]) / 
                                  (RealEigValsScaled[i+1] - RealEigValsScaled[i]);
  }

  RealMin = *min_element(std::begin(RealEigValsScaled), std::end(RealEigValsScaled));
  ImagMax = *max_element(std::begin(ImagEigValsScaled), std::end(ImagEigValsScaled));

  NumConstr = NumEigVals + ConsOrder - 1; // -1 Since one consistency order comes for free

  OddDegree = Degree % 2;
  NumRoots  = NumStages / 2; // Note: Integer division is here desired
  
  NumUnknowns = 2 * NumRoots + 1;

  read_x0("./Real_Optimized_" + std::to_string(NumStages) + ".txt", xy0, NumUnknowns);
  std::cout << "Initial values are: " << std::endl;
  for(size_t i = 0; i < NumUnknowns; i++)
    std::cout << xy0[i] << std::endl;
  std::cout << std::endl;

  UseHull = false;
  std::cout << "Use hull for interpolation? " << UseHull << std::endl << std::endl;
}

Roots_RealImag::Roots_RealImag(
   const int NumStages_,
   const int ConsOrder_,
   const int NumStagesRef_, 
   const Number dtRef_,
   const std::string EigValFileName,
   const std::string HullPointPath
) : NumStages(NumStages_), Degree(NumStages_), ConsOrder(ConsOrder_), NumStagesRef(NumStagesRef_), 
    dtRef(dtRef_), dtExp((dtRef / NumStagesRef_) *  NumStages_)
{
  std::cout << std::endl << "Optimize a " << Degree << " degree stability polynomial" << std::endl
            << std::endl << "Optimize roots of the " << Degree - 1 << " \"lower\" degree polynomial" 
            << std::endl << std::endl;

  std::cout << "The expected maximal stable timestep is: " << dtExp  << std::endl << std::endl;

  NumEigVals = -1;
  read_EigVals(EigValFileName, NumEigVals, RealEigValsScaled, ImagEigValsScaled);

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEigValsScaled[i] *= dtExp;
    ImagEigValsScaled[i] *= dtExp;
  }

  RealMin = *min_element(std::begin(RealEigValsScaled), std::end(RealEigValsScaled));
  ImagMax = *max_element(std::begin(ImagEigValsScaled), std::end(ImagEigValsScaled));

  NumConstr = NumEigVals + ConsOrder - 1;  // -1 Since one consistency order comes for free

  OddDegree = Degree % 2;
  NumRoots  = NumStages / 2; // Note: Integer division is here desired
  
  NumUnknowns = 2 * NumRoots + 1;

  read_x0("./Real_Optimized_" + std::to_string(NumStages) + ".txt", xy0, NumUnknowns);
  std::cout << "Initial values are: " << std::endl;
  for(size_t i = 0; i < NumUnknowns; i++)
    std::cout << xy0[i] << std::endl;
  std::cout << std::endl;

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
}

// destructor
Roots_RealImag::~Roots_RealImag()
{}

// [TNLP_get_nlp_info]
// returns the size of the problem
bool Roots_RealImag::get_nlp_info(
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

   // Jacobian will in general be dense
   nnz_jac_g = NumConstr * NumUnknowns;

   // The Hessian is also dense and has NumUnknowns*NumUnknowns total nonzeros, but we
   // only need the lower left corner (since it is symmetric)
   nnz_h_lag = (NumUnknowns * (NumUnknowns + 1))/2;

   // use the C style indexing (0-based)
   index_style = TNLP::C_STYLE;

   return true;
}
// [TNLP_get_nlp_info]

// [TNLP_get_bounds_info]
// returns the variable bounds
bool Roots_RealImag::get_bounds_info(
   Index   n,
   Number* xy_l,
   Number* xy_u,
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

   // Real part
   // Place relatively strict bounds (trust previous result)
   // TODO: Normalize this by time, otherwise the margin becomes for larger times bigger, although more roots are used
   const Number RealPartMargin = 0.02 * std::abs(RealMin);
   std::cout << "RealPartMargin: " << RealPartMargin << std::endl << std::endl;
   for( Index i = 0; i < NumRoots; i++ ) {
      xy_l[i] = xy0[i] - RealPartMargin;
      
      // Division-by-zero guard
      xy_u[i] = std::min(xy0[i] + RealPartMargin, -1e-9);
   }

   // Imag part
   const Number IntPolErrorTol = 0.02 * ImagMax;
   std::cout << "IntPolErrorTol: " << IntPolErrorTol << std::endl << std::endl;
   for( Index i = NumRoots; i < NumUnknowns; i++ ) {
      xy_l[i] = std::max(-IntPolErrorTol, 0.); // Restrict to second quadrant
      xy_u[i] =  IntPolErrorTol;
   }

   // Timestep
   xy_l[NumUnknowns-1] = 0.;
   //xy_l[NumUnknowns-1] = xy0[NumRoots]; // Try to be at least as good as real-only solution
   
   xy_u[NumUnknowns-1] = 1.1 * dtExp; // Allow for some increase in maximum timestep

   for( Index i = 0; i < NumEigVals; i++ ) {
     g_l[i] = 0.; // Automatically ensured by abs, seems to benefit optimization
     //g_l[i] = -2e19; // Two times ipopt default value for treating number as inf - seems to harm here
     
     g_u[i] = 1.; // Actual stability bound
   }

   for ( Index i = NumEigVals; i < NumEigVals + ConsOrder - 1; i++ ) {
      g_l[i] = 0.;
      g_u[i] = 0.;
   }

   return true;
}
// [TNLP_get_bounds_info]

// [TNLP_get_starting_point]
// returns the initial point for the problem
bool Roots_RealImag::get_starting_point(
   Index   n,
   bool    init_xy,
   Number* xy,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
)
{
   // Here, we assume we only have starting values for x, if you code
   // your own NLP, you can provide starting values for the dual variables
   // if you wish
   assert(init_xy == true);

   // initialize to the given starting point
   for( Index i = 0; i < NumRoots; i++ ) {
      xy[i] = xy0[i];
   }

   for( Index i = NumRoots; i < NumUnknowns; i++ ) {
      xy[i] = 0.;
   }
   
   // TODO: Start for "simple" problems with maximum timestep here?
   xy[NumUnknowns - 1] = xy0[NumRoots];
   //xy[NumUnknowns - 1] = dtExp;

   if(init_z) {
     for(Index i = 0; i < NumUnknowns; i++) {
        z_L[i] = 0.;
        z_U[i] = 0.;
     }
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
bool Roots_RealImag::eval_f(
   Index         n,
   const Number* xy,
   bool          new_x,
   Number&       obj_value
)
{
   //assert(n == NumUnknowns); // Removed for performance

   obj_value = -xy[n-1]; // This is the (negated) timestep

   return true;
}
// [TNLP_eval_f]

// [TNLP_eval_grad_f]
// return the gradient of the objective function grad_{x} f(x)
bool Roots_RealImag::eval_grad_f(
   Index         n,
   const Number* xy,
   bool          new_x,
   Number*       grad_f
)
{
   //assert(n == NumUnknowns); // Removed for performances

   for(size_t i = 0; i < n-1; i++)
     grad_f[i] = 0.;
   
   grad_f[n-1] = -1.;

   return true;
}
// [TNLP_eval_grad_f]

// [TNLP_eval_g]
// return the value of the constraints: g(x)
bool Roots_RealImag::eval_g(
   Index         n,
   const Number* xy,
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
         StabConstr_RealImag(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, HullRealScaled, HullImagScaled, dtExp);
      else   
         StabConstr_RealImag(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, dtExp);

      if(ConsOrder >= 2) {
         if(UseHull)
            SecOrder(xy, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
         else
            SecOrder(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);

         if(ConsOrder >= 3) {
            if(UseHull)
               ThirdOrder(xy, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
            else
               ThirdOrder(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
         }
      }
   }
   else {
      i_min = 0;
      for(size_t i = 1; i < NumRoots; i++) {
         if(xy[i] < xy[i_min])
            i_min = i;
      }

      if(UseHull)
         StabConstr_RealImag(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                             HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
      else   
         StabConstr_RealImag(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                             ImagDiff_over_RealDiff, dtExp, i_min);

      if(ConsOrder >= 2) {
         if(UseHull)
            SecOrder(xy, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                     ImagDiff_over_RealDiff, i_min);
         else
            SecOrder(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                     ImagDiff_over_RealDiff, i_min);

         if(ConsOrder >= 3) {
            if(UseHull)
               ThirdOrder(xy, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                          ImagDiff_over_RealDiff, i_min);
            else
               ThirdOrder(xy, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                          ImagDiff_over_RealDiff, i_min);
         }
      }
   }
  
   return true;
}
// [TNLP_eval_g]

// [TNLP_eval_jac_g]
// return the structure or values of the Jacobian
bool Roots_RealImag::eval_jac_g(
   Index         n,
   const Number* xy,
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

   if( values == NULL )
   {
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

      std::vector<DCO_T> xy_dco(NumUnknowns);
      for(size_t i = 0; i < NumUnknowns; i++)
         dco::value(xy_dco[i]) = xy[i];

      DCO_M::global_tape = DCO_TT::create(); // Prepare / "touch" tape
      DCO_M::global_tape->register_variable(xy_dco);

      std::vector<DCO_T> g(m);
      if(OddDegree) {
         if(UseHull)
            StabConstr_RealImag(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, HullRealScaled, HullImagScaled, dtExp);
         else
            StabConstr_RealImag(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, dtExp);
      }
      else {
         i_min = 0;
         for(size_t i = 1; i < NumRoots; i++) {
            if(xy[i] < xy[i_min])
               i_min = i;
         }

         if(UseHull)
            StabConstr_RealImag(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                                HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
         else
            StabConstr_RealImag(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                                ImagDiff_over_RealDiff, dtExp, i_min);
      }

      DCO_M::global_tape->register_output_variable(g); // Record active output

      size_t ind = 0;
      for(size_t i = 0; i < NumEigVals; i++) {
         dco::derivative(g)[i] = 1.; // Seed first component

         DCO_M::global_tape->interpret_adjoint(); // Interpret (stored) tape

         // Harvest
         for (size_t j = 0; j < NumUnknowns; j++) {
            values[ind] = dco::derivative(xy_dco[j]);
            ind++;
         }

         // Reset adjoints
         DCO_M::global_tape->zero_adjoints();
      }

      if(ConsOrder >= 2) {
         if(OddDegree) {
            if(UseHull)
               SecOrder(xy_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
            else            
               SecOrder(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
         }
         else {
            if(UseHull)
               SecOrder(xy_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                        ImagDiff_over_RealDiff, i_min);
            else            
               SecOrder(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                        ImagDiff_over_RealDiff, i_min);
         }
         
         dco::derivative(g)[NumEigVals] = 1.; // Seed component

         DCO_M::global_tape->interpret_adjoint(); // Interpret (stored) tape

         // Harvest
         for (size_t j = 0; j < NumUnknowns; j++) {
            values[ind] = dco::derivative(xy_dco[j]);
            ind++;
         }

         if(ConsOrder >= 3) {
            // Reset adjoints
            DCO_M::global_tape->zero_adjoints();
            
            if(OddDegree) {
               if(UseHull)
                  ThirdOrder(xy_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
               else            
                  ThirdOrder(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
            }
            else {
               if(UseHull)
                  ThirdOrder(xy_dco, g, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                             ImagDiff_over_RealDiff, i_min);
               else            
                  ThirdOrder(xy_dco, g, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                             ImagDiff_over_RealDiff, i_min);
            }
            
            dco::derivative(g)[NumEigVals+1] = 1.; // Seed component

            DCO_M::global_tape->interpret_adjoint(); // Interpret (stored) tape

            // Harvest
            for (size_t j = 0; j < NumUnknowns; j++) {
               values[ind] = dco::derivative(xy_dco[j]);
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
bool Roots_RealImag::eval_h(
   Index         n,
   const Number* xy,
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
   // Removed for efficiency
   /*
   assert(n == NumUnknowns);
   assert(m == NumEigVals+ConsOrder-1);
   */

   if( values == NULL )
   {
      // return the structure. This is a symmetric matrix, fill the lower left
      // triangle only.

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
   else
   {
      // return the values. This is a symmetric matrix, fill the lower left
      // triangle only

      // fill the objective portion
      // Not existent for this problem (feasibility problem)
      // Constraints
        
      using DCO_BM  = typename dco::ga1s<Number>; // adjoint (base) mode
      using DCO_BT  = typename DCO_BM::type; // adjoint (base) type
      using DCO_BTT = typename DCO_BM::tape_t; // base tape type
      using DCO_M   = typename dco::ga1s<DCO_BT>; // adjoint of base mode
      using DCO_T   = typename DCO_M::type; // adjoint of base type
      using DCO_TT  = typename DCO_M::tape_t; // tape type

      std::vector<DCO_T> xy_dco(NumUnknowns);

      DCO_BM::global_tape = DCO_BTT::create(); // "touch" base tape
      DCO_M::global_tape  = DCO_TT::create(); // "touch" tape

      DCO_T g; // Scalar output

      for (size_t i = 0; i < NumUnknowns; i++) {
         DCO_BM::global_tape->register_variable(dco::value(xy_dco[i]) ); // record active input
         DCO_BM::global_tape->register_variable(dco::derivative(xy_dco[i]) ); // record active input

         dco::value(dco::value(xy_dco[i]) ) = xy[i];

         DCO_M::global_tape->register_variable(xy_dco[i]);
      }

      if(OddDegree) {
         if(UseHull)
            g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, 
                                      HullRealScaled, HullImagScaled, dtExp);
         else
            g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, dtExp);
      }
      else {
         i_min = 0;
         for(size_t i = 1; i < NumRoots; i++) {
            if(xy[i] < xy[i_min])
               i_min = i;
         }

         if(UseHull)
            g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, 
                                      HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
         else
            g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 0, 
                                      ImagDiff_over_RealDiff, dtExp, i_min);
      }

      dco::value(dco::derivative(g) ) = 1.; // Seed
      DCO_M::global_tape->interpret_adjoint();

      int ind = 0;
      for(size_t i = 0; i < NumUnknowns; i++) {
         dco::derivative(dco::derivative(xy_dco[i]) ) = 1.; // Seed

         DCO_BM::global_tape->interpret_adjoint();
         for(size_t j = 0; j <= i; j++) {
            values[ind] = lambda[0] * dco::derivative(dco::value(xy_dco[j]) );
            ind++;
         }

         //dco::derivative(dco::derivative(xy_dco[i]) ) = 0.; // Unseed

         DCO_BM::global_tape->zero_adjoints();
      }

      // Clear tape only, do not remove
      DCO_M::global_tape->reset();
      DCO_BM::global_tape->reset();

      /// Remaining Eigenvalues ///
      for(size_t i = 1; i < NumEigVals; i++) {
         for (size_t j = 0; j < NumUnknowns; j++) {
            DCO_BM::global_tape->register_variable(dco::value(xy_dco[j]) ); // record active input
            DCO_BM::global_tape->register_variable(dco::derivative(xy_dco[j]) ); // record active input

            dco::value(dco::value(xy_dco[j]) ) = xy[j];

            DCO_M::global_tape->register_variable(xy_dco[j]);
         }

         if(OddDegree) {
            if(UseHull)
               g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, 
                                         HullRealScaled, HullImagScaled, dtExp);
            else
               g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, dtExp);
         }
         else {
            if(UseHull)
               g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, 
                                         HullRealScaled, HullImagScaled, 
                                         ImagDiff_over_RealDiff, dtExp, i_min);
            else
               g = StabConstr_RealImag_i(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, i, 
                                         ImagDiff_over_RealDiff, dtExp, i_min);
         }

         dco::value(dco::derivative(g) ) = 1.; // Seed
         DCO_M::global_tape->interpret_adjoint();

         ind = 0;
         for(size_t k = 0; k < NumUnknowns; k++) {
            dco::derivative(dco::derivative(xy_dco[k]) ) = 1.; // Seed

            DCO_BM::global_tape->interpret_adjoint();
            for(size_t j = 0; j <= k; j++) {
               values[ind] += lambda[i] * dco::derivative(dco::value(xy_dco[j]) );
               ind++;
            }

            //dco::derivative(dco::derivative(xy_dco[k]) ) = 0.; // Unseed

            DCO_BM::global_tape->zero_adjoints();
         }
         // Clear tape only, do not remove
         DCO_M::global_tape->reset();
         DCO_BM::global_tape->reset();
      }

      /// Stability constraints ///

      if(ConsOrder >= 2) {
         for (size_t j = 0; j < NumUnknowns; j++) {
            DCO_BM::global_tape->register_variable(dco::value(xy_dco[j]) ); // record active input
            DCO_BM::global_tape->register_variable(dco::derivative(xy_dco[j]) ); // record active input

            dco::value(dco::value(xy_dco[j]) ) = xy[j];

            DCO_M::global_tape->register_variable(xy_dco[j]);
         }

         if(OddDegree) {
            if(UseHull)
               g = SecOrder(xy_dco, NumRoots, HullRealScaled, HullImagScaled);
            else
               g = SecOrder(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled);
         }
         else {
            if(UseHull)
               g = SecOrder(xy_dco, NumRoots, HullRealScaled, HullImagScaled, 
                            ImagDiff_over_RealDiff, i_min);
            else
               g = SecOrder(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 
                            ImagDiff_over_RealDiff, i_min);
         }

         dco::value(dco::derivative(g) ) = 1.; // Seed
         DCO_M::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint

         ind = 0;
         for(size_t k = 0; k < NumUnknowns; k++) {
            dco::derivative(dco::derivative(xy_dco[k]) ) = 1.; // Seed

            DCO_BM::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint
            for(size_t j = 0; j <= k; j++) {
               values[ind] += lambda[NumEigVals] * dco::derivative(dco::value(xy_dco[j]) );
               ind++;
            }

            DCO_BM::global_tape->zero_adjoints(); // Unseed
         }

         if(ConsOrder >= 3) {
            for (size_t j = 0; j < NumUnknowns; j++) {
               DCO_BM::global_tape->register_variable(dco::value(xy_dco[j]) ); // record active input
               DCO_BM::global_tape->register_variable(dco::derivative(xy_dco[j]) ); // record active input

               dco::value(dco::value(xy_dco[j]) ) = xy[j];

               DCO_M::global_tape->register_variable(xy_dco[j]);
            }

            if(OddDegree) {
               if(UseHull)
                  g = ThirdOrder(xy_dco, NumRoots, HullRealScaled, HullImagScaled);
               else
                  g = ThirdOrder(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled);
            }
            else {
               if(UseHull)
                  g = ThirdOrder(xy_dco, NumRoots, HullRealScaled, HullImagScaled, 
                                 ImagDiff_over_RealDiff, i_min);
               else
                  g = ThirdOrder(xy_dco, NumRoots, RealEigValsScaled, ImagEigValsScaled, 
                                 ImagDiff_over_RealDiff, i_min);
            }

            dco::value(dco::derivative(g) ) = 1.; // Seed
            DCO_M::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint

            ind = 0;
            for(size_t k = 0; k < NumUnknowns; k++) {
               dco::derivative(dco::derivative(xy_dco[k]) ) = 1.; // Seed

               DCO_BM::global_tape->interpret_adjoint(); // Back-propagate from output/adjoint
               for(size_t j = 0; j <= k; j++) {
                  values[ind] += lambda[NumEigVals+1] * dco::derivative(dco::value(xy_dco[j]) );
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
bool Roots_RealImag::intermediate_callback(
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
   // Terminate at first feasible point.
   /*
   if(mode == RegularMode) { // Otherwise 'ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX)' returns erronous 0!
      // Terminate at first feasible point.
      if(ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX) <= std::numeric_limits<Number>::epsilon())
         return false;
      else
         return true;
   }
   */
   return true;
}
// [TNLP_intermediate_callback]

// [TNLP_finalize_solution]
void Roots_RealImag::finalize_solution(
   SolverReturn               status,
   Index                      n,
   const Number*              xy,
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
   // here is where we would store the solution to variables, or write to a file, etc
   // so we could use the solution.
   std::cout.precision(std::numeric_limits<Number>::max_digits10);
   std::cout << std::endl << std::endl << std::endl << "### RESULTS ###" << std::endl;

   std::vector<Number> Reals(NumRoots);
   std::vector<Number> Imags(NumRoots);

   // For this example, we write the solution to the console
   std::cout << std::endl << std::endl << "Real part of optimized roots:" << std::endl;
   for( Index i = 0; i < n/2; i++ ) {
      std::cout << "xy[" << i << "] = " << xy[i] << std::endl;
      Reals[i] = xy[i];
   }

   std::cout << std::endl << std::endl << "Corrections of real part of optimized roots:" << std::endl;
   for( Index i = 0; i < n/2; i++ ) {
      std::cout << "xy - xy0 [" << i << "] = " << xy[i] - xy0[i] << std::endl;
   }

   std::cout << std::endl << std::endl << "Imaginary part of optimized roots:" << std::endl;
   for( Index i = 0; i < n/2; i++ ) {
      if(UseHull) {
         Imags[i] = Lin_IntPol(xy[i], HullRealScaled, HullImagScaled) + xy[i+n/2];
         std::cout << "y[" << i << "] = " << Imags[i] << std::endl;
      }
      else {
         Imags[i] = Lin_IntPol(xy[i], RealEigValsScaled, ImagEigValsScaled) + xy[i+n/2];
         std::cout << "y[" << i << "] = " << Imags[i] << std::endl;
      }
   }

   std::cout << std::endl << std::endl << "Corrections of imaginary part of optimized roots:" << std::endl;
   for( Index i = n/2; i < n - 1; i++ ) {
      std::cout << "xy[" << i << "] = " << xy[i] << std::endl;
   }

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

   std::cout << std::endl << std::endl << "Optimal timestep: " << xy[n-1]
             << std::endl << "This corresponds to " << xy[n-1] / dtExp << " Efficiency" 
             << std::endl << "and a reference timestep of: " << dtExp 
             << std::endl << std::endl;

   Number Constr[NumConstr];
   if(OddDegree)
      if(UseHull)
         StabConstr_RealImag(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, HullRealScaled, HullImagScaled, dtExp);
      else 
         StabConstr_RealImag(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, dtExp);
   else {
      i_min = 0;
      for(size_t i = 1; i < NumRoots; i++) {
         if(xy[i] < xy[i_min])
            i_min = i;
      }

      if(UseHull)
         StabConstr_RealImag(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                             HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff, dtExp, i_min);
      else 
         StabConstr_RealImag(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                             ImagDiff_over_RealDiff, dtExp, i_min);
   }

   std::cout << std::endl << "Final value of violating eigenvalue constraints:" << std::endl;
   for( Index i = 0; i < NumEigVals; i++ )
   {
      if(Constr[i] >= 1.)
         std::cout << "g(" << i << ") = " << Constr[i] << std::endl;
   }

   if(ConsOrder >= 2) {   
      if(OddDegree) {
         if(UseHull)
            SecOrder(xy, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
         else
            SecOrder(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
         }
      else {
         if(UseHull)
            SecOrder(xy, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                     ImagDiff_over_RealDiff, i_min);
         else
            SecOrder(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                     ImagDiff_over_RealDiff, i_min);
      }
      std::cout << std::endl << "Final value of 2nd order constraint: " << Constr[NumEigVals] << std::endl;

      if(ConsOrder >= 3) {   
         if(OddDegree) {
            if(UseHull)
               ThirdOrder(xy, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled);
            else
               ThirdOrder(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled);
            }
         else {
            if(UseHull)
               ThirdOrder(xy, Constr, NumRoots, NumEigVals, HullRealScaled, HullImagScaled, 
                          ImagDiff_over_RealDiff, i_min);
            else
               ThirdOrder(xy, Constr, NumRoots, NumEigVals, RealEigValsScaled, ImagEigValsScaled, 
                          ImagDiff_over_RealDiff, i_min);
         }
         std::cout << "Final value of 3rd order constraint: " << Constr[NumEigVals+1] << std::endl;
      }
   }

   std::ofstream RealImagOptFile("./RealImag_Optimized_" + std::to_string(NumStages) + ".txt");
   for(size_t i = 0; i < n/2; i++) {
      std::stringstream StringStr; // On purpose within loop (automatic reset)
      StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
      StringStr << xy[i];
      StringStr << " + ";
      if(UseHull)
         StringStr << Lin_IntPol(xy[i], HullRealScaled, HullImagScaled) + xy[i+n/2];
      else
         StringStr << Lin_IntPol(xy[i], RealEigValsScaled, ImagEigValsScaled) + xy[i+n/2];
      StringStr << "i";
      RealImagOptFile << StringStr.str();
      if(i != n/2 - 1)
         RealImagOptFile << "\n";
   }
   RealImagOptFile.close();

   int i_min = 0;
   if(!OddDegree)
      for(size_t i = 1; i < NumRoots; i++)
         if(xy[i] < xy[i_min])
            i_min = i;

   std::ofstream PureRealFile("./PureReal" + std::to_string(NumStages) + ".txt");
   std::stringstream StringStr; // On purpose within loop (automatic reset)
   StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
   StringStr << xy[i_min];
   PureRealFile << StringStr.str();
   PureRealFile.close();

   std::ofstream TrueComplexFile("./TrueComplex" + std::to_string(NumStages) + ".txt");
   for(size_t j = 0; j < i_min; j++) {
      std::stringstream StringStr; // On purpose within loop (automatic reset)
      StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
      StringStr << xy[j];
      StringStr << "+";
      if(UseHull)
         StringStr << Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j+n/2];
      else
         StringStr << Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j+n/2];
      StringStr << "im";
      TrueComplexFile << StringStr.str();
      if(j != i_min - 1)
         TrueComplexFile << "\n";
   }

   for(size_t j = i_min + 1; j < NumRoots; j++) {
      std::stringstream StringStr; // On purpose within loop (automatic reset)
      StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
      StringStr << xy[j];
      StringStr << "+";
      if(UseHull)
         StringStr << Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j+n/2];
      else
         StringStr << Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j+n/2];
      StringStr << "im";
      TrueComplexFile << StringStr.str();
      if(j != NumRoots - 1)
         TrueComplexFile << "\n";
    }
    TrueComplexFile.close();

   // Make P-ERK ready
   compute_a_coeffs(NumStages, NumStages, OddDegree, ConsOrder, Reals, Imags, 
                    NumEigVals, RealEigValsScaled, ImagEigValsScaled, xy[n-1], dtExp);
}
// [TNLP_finalize_solution]
