// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Daniel Doehring                  RWTH Aachen University 2022-11-20


#ifndef __RKCOEFFS_HPP__
#define __RKCOEFFS_HPP__

#include <vector>
#include <complex>
#include <cassert>
#include <algorithm> // For reversing a vector

// For check of stability property
//#include <eigen3/Eigen/Sparse>
// Required for determinant - unfortunately, eigen does not offer sparse matrix determinant computation
#include <eigen3/Eigen/Dense>

// Multi-precision
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp> // For quad, oct precision

#include <fstream>
#include <iostream>

//using MP_Real = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<71>>;
//using MP_Real = boost::multiprecision::cpp_dec_float_100;

// 16 Digits (for sanity check only)
//using MP_Real = double;

// 18 Digits
//using MP_Real = long double;

// 33 Digits:
using MP_Real = boost::multiprecision::cpp_bin_float_quad;

// 71 Digits:
//using MP_Real = boost::multiprecision::cpp_bin_float_oct;


std::vector<MP_Real> ComputeCoeffs(const bool OddDegree, const int ConsOrder, const int NumStages,
                                   const std::vector<double>& RootsReal, const std::vector<double>& RootsImag) {

  size_t NumRoots = 2 * RootsReal.size(); // True for Base Poly Odd Degree => Lower Degree Polynom Even Degree
  if(!OddDegree) // Base Poly Even Degree => Lower Degree Polynom Odd Degree
    NumRoots -= 1;

  std::vector<std::complex<MP_Real>> Roots = std::vector<std::complex<MP_Real>>(NumRoots);
  if(OddDegree) {
    for(size_t i = 0; i < RootsReal.size(); i++) {
      Roots[i] = std::complex<MP_Real>(RootsReal[i], RootsImag[i]);
      Roots[i+RootsReal.size()] = std::complex<MP_Real>(RootsReal[i], -RootsImag[i]);
    }
  }
  else { // Base Poly Even Degree => Lower Degree Polynom Odd Degree
    size_t i_min = 0;
    for(size_t i = 1; i < RootsReal.size(); i++)
      if(RootsReal[i] < RootsReal[i_min])
        i_min = i;

    Roots[0] = std::complex<MP_Real>(RootsReal[i_min], 0.);

    for(size_t i = 0; i < i_min; i++) {
      Roots[2*i+1] = std::complex<MP_Real>(RootsReal[i],  RootsImag[i]);
      Roots[2*i+2] = std::complex<MP_Real>(RootsReal[i], -RootsImag[i]);
    }

    for(size_t i = i_min + 1; i < RootsReal.size(); i++) {
      Roots[2*i-1] = std::complex<MP_Real>(RootsReal[i],  RootsImag[i]);
      Roots[2*i]   = std::complex<MP_Real>(RootsReal[i], -RootsImag[i]);
    }
  }

  // Initialize (lower-degree) polynomial with 1
  std::vector<std::complex<MP_Real>> MonCoeffsComplex{std::complex<MP_Real>(1., 0.)};
  for(size_t i = 0; i < NumRoots; i++) {
    std::vector<std::complex<MP_Real>> temp{MonCoeffsComplex.begin(), MonCoeffsComplex.end()};

    for(auto & item : temp) {
      item /= -Roots[i];
      // Stabilized version:
      //item /= -Roots[i] / static_cast<MP_Real>(NumStageEvals);
    }

    MonCoeffsComplex.push_back(std::complex<MP_Real>(0., 0.));
    for(size_t j = 1 ; j < MonCoeffsComplex.size(); j++) {
      MonCoeffsComplex[j] += temp[j-1];
    }

    // Truncate to real after each complex-conjugated pair (Note: Works only for even-degree!)
    if(i % 2 == 0 && !OddDegree)
      for(size_t j = 1 ; j < MonCoeffsComplex.size(); j++) {
        MonCoeffsComplex[j] = std::complex<MP_Real>(MonCoeffsComplex[j].real(), 0.0);
      }
  }

  /*
  std::cout << std::endl << "Monomial coefficients are: " << std::endl << std::endl;
  std::cout << std::setprecision(std::numeric_limits<std::complex<MP_Real>::value_type>::digits10);
  for(size_t i = 0; i < NumRoots; i++) {
    std::cout << MonCoeffsComplex[i].real() << std::endl << MonCoeffsComplex[i].imag() << std::endl << std::endl;
  }
  */

  // 'MonCoeffs' are the monomial coefficients of the lower degree polynomial
  std::vector<MP_Real> MonCoeffs = std::vector<MP_Real>(NumStages);
  for(size_t i = 0; i < NumStages; i++)
    MonCoeffs[i] = static_cast<MP_Real>(real(MonCoeffsComplex[i]));

  return MonCoeffs;
}

// Compute SE Factors based on observation
std::vector<MP_Real> Compute_SE_Factors(const int NumStages, const int NumStageEvals, const int ConsOrder) {
  std::vector<MP_Real> c(NumStages);

  for(size_t i = 0; i < NumStages; i++)
    c[i] = i / (static_cast<MP_Real>(2.0) * (NumStages - 1));

  std::vector<MP_Real> SE_Factors(NumStageEvals - ConsOrder);
  for(size_t i = 0; i < NumStageEvals - ConsOrder; i++) {
    SE_Factors[i] = c[NumStages - i - 2];
  }

  std::ofstream SE_File("SE_Factors" + std::to_string(NumStages) + "_" + std::to_string(NumStageEvals) + ".txt");
  for(size_t i = 0; i < NumStageEvals - ConsOrder; i++) {
    std::stringstream StringStr; // On purpose within loop (automatic reset)
    // Multi-precision
    StringStr << std::setprecision(std::numeric_limits<MP_Real>::max_digits10);
    StringStr << SE_Factors[i];

    SE_File << StringStr.str();
    if(i != NumStageEvals - ConsOrder - 1)
      SE_File << "\n";
  }
  SE_File.close();

  return SE_Factors;
}


// Chose float_type sinze even long long int runs out of "space" at ~22 digits
template<typename float_type>
inline float_type factorial(const size_t n, const std::vector<float_type>& TellCppWhatTypeToUse) {
  float_type res = 1.0;
  for(int i = 1; i <= n; i++)
    res *= i;
  return res;
}

template<typename float_type>
void CheckStability(const int NumStages, const int NumStageEvals, const int ConsOrder,
                    const std::vector<float_type>& MonCoeffs,
                    const std::vector<float_type>& a, const std::vector<float_type>& SE_Factors,
                    const int NumEigVals,
                    const std::vector<double>& RealEigValsScaled, const std::vector<double>& ImagEigValsScaled,
                    const double dtScaling) {

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10);
  std::complex<float_type> z, z_power, StabPnom;
  float_type AbsDet;
  for(size_t i = 0; i < NumEigVals; i++) {
    z = std::complex<float_type>(RealEigValsScaled[i]*dtScaling, ImagEigValsScaled[i]*dtScaling);

    StabPnom = std::complex<float_type>(1., 0.) + z;
    z_power = z;
    for (size_t i = 2; i <= ConsOrder; i++) {
      z_power *= z;
      StabPnom += z_power * MonCoeffs[i-1];
    }

    for (size_t i = ConsOrder + 1; i <= NumStageEvals; i++) {
      z_power *= z * a[i - (ConsOrder + 1)];
      StabPnom += z_power * SE_Factors[i-(ConsOrder + 1)];
      /*
      std::cout << z << std::endl << z_power << std::endl << a[i - (ConsOrder + 1)] << std::endl 
                << SE_Factors[i-(ConsOrder + 1)] << std::endl << std::endl;
      */
    }

    //AbsDet = std::abs(AbsDet); // Does not compile => Manual implementation
    AbsDet = sqrt(real(StabPnom)*real(StabPnom) + imag(StabPnom)*imag(StabPnom)); // Seems to be the right thing
    if(AbsDet > static_cast<float_type>(1.))
      std::cout << i <<"'th eigenvalue violates constraint with stability polynomial value: "
                << std::endl << AbsDet << std::endl << std::endl;
  }
}

template<typename float_type>
void CheckStability(const int NumStages, const int NumStageEvals, const int ConsOrder,
                    const std::vector<float_type>& a, const std::vector<float_type>& SE_Factors,
                    const int NumEigVals,
                    const std::vector<double>& RealEigValsScaled, const std::vector<double>& ImagEigValsScaled,
                    const double dtScaling) {

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10);
  std::complex<float_type> z, z_power, StabPnom;
  float_type AbsDet;
  for(size_t i = 0; i < NumEigVals; i++) {
    z = std::complex<float_type>(RealEigValsScaled[i]*dtScaling, ImagEigValsScaled[i]*dtScaling);

    StabPnom = std::complex<float_type>(1., 0.) + z;
    z_power = z;
    for (size_t i = 2; i <= ConsOrder; i++) {
      z_power *= z;
      StabPnom += z_power / factorial(i, a);
    }

    for (size_t i = ConsOrder + 1; i <= NumStageEvals; i++) {
      z_power *= z * a[i - (ConsOrder + 1)];
      StabPnom += z_power * SE_Factors[i-(ConsOrder + 1)];
      /*
      std::cout << z << std::endl << z_power << std::endl << a[i - (ConsOrder + 1)] << std::endl 
                << SE_Factors[i-(ConsOrder + 1)] << std::endl << std::endl;
      */
    }

    //AbsDet = std::abs(AbsDet); // Does not compile => Manual implementation
    AbsDet = sqrt(real(StabPnom)*real(StabPnom) + imag(StabPnom)*imag(StabPnom)); // Seems to be the right thing
    if(AbsDet > static_cast<float_type>(1.))
      std::cout << i <<"'th eigenvalue violates constraint with stability polynomial value: "
                << std::endl << AbsDet << std::endl << std::endl;
  }
}

template<typename float_type>
float_type MaxAbsPnom(const int NumStages, const int NumStageEvals, const int ConsOrder,
                      const std::vector<float_type>& a, const std::vector<float_type>& SE_Factors,
                      const int NumEigVals,
                      const std::vector<double>& RealEigValsScaled, const std::vector<double>& ImagEigValsScaled,
                      const float_type dtScaling) {
  std::complex<float_type> z, z_power, StabPnom;
  float_type AbsDet, AbsDetMax = 0.;
  for(size_t i = 0; i < NumEigVals; i++) {
    z = std::complex<float_type>(RealEigValsScaled[i] * dtScaling, ImagEigValsScaled[i] * dtScaling);

    StabPnom = std::complex<float_type>(1., 0.) + z;
    z_power = z;
    for (size_t i = 2; i <= ConsOrder; i++) {
      z_power *= z;
      StabPnom += z_power / factorial(i, a);
    }

    for (size_t i = ConsOrder + 1; i <= NumStageEvals; i++) {
      z_power *= z * a[i - (ConsOrder + 1)];

      StabPnom += z_power * SE_Factors[i-(ConsOrder + 1)];
    }

    //AbsDet = std::abs(AbsDet); // Does not compile => Manual implementation
    AbsDet = sqrt(real(StabPnom)*real(StabPnom) + imag(StabPnom)*imag(StabPnom)); // Seems to be the right thing
    if(AbsDet > AbsDetMax)
      AbsDetMax = AbsDet;
  }

  return AbsDetMax;
}

template<typename float_type>
float_type MaxAbsPnom(const int NumStages, const int NumStageEvals, const int ConsOrder,
                      const std::vector<float_type>& MonCoeffs,
                      const std::vector<float_type>& a, const std::vector<float_type>& SE_Factors,
                      const int NumEigVals,
                      const std::vector<double>& RealEigValsScaled, const std::vector<double>& ImagEigValsScaled,
                      const float_type dtScaling) {
  std::complex<float_type> z, z_power, StabPnom;
  float_type AbsDet, AbsDetMax = 0.;
  for(size_t i = 0; i < NumEigVals; i++) {
    z = std::complex<float_type>(RealEigValsScaled[i] * dtScaling, ImagEigValsScaled[i] * dtScaling);

    StabPnom = std::complex<float_type>(1., 0.) + z;
    z_power = z;
    for (size_t i = 2; i <= ConsOrder; i++) {
      z_power *= z;
      StabPnom += z_power * MonCoeffs[i-1];
    }

    for (size_t i = ConsOrder + 1; i <= NumStageEvals; i++) {
      z_power *= z * a[i - (ConsOrder + 1)];

      StabPnom += z_power * SE_Factors[i-(ConsOrder + 1)];
    }

    //AbsDet = std::abs(AbsDet); // Does not compile => Manual implementation
    AbsDet = sqrt(real(StabPnom)*real(StabPnom) + imag(StabPnom)*imag(StabPnom)); // Seems to be the right thing
    if(AbsDet > AbsDetMax)
      AbsDetMax = AbsDet;
  }

  return AbsDetMax;
}

/* NOTE: Seems unreliable, i.e., assumption of monotonicity in timestep as taken in 
   https://msp.org/camcos/2012/7-2/p04.xhtml 
   is not observed in practice!
*/ 
template<typename float_type>
float_type FindMaxTimeStep(const int NumStages, const int NumStageEvals, const int ConsOrder,
                           const std::vector<float_type>& a, const std::vector<float_type>& SE_Factors,
                           const int NumEigVals,
                           const std::vector<double>& RealEigValsScaled, const std::vector<double>& ImagEigValsScaled,
                           const double dtScalingOpt) {
  
  //float_type dtScalMax = 1.; // Suffering from non-monotonicity
  float_type dtScalMax = 2. * dtScalingOpt;
  float_type dtScalMin = 0.;
  const double dtScalEps = 1e-12;

  float_type Violation, dtScaling;
  while(dtScalMax - dtScalMin > dtScalEps) {
    dtScaling = 0.5 * (dtScalMax + dtScalMin);
    Violation = MaxAbsPnom(NumStages, NumStageEvals, ConsOrder, a, SE_Factors, NumEigVals, 
                           RealEigValsScaled, ImagEigValsScaled, dtScaling);

    if(Violation < static_cast<float_type>(1.0))
      dtScalMin = dtScaling;
    else
      dtScalMax = dtScaling;
  } 

  return dtScalMin; // Use conservative value (minimum)
}

/* NOTE: Seems unreliable, i.e., assumption of monotonicity in timestep as taken in 
   https://msp.org/camcos/2012/7-2/p04.xhtml 
   is not observed in practice!
*/ 
template<typename float_type>
float_type FindMaxTimeStep(const int NumStages, const int NumStageEvals, const int ConsOrder,
                           const std::vector<float_type>& MonCoeffs,
                           const std::vector<float_type>& a, const std::vector<float_type>& SE_Factors,
                           const int NumEigVals,
                           const std::vector<double>& RealEigValsScaled, const std::vector<double>& ImagEigValsScaled,
                           const double dtScalingOpt) {
  
  //float_type dtScalMax = 1.0; // Suffering from non-monotonicity
  float_type dtScalMax = 2. * dtScalingOpt;
  float_type dtScalMin = 0.;
  const double dtScalEps = 1e-12;

  float_type Violation, dtScaling;
  while(dtScalMax - dtScalMin > dtScalEps) {
    dtScaling = 0.5 * (dtScalMax + dtScalMin);
    Violation = MaxAbsPnom(NumStages, NumStageEvals, ConsOrder, MonCoeffs, a, SE_Factors, 
                           NumEigVals, RealEigValsScaled, ImagEigValsScaled, dtScaling);

    if(Violation < static_cast<float_type>(1.0))
      dtScalMin = dtScaling;
    else
      dtScalMax = dtScaling;
  } 

  return dtScalMin; // Use conservative value (minimum)
}

void compute_a_coeffs(const int NumStages, const int NumStageEvals, const bool OddDegree, const int ConsOrder,
                      const std::vector<double>& RootsReal, const std::vector<double>& RootsImag,
                      const int NumEigVals,
                      const std::vector<double>& RealEigValsScaled, const std::vector<double>& ImagEigValsScaled,
                      const double dt, const double dtExp) {

  std::cout << std::endl << std::endl << std::endl << "### Runge-Kutta Coefficient Computing ###" << std::endl;

  std::cout << std::endl << "Multi-precision type has " 
            << std::numeric_limits<MP_Real>::digits10 << " (significant) digits" << std::endl;

  const std::vector<MP_Real> MonCoeffs = ComputeCoeffs(OddDegree, ConsOrder, NumStages, RootsReal, RootsImag);

  std::ofstream MC_file("gamma_" + std::to_string(NumStages) + ".txt");
  for(size_t i = 0; i < NumStages; i++) {
    std::stringstream StringStr; // On purpose within loop (automatic reset)
    // Double precision
    //StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
    //StringStr << a_DBL[i];

    // Multi-precision
    StringStr << std::setprecision(std::numeric_limits<MP_Real>::max_digits10);
    StringStr << MonCoeffs[i];

    MC_file << StringStr.str();
    if(i != NumStages - 1)
      MC_file << "\n";
  }
  MC_file.close();

  
  const std::vector<MP_Real> SE_Factors = Compute_SE_Factors(NumStages, NumStageEvals, ConsOrder);

  std::vector<MP_Real> a_MP = std::vector<MP_Real>(MonCoeffs.begin() + ConsOrder, MonCoeffs.end());
  for(size_t i = 0; i < NumStageEvals - ConsOrder; i++) {
    a_MP[i] /= SE_Factors[i];
    // Stabilized version:
    //a_MP[i] /= SE_Factors[i] * std::pow(NumStageEvals, ConsOrder);
    for(size_t j = 0; j < i; j++) {
      a_MP[i] /= a_MP[j];
      // Stabilized version:
      //a_MP[i] /= a_MP[j] * NumStageEvals;
    }
  }

  /*
  // Sanity check: Re-compute gamma from a
  std::cout << "Sanity check: Re-compute gamma from a" << std::endl << std::endl;
  MP_Real temp;
  for(size_t i = 0; i < NumStageEvals - ConsOrder; i++) {
    // Stabilized version
    //temp = a_MP[i] * SE_Factors[i] * std::pow(NumStageEvals, ConsOrder);

    temp = a_MP[i] * SE_Factors[i];
    for(size_t j = 0; j < i; j++) {
      // Stabilized version
      //temp *= a_MP[j] * NumStageEvals;

      temp *= a_MP[j];
    }
    std::cout << temp << std::endl;
  }
  */

  std::cout << std::endl << "Checking stab. constr. computed from a-coeffs in Multiprecision " << std::endl << std::endl;
  CheckStability(NumStages, NumStageEvals, ConsOrder, MonCoeffs, a_MP, SE_Factors, NumEigVals, 
                 RealEigValsScaled, ImagEigValsScaled, dt/dtExp);

  const MP_Real dtScaling_MP = FindMaxTimeStep(NumStages, NumStageEvals, ConsOrder, MonCoeffs, a_MP, SE_Factors, 
                                               NumEigVals, RealEigValsScaled, ImagEigValsScaled, dt/dtExp);

  std::cout << std::endl << "In theory optimal timestep degenerates to: " << dtScaling_MP 
            << " of desired value, i.e., " << dtExp * dtScaling_MP  << std::endl << std::endl;

  // Conversion from MP to double
  std::vector<double> a_DBL(a_MP.begin(), a_MP.end());
  std::vector<double> SE_Factors_DBL(SE_Factors.begin(), SE_Factors.end());

  std::cout << std::endl << "Checking stab. constr. computed from a-coeffs in double " << std::endl << std::endl;
  CheckStability(NumStages, NumStageEvals, ConsOrder, a_DBL, SE_Factors_DBL, NumEigVals, 
                 RealEigValsScaled, ImagEigValsScaled, dt/dtExp);
  
  const double dtScaling = FindMaxTimeStep(NumStages, NumStageEvals, ConsOrder, a_DBL, SE_Factors_DBL, NumEigVals, 
                                           RealEigValsScaled, ImagEigValsScaled, dt/dtExp);
  std::cout << std::endl << "In theory optimal timestep degenerates to: " << dtScaling 
            << " of desired value, i.e., " << dtExp * dtScaling  << std::endl;
  

  // Need to flip coefficients since we compute the higher ones first
  std::reverse(a_MP.begin(), a_MP.end());
  std::reverse(a_DBL.begin(), a_DBL.end());

  std::ofstream a_file("a_Unknown" + std::to_string(NumStageEvals) + ".txt");
  for(size_t i = 0; i < NumStageEvals - ConsOrder; i++) {
    std::stringstream StringStr; // On purpose within loop (automatic reset)
    // Double precision
    //StringStr << std::setprecision(std::numeric_limits<Number>::max_digits10);
    //StringStr << a_DBL[i];

    // Multi-precision
    StringStr << std::setprecision(std::numeric_limits<MP_Real>::max_digits10);
    StringStr << a_MP[i];

    a_file << StringStr.str();
    if(i != NumStageEvals - ConsOrder - 1)
      a_file << "\n";
  }
  a_file.close();
}

#endif
