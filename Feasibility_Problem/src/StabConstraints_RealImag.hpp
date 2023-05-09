// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Daniel Doehring                  RWTH Aachen University 2022-11-20

#ifndef __STABCONSTRAINTSREALIMAG_HPP__
#define __STABCONSTRAINTSREALIMAG_HPP__

#include "Interpolation.hpp"

#include <complex>
#include <vector>

// For Odd Base Polynom => Even Lower Degree Polynomial
template <typename T>
void StabConstr_RealImag(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                         const std::vector<T>& RealEigValsScaled, const std::vector<T>& ImagEigValsScaled,
                         const std::vector<T>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    b = Lin_IntPol(xy[0], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[0 + NumRoots];
    Radius = xy[0]*xy[0] + b*b;
    
    Real = (xy[0] * (xy[0] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag = (b*RealEV - xy[0]*ImagEV) / Radius;
    Prod = std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);

    for(size_t j = 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template <typename T>
void StabConstr_RealImag(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                         const std::vector<T>& RealEigValsScaled, const std::vector<T>& ImagEigValsScaled,
                         const size_t i_min, const std::vector<T>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    Prod = std::complex<T>(1. - RealEV / xy[i_min], -ImagEV / xy[i_min]);

    for(size_t j = 0; j < i_min; j++) {
      b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    for(size_t j = i_min + 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template <typename T>
void StabConstr_RealImag(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                         const std::vector<T>& RealEigValsScaled, const std::vector<T>& ImagEigValsScaled,
                         const std::vector<T>& HullRealScaled, const std::vector<T>& HullImagScaled,
                         const std::vector<T>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    b = Lin_IntPol(xy[0], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[0 + NumRoots];
    Radius = xy[0]*xy[0] + b*b;
    
    Real = (xy[0] * (xy[0] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag = (b*RealEV - xy[0]*ImagEV) / Radius;
    Prod = std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);

    for(size_t j = 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template <typename T>
void StabConstr_RealImag(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                         const std::vector<T>& RealEigValsScaled, const std::vector<T>& ImagEigValsScaled,
                         const std::vector<T>& HullRealScaled, const std::vector<T>& HullImagScaled,
                         const size_t i_min, const std::vector<T>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    Prod = std::complex<T>(1. - RealEV / xy[i_min], -ImagEV / xy[i_min]);

    for(size_t j = 0; j < i_min; j++) {
      b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    for(size_t j = i_min + 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}


// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void StabConstr_RealImag(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                         const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                         const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    b = Lin_IntPol(xy[0], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[0 + NumRoots];
    Radius = xy[0]*xy[0] + b*b;
    
    Real = (xy[0] * (xy[0] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag = (b*RealEV - xy[0]*ImagEV) / Radius;
    Prod = std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);

    for(size_t j = 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void StabConstr_RealImag(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                         const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                         const size_t i_min, const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    Prod = std::complex<T>(1. - RealEV / xy[i_min], -ImagEV / xy[i_min]);

    for(size_t j = 0; j < i_min; j++) {
      b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    for(size_t j = i_min + 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void StabConstr_RealImag(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                         const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                         const std::vector<PT>& HullRealScaled, const std::vector<PT>& HullImagScaled,
                         const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    b = Lin_IntPol(xy[0], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[0 + NumRoots];
    Radius = xy[0]*xy[0] + b*b;
    
    Real = (xy[0] * (xy[0] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag = (b*RealEV - xy[0]*ImagEV) / Radius;
    Prod = std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);

    for(size_t j = 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void StabConstr_RealImag(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                         const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                         const std::vector<PT>& HullRealScaled, const std::vector<PT>& HullImagScaled,
                         const size_t i_min, const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T RealEV, ImagEV, b, Radius, Real, Imag;
  std::complex<T> Prod;

  for(size_t i = 0; i < NumEigVals; i++) {
    RealEV = RealEigValsScaled[i];
    ImagEV = ImagEigValsScaled[i];

    Prod = std::complex<T>(1. - RealEV / xy[i_min], -ImagEV / xy[i_min]);

    for(size_t j = 0; j < i_min; j++) {
      b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    for(size_t j = i_min + 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
      Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
      Prod *= std::complex<T>(Real, Imag);

      Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
    }

    // From lower-degree to actual stability polynomial
    Prod *= std::complex<T>(RealEV, ImagEV);
    Prod += 1.;

    g[i] = std::abs(Prod);
  }
}


// Constraint for Hessian: Compute one constraint at a time
// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T StabConstr_RealImag_i(const std::vector<T>& xy, const int NumRoots,
                        const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                        const int EigValInd, const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T b, Radius, Real, Imag;
  std::complex<T> Prod;

  const T RealEV = RealEigValsScaled[EigValInd];
  const T ImagEV = ImagEigValsScaled[EigValInd];

  b = Lin_IntPol(xy[0], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[0 + NumRoots];
  Radius = xy[0]*xy[0] + b*b;
  
  Real = (xy[0] * (xy[0] - RealEV) + b * (b - ImagEV)) / Radius;
  Imag = (b*RealEV - xy[0]*ImagEV) / Radius;
  Prod = std::complex<T>(Real, Imag);

  Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);

  for(size_t j = 1; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
    Prod *= std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
  }

  // From lower-degree to actual stability polynomial
  Prod *= std::complex<T>(RealEV, ImagEV);
  Prod += 1.;

  return std::abs(Prod);
}

// Constraint for Hessian: Compute one constraint at a time
// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T StabConstr_RealImag_i(const std::vector<T>& xy, const int NumRoots,
                        const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                        const int EigValInd, const size_t i_min, const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T b, Radius, Real, Imag;
  std::complex<T> Prod;

  const T RealEV = RealEigValsScaled[EigValInd];
  const T ImagEV = ImagEigValsScaled[EigValInd];

  Prod = std::complex<T>(1. - RealEV / xy[i_min], -ImagEV / xy[i_min]);

  for(size_t j = 0; j < i_min; j++) {
    b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
    Prod *= std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], RealEigValsScaled, ImagEigValsScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
    Prod *= std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
  }

  // From lower-degree to actual stability polynomial
  Prod *= std::complex<T>(RealEV, ImagEV);
  Prod += 1.;

  return std::abs(Prod);
}

// Constraint for Hessian: Compute one constraint at a time
// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T StabConstr_RealImag_i(const std::vector<T>& xy, const int NumRoots,
                        const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                        const int EigValInd,
                        const std::vector<PT>& HullRealScaled, const std::vector<PT>& HullImagScaled,
                        const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T b, Radius, Real, Imag;
  std::complex<T> Prod;

  const T RealEV = RealEigValsScaled[EigValInd];
  const T ImagEV = ImagEigValsScaled[EigValInd];

  b = Lin_IntPol(xy[0], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[0 + NumRoots];
  Radius = xy[0]*xy[0] + b*b;
  
  Real = (xy[0] * (xy[0] - RealEV) + b * (b - ImagEV)) / Radius;
  Imag = (b*RealEV - xy[0]*ImagEV) / Radius;
  Prod = std::complex<T>(Real, Imag);

  Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);

  for(size_t j = 1; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
    Prod *= std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
  }

  // From lower-degree to actual stability polynomial
  Prod *= std::complex<T>(RealEV, ImagEV);
  Prod += 1.;

  return std::abs(Prod);
}

// Constraint for Hessian: Compute one constraint at a time
// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T StabConstr_RealImag_i(const std::vector<T>& xy, const int NumRoots,
                        const std::vector<PT>& RealEigValsScaled, const std::vector<PT>& ImagEigValsScaled,
                        const int EigValInd,
                        const std::vector<PT>& HullRealScaled, const std::vector<PT>& HullImagScaled,
                        const size_t i_min, const std::vector<PT>& ImagDiff_over_RealDiff)
{
  T b, Radius, Real, Imag;
  std::complex<T> Prod;

  const T RealEV = RealEigValsScaled[EigValInd];
  const T ImagEV = ImagEigValsScaled[EigValInd];

  Prod = std::complex<T>(1. - RealEV / xy[i_min], -ImagEV / xy[i_min]);

  for(size_t j = 0; j < i_min; j++) {
    b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
    Prod *= std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], HullRealScaled, HullImagScaled, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    Real  = (xy[j]*(xy[j] - RealEV) + b * (b - ImagEV)) / Radius;
    Imag  = (b*RealEV - xy[j]*ImagEV) / Radius;
    Prod *= std::complex<T>(Real, Imag);

    Prod *= std::complex<T>(Real + 2.*b*ImagEV/Radius, Imag - 2.*b*RealEV/Radius);
  }

  // From lower-degree to actual stability polynomial
  Prod *= std::complex<T>(RealEV, ImagEV);
  Prod += 1.;

  return std::abs(Prod);
}

#endif
