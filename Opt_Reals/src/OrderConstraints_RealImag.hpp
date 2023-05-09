// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Daniel Doehring                  RWTH Aachen University 2022-11-20

#ifndef __ORDERCONSTRAINTS_REALIMAG_HPP__
#define __ORDERCONSTRAINTS_REALIMAG_HPP__

#include <vector>

#include "Interpolation.hpp"

// NOTE: The constraints act on the pseudo/lower-degree polynomial!

// For Odd Base Polynom => Even Lower Degree Polynomial
template <typename T>
void SecOrder(const T* xy, T* g, const int NumRoots, const int NumEigVals,
              const std::vector<T>& RealRange, const std::vector<T>& ImagRange)
{
  g[NumEigVals] = 0.25;

  T b, Radius;

  for(size_t j = 0; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    g[NumEigVals] += xy[j]/Radius;
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template <typename T>
void SecOrder(const T* xy, T* g, const int NumRoots, const int NumEigVals,
              const std::vector<T>& RealRange, const std::vector<T>& ImagRange,
              const std::vector<T>& ImagDiff_over_RealDiff, const size_t i_min)
{
  g[NumEigVals] = 0.25;

  T b, Radius;

    g[NumEigVals] += 0.5/xy[i_min];

    for(size_t j = 0; j < i_min; j++) {
      b = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      g[NumEigVals] += xy[j]/Radius;
    }

    for(size_t j = i_min + 1; j < NumRoots; j++) {
      b = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
      Radius = xy[j]*xy[j] + b*b;

      g[NumEigVals] += xy[j]/Radius;
    }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void SecOrder(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
              const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange)
{
  g[NumEigVals] = 0.25;

  T b, Radius;

  for(size_t j = 0; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    g[NumEigVals] += xy[j]/Radius;
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void SecOrder(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
              const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange, 
              const std::vector<PT>& ImagDiff_over_RealDiff, const size_t i_min)
{
  g[NumEigVals] = 0.25;

  T b, Radius;
    g[NumEigVals] += 0.5/xy[i_min];

  for(size_t j = 0; j < i_min; j++) {
    b = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    g[NumEigVals] += xy[j]/Radius;
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    g[NumEigVals] += xy[j]/Radius;
  }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T SecOrder(const std::vector<T>& xy, const int NumRoots,
           const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange)
{
  T g = 0.25;

  T b, Radius;
  for(size_t j = 0; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    g += xy[j]/Radius;
  }

  return g;
}

// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T SecOrder(const std::vector<T>& xy, const int NumRoots,
           const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange,
           const std::vector<PT>& ImagDiff_over_RealDiff, const size_t i_min)
{
  T g = 0.25;

  T b, Radius;

  g += 0.5/xy[i_min];

  for(size_t j = 0; j < i_min; j++) {
    b = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    g += xy[j]/Radius;
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius = xy[j]*xy[j] + b*b;

    g += xy[j]/Radius;
  }

  return g;
}

#endif