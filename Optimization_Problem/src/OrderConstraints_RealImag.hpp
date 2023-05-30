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

/*

For second order, we demand: 0.5 = - sum_i^S 1/r_i.
Since we have that for complex-conjugated roots 1/r + 1/r* = 2 * Re(r)/Radius(r)^2 we can simplify to 
0.5 = - sum_i^(S/2) 2 * Re(r_i)/Radius(r_i)^2
=> 0.25 = - sum_i^(S/2) Re(r_i)/Radius(r_i)^2

For third order, we demand 1/6 = sum_{i,j i != j}^S 1/(r_i r_j)
For complex conjugated roots we have 
1/(r_i r_j) + 1/(r_i r_j*) + 1/(r_i* r_j) + 1/(r_i* r_j*) 
= 4 * Re(r_i) * Re(r_j) / (Radius(r_i)^2 + Radius(r_j)^2)
=> 1.0/24.0 = sum_{i,j i != j}^(S/2) Re(r_i) * Re(r_j) / (Radius(r_i)^2 + Radius(r_j)^2)

*/

/// SECOND ORDER ///

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


/// THIRD ORDER ///


// For Odd Base Polynom => Even Lower Degree Polynomial
template <typename T>
void ThirdOrder(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                const std::vector<T>& RealRange, const std::vector<T>& ImagRange)
{
  g[NumEigVals+1] = 1.0/24.0;

  T b1, b2, Radius1, Radius2;

  for(size_t j = 0; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.25 / Radius1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template <typename T>
void ThirdOrder(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                const std::vector<T>& RealRange, const std::vector<T>& ImagRange,
                const std::vector<T>& ImagDiff_over_RealDiff, const size_t i_min)
{
  g[NumEigVals+1] = 1.0/24.0;

  T b1, b2, Radius1, Radius2;

  // Handle combination: Real with complex conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.5 * xy[j]/(xy[i_min] * Radius1);
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.5 * xy[j]/(xy[i_min] * Radius1);
  }

  // Now: Complex Conjugated <-> Complex Conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.25 / Radius1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.25 / Radius1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void ThirdOrder(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange)
{
  g[NumEigVals+1] = 1.0/24.0;

  T b1, b2, Radius1, Radius2;

  for(size_t j = 0; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.25 / Radius1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void ThirdOrder(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange, 
                const std::vector<PT>& ImagDiff_over_RealDiff, const size_t i_min)
{
  g[NumEigVals+1] = 1.0/24.0;

  T b1, b2, Radius1, Radius2;

  // Handle combination: Real with complex conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff)  + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.5 * xy[j]/(xy[i_min] * Radius1);
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.5 * xy[j]/(xy[i_min] * Radius1);
  }

  // Now: Complex Conjugated <-> Complex Conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.25 / Radius1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g[NumEigVals+1] -= 0.25 / Radius1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+1] -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T ThirdOrder(const std::vector<T>& xy, const int NumRoots,
             const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange)
{
  T g = 1.0/24.0;

  T b1, b2, Radius1, Radius2;

  for(size_t j = 0; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g -= 0.25 / Radius1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }

  return g;
}

// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T ThirdOrder(const std::vector<T>& xy, const int NumRoots,
             const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange,
             const std::vector<PT>& ImagDiff_over_RealDiff, const size_t i_min)
{
  T g = 1.0/24.0;

  T b1, b2, Radius1, Radius2;

  // Handle combination: Real with complex conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g -= 0.5 * xy[j]/(xy[i_min] * Radius1);
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g -= 0.5 * xy[j]/(xy[i_min] * Radius1);
  }

  // Now: Complex Conjugated <-> Complex Conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g -= 0.25 / Radius1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g -= xy[j] * xy[k] / (Radius1 * Radius2);
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    g -= 0.25 / Radius1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }

  return g;
}

#endif