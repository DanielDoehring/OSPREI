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
1/(r_i r_i*) + 1/(r_j r_j*) + 1/(r_i r_j) + 1/(r_i r_j*) + 1/(r_i* r_j) + 1/(r_i* r_j*)
= 1/Radius(r_i) + 1/Radius(r_j) + 4 * Re(r_i) * Re(r_j) / (Radius(r_i)^2 + Radius(r_j)^2)
=> 1.0/24.0 = sum_{i}^(S/2) 0.25/Radius(r_i) + sum_{j; i != j}^(S/2) Re(r_i) * Re(r_j) / (Radius(r_i)^2 + Radius(r_j)^2)

If r_i is real, this simplifies to 
1/(r_i r_j) + 1/(r_i r_j*) = 2 * Re(r_j) / ( r_i Radius(r_j)^2 ) 


For fourth order, we demand 1/24 = sum_{i,j,k i != j != k}^S 1/(r_i r_j r_k)
For complex conjugated roots we have
  1/(r_i r_i* r_j) + 1/(r_i r_i* r_j*) + 1/(r_i r_i* r_k) + 1/(r_i r_i* r_k*)
+ 1/(r_i r_j r_j*) + 1/(r_i r_j r_k) + 1/(r_i r_j r_k*) 
+ 1/(r_i r_j* r_k) + 1/(r_i r_j* r_k*) + 1/(r_i r_k r_k*)
+ 1/(r_i* r_j r_j*) + 1/(r_i* r_j r_k) + 1/(r_i* r_j r_k*
+ 1/(r_i* r_j* r_k) + 1/(r_i* r_j* r_k*) + 1/(r_i* r_k r_k*)

+ 1/(r_j r_j* r_k) + 1/(r_j r_j* r_k*) + 1/(r_j r_k r_k*) + 1/(r_j* r_k r_k*) 

The first 5 lines can be rewritten as
2 * [Re(r_i) * ( Re(r_j)^2 + 4Re(r_j)Re(r_k) + Im(r_j)^2 + Re(r_k)^2 + Im(r_k)^2) + 
     Re(r_k) * ( Re(r_j) * ( Re(r_j) + Re(r_k) ) + Im(r_j)^2 ) +
     Re(r_j) * Im(r_k)^2] 
/ ( Radius(r_i)^2 * Radius(r_j)^2 * Radius(r_k)^2 )

and this is what will be implemented, similar to the third order constraint.

If r_i is real, this simplifies to 
  1/(r_i r_j r_j*) + 1/(r_i r_j r_k) + 1/(r_i r_j r_k*)
+ 1/(r_i r_j* r_k) + 1/(r_i r_j* r_k*)
= [4 Re(r_j) Re(r_k) + Re(r_k^2) + Im(r_k)^2] 
  / [r_i Radius(r_j)^2 Radius(r_k)^2]
  
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

  // Single purely real root
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

  // Single purely real root
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

  // Single purely real root
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
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

    // 1/(r_j r_j*) = 1/Radius(r_j)
    g -= 0.25 / Radius1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g -= xy[j] * xy[k] / (Radius1 * Radius2);
    }
  }

  return g;
}


/// FOURTH ORDER ///


// For Odd Base Polynom => Even Lower Degree Polynomial
template <typename T>
void FourthOrder(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                 const std::vector<T>& RealRange, const std::vector<T>& ImagRange)
{
  g[NumEigVals+2] = 1.0/48.0;

  T b1, b2, b3, Radius1, Radius2, Radius3;

  for(size_t j = 0; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }
}

// For Even Base Polynomial => Odd Lower degree polynomial
template <typename T>
void FourthOrder(const T* xy, T* g, const int NumRoots, const int NumEigVals,
                 const std::vector<T>& RealRange, const std::vector<T>& ImagRange, 
                 const std::vector<T>& ImagDiff_over_RealDiff, const size_t i_min)
{
  g[NumEigVals+2] = 1.0/48.0;

  T b1, b2, b3, Radius1, Radius2, Radius3;

  // Handle combination: Real with complex conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff)  + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+2] += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff)  + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+2] += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+2] += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }
  }

  // Now: Complex Conjugated <-> Complex Conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = l+1; l < i_min; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }

      for(size_t l = i_min+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void FourthOrder(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                 const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange)
{
  g[NumEigVals+2] = 1.0/48.0;

  T b1, b2, b3, Radius1, Radius2, Radius3;

  for(size_t j = 0; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }
}

// For Even Base Polynom => Odd Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
void FourthOrder(const std::vector<T>& xy, std::vector<T>& g, const int NumRoots, const int NumEigVals,
                 const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange, 
                 const std::vector<PT>& ImagDiff_over_RealDiff, const size_t i_min)
{
  g[NumEigVals+2] = 1.0/48.0;

  T b1, b2, b3, Radius1, Radius2, Radius3;

  // Handle combination: Real with complex conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+2] += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+2] += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g[NumEigVals+2] += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }
  }

  // Now: Complex Conjugated <-> Complex Conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = l+1; l < i_min; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }

      for(size_t l = i_min+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g[NumEigVals+2] += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
                            xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
                            xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }
}

// For Odd Base Polynom => Even Lower Degree Polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T FourthOrder(const std::vector<T>& xy, const int NumRoots,
              const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange)
{
  T g = 1.0/48.0;

  T b1, b2, b3, Radius1, Radius2, Radius3;

  for(size_t j = 0; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
              xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
              xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }

  return g;
}

// For Even Base Polynomial => Odd Lower degree polynomial
template<typename T, typename PT> // T: for dco types PT: Passive Type: For usual real types
T FourthOrder(const std::vector<T>& xy, const int NumRoots,
              const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange, 
              const std::vector<PT>& ImagDiff_over_RealDiff, const size_t i_min)
{
  T g = 1.0/48.0;

  T b1, b2, b3, Radius1, Radius2, Radius3;

  // Handle combination: Real with complex conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      g += 0.5 * (4.0*xy[j]*xy[k] + xy[k]*xy[k] + b2*b2)/(xy[i_min] * Radius1 * Radius2);
    }
  }

  // Now: Complex Conjugated <-> Complex Conjugated
  for(size_t j = 0; j < i_min; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < i_min; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = l+1; l < i_min; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
              xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
              xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }

      for(size_t l = i_min+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
              xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
              xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }

    for(size_t k = i_min+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
              xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
              xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }

  for(size_t j = i_min + 1; j < NumRoots; j++) {
    b1 = Lin_IntPol(xy[j], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[j + NumRoots];
    Radius1 = xy[j]*xy[j] + b1*b1;

    for(size_t k = j+1; k < NumRoots; k++) {
      b2 = Lin_IntPol(xy[k], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[k + NumRoots];
      Radius2 = xy[k]*xy[k] + b2*b2;

      for(size_t l = k+1; l < NumRoots; l++) {
        b3 = Lin_IntPol(xy[l], RealRange, ImagRange, ImagDiff_over_RealDiff) + xy[l + NumRoots];
        Radius3 = xy[l]*xy[l] + b3*b3;

        g += (xy[j] * (b2*b2 + 4 * xy[j] * xy[k] + xy[l]*xy[l] + b3*b3) + 
              xy[l] * (xy[k] * (xy[k] + xy[l]) + b2*b2) +
              xy[k] * b3*b3) / (Radius1 * Radius2 * Radius3);
      }
    }
  }

  return g;
}

#endif