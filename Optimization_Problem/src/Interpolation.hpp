// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Daniel Doehring                  RWTH Aachen University 2022-11-20

#ifndef __INTERPOLATION_HPP__
#define __INTERPOLATION_HPP__

#include <vector>
#include <algorithm>

template<typename T>
inline T Lin_IntPol(const T Real, const std::vector<T>& RealRange, const std::vector<T>& ImagRange)
{
  if(Real <= RealRange[0]){ // Catch case for which interpolation doesn't make sense
    return ImagRange[0];
  }
  else if(Real > 0.) {
    // Alternative: Mirror the imaginary part by flipping the sign of the real part
    // Needed for fourth order stability polynomial
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), -Real) - RealRange.begin();

    return ImagRange[i-1] + (ImagRange[i] - ImagRange[i-1]) / (RealRange[i] - RealRange[i-1]) *
                            (-Real - RealRange[i-1]);
  }
  else {
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), Real) - RealRange.begin();

    return ImagRange[i-1] + (ImagRange[i] - ImagRange[i-1]) / (RealRange[i] - RealRange[i-1]) *
                            (Real - RealRange[i-1]) ;
  }
}

// T: for dco types PT: Passive Type: For usual real types
template<typename T, typename PT>
inline T Lin_IntPol(const T& Real, const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange)
{
  if(Real <= RealRange[0]){ // Catch case for which interpolation doesn't make sense
    return ImagRange[0];
  }
  else if(Real > 0.) {
    // Alternative: Mirror the imaginary part by flipping the sign of the real part
    // Needed for fourth order stability polynomial
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), -Real) - RealRange.begin();

    return ImagRange[i-1] + (ImagRange[i] - ImagRange[i-1]) / (RealRange[i] - RealRange[i-1]) *
                            (-Real - RealRange[i-1]);
  }
  else {
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), Real) - RealRange.begin();
    
    return ImagRange[i-1] + (ImagRange[i] - ImagRange[i-1]) / (RealRange[i] - RealRange[i-1]) *
                            (Real - RealRange[i-1]) ;
  }
}

/// Pre-computed (I[i] - I[i-1]) / (R[i] - R[i-1])

template<typename T>
inline T Lin_IntPol(const T Real, const std::vector<T>& RealRange, const std::vector<T>& ImagRange,
                    const std::vector<T>& ImagDiff_over_RealDiff)
{
  if(Real <= RealRange[0]){ // Catch case for which interpolation doesn't make sense
    return ImagRange[0];
  }
  else if(Real > 0.) {  
    // Alternative: Mirror the imaginary part by flipping the sign of the real part
    // Needed for fourth order stability polynomial
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), -Real) - RealRange.begin();

    return ImagRange[i-1] + (-Real - RealRange[i-1]) * ImagDiff_over_RealDiff[i-1];
  }
  else {
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), Real) - RealRange.begin();

    return ImagRange[i-1] + (Real - RealRange[i-1]) * ImagDiff_over_RealDiff[i-1];
  }
}

// T: for dco types PT: Passive Type: For usual real types
template<typename T, typename PT>
inline T Lin_IntPol(const T& Real, const std::vector<PT>& RealRange, const std::vector<PT>& ImagRange,
                    const std::vector<PT>& ImagDiff_over_RealDiff)
{
  if(Real <= RealRange[0]){ // Catch case for which interpolation doesn't make sense
    return ImagRange[0];
  }
  else if(Real > 0.) {  
    // Alternative: Mirror the imaginary part by flipping the sign of the real part
    // Needed for fourth order stability polynomial
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), -Real) - RealRange.begin();

    return ImagRange[i-1] + (-Real - RealRange[i-1]) * ImagDiff_over_RealDiff[i-1];
  }
  else {
    const size_t i = std::upper_bound(RealRange.begin(), RealRange.end(), Real) - RealRange.begin();
    
    return ImagRange[i-1] + (Real - RealRange[i-1]) * ImagDiff_over_RealDiff[i-1];
  }
}

#endif
