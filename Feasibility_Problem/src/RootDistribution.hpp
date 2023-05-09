// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Daniel Doehring                  RWTH Aachen University 2022-11-20

#ifndef __ROOTDISTRIBUTION_HPP__
#define __ROOTDISTRIBUTION_HPP__

#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

template <typename T>
inline T Dist(const T Re1, const T Im1, const T Re2, const T Im2) {
  return sqrt((Re1 - Re2) * (Re1 - Re2) + (Im1 - Im2) * (Im1 - Im2));
}

template <typename T>
std::vector<T> InitialRootDistr(const int NumEigVals, const int NumUnknowns, const int NumStages,
                                const std::vector<T>& RealEigValsScaled, const std::vector<T>& ImagEigValsScaled) {
  
  std::vector<T> ArcLengths(NumEigVals);
  ArcLengths[0] = 0.;
  for(size_t i = 1; i < NumEigVals; i++)
    ArcLengths[i] = ArcLengths[i-1] + Dist(RealEigValsScaled[i],   ImagEigValsScaled[i],
                                           RealEigValsScaled[i-1], ImagEigValsScaled[i-1]);

  // Reasonable approximation for initial distribution of eigenvalues
  // NOTE: Distribution with "additional root" in (0, 0) => Do NOT subtract 1 from NumUnknowns
  const T ArcLengthAVG = ArcLengths[NumEigVals - 1] / (NumUnknowns - 0);
  
  std::vector<T> x0(NumUnknowns);
  // CAVEAT: Hard-coded to even stability polynoms only!
  assert(NumStages % 2 == 0);
  x0[0] = RealEigValsScaled[0]; // This follows from the even stability polynom assumption
  T span;
  size_t ind;
  for(size_t i = 1; i < NumUnknowns; i++) {
    span = i*ArcLengthAVG;
    ind = std::upper_bound(ArcLengths.begin(), ArcLengths.end(), span) - ArcLengths.begin();

    x0[i] = RealEigValsScaled[ind-1] + (RealEigValsScaled[ind] - RealEigValsScaled[ind-1]) *
            (span - ArcLengths[ind - 1]) / (ArcLengths[ind] - ArcLengths[ind-1]);
  }

  return x0;
}

template <typename T>
std::vector<T> InitialRootDistr(const int NumHullPoints, const int NumUnknowns, const int NumStages,
                                const T RealEigValMin, const std::vector<T>& HullRealScaled, const std::vector<T>& HullImagScaled) {
  
  std::vector<T> ArcLengths(NumHullPoints);
  ArcLengths[0] = 0.;
  for(size_t i = 1; i < NumHullPoints; i++)
    ArcLengths[i] = ArcLengths[i-1] + Dist(HullRealScaled[i], HullImagScaled[i], 
                                           HullRealScaled[i-1], HullImagScaled[i-1]);

  // Reasonable approximation for initial distribution of eigenvalues
  // NOTE: Distribution with "additional root" in (0, 0) => Do NOT subtract 1 from NumUnknowns
  const T ArcLengthAVG = ArcLengths[NumHullPoints-1] / (NumUnknowns - 0);
  
  std::vector<T> x0(NumUnknowns);
  // CAVEAT: Hard-coded to even stability polynoms only!
  assert(NumStages % 2 == 0);
  x0[0] = RealEigValMin; // This follows from the even stability polynom assumption
  T span;
  size_t ind;
  for(size_t i = 1; i < NumUnknowns; i++) {
    span = i*ArcLengthAVG;
    ind = std::upper_bound(ArcLengths.begin(), ArcLengths.end(), span) - ArcLengths.begin();

    x0[i] = HullRealScaled[ind-1] + (HullRealScaled[ind] - HullRealScaled[ind-1]) *
            (span - ArcLengths[ind - 1]) / (ArcLengths[ind] - ArcLengths[ind-1]);
  }

  return x0;
}

template <typename T>
std::vector<T> InitialRootDistr(const int NumCurvePoints, const int NumUnknowns, const int NumStages,
                                const std::vector<T>& CurveRealScaled, const std::vector<T>& CurveImagScaled,
                                const int NumPEHalfStages, const std::vector<T>& Real_PE_HalfStagesScaled) {
  
  // CAVEAT: Hard-coded to even stability polynoms only!
  assert(NumStages % 2 == 0);                                  

  std::vector<T> ArcLengths(NumCurvePoints);
  ArcLengths[0] = 0.;
  for(size_t i = 1; i < NumCurvePoints; i++)
    ArcLengths[i] = ArcLengths[i-1] + Dist(CurveRealScaled[i], CurveImagScaled[i], 
                                           CurveRealScaled[i-1], CurveImagScaled[i-1]);

  std::vector<T> ArcLengthsPEHalfStages(NumPEHalfStages);
  size_t ind;
  ArcLengthsPEHalfStages[0] = 0.; // First value is pure real, left end of spectrum
  for(size_t i = 1; i < NumPEHalfStages; i++) {
    ind = std::upper_bound(CurveRealScaled.begin(), CurveRealScaled.end(), Real_PE_HalfStagesScaled[i]) 
         - CurveRealScaled.begin();

    ArcLengthsPEHalfStages[i] = ArcLengths[ind-1] + (ArcLengths[ind] - ArcLengths[ind-1]) *
                                 (Real_PE_HalfStagesScaled[i] - CurveRealScaled[ind - 1]) 
                               / (CurveRealScaled[ind] - CurveRealScaled[ind-1]);
  }

  std::vector<T> RealPartsPEIntermediates(NumPEHalfStages - 1);
  T ArcLenghtIntermediate;
  for(size_t i = 0; i < NumPEHalfStages - 1; i++) {
    ArcLenghtIntermediate = ArcLengthsPEHalfStages[i] + 0.5 * (ArcLengthsPEHalfStages[i+1] - ArcLengthsPEHalfStages[i]);

    ind = std::upper_bound(ArcLengths.begin(), ArcLengths.end(), ArcLenghtIntermediate) - ArcLengths.begin();

    RealPartsPEIntermediates[i] = CurveRealScaled[ind-1] + (CurveRealScaled[ind] - CurveRealScaled[ind-1]) *
            (ArcLenghtIntermediate - ArcLengths[ind - 1]) / (ArcLengths[ind] - ArcLengths[ind-1]);
  }

  std::vector<T> x0(NumUnknowns);
  x0[0] = Real_PE_HalfStagesScaled[0];
  for(size_t i = 0; i < NumPEHalfStages - 1; i++) {
    x0[2*i+1] = RealPartsPEIntermediates[i];
    x0[2*i+2] = Real_PE_HalfStagesScaled[i+1];
  }
  
  T LastArcLength = ((NumUnknowns - 1)*ArcLengths[NumCurvePoints-1]) / NumUnknowns;
  ind = std::upper_bound(ArcLengths.begin(), ArcLengths.end(), LastArcLength) - ArcLengths.begin();

  x0[NumUnknowns - 1] = CurveRealScaled[ind-1] + (CurveRealScaled[ind] - CurveRealScaled[ind-1]) *
                        (LastArcLength - ArcLengths[ind - 1]) / (ArcLengths[ind] - ArcLengths[ind-1]);

  return x0;
}

#endif
