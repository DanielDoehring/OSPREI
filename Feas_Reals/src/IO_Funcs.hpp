// This code is published under the Eclipse Public License.
// If you have not obtained a copy of the EPL, you can find the license (version 2.0) at
// https://www.eclipse.org/legal/epl-2.0/
//
// Authors:  Daniel Doehring                  RWTH Aachen University 2022-11-20

#ifndef __IOFUNCS_HPP__
#define __IOFUNCS_HPP__

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <sstream>

template <typename T>
void read_EigVals(const std::string EigValFileName, int& NumEigVals,
                  std::vector<T>& RealEigVals, std::vector<T>& ImagEigVals) {
  int i = 0;
  std::ifstream EigValsFile(EigValFileName);
  assert(EigValsFile);
  std::string line;
  while (std::getline(EigValsFile, line))
    i++;
  EigValsFile.close();

  NumEigVals = i; // This is the number cones (and also number of inequalities)
  std::cout << "Number of Eigenvalues is: " << NumEigVals << std::endl << std::endl;

  RealEigVals.resize(NumEigVals);
  ImagEigVals.resize(NumEigVals);
  T real, imag;

  char plus;
  i = 0;
  EigValsFile.open(EigValFileName);
  while (std::getline(EigValsFile, line)) {
    std::stringstream stream(line);
    assert(i < NumEigVals);
    while (stream >> real >> plus >> imag) {
      if(real > 0)
        std::cout << "CARE: Eigenvalue with positive real part, should be removed!" << std::endl;
      RealEigVals[i] = real;

      if(imag < 0)
        std::cout << "CARE: Eigenvalue with negative imag part, should be removed!" << std::endl;
      ImagEigVals[i] = imag;
    }
    i++;
  }
  EigValsFile.close();
}

template <typename T>
void read_x0(const std::string x0_FileName, std::vector<T>& x0, const int NumUnknowns) {
  int i = 0;
  std::ifstream x0_File(x0_FileName);
  assert(x0_File);
  std::string line;
  T initial_pos;

  x0.resize(NumUnknowns); // Fill x0 with zeros (have some values for imaginary corrections)
  while (std::getline(x0_File, line)) {
    std::stringstream stream(line);
    assert(i < NumUnknowns);
    while (stream >> initial_pos) {
      x0[i] = initial_pos;
    }
    i++;
  }
  x0_File.close();
}

template <typename T>
void read_Hull(const std::string Hull_FileName, std::vector<T>& hull) {
  std::ifstream Hull_File(Hull_FileName);
  assert(Hull_File);
  std::string line;
  T point;

  while (std::getline(Hull_File, line)) {
    std::stringstream stream(line);
    while (stream >> point) {
      hull.push_back(point);
    }
  }
  Hull_File.close();
}

template <typename T>
void read_PE(const std::string PEFileName, const int NumPE, std::vector<T>& Real_PE_HalfStages) {
  std::ifstream PEFile(PEFileName);
  assert(PEFile);

  Real_PE_HalfStages.resize(NumPE);
  T real, imag;
  char plus;
  int i = 0;
  std::string line;
  while (std::getline(PEFile, line)) {
    std::stringstream stream(line);
    assert(i < NumPE);
    while (stream >> real >> plus >> imag) {
      Real_PE_HalfStages[i] = real;
    }
    i++;
  }
  PEFile.close();
}

#endif
