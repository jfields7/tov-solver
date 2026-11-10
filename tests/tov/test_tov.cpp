//========================================================================================
// TOVSolver library
// Copyright(C) 2026 Jacob Fields <fields@ias.edu>
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file test_tov.cpp
//! \brief Tests for the TOV solver.

#include <iostream>
#include <cmath>

#include <ts_tov.hpp>
#include <ts_barotropic.hpp>
#include <ts_polytrope.hpp>

#include <testing.hpp>

bool TestBuild() {
  ts::TOV<> tov;

  if (tov.IsInitialized()) {
    std::cout << "TOV is unexpectedly initialized.\n";
    return false;
  }

  double rhoc = tov.GetCentralDensity();
  if (!ts_math::isnan(rhoc)) {
    std::cout << "Unexpectedly returned non-NaN density without being initialized.\n";
    return false;
  }

  double dr = tov.GetSpacing();
  if (!ts_math::isnan(dr)) {
    std::cout << "Unexpectedly returned non-NaN spacing without being initialized.\n";
    return false;
  }

  unsigned int npoints = tov.GetNPoints();
  if (npoints > 0) {
    std::cout << "Unexpectedly returned non-zero points without being initialized.\n";
    return false;
  }

  return true;
}

/*bool TestPolytrope() {
  ts::TOV<> tov;
  ts::BarotropicEOS<ts::Polytrope<>> eos{2.0, 100.0};
  
  const double rhoc = 1.28e-3;
  const unsigned int npoints = 10000;
  const double dr = 1e-3;

  tov.ComputeSolution(eos, rhoc, npoints, dr);
  return true;
}*/

int main(int argc, char *argv[]) {
  UnitTests tester{"EOS"};

  tester.RunTest(&TestBuild, "Build TOV object");

  //tester.RunTest(&TestPolytrope, "Solve TOV with polytropic EOS");

  tester.PrintSummary();

  if (tester.tests_failed.size() > 0) {
    return -1;
  }
  return 0;
}
