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

  unsigned int npoints = tov.GetNPoints();
  if (npoints > 0) {
    std::cout << "Unexpectedly returned non-zero points without being initialized.\n";
    return false;
  }

  double mass = tov.GetMass();

  if (!ts_math::isnan(mass)) {
    std::cout << "Unexpectedly returned non-NaN mass without being initialized.\n";
    return false;
  }

  double radius = tov.GetRadius();

  if (!ts_math::isnan(radius)) {
    std::cout << "Unexpectedly returned non-NaN radius without being initialized.\n";
    return false;
  }

  double iso_radius = tov.GetIsotropicRadius();

  if (!ts_math::isnan(iso_radius)) {
    std::cout <<
      "Unexpectedly returned non-NaN isotropic radius without being initialized.\n";
    return false;
  }

  return true;
}

bool TestPolytrope() {
  ts::TOV<> tov;
  ts::BarotropicEOS<ts::Polytrope<>> eos{2.0, 100.0};
  
  const double rhoc = 1.28e-3;
  const unsigned int npoints = 1280;

  tov.ComputeSolution(eos, rhoc, npoints);

  double mass = tov.GetMass();
  double rad = tov.GetRadius();
  double iso_rad = tov.GetIsotropicRadius();

  if (ts_math::isnan(mass) || ts_math::isnan(rad) || ts_math::isnan(iso_rad)) {
    std::cout << "TOV has NaNs after solution.\n";
    return false;
  }

  double err = GetError(1.4, mass);

  if (err > 1.e-3) {
    std::cout << "TOV solution did not return the correct mass.\n";
    PrintError(1.4, mass);

    return false;
  }

  err = GetError(9.59, rad);

  if (err > 1e-3) {
    std::cout << "TOV solution did not return the correct radius.\n";
    PrintError(9.59, rad);
    
    return false;
  }

  err = GetError(8.1335, iso_rad);

  if (err > 1e-3) {
    std::cout <<"TOV solution did not return the correct isotropic radius.\n";
    PrintError(8.1335, iso_rad);

    return false;
  }

  return true;
}

int main(int argc, char *argv[]) {
  UnitTests tester{"EOS"};

  tester.RunTest(&TestBuild, "Build TOV object");

  tester.RunTest(&TestPolytrope, "Solve TOV with polytropic EOS");

  tester.PrintSummary();

  if (tester.tests_failed.size() > 0) {
    return -1;
  }
  return 0;
}
