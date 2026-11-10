//========================================================================================
// TOVSolver library
// Copyright(C) 2026 Jacob Fields <fields@ias.edu>
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file test_polytrope.cpp
//! \brief Unit tests for the Polytrope class
#include <iostream>

#include <ts_barotropic.hpp>
#include <ts_polytrope.hpp>

#include <testing.hpp>

bool TestBuild() {
  ts::BarotropicEOS<ts::Polytrope<>> eos{2.0, 100.0};

  if (eos.GetGamma() != 2.0) {
    return false;
  }

  if (eos.GetKappa() != 100.0) {
    return false;
  }

  return true;
}

bool TestPressure() {
  ts::BarotropicEOS<ts::Polytrope<>> eos{2.0, 100.0};

  // First compute this all for density
  double rho = 1.28e-3;

  double p = eos.GetPressure(rho);

  double rho_calc = eos.GetRhoFromPressure(p);

  double err = GetError(rho, rho_calc);

  if (err > 1e-15) {
    std::cout << "Pressure/density relationship is not invertible.\n";
    PrintError(rho, rho_calc);
    return false;
  }

  // Now perform the same calculation for pseudo-enthalpy.
  double lnh = 0.1;

  p = eos.GetPressureFromLogEnthalpy(lnh);

  double lnh_calc = eos.GetLogEnthalpyFromPressure(p);

  err = GetError(lnh, lnh_calc);

  if (err > 1e-15) {
    std::cout << "Pressure/enthalpy relationship is not invertible.\n";
    PrintError(lnh, lnh_calc);
    return false;
  }

  return true;
}

bool TestEnergy() {
  ts::BarotropicEOS<ts::Polytrope<>> eos{2.0, 100.0};

  // First compute this all for density
  double rho = 1.28e-3;

  double e = eos.GetEnergy(rho);

  double rho_calc = eos.GetRhoFromEnergy(e);

  double err = GetError(rho, rho_calc);

  if (err > 1e-15) {
    std::cout << "Energy/density relationship is not invertible.\n";
    PrintError(rho, rho_calc);
    return false;
  }

  // Now compute this all for the pseudo-enthalpy.
  double lnh = 0.1;

  e = eos.GetEnergyFromLogEnthalpy(lnh);

  double lnh_calc = eos.GetLogEnthalpyFromEnergy(e);

  err = GetError(lnh, lnh_calc);

  if (err > 1e-15) {
    std::cout << "Energy/enthalpy relationship is not invertible.\n";
    PrintError(lnh, lnh_calc);
    return false;
  }

  return true;
}

bool TestEnthalpy() {
  ts::BarotropicEOS<ts::Polytrope<>> eos{2.0, 100.0};

  double rho = 1.28e-3;

  double lnh = eos.GetLogEnthalpy(rho);

  double rho_calc = eos.GetRhoFromLogEnthalpy(lnh);

  double err = GetError(rho, rho_calc);

  if (err > 1e-15) {
    std::cout << "Log pseudo-enthalpy/density relationship is not invertible.\n";
    PrintError(rho, rho_calc);
    return false;
  }

  return true;
}

int main(int argc, char* argv[]) {
  UnitTests tester{"Polytrope"};

  tester.RunTest(&TestBuild, "Build test");

  tester.RunTest(&TestPressure, "Pressure test");

  tester.RunTest(&TestEnergy, "Energy test");

  tester.RunTest(&TestEnthalpy, "Log pseudo-enthalpy test");

  tester.PrintSummary();

  if (tester.tests_failed.size() > 0) {
    return -1;
  }

  return 0;
}
