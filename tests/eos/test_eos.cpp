//========================================================================================
// TOVSolver library
// Copyright(C) 2026 Jacob Fields <fields@ias.edu>
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file test_eos.cpp
//! \brief Tests for constructing an EOS template. Note that these tests do not check the
//  validity of the EOS; they check that the forwarding functions work correctly and
//  enforce a strict interface.

#include <iostream>
#include <cmath>

#include <ts_barotropic.hpp>

#include <testing.hpp>

class DummyEOS {
 private:
 protected:
  // Forward functions
  double Pressure(double rho) {
    return rho;
  }
  double Energy(double rho) {
    return 2.0*rho;
  }
  double LogEnthalpy(double rho) {
    // Note that this is not strictly physical! This is only valid for testing.
    return std::log(3.0*rho);
  }

  // Inverse functions
  double RhoFromPressure(double P) {
    return P;
  }
  double RhoFromEnergy(double e) {
    return e/2.0;
  }
  double RhoFromLogEnthalpy(double lnh) {
    return std::exp(lnh)/3.0;
  }

  // Enthalpy functions
  double PressureFromLogEnthalpy(double lnh) {
    return std::exp(lnh)/3.0;
  }
  double EnergyFromLogEnthalpy(double lnh) {
    return 2.0*std::exp(lnh)/3.0;
  }

  // Enthalpy inverses
  double LogEnthalpyFromPressure(double P) {
    return std::log(3.0*P);
  }
  double LogEnthalpyFromEnergy(double e) {
    return std::log(3.0*e/2.0);
  }
 public:
};

bool TestBuildEOS() {
  [[maybe_unused]]
  ts::BarotropicEOS<DummyEOS> eos{};

  return true;
}

bool TestForwardFunctions() {
  ts::BarotropicEOS<DummyEOS> eos{};

  double rho = 1.0;

  double P = eos.GetPressure(rho);

  if (P != rho) {
    return false;
  }

  double e = eos.GetEnergy(rho);

  if (e != 2.0*rho) {
    return false;
  }

  double lnh = eos.GetLogEnthalpy(rho);

  if (lnh != ts_math::log(3.0*rho)) {
    return false;
  }

  return true;
}

bool TestInverseFunctions() {
  ts::BarotropicEOS<DummyEOS> eos{};

  double P = 1.0;

  double rho = eos.GetRhoFromPressure(P);

  if (rho != P) {
    return false;
  }

  double e = 2.0;

  rho = eos.GetRhoFromEnergy(e);

  if (rho != e/2.0) {
    return false;
  }

  double lnh = 3.5;
  
  rho = eos.GetRhoFromLogEnthalpy(lnh);

  if (rho != ts_math::exp(lnh)/3.0) {
    return false;
  }

  return true;
}

bool TestEnthalpyFunctions() {
  ts::BarotropicEOS<DummyEOS> eos{};

  double lnh = 1.0;

  double P = eos.GetPressureFromLogEnthalpy(lnh);

  if (P != ts_math::exp(lnh)/3.0) {
    return false;
  }

  double e = eos.GetEnergyFromLogEnthalpy(lnh);

  if (e != 2.0*ts_math::exp(lnh)/3.0) {
    return false;
  }

  lnh = eos.GetLogEnthalpyFromPressure(P);

  if (lnh != ts_math::log(3.0*P)) {
    return false;
  }

  lnh = eos.GetLogEnthalpyFromEnergy(e);

  if (lnh != ts_math::log(3.0*e/2.0)) {
    return false;
  }

  return true;
}

int main(int argc, char *argv[]) {
  UnitTests tester{"EOS"};

  tester.RunTest(&TestBuildEOS, "Construct new EOS");

  tester.RunTest(&TestForwardFunctions, "Compute variables from rest-mass density");

  tester.RunTest(&TestInverseFunctions, "Compute rest-mass density from other variables");

  tester.RunTest(&TestEnthalpyFunctions,
                 "Compute variables from log specific pseudo-enthalpy");

  tester.PrintSummary();

  if (tester.tests_failed.size() > 0) {
    return -1;
  }

  return 0;
}
