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

bool TestPressure() {
  ts::BarotropicEOS<ts::Polytrope> eos{};

  return true;
}

int main(int argc, char* argv[]) {
  UnitTests tester{"Polytrope"};

  tester.RunTest(&TestPressure, "Pressure test");

  tester.PrintSummary();

  if (tester.tests_failed.size() > 0) {
    return -1;
  }

  return 0;
}
