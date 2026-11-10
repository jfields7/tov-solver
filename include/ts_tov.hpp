#ifndef TS_TOV_HPP_
#define TS_TOV_HPP_
//========================================================================================
// TOVSolver library
// Copyright(C) 2026 Jacob Fields <fields@ias.edu>
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ts_tov.hpp
//! \brief The TOV solution and solver.

#include <ts_types.hpp>
#include <ts_barotropic.hpp>

namespace ts {

template<class T = double>
class TOV {
 private:
  bool initialized;
  T rhoc;
  unsigned int npoints;
  T dr;
 public:
  TOV() : initialized(false), rhoc(ts_math::quiet_NaN<T>), npoints(0),
          dr(ts_math::quiet_NaN<T>) {
  }

  template<class EOSPolicy>
  void ComputeSolution(BarotropicEOS<EOSPolicy>& eos, T rhoc_, unsigned int npoints_,
                       T dr_) {
    rhoc = rhoc_;
    npoints = npoints_;
  }

  TS_INLINE
  bool IsInitialized() const {
    return initialized;
  }

  TS_INLINE
  T GetCentralDensity() const {
    return rhoc;
  }

  TS_INLINE
  T GetSpacing() const {
    return dr;
  }

  TS_INLINE
  unsigned int GetNPoints() const {
    return npoints;
  }
};

}; // namespace ts

#endif // TS_TOV_HPP_
