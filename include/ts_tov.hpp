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
  T lnhc;

  T mass;
  T radius;
  T iso_radius;

  void SetToNaN() {
    initialized = false;
    rhoc = ts_math::quiet_NaN<T>;
    npoints = 0;

    mass = ts_math::quiet_NaN<T>;
    radius = ts_math::quiet_NaN<T>;
    iso_radius = ts_math::quiet_NaN<T>;
  }
 public:
  TOV() {
    SetToNaN();
  }

  template<class EOSPolicy>
  void ComputeSolution(BarotropicEOS<EOSPolicy>& eos, T rhoc_, unsigned int npoints_) {
    rhoc = rhoc_;
    npoints = npoints_;

    // Compute the central pseudo-enthalpy, which we use as the starting value for our
    // integration.
    lnhc = eos.GetLogEnthalpy(rhoc);

    // Initial conditions
    T u = T(0.);
    T v = T(0.);
    T q = T(0.);
    // For simplicity, compute the solution with an explicit midpoint method. We'll fix
    // this later.
    T dlnh = lnhc/(npoints-1);
    T one_half = T(1.0/2.0);
    for (unsigned int i = 0; i < npoints; i++) {
      T lnh = lnhc - i*dlnh;
      // First stage
      T du, dv, dq, u1, v1, q1;
      ComputeRHS(eos, lnh, u, v, q, du, dv, dq);

      u1 = u - du*one_half*dlnh;
      v1 = v - dv*one_half*dlnh;
      q1 = q - dq*one_half*dlnh;

      // Second stage
      lnh = lnhc - (i+one_half)*dlnh;

      ComputeRHS(eos, lnh, u1, v1, q1, du, dv, dq);

      u = u - du*dlnh;
      v = v - dv*dlnh;
      q = q - dq*dlnh;
    }

    if (ts_math::isnan(u) || ts_math::isnan(v) || u < T(0) || v < T(0)) {
      SetToNaN();
      return;
    }

    radius = ts_math::sqrt(u);
    mass = v*radius;
    // We need to match the isotropic radius to the Schwarzschild solution.
    iso_radius = 0.5*(ts_math::sqrt(radius*(radius - 2*mass)) + radius - mass);
    initialized = true;
  }

  template<class EOSPolicy>
  void ComputeRHS(BarotropicEOS<EOSPolicy>& eos, T lnh, T u, T v, T q, T& du, T& dv,
                  T& dq) const {
    T P = eos.GetPressureFromLogEnthalpy(lnh);
    T e = eos.GetEnergyFromLogEnthalpy(lnh);

    // A regularized solution to use at the center.
    if (lnh == lnhc) {
      T D = e + 3*P;
      du = -T(3.0)/(2*ts_math::pi_v<T>*D);
      dv = -2*e/D;
      dq = - e/D;
    } else {
      du = -T(2.0)*(u*(T(1.) - 2*v))/(v + 4*ts_math::pi_v<T>*u*P);
      dv = (2*ts_math::pi_v<T>*e - v/(2*u))*du;
      dq = ((T(1.0) - 2*v) - ts_math::sqrt(T(1.0) - 2*v))/(v + 4*ts_math::pi_v<T>*u*P);
    }
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
  unsigned int GetNPoints() const {
    return npoints;
  }

  TS_INLINE
  T GetMass() const {
    return mass;
  }

  TS_INLINE
  T GetRadius() const {
    return radius;
  }

  TS_INLINE
  T GetIsotropicRadius() const {
    return iso_radius;
  }
};

}; // namespace ts

#endif // TS_TOV_HPP_
