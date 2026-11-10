#ifndef TS_POLYTROPE_HPP_
#define TS_POLYTROPE_HPP_
//========================================================================================
// TOVSolver library
// Copyright(C) 2026 Jacob Fields <fields@ias.edu>
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ts_polytrope.hpp
//! \brief EOS policy for a fixed polytropic EOS

#include <ts_types.hpp>

namespace ts {

//! \class Polytrope
//! \brief An EOS policy that represents the polytrope $P = K*\rho^\gamma$
template<class T = double>
class Polytrope {
 private:
  const T gamma;
  const T kappa;
  const T igamma;
  const T igammam1;
 protected:
  Polytrope(T gamma_, T kappa_) : gamma(gamma_), kappa(kappa_), igamma(1.0/gamma_),
                                  igammam1(1.0/(gamma_ - 1.0)) {

  }

  // Pressure functions
  TS_INLINE
  T Pressure(T rho) const {
    return kappa*ts_math::pow(rho, gamma);
  }
  TS_INLINE
  T RhoFromPressure(T P) const {
    return ts_math::pow(P/kappa, igamma);
  }
  TS_INLINE
  T PressureFromLogEnthalpy(T lnh) const {
    T h = ts_math::exp(lnh);

    return ts_math::pow(ts_math::pow(kappa,-igamma)*(h-1.)*(gamma - 1.)*igamma,
                        gamma*igammam1);
  }
  TS_INLINE
  T LogEnthalpyFromPressure(T P) const {
    T rho = RhoFromPressure(P);

    return ts_math::log(1.0 + gamma*igammam1*(P/rho));
  }

  // Energy functions
  TS_INLINE
  T Energy(T rho) const {
    return rho + Pressure(rho)*igammam1;
  }
  TS_INLINE
  T RhoFromEnergy(T e) const {
    // Note that this relationship is not generally analytically invertible. We therefore
    // employ a bracketed Newton-Raphson solve to find it.
    T ub = e;
    T lb = 0.0;
    T rho = 0.5*(lb + ub);
    const int N = 10;
    const T tol = ts_math::epsilon_v<T>;
    for (int i = 0; i < N; i++) {
      // Compute the candidate solution and update the bounds
      T f = e - (rho + kappa*ts_math::pow(rho, gamma)*igammam1);
      if (f > 0) {
        lb = rho;
      } else {
        ub = rho;
      }
      if (ts_math::fabs(f) < tol) {
        break;
      }
      T df = (gamma - 1.0) + (f - e)*gamma/rho;

      rho = rho - f/df;
      if (rho > ub || rho < lb) {
        // Update with bisection if the root is diverging.
        rho = 0.5*(lb + ub);
      }
    }
    return rho;
  }
  TS_INLINE
  T EnergyFromLogEnthalpy(T lnh) const {
    const T h = ts_math::exp(lnh);
    const T rho = ts_math::pow((h - 1.0)*(gamma - 1.0)*igamma/kappa, igammam1);
    return rho*h - Pressure(rho);
  }
  TS_INLINE
  T LogEnthalpyFromEnergy(T e) const {
    const T rho = RhoFromEnergy(e);
    return ts_math::log(1.0 + gamma*igammam1*kappa*ts_math::pow(rho, gamma - 1.0));
  }

  // Enthalpy functions
  TS_INLINE
  T LogEnthalpy(T rho) const {
    return ts_math::log(1.0 + gamma*igammam1*kappa*ts_math::pow(rho, gamma - 1.0));
  }
  TS_INLINE
  T RhoFromLogEnthalpy(T lnh) const {
    const T h = ts_math::exp(lnh);
    return ts_math::pow((h - 1.0)*(gamma - 1.0)*igamma/kappa, igammam1);
  }
 public:
  TS_INLINE
  T GetGamma() {
    return gamma;
  }

  TS_INLINE
  T GetKappa() {
    return kappa;
  }
};

}; // namespace ts

#endif // TS_POLYTROPE_HPP_
