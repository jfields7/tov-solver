#ifndef TS_EOS_HPP_
#define TS_EOS_HPP_
//========================================================================================
// TOVSolver library
// Copyright(C) 2026 Jacob Fields <fields@ias.edu>
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ts_barotropic.hpp
//! \brief A policy template for handling a barotropic equation of state.

#include <ts_types.hpp>

namespace ts {

#define FORWARD_FUNC(a, b) \
  using arg = decltype((a)); \
  static_assert(std::is_floating_point_v<std::decay_t<arg>>, \
    #b " must accept a floating-point type."); \
  using res = decltype((b((a)))); \
  static_assert(std::is_floating_point_v<std::decay_t<res>>, \
    #b " must return a floating-point type."); \
  return b(a);


template<class EOSPolicy>
class BarotropicEOS : public EOSPolicy {
 private:
 public:
  // Forward functions

  //! \func GetPressure
  //! \brief Compute the pressure as a function of rest-mass density
  template<class T>
  TS_INLINE
  decltype(auto) GetPressure(T rho) {
    FORWARD_FUNC(rho, EOSPolicy::Pressure);
  }

  //! \func GetEnergy
  //! \brief Compute the total energy density as a function of rest-mass density
  template<class T>
  TS_INLINE
  decltype(auto) GetEnergy(T rho) {
    FORWARD_FUNC(rho, EOSPolicy::Energy);
  }

  //! \func GetLogEnthalpy
  //! \brief Compute the log specific pseudo-enthalpy as a function of rest-mass density
  template<class T>
  TS_INLINE
  decltype(auto) GetLogEnthalpy(T rho) {
    FORWARD_FUNC(rho, EOSPolicy::LogEnthalpy);
  }

  // Inverse functions

  //! \func GetRhoFromPressure
  //! \brief Compute the rest-mass density as a function of pressure
  template<class T>
  TS_INLINE
  decltype(auto) GetRhoFromPressure(T P) {
    FORWARD_FUNC(P, EOSPolicy::RhoFromPressure);
  }

  //! \func GetRhoFromEnergy
  //! \brief Compute the rest-mass density as a function of energy
  template<class T>
  TS_INLINE
  decltype(auto) GetRhoFromEnergy(T e) {
    FORWARD_FUNC(e, EOSPolicy::RhoFromEnergy);
  }

  //! \func GetRhoFromLogEnthalpy
  //! \brief Compute the rest-mass density as a function of log specific pseudo-enthalpy.
  template<class T>
  TS_INLINE
  decltype(auto) GetRhoFromLogEnthalpy(T lnh) {
    FORWARD_FUNC(lnh, EOSPolicy::RhoFromLogEnthalpy);
  }

  // Enthalpy functions

  //! \func GetPressureFromLogEnthalpy
  //! \brief Compute the pressure as a function of log specific pseudo-enthalpy.
  template<class T>
  TS_INLINE
  decltype(auto) GetPressureFromLogEnthalpy(T lnh) {
    FORWARD_FUNC(lnh, EOSPolicy::PressureFromLogEnthalpy);
  }
};

#undef FORWARD_FUNC

}; // namespace ts

#endif // TS_EOS_HPP_
