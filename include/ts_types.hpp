#ifndef TS_TYPES_HPP_
#define TS_TYPES_HPP_
//========================================================================================
// TOVSolver library
// Copyright(C) 2026 Jacob Fields <fields@ias.edu>
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ts_types.hpp
//! \brief A collection of types and macros used by TOVSolver

#if __has_include(<Kokkos_Core.hpp>)
#include <Kokkos_Core.hpp>
// Inline functions
#define TS_INLINE KOKKOS_INLINE_FUNCTION
// Math libraries
namespace ts_math = Kokkos;

#else
// Inline functions
#define TS_INLINE inline
// Math libraries
#include <cmath>
#include <limits>
namespace ts_math = std;

namespace std {

template<class T>
constexpr T epsilon_v = std::numeric_limits<T>::epsilon();

} // namespace ts_math

#endif

#endif // TS_TYPES_HPP_
