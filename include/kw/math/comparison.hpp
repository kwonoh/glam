// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef KW_MATH_COMPARISON_HPP
#define KW_MATH_COMPARISON_HPP

#include <boost/math/special_functions/relative_difference.hpp>

namespace kw {

template <typename T>
bool
epsilon_equal(T const& a, T const& b) {
  return boost::math::epsilon_difference(a, b) <= T(1.0);
}

}  // namespace kw

#endif  // KW_MATH_COMPARISON_HPP
