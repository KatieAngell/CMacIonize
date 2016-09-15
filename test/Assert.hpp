/**
 * @file Assert.hpp
 *
 * @brief Custom assert macros
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ASSERT_HPP
#define ASSERT_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>

#define assert_values_equal(a, b)                                              \
  if (abs(a - b) > 1.e-4 && abs(a - b) > 1.e-4 * abs(a + b)) {                 \
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;                     \
    std::cerr << "Assertion failed: " << a << " != " << b << std::endl;        \
    abort();                                                                   \
  }

#endif // ASSERT_HPP