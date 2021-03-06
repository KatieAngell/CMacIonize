/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file testWeightedSpectrumTracker.cpp
 *
 * @brief Unit test for the WeightedSpectrumTracker class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "RandomGenerator.hpp"
#include "WeightedSpectrumTracker.hpp"

/**
 * @brief Unit test for the WeightedSpectrumTracker class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// direction along main axis
  {
    const CoordinateVector<> n(1., 0., 0.);
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    assert_condition(area == 1.);
  }
  {
    const CoordinateVector<> n(0., 1., 0.);
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    assert_condition(area == 1.);
  }
  {
    const CoordinateVector<> n(0., 0., 1.);
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    assert_condition(area == 1.);
  }

  {
    const CoordinateVector<> n(-1., 0., 0.);
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    assert_condition(area == 1.);
  }
  {
    const CoordinateVector<> n(0., -1., 0.);
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    assert_condition(area == 1.);
  }
  {
    const CoordinateVector<> n(0., 0., -1.);
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    assert_condition(area == 1.);
  }

  /// direction along planar diagonal
  {
    CoordinateVector<> n(1., 1., 0.);
    n /= n.norm();
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    assert_condition(area == std::sqrt(2.));
  }

  /// direction along cube diagonal
  {
    CoordinateVector<> n(1., 1., 1.);
    n /= n.norm();
    const double area = WeightedSpectrumTracker::get_projected_area(n);
    cmac_warning("area: %g", area);
    assert_values_equal_rel(area, std::sqrt(3), 1.e-16);
  }

  /// a bunch of random directions to make sure the result is always valid
  {
    RandomGenerator rg(42);
    for (uint_fast32_t i = 0; i < 1e6; ++i) {
      // draw two pseudo random numbers
      const double cost = 2. * rg.get_uniform_random_double() - 1.;
      const double phi = 2. * M_PI * rg.get_uniform_random_double();

      // now use them to get all directional angles
      const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
      const double cosp = std::cos(phi);
      const double sinp = std::sin(phi);

      // set the direction...
      CoordinateVector<> direction(sint * cosp, sint * sinp, cost);
      direction /= direction.norm();
      const double area =
          WeightedSpectrumTracker::get_projected_area(direction);
      assert_condition(area > 0.);
      assert_condition(area == area);
    }
  }

  return 0;
}
