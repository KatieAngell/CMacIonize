/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testVoronoiDensityGrid.cpp
 *
 * @brief Unit test for the VoronoiDensityGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "Assert.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "VoronoiDensityGrid.hpp"
#include "VoronoiDensityGridPositions.hpp"

/**
 * @brief Test implementation of VoronoiDensityGridPositions.
 */
class TestVoronoiDensityGridPositions : public VoronoiDensityGridPositions {
public:
  /**
   * @brief Get the number of positions that is generated by this object.
   *
   * @return 100.
   */
  virtual unsigned int get_number_of_positions() const { return 100; }

  /**
   * @brief Get a uniform random position.
   *
   * @return Uniform random position (in m).
   */
  virtual CoordinateVector<> get_position() const {
    return Utilities::random_position();
  }
};

/**
 * @brief Unit test for the VoronoiDensityGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  TestVoronoiDensityGridPositions test_positions;
  HomogeneousDensityFunction density_function(1., 2000.);
  Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  VoronoiDensityGrid grid(test_positions, density_function, box, false, false,
                          nullptr);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  assert_values_equal(1., grid.get_total_hydrogen_number());
  assert_values_equal(2000., grid.get_average_temperature());

  return 0;
}
