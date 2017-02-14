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
 * @file testHydroIntegrator.cpp
 *
 * @brief Unit test for the HydroIntegrator class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "CartesianDensityGrid.hpp"
#include "HomogeneousDensityFunction.hpp"
#include "HydroIntegrator.hpp"

/**
 * @brief Unit test for the HydroIntegrator class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {
  HydroIntegrator integrator(5. / 3.);

  Box box(CoordinateVector<>(0.), CoordinateVector<>(1.));
  CoordinateVector< int > ncell(100, 1, 1);
  HomogeneousDensityFunction density_function(1.);
  CoordinateVector< bool > periodic(true);
  CartesianDensityGrid grid(box, ncell, density_function, periodic, true);
  std::pair< unsigned long, unsigned long > block =
      std::make_pair(0, grid.get_number_of_cells());
  grid.initialize(block);

  integrator.do_hydro_step(grid, 0.1);

  return 0;
}
