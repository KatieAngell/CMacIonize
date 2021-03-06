/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testDensitySubGrid.cpp
 *
 * @brief Unit test for the DensitySubGrid class.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */

#include "Assert.hpp"
#include "DensitySubGrid.hpp"
#include "IonizationStateCalculator.hpp"
#include "RandomGenerator.hpp"

#include <fstream>

/**
 * @brief Unit test for the DensitySubGrid class.
 *
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  const double box[6] = {-1.543e17, -1.543e17, -1.543e17,
                         3.086e17,  3.086e17,  3.086e17};
  const CoordinateVector< int_fast32_t > ncell(16, 16, 16);
  DensitySubGrid grid(box, ncell);
  RandomGenerator random_generator(42);

  for (auto cellit = grid.begin(); cellit != grid.end(); ++cellit) {
    cellit.get_ionization_variables().set_number_density(1.e8);
    cellit.get_ionization_variables().set_ionic_fraction(ION_H_n, 1.e-6);
  }

  for (uint_fast32_t iloop = 0; iloop < 10; ++iloop) {
    for (uint_fast32_t i = 0; i < 1e5; ++i) {
      PhotonPacket photon;

      photon.set_energy(3.288e15);
      for (int_fast32_t i = 0; i < NUMBER_OF_IONNAMES; ++i) {
        photon.set_photoionization_cross_section(i, 0.);
      }

      const double cost =
          2. * random_generator.get_uniform_random_double() - 1.;
      const double phi =
          2. * M_PI * random_generator.get_uniform_random_double();
      const double sint = std::sqrt(std::max(1. - cost * cost, 0.));
      const double cosp = std::cos(phi);
      const double sinp = std::sin(phi);
      const CoordinateVector<> d(sint * cosp, sint * sinp, cost);

      const double tau =
          -std::log(random_generator.get_uniform_random_double());

      photon.set_position(CoordinateVector<>(0.));
      photon.set_direction(d);
      photon.set_photoionization_cross_section(ION_H_n, 6.3e-22);
      photon.set_weight(1.);
      photon.set_target_optical_depth(tau);

      grid.interact(photon, TRAVELDIRECTION_INSIDE);
    }

    for (auto cellit = grid.begin(); cellit != grid.end(); ++cellit) {
      const double jH =
          (1.e49 / (1.e5 * cellit.get_volume())) *
          cellit.get_ionization_variables().get_mean_intensity(ION_H_n);
      const double alphaH = 4.e-19;
      const double nH = cellit.get_ionization_variables().get_number_density();
      const double xH =
          IonizationStateCalculator::compute_ionization_state_hydrogen(alphaH,
                                                                       jH, nH);
      cellit.get_ionization_variables().set_ionic_fraction(ION_H_n, xH);
      cellit.get_ionization_variables().reset_mean_intensities();
    }
  }

  /// write a restart file
  {
    RestartWriter writer("test_densitysubgrid.restart");
    grid.write_restart_file(writer);
  }

  /// read the restart file and check that both grids are the same
  {
    RestartReader reader("test_densitysubgrid.restart");
    DensitySubGrid grid2(reader);
    assert_condition(grid.get_number_of_cells() == grid2.get_number_of_cells());
    auto it = grid.begin();
    auto it2 = grid2.begin();
    while (it != grid.end() && it2 != grid2.end()) {
      assert_condition(
          it.get_ionization_variables().get_ionic_fraction(ION_H_n) ==
          it2.get_ionization_variables().get_ionic_fraction(ION_H_n));
      ++it;
      ++it2;
    }
  }

  std::ofstream ofile("testDensitySubGrid_output.txt");
  ofile << "# x (m)\ty (m)\tz (m)\txH\tnH (m^-3)\n";
  grid.print_intensities(ofile);

  return 0;
}
