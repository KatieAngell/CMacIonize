/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file HomogeneousDensityFunction.hpp
 *
 * @brief Homogeneous DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PLUMMERDENSITYFUNCTION_HPP
#define PLUMMERDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "cmath"
/**
 * @brief DensityFunction that returns a constant value for all coordinates,
 * corresponding to a homogeneous density field.
 */
class PlummerDensityFunction : public DensityFunction {
private:
  /*! @brief Single density value for the entire box (in m^-3). */
  const double _density;

  /*! @brief Single temperature value for the entire box (in K). */
  const double _temperature;

  /*! @brief Initial hydrogen neutral fraction for the entire box. */
  const double _neutral_fraction_H;

public:
  /**
   * @brief Constructor.
   *
   * @param density Single density value for the entire box (in m^-3).
   * @param temperature Single temperature value for the entire box (in K).
   * @param neutral_fraction_H Single hydrogen neutral fraction value for the
   * entire box.
   * @param log Log to write logging information to.
   */
  PlummerDensityFunction(double density = 1., double temperature = 8000.,
                             double neutral_fraction_H = 1.e-6,
                             Log *log = nullptr)
      : _density(density), _temperature(temperature),
        _neutral_fraction_H(neutral_fraction_H) {

    if (log) {
      log->write_status(
          "Created HomogeneousDensityFunction with constant density ", _density,
          " m^-3 and constant temperature ", _temperature, " K.");
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - density: Constant number density value (default: 100. cm^-3)
   *  - temperature: Constant initial temperature value (default: 8000. K)
   *  - neutral fraction H: Constant initial neutral fraction value
   *    (default: 1.e-6)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   */
  PlummerDensityFunction(ParameterFile &params, Log *log = nullptr)
      : PlummerDensityFunction(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:density", "100. cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "8000. K"),
            params.get_value< double >("DensityFunction:neutral fraction H",
                                       1.e-6),
            log) {}

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {
    DensityValues values;

	double a = 1e13;
    double M = 2e30;

    auto midpoint = cell.get_cell_midpoint();
    auto radius = midpoint.norm();
    if (radius == 0) {
      // Take density 10% of cell width from r=0
      radius = pow(cell.get_volume(), 1.0 / 3) / 10;
    }
    auto plummer =        
        3*M/(4*M_PI*pow(a,3)) * pow(1 + pow(radius/a, 2), -2.5)/PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS);
    // auto hernquist = 8.0 / radius * pow(1 + radius, -3)*cell.get_volume();
    values.set_number_density(plummer);
    values.set_temperature(8000);
    return values;
  }
};

#endif // HOMOGENEOUSDENSITYFUNCTION_HPP
