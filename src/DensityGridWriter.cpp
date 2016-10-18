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
 * @file DensityGridWriter.cpp
 *
 * @brief DensityGridWriter implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "DensityGridWriter.hpp"
#include "Box.hpp"
#include "CoordinateVector.hpp"
#include "DensityGrid.hpp"
#include "DensityValues.hpp"
#include "HDF5Tools.hpp"
#include <vector>

/**
 * @brief Constructor.
 *
 * @param name Name of the file to write.
 * @param grid DensityGrid containing the data to write.
 */
DensityGridWriter::DensityGridWriter(std::string name, DensityGrid &grid)
    : _name(name), _grid(grid) {
  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();
}

/**
 * @brief Write the file.
 */
void DensityGridWriter::write() {
  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(_name, HDF5Tools::HDF5FILEMODE_WRITE);

  // write header
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  Box box = _grid.get_box();
  CoordinateVector<> boxsize = box.get_sides() - box.get_anchor();
  HDF5Tools::write_attribute< CoordinateVector<> >(group, "BoxSize", boxsize);
  int dimension = 3;
  HDF5Tools::write_attribute< int >(group, "Dimension", dimension);
  std::vector< unsigned int > flag_entropy(6, 0);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "Flag_Entropy_ICs", flag_entropy);
  std::vector< double > masstable(6, 0.);
  HDF5Tools::write_attribute< std::vector< double > >(group, "MassTable",
                                                      masstable);
  int numfiles = 1;
  HDF5Tools::write_attribute< int >(group, "NumFilesPerSnapshot", numfiles);
  std::vector< unsigned int > numpart(6, 0);
  numpart[0] = _grid.get_number_of_cells();
  std::vector< unsigned int > numpart_high(6, 0);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_ThisFile", numpart);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_Total", numpart);
  HDF5Tools::write_attribute< std::vector< unsigned int > >(
      group, "NumPart_Total_HighWord", numpart_high);
  double time = 0.;
  HDF5Tools::write_attribute< double >(group, "Time", time);
  HDF5Tools::close_group(group);

  // write units
  group = HDF5Tools::create_group(file, "Units");
  double unit_value = 1;
  HDF5Tools::write_attribute< double >(group, "Unit current in cgs (U_I)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit length in cgs (U_L)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit mass in cgs (U_M)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit temperature in cgs (U_T)",
                                       unit_value);
  HDF5Tools::write_attribute< double >(group, "Unit time in cgs (U_t)",
                                       unit_value);
  HDF5Tools::close_group(group);

  // write particles
  group = HDF5Tools::create_group(file, "PartType0");
  std::vector< CoordinateVector<> > coords(numpart[0]);
  std::vector< double > ntot(numpart[0]);
  std::vector< double > nfracH(numpart[0]);
  std::vector< double > nfracHe(numpart[0]);
  std::vector< double > temperature(numpart[0]);
  unsigned int index = 0;
  for (auto it = _grid.begin(); it != _grid.end(); ++it) {
    Box cellbox = it.get_cell();
    DensityValues cellvals = it.get_values();
    coords[index] =
        cellbox.get_anchor() + 0.5 * cellbox.get_sides() - box.get_anchor();
    ntot[index] = cellvals.get_total_density();
    nfracH[index] = cellvals.get_neutral_fraction_H();
    nfracHe[index] = cellvals.get_neutral_fraction_He();
    temperature[index] = cellvals.get_temperature();
    ++index;
  }
  HDF5Tools::write_dataset< CoordinateVector<> >(group, "Coordinates", coords);
  HDF5Tools::write_dataset< double >(group, "NumberDensity", ntot);
  HDF5Tools::write_dataset< double >(group, "NeutralFractionH", nfracH);
  HDF5Tools::write_dataset< double >(group, "NeutralFractionHe", nfracHe);
  HDF5Tools::write_dataset< double >(group, "Temperature", temperature);
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}