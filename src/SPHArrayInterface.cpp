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
 * @file SPHArrayInterface.cpp
 *
 * @brief SPHArrayInterface implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "SPHArrayInterface.hpp"
#include "DensityGrid.hpp"
#include "DensityGridTraversalJobMarket.hpp"
#include "Lock.hpp"
#include <cfloat>

/**
 * @brief Constructor.
 *
 * Non-periodic version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 * @param mapping_type Type of density mapping to use.
 */
SPHArrayInterface::SPHArrayInterface(const double unit_length_in_SI,
                                     const double unit_mass_in_SI,
                                     const std::string mapping_type)
    : DensityGridWriter("", false, DensityGridWriterFields(false), nullptr),
      _unit_length_in_SI(unit_length_in_SI), _unit_mass_in_SI(unit_mass_in_SI),
      _is_periodic(false), _mapping_type(get_mapping_type(mapping_type)),
      _octree(nullptr) {

  if (_mapping_type == SPHARRAY_MAPPING_PETKOVA) {
    gridding();
  }
  _time_log.output("time-log-file.txt", true);
}

/**
 * @brief Constructor.
 *
 * Periodic double precision version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 * @param box_anchor Coordinates of the left front bottom corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 * @param mapping_type Type of density mapping to use.
 */
SPHArrayInterface::SPHArrayInterface(const double unit_length_in_SI,
                                     const double unit_mass_in_SI,
                                     const double *box_anchor,
                                     const double *box_sides,
                                     const std::string mapping_type)
    : DensityGridWriter("", false, DensityGridWriterFields(false), nullptr),
      _unit_length_in_SI(unit_length_in_SI), _unit_mass_in_SI(unit_mass_in_SI),
      _is_periodic(true), _mapping_type(get_mapping_type(mapping_type)),
      _octree(nullptr) {

  _box.get_anchor()[0] = box_anchor[0] * _unit_length_in_SI;
  _box.get_anchor()[1] = box_anchor[1] * _unit_length_in_SI;
  _box.get_anchor()[2] = box_anchor[2] * _unit_length_in_SI;
  _box.get_sides()[0] = box_sides[0] * _unit_length_in_SI;
  _box.get_sides()[1] = box_sides[1] * _unit_length_in_SI;
  _box.get_sides()[2] = box_sides[2] * _unit_length_in_SI;

  if (_mapping_type == SPHARRAY_MAPPING_PETKOVA) {
    gridding();
  }
  _time_log.output("time-log-file.txt", true);
}

/**
 * @brief Constructor.
 *
 * Periodic single precision version.
 *
 * @param unit_length_in_SI Length unit used in the input arrays (in m).
 * @param unit_mass_in_SI Mass unit used in the input arrays (in kg).
 * @param box_anchor Coordinates of the left front bottom corner of the
 * simulation box (in the given length unit).
 * @param box_sides Side lengths of the simulation box (in the given length
 * unit).
 * @param mapping_type Type of density mapping to use.
 */
SPHArrayInterface::SPHArrayInterface(const double unit_length_in_SI,
                                     const double unit_mass_in_SI,
                                     const float *box_anchor,
                                     const float *box_sides,
                                     const std::string mapping_type)
    : DensityGridWriter("", false, DensityGridWriterFields(false), nullptr),
      _unit_length_in_SI(unit_length_in_SI), _unit_mass_in_SI(unit_mass_in_SI),
      _is_periodic(true), _mapping_type(get_mapping_type(mapping_type)),
      _octree(nullptr) {

  _box.get_anchor()[0] = box_anchor[0] * _unit_length_in_SI;
  _box.get_anchor()[1] = box_anchor[1] * _unit_length_in_SI;
  _box.get_anchor()[2] = box_anchor[2] * _unit_length_in_SI;
  _box.get_sides()[0] = box_sides[0] * _unit_length_in_SI;
  _box.get_sides()[1] = box_sides[1] * _unit_length_in_SI;
  _box.get_sides()[2] = box_sides[2] * _unit_length_in_SI;

  if (_mapping_type == SPHARRAY_MAPPING_PETKOVA) {
    gridding();
  }
  _time_log.output("time-log-file.txt", true);
}

/**
 * @brief Destructor.
 *
 * Frees up memory used by the internal Octree.
 */
SPHArrayInterface::~SPHArrayInterface() {
  delete _octree;
  _time_log.output("time-log-file.txt", true);
}

/**
 * @brief Reset the internal data values.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayInterface::reset(const double *x, const double *y, const double *z,
                              const double *h, const double *m,
                              const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  _neutral_fractions.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Reset the internal data values.
 *
 * This version uses single precision smoothing length and masses.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayInterface::reset(const double *x, const double *y, const double *z,
                              const float *h, const float *m,
                              const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  _neutral_fractions.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Reset the internal data values.
 *
 * This version uses single precision coordinates, smoothing lengths and masses.
 *
 * @param x Array containing x coordinates (in the given length unit).
 * @param y Array containing y coordinates (in the given length unit).
 * @param z Array containing z coordinates (in the given length unit).
 * @param h Array containing smoothing lengths (in the given length unit).
 * @param m Array containing masses (in the given mass unit).
 * @param npart Number of elements in each of the arrays.
 */
void SPHArrayInterface::reset(const float *x, const float *y, const float *z,
                              const float *h, const float *m,
                              const size_t npart) {

  delete _octree;

  _positions.resize(npart);
  _smoothing_lengths.resize(npart, 0.);
  _masses.resize(npart, 0.);
  _neutral_fractions.resize(npart, 0.);
  for (size_t i = 0; i < npart; ++i) {
    _positions[i][0] = x[i] * _unit_length_in_SI;
    _positions[i][1] = y[i] * _unit_length_in_SI;
    _positions[i][2] = z[i] * _unit_length_in_SI;
    _smoothing_lengths[i] = h[i] * _unit_length_in_SI;
    _masses[i] = m[i] * _unit_mass_in_SI;
  }

  if (!_is_periodic) {
    CoordinateVector<> minpos(DBL_MAX);
    CoordinateVector<> maxpos(-DBL_MAX);
    for (size_t i = 0; i < npart; ++i) {
      minpos = CoordinateVector<>::min(minpos, _positions[i]);
      maxpos = CoordinateVector<>::max(maxpos, _positions[i]);
    }

    maxpos -= minpos;
    _box.get_anchor() = minpos - 0.005 * maxpos;
    _box.get_sides() = 1.01 * maxpos;
  }
}

/**
 * @brief Get a pointer to the internal Octree.
 *
 * @return Pointer to the internal Octree.
 */
Octree *SPHArrayInterface::get_octree() { return _octree; }

/**
 * @brief Initialize the internal Octree.
 */
void SPHArrayInterface::initialize() {
  _octree = new Octree(_positions, _box, _is_periodic);
  _octree->set_auxiliaries(_smoothing_lengths, Octree::max< double >);
  //_dens_map = new DensityMapping();
}

/**
 * @brief Initialize the pre-computed array of density values.
 */
void SPHArrayInterface::gridding() {
  double phi, r0, R_0, mu0;
  int_fast32_t i, j, k;

  _time_log.start("Gridding");

  const double h = 1.0;
  const int_fast32_t n = 150;
  const int_fast32_t nr1 = 50;
  const int_fast32_t nr2 = 199;
  const double rl = 0.1;
  const double mul = 0.98;
  const double cphil = 0.98;

  _density_values.resize(nr1 + nr2 + 2);
  for (i = 0; i <= nr1 + nr2 + 1; ++i) {
    _density_values[i].resize(2 * n + 1);
    for (j = 0; j <= 2 * n; ++j) {
      _density_values[i][j].resize(2 * n + 1, 0.);
    }
  }

  for (i = 0; i < nr1; ++i) {
    r0 = (rl / nr1) * (i + 1) * h;
    for (j = 0; j < n; ++j) {
      mu0 = (mul / (n - 1)) * j;
      R_0 = r0 * (std::sqrt(1.0 - mu0 * mu0)) / mu0;
      for (k = 0; k < n; ++k) {
        const double cosphi = (cphil / (n - 1)) * k;
        phi = std::acos(cosphi);
        _density_values[i + 1][j][k] = full_integral(phi, cosphi, r0, R_0, h);
      }

      for (k = 1; k <= n; ++k) {
        const double cosphi = cphil + ((1.0 - cphil) / n) * k;
        phi = std::acos(cosphi);
        _density_values[i + 1][j][n + k - 1] =
            full_integral(phi, cosphi, r0, R_0, h);
      }
    }

    for (j = 1; j <= n; ++j) {
      mu0 = mul + ((1.0 - mul) / n) * j;
      R_0 = r0 * (std::sqrt(1.0 - mu0 * mu0)) / mu0;
      for (k = 0; k < n; ++k) {
        const double cosphi = (cphil / (n - 1)) * k;
        phi = std::acos(cosphi);
        _density_values[i + 1][n + j - 1][k] =
            full_integral(phi, cosphi, r0, R_0, h);
      }

      for (k = 1; k <= n; ++k) {
        const double cosphi = cphil + ((1.0 - cphil) / n) * k;
        phi = std::acos(cosphi);
        _density_values[i + 1][n + j - 1][n + k - 1] =
            full_integral(phi, cosphi, r0, R_0, h);
      }
    }
  }

  for (i = 1; i <= nr2; ++i) {
    r0 = rl + ((2.0 - rl) / nr2) * i * h;
    for (j = 0; j < n; ++j) {
      mu0 = (mul / (n - 1)) * j;
      R_0 = r0 * (std::sqrt(1.0 - mu0 * mu0)) / mu0;

      for (k = 0; k < n; ++k) {
        const double cosphi = (cphil / (n - 1)) * k;
        phi = std::acos(cosphi);
        _density_values[nr1 + i][j][k] = full_integral(phi, cosphi, r0, R_0, h);
      }

      for (k = 1; k <= n; ++k) {
        const double cosphi = cphil + ((1.0 - cphil) / n) * k;
        phi = std::acos(cosphi);
        _density_values[nr1 + i][j][n + k - 1] =
            full_integral(phi, cosphi, r0, R_0, h);
      }
    }

    for (j = 1; j <= n; ++j) {
      mu0 = mul + ((1.0 - mul) / n) * j;
      R_0 = r0 * (std::sqrt(1.0 - mu0 * mu0)) / mu0;

      for (k = 0; k < n; ++k) {
        const double cosphi = (cphil / (n - 1)) * k;
        phi = std::acos(cosphi);
        _density_values[nr1 + i][n + j - 1][k] =
            full_integral(phi, cosphi, r0, R_0, h);
      }

      for (k = 1; k <= n; ++k) {
        const double cosphi = cphil + ((1.0 - cphil) / n) * k;
        phi = std::acos(cosphi);
        _density_values[nr1 + i][n + j - 1][n + k - 1] =
            full_integral(phi, cosphi, r0, R_0, h);
      }
    }
  }

  i = nr1 + nr2;
  for (j = 0; j < 2 * n; ++j) {
    for (k = 0; k < 2 * n; ++k) {
      _density_values[i + 1][j][k] = _density_values[i][j][k];
    }
  }

  _time_log.end("Gridding");
}

/**
 * @brief Function that gives the 3D integral of the kernel of
 * a particle for a given vertex of a cell face, interpolated from a
 * pre-computed grid.
 *
 * @param phi Azimuthal angle of the vertex.
 * @param cosphi Cosine of the azimuthal angle.
 * @param r0_old Distance from the particle to the face of the cell.
 * @param R_0_old Distance from the orthogonal projection of the particle
 * onto the face of the cell to a side of the face (containing the vertex).
 * @param h_old The kernel smoothing length of the particle.
 * @return The integral of the kernel for the given vertex.
 */

double SPHArrayInterface::gridded_integral(const double phi,
                                           const double cosphi,
                                           const double r0_old,
                                           const double R_0_old,
                                           const double h_old) const {

  double r0, R_0, cphi, mu0, h;
  int_fast32_t i, j, k;
  double fx1, fx2, fx3, fx4, fy1, fy2, fz, frac;
  const int_fast32_t nr01 = 50;
  const int_fast32_t nr02 = 199;
  const int_fast32_t nR_01 = 150;
  const int_fast32_t nR_02 = 150;
  const int_fast32_t nphi1 = 150;
  const int_fast32_t nphi2 = 150;
  const double rl = 0.1;
  const double mul = 0.98;
  const double cphil = 0.98;

  const double h_old_inverse = 1. / h_old;
  h = 1.0;
  r0 = r0_old * h_old_inverse;
  R_0 = R_0_old * h_old_inverse;
  cphi = cosphi;

  cmac_assert_message(r0 >= 0., "r0 < 0: %g %g", r0, r0_old);
  cmac_assert_message(R_0 >= 0., "R_0 < 0: %g %g", R_0, R_0_old);

  if (r0 == 0.0) {
    return 0.0;
  }
  if (R_0 == 0.0) {
    return 0.0;
  }
  if (phi == 0.0) {
    return 0.0;
  }

  mu0 = r0 / std::sqrt(r0 * r0 + R_0 * R_0);

  if (r0 > 2.0)
    r0 = 2.0;

  if (r0 < rl) {
    i = static_cast< int_fast32_t >(r0 * nr01 / rl);
  } else {
    i = nr01 + static_cast< int_fast32_t >((r0 - rl) * nr02 / ((2.0 - rl) * h));
  }

  if (mu0 < mul) {
    j = static_cast< int_fast32_t >(mu0 * (nR_01 - 1) / mul);
  } else {
    j = nR_01 - 1 +
        static_cast< int_fast32_t >((mu0 - mul) * nR_02 / ((1.0 - mul)));
  }

  if (cphi < cphil) {
    k = static_cast< int_fast32_t >(cphi * (nphi1 - 1) / cphil);
  } else {
    k = nphi1 - 1 +
        static_cast< int_fast32_t >((cphi - cphil) * nphi2 / ((1.0 - cphil)));
  }

  if (i < 0 || i > nr01 + nr02 || j < 0 || j > nR_01 + nR_02 - 1 || k < 0 ||
      k > nphi1 + nphi2) {
    cmac_warning("i, j, k: %" PRIiFAST32 " %" PRIiFAST32 " %" PRIiFAST32 " \n",
                 i, j, k);
  }

  if (r0 < rl) {
    frac = (r0 * nr01 - rl * i * h) / (rl * h);
  } else {
    frac = (r0 * nr02 - rl * nr02 - (2.0 - rl) * (i - nr01) * h) /
           ((2.0 - rl) * h);
  }
  fx1 = frac * SPHArrayInterface::get_gridded_density_value(i + 1, j, k) +
        (1. - frac) * SPHArrayInterface::get_gridded_density_value(i, j, k);
  fx2 = frac * SPHArrayInterface::get_gridded_density_value(i + 1, j + 1, k) +
        (1. - frac) * SPHArrayInterface::get_gridded_density_value(i, j + 1, k);
  fx3 = frac * SPHArrayInterface::get_gridded_density_value(i + 1, j, k + 1) +
        (1. - frac) * SPHArrayInterface::get_gridded_density_value(i, j, k + 1);
  fx4 =
      frac * SPHArrayInterface::get_gridded_density_value(i + 1, j + 1, k + 1) +
      (1. - frac) *
          SPHArrayInterface::get_gridded_density_value(i, j + 1, k + 1);

  if (mu0 < mul) {
    frac = (mu0 * (nR_01 - 1) - mul * j) / mul;
  } else {
    frac = (mu0 * nR_02 - mul * nR_02 - (1.0 - mul) * (j - nR_01 + 1)) /
           (1.0 - mul);
  }
  fy1 = frac * fx2 + (1. - frac) * fx1;
  fy2 = frac * fx4 + (1. - frac) * fx3;

  if (cphi < cphil) {
    frac = (cphi * (nphi1 - 1) - cphil * k) / cphil;
  } else {
    frac = (cphi * nphi2 - cphil * nphi2 - (1.0 - cphil) * (k - nphi1 + 1)) /
           (1.0 - cphil);
  }
  fz = frac * fy2 + (1 - frac) * fy1;

  if ((j == (nR_01 + nR_02 - 1)) || (k == (nphi1 + nphi2 - 1))) {
    fz = 0.0;
  }

  return fz;
}

/**
 * @brief Function that gives the 3D integral of the kernel of
 * a particle for a given vertex of a cell face.
 *
 * @param phi Azimuthal angle of the vertex.
 * @param cosphi Cosine of the azimuthal angle.
 * @param r0 Distance from the particle to the face of the cell.
 * @param R_0 Distance from the orthogonal projection of the particle
 * onto the face of the cell to a side of the face (containing the vertex).
 * @param h The kernel smoothing length of the particle.
 * @return The integral of the kernel for the given vertex.
 */

double SPHArrayInterface::full_integral(const double phi, const double cosphi,
                                        const double r0, const double R_0,
                                        const double h) {

  double B1, B2, B3, mu, a, logs, u;
  double full_int;
  double a2, cosp, cosp2, tanp;
  double r2, R, linedist2, phi1, phi2;
  double I0, I1, I_1, I_2, I_3, I_4, I_5;
  double D2, D3;

  B1 = 0.0;
  B2 = 0.0;
  B3 = 0.0;

  if (r0 == 0.0) {
    return 0.0;
  }
  if (R_0 == 0.0) {
    return 0.0;
  }
  if (phi == 0.0) {
    return 0.0;
  }

  const double h_inv = 1. / h;
  const double r0_inv = 1. / r0;
  const double h2 = h * h;
  const double r02 = r0 * r0;
  const double R_02 = R_0 * R_0;
  const double r03 = r02 * r0;
  const double r0h = r0 * h_inv;
  const double r0h2 = r0h * r0h;
  const double r0h3 = r0h2 * r0h;
  const double r0h_2 = h2 * r0_inv * r0_inv;
  const double r0h_3 = r0h_2 * h * r0_inv;

  // Setting up the B1, B2, B3 constants of integration.

  if (r0 >= 2.0 * h) {
    B3 = 0.25 * h2 * h;
  } else if (r0 > h) {
    B3 = 0.25 * r03 *
         (-4. / 3. + r0h - 0.3 * r0h2 + 1. / 30. * r0h3 - 1. / 15. * r0h_3 +
          8. / 5. * r0h_2);
    B2 = 0.25 * r03 *
         (-4. / 3. + r0h - 0.3 * r0h2 + 1. / 30. * r0h3 - 1. / 15. * r0h_3);
  } else {
    B3 = 0.25 * r03 * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3 + 7. / 5. * r0h_2);
    B2 = 0.25 * r03 * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3 - 1. / 5. * r0h_2);
    B1 = 0.25 * r03 * (-2. / 3. + 0.3 * r0h2 - 0.1 * r0h3);
  }

  a = R_0 * r0_inv;
  a2 = a * a;

  linedist2 = r02 + R_02;
  R = R_0 / cosphi;
  r2 = r02 + R * R;

  full_int = 0.0;
  D2 = 0.0;
  D3 = 0.0;

  if (linedist2 <= h2) {
    ////// phi1 business /////
    cosp = R_0 / std::sqrt(h2 - r02);
    const double sinphi1 = std::sqrt((1. + cosp) * (1. - cosp));
    phi1 = acos(cosp);

    cosp2 = cosp * cosp;
    mu = cosp / std::sqrt(a2 + cosp2);

    tanp = sinphi1 / cosp;

    I0 = phi1;
    I_2 = phi1 + a2 * tanp;
    I_4 = phi1 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sinphi1 * std::sqrt((1. - mu) * (1. + mu));
    const double u2 = u * u;
    const double u3 = u2 * u;
    logs = std::log((1. + u) / (1. - u));
    I1 = std::atan(u / a);

    I_1 = 0.5 * a * logs + I1;
    I_3 = I_1 + 0.25 * a * (1. + a2) * (2 * u / (1. - u2) + logs);
    I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. *
                    ((10 * u - 6 * u3) / ((1. - u2) * (1. - u2)) + 3. * logs);

    D2 = -1. / 6. * I_2 + 0.25 * r0h * I_3 - 0.15 * r0h2 * I_4 +
         1. / 30. * r0h3 * I_5 - 1. / 60. * r0h_3 * I1 + (B1 - B2) / r03 * I0;

    ////// phi2 business /////
    cosp = R_0 / std::sqrt(4.0 * h * h - r0 * r0);
    const double sinphi2 = std::sqrt((1. - cosp) * (1. + cosp));
    phi2 = std::acos(cosp);

    cosp2 = cosp * cosp;
    mu = cosp / std::sqrt(a2 + cosp2);

    tanp = sinphi2 / cosp;

    I0 = phi2;
    I_2 = phi2 + a2 * tanp;
    I_4 = phi2 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sinphi2 * std::sqrt((1. - mu) * (1. + mu));
    const double uu2 = u * u;
    const double uu3 = uu2 * u;
    logs = std::log((1. + u) / (1. - u));
    I1 = std::atan(u / a);

    I_1 = 0.5 * a * logs + I1;
    I_3 = I_1 + 0.25 * a * (1. + a2) * (2 * u / (1. - uu2) + logs);
    I_5 = I_3 +
          a * (1. + a2) * (1. + a2) / 16. *
              ((10. * u - 6. * uu3) / ((1. - uu2) * (1. - uu2)) + 3. * logs);

    D3 = 1. / 3. * I_2 - 0.25 * r0h * I_3 + 3. / 40. * r0h2 * I_4 -
         1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 +
         D2;
  } else if (linedist2 <= 4.0 * h2) {
    ////// phi2 business /////
    cosp = R_0 / std::sqrt(4.0 * h2 - r02);
    const double sinphi2 = std::sqrt((1. - cosp) * (1. + cosp));
    phi2 = std::acos(cosp);

    cosp2 = cosp * cosp;
    mu = cosp / std::sqrt(a2 + cosp2);

    tanp = sinphi2 / cosp;

    I0 = phi2;
    I_2 = phi2 + a2 * tanp;
    I_4 = phi2 + 2. * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    u = sinphi2 * std::sqrt((1. - mu) * (1. + mu));
    const double u2 = u * u;
    const double u3 = u2 * u;
    logs = std::log((1. + u) / (1. - u));
    I1 = std::atan(u / a);

    I_1 = 0.5 * a * logs + I1;
    I_3 = I_1 + 0.25 * a * (1. + a2) * (2. * u / (1. - u2) + logs);
    I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. *
                    ((10. * u - 6. * u3) / ((1. - u2) * (1. - u2)) + 3. * logs);

    D3 = 1. / 3. * I_2 - 0.25 * r0h * I_3 + 3. / 40. * r0h2 * I_4 -
         1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 +
         D2;
  }

  //////////////////////////////////
  // Calculating I_n expressions. //
  //////////////////////////////////

  cosp = cosphi;
  const double sinphi = std::sqrt((1. - cosp) * (1. + cosp));
  cosp2 = cosp * cosp;
  mu = cosp / std::sqrt(a2 + cosp2);

  tanp = sinphi / cosphi;

  I0 = phi;
  I_2 = phi + a2 * tanp;
  I_4 = phi + 2. * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

  u = sinphi * std::sqrt((1. - mu) * (1. + mu));
  const double u2 = u * u;
  logs = std::log((1. + u) / (1. - u));
  I1 = std::atan(u / a);

  I_1 = 0.5 * a * logs + I1;
  I_3 = I_1 + 0.25 * a * (1. + a2) * (2 * u / (1. - u2) + logs);
  I_5 =
      I_3 + a * (1. + a2) * (1. + a2) / 16. *
                ((10. * u - 6. * u2 * u) / ((1. - u2) * (1. - u2)) + 3. * logs);

  // Calculating the integral expression.

  if (r2 < h2) {
    full_int = M_1_PI * r0h3 *
               (1. / 6. * I_2 - 3. / 40. * r0h2 * I_4 + 1. / 40. * r0h3 * I_5 +
                B1 / r03 * I0);
  } else if (r2 < 4.0 * h2) {
    full_int = M_1_PI * r0h3 *
               (0.25 * (4. / 3. * I_2 - r0h * I_3 + 0.3 * r0h2 * I_4 -
                        1. / 30. * r0h3 * I_5 + 1. / 15. * r0h_3 * I1) +
                B2 / r03 * I0 + D2);
  } else {
    full_int = M_1_PI * r0h3 * (-0.25 * r0h_3 * I1 + B3 / r03 * I0 + D3);
  }

  return full_int;
}

/**
 * @brief Function that calculates the mass contribution of a particle
 * towards the total mass of a cell.
 *
 * @param cell Geometrical information about the cell.
 * @param particle The particle position.
 * @param h The kernel smoothing length of the particle.
 * @return The mass contribution of the particle to the cell.
 */
double SPHArrayInterface::mass_contribution(const Cell &cell,
                                            const CoordinateVector<> particle,
                                            const double h) const {

  double M, Msum;

  Msum = 0.;
  M = 0.;

  std::vector< Face > face_vector = cell.get_faces();

  // Loop over each face of a cell.
  for (size_t i = 0; i < face_vector.size(); i++) {

    CoordinateVector<> vert_position1;
    CoordinateVector<> projected_particle;
    double r0 = 0.;
    double ar0 = 0.;
    double s2 = 0.;

    // Loop over the vertices of each face.
    for (Face::Vertices j = face_vector[i].first_vertex();
         j != face_vector[i].last_vertex(); ++j) {
      if (j == face_vector[i].first_vertex()) { // Calculating the distance from
                                                // particle to each face of a
                                                // cell
        // http://mathinsight.org/distance_point_plane
        // http://mathinsight.org/forming_planes
        Face::Vertices j_twin = j;
        vert_position1 = j_twin.get_position();
        const CoordinateVector<> vert_position2 = (++j_twin).get_position();
        const CoordinateVector<> vert_position3 = (++j_twin).get_position();

        const double A = (vert_position2[1] - vert_position1[1]) *
                             (vert_position3[2] - vert_position1[2]) -
                         (vert_position3[1] - vert_position1[1]) *
                             (vert_position2[2] - vert_position1[2]);
        const double B = (vert_position2[2] - vert_position1[2]) *
                             (vert_position3[0] - vert_position1[0]) -
                         (vert_position3[2] - vert_position1[2]) *
                             (vert_position2[0] - vert_position1[0]);
        const double C = (vert_position2[0] - vert_position1[0]) *
                             (vert_position3[1] - vert_position1[1]) -
                         (vert_position3[0] - vert_position1[0]) *
                             (vert_position2[1] - vert_position1[1]);
        const double D = -A * vert_position1[0] - B * vert_position1[1] -
                         C * vert_position1[2];

        const double norm = std::sqrt(A * A + B * B + C * C);
        const double inverse_norm = 1. / norm;
        r0 = (A * particle[0] + B * particle[1] + C * particle[2] + D) *
             inverse_norm;
        ar0 = std::fabs(r0);

        // Calculate of the orthogonal projection of the particle position onto
        // the face.
        projected_particle[0] = particle[0] - r0 * A * inverse_norm;
        projected_particle[1] = particle[1] - r0 * B * inverse_norm;
        projected_particle[2] = particle[2] - r0 * C * inverse_norm;

        // s2 contains information about the orientation of the face vertices.
        s2 = vert_position1[0] * (vert_position2[1] * vert_position3[2] -
                                  vert_position2[2] * vert_position3[1]) +
             vert_position1[1] * (vert_position2[2] * vert_position3[0] -
                                  vert_position2[0] * vert_position3[2]) +
             vert_position1[2] * (vert_position2[0] * vert_position3[1] -
                                  vert_position2[1] * vert_position3[0]);
      }

      Face::Vertices j_twin = j;
      const CoordinateVector<> vert_position2 = j_twin.get_position();
      CoordinateVector<> vert_position3;

      if (++j_twin == face_vector[i].last_vertex()) {
        vert_position3 = vert_position1;
      } else {
        vert_position3 = j_twin.get_position();
      }

      const double r23 = (vert_position2 - vert_position3).norm();
      const double r12 = (projected_particle - vert_position2).norm();
      const double r13 = (projected_particle - vert_position3).norm();
      const double cosa = ((vert_position3[0] - vert_position2[0]) *
                               (projected_particle[0] - vert_position2[0]) +
                           (vert_position3[1] - vert_position2[1]) *
                               (projected_particle[1] - vert_position2[1]) +
                           (vert_position3[2] - vert_position2[2]) *
                               (projected_particle[2] - vert_position2[2])) /
                          (r12 * r23);

      double cosphi1 = 0.;
      double phi1 = 0.;
      double phi2 = 0.;

      if (std::fabs(cosa) < 1.0) {
        cosphi1 = std::sqrt((1. - cosa) * (1. + cosa));
      } else {
        if (std::fabs(cosa) - 1.0 < 0.00001) {
          cosphi1 = 0.0;
        } else {
          cmac_warning("Error: cosa > 1: %g\n", cosa);
        }
      }
      const double R_0 = r12 * cosphi1;
      const double cosphi2 = R_0 / r13;

      const double s1 =
          projected_particle[0] * (vert_position2[1] * vert_position3[2] -
                                   vert_position2[2] * vert_position3[1]) +
          projected_particle[1] * (vert_position2[2] * vert_position3[0] -
                                   vert_position2[0] * vert_position3[2]) +
          projected_particle[2] * (vert_position2[0] * vert_position3[1] -
                                   vert_position2[1] * vert_position3[0]);

      if (R_0 < r12) {
        phi1 = std::acos(cosphi1);
      } else {
        if ((R_0 - r12) < 0.00001 * h) {
          phi1 = 0.0;
        } else {
          cmac_warning("Error: R0 > r12: %g\n", R_0 - r12);
        }
      }
      if (R_0 < r13) {
        phi2 = std::acos(cosphi2);
      } else {
        if ((R_0 - r13) < 0.00001 * h) {
          phi2 = 0.0;
        } else {
          cmac_warning("Error: R0 > r13: %g\n", R_0 - r13);
        }
      }

      // Find out if the vertex integral will contribute positively or
      // negatively to the cell mass.
      if (s1 * s2 * r0 <= 0) {
        M = -1.;
      } else {
        M = 1.;
      }

      cmac_assert_message(ar0 >= 0., "Wrong from mass_contribution: r0 = %g",
                          ar0);
      cmac_assert_message(R_0 >= 0., "Wrong from mass_contribution: R_0 = %g",
                          R_0);

      // Calculate the vertex integral.
      const bool is_pre_computed = true;
      const double sinphi1 = std::sqrt((1. - cosphi1) * (1. + cosphi1));
      const double sinphi2 = std::sqrt((1. - cosphi2) * (1. + cosphi2));
      if ((r12 * sinphi1 >= r23) || (r13 * sinphi2 >= r23)) {
        if (phi1 >= phi2) {
          if (is_pre_computed) {
            M = M * (gridded_integral(phi1, cosphi1, ar0, R_0, h) -
                     gridded_integral(phi2, cosphi2, ar0, R_0, h));
          } else {
            M = M * (full_integral(phi1, cosphi1, ar0, R_0, h) -
                     full_integral(phi2, cosphi2, ar0, R_0, h));
          }
        } else {
          if (is_pre_computed) {
            M = M * (gridded_integral(phi2, cosphi2, ar0, R_0, h) -
                     gridded_integral(phi1, cosphi1, ar0, R_0, h));
          } else {
            M = M * (full_integral(phi2, cosphi2, ar0, R_0, h) -
                     full_integral(phi1, cosphi1, ar0, R_0, h));
          }
        }
      } else {
        if (is_pre_computed) {
          M = M * (gridded_integral(phi1, cosphi1, ar0, R_0, h) +
                   gridded_integral(phi2, cosphi2, ar0, R_0, h));
        } else {
          M = M * (full_integral(phi1, cosphi1, ar0, R_0, h) +
                   full_integral(phi2, cosphi2, ar0, R_0, h));
        }
      }
      Msum = Msum + M;
    }
  }

  // Ensure there is no negative mass
  Msum = std::max(Msum, 1.e-6);

  return Msum;
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
DensityValues SPHArrayInterface::operator()(const Cell &cell) {

  // time_log.start("Density_mapping");

  DensityValues values;

  const CoordinateVector<> position = cell.get_cell_midpoint();
  double density = 0.;

  if (_mapping_type == SPHARRAY_MAPPING_M_OVER_V) {
    density = _masses[0] / cell.get_volume();
  } else if (_mapping_type == SPHARRAY_MAPPING_CENTROID) {
    const std::vector< uint_fast32_t > ngbs = _octree->get_ngbs(position);
    const size_t numngbs = ngbs.size();
    for (size_t i = 0; i < numngbs; ++i) {
      const size_t index = ngbs[i];
      double r;
      if (!_box.get_sides().x()) {
        r = (position - _positions[index]).norm();
      } else {
        r = _box.periodic_distance(position, _positions[index]).norm();
      }
      const double h = _smoothing_lengths[index];
      const double u = r / h;
      const double m = _masses[index];
      const double splineval = m * CubicSplineKernel::kernel_evaluate(u, h);
      density += splineval;
    }
  } else if (_mapping_type == SPHARRAY_MAPPING_PETKOVA) {
    CoordinateVector<> position = cell.get_cell_midpoint();

    // Find the vertex that is furthest away from the cell midpoint.
    std::vector< Face > face_vector = cell.get_faces();
    double radius = 0.0;
    for (unsigned int i = 0; i < face_vector.size(); i++) {
      for (Face::Vertices j = face_vector[i].first_vertex();
           j != face_vector[i].last_vertex(); ++j) {
        double distance = (j.get_position() - position).norm();
        if (distance > radius)
          radius = distance;
      }
    }

    // Find the neighbours that are contained inside of a sphere of centre the
    // cell midpoint
    // and radius given by the distance to the furthest vertex.
    std::vector< uint_fast32_t > ngbs =
        _octree->get_ngbs_sphere(position, radius);
    const unsigned int numngbs = ngbs.size();

    // Loop over all the neighbouring particles and calculate their mass
    // contributions.
    for (unsigned int i = 0; i < numngbs; i++) {
      const unsigned int index = ngbs[i];
      const double h = _smoothing_lengths[index] / 2.0;
      const CoordinateVector<> particle = _positions[index];
      if (h < 0)
        cmac_warning("h < 0: %g, %u", h, index);
      density += SPHArrayInterface::mass_contribution(cell, particle, h) *
                 _masses[index];
    }

    // Divide the cell mass by the cell volume to get density.
    density = density / cell.get_volume();
  }

  // Ensure that the density > 0
  if (density <= 0.0)
    density = _masses[0] / cell.get_volume() * 1e-6;

  values.set_number_density(density / 1.6737236e-27);
  values.set_temperature(8000.);
  values.set_ionic_fraction(ION_H_n, 1.e-6);
#ifdef HAS_HELIUM
  values.set_ionic_fraction(ION_He_n, 1.e-6);
#endif

  return values;
}

/**
 * @brief Fill the given array with the remapped neutral fractions.
 *
 * Double precision version.
 *
 * @param nH Array to fill.
 */
void SPHArrayInterface::fill_array(double *nH) {
  for (size_t i = 0; i < _neutral_fractions.size(); ++i) {
    nH[i] = _neutral_fractions[i];
  }
}

/**
 * @brief Fill the given array with the remapped neutral fractions.
 *
 * Single precision version.
 *
 * @param nH Array to fill.
 */
void SPHArrayInterface::fill_array(float *nH) {
  for (size_t i = 0; i < _neutral_fractions.size(); ++i) {
    nH[i] = _neutral_fractions[i];
  }
}

/**
 * @brief Write a snapshot.
 *
 * @param grid DensityGrid to write.
 * @param iteration Iteration number to use in the snapshot file name(s).
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 * @param hydro_units Internal unit system for the hydro.
 */
void SPHArrayInterface::write(DensityGrid &grid, uint_fast32_t iteration,
                              ParameterFile &params, double time,
                              const InternalHydroUnits *hydro_units) {

  _time_log.start("Inverse_mapping");

  for (unsigned int i = 0; i < _neutral_fractions.size(); ++i) {
    _neutral_fractions[i] = 1.0;
  }

  std::vector< Lock > locks(_neutral_fractions.size());

  std::pair< cellsize_t, cellsize_t > block =
      std::make_pair(0, grid.get_number_of_cells());
  WorkDistributor< DensityGridTraversalJobMarket< InverseMappingFunction >,
                   DensityGridTraversalJob< InverseMappingFunction > >
      workers;

  InverseMappingFunction do_calculation(*this, locks);

  DensityGridTraversalJobMarket< InverseMappingFunction > jobs(
      grid, do_calculation, block);
  workers.do_in_parallel(jobs);

  _time_log.end("Inverse_mapping");
  _time_log.output("time-log-file.txt", true);
}

/**
 * @brief Write a snapshot for a split grid.
 *
 * @param grid_creator Grid.
 * @param counter Counter value to add to the snapshot file name.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void SPHArrayInterface::write(
    DensitySubGridCreator< DensitySubGrid > &grid_creator,
    const uint_fast32_t counter, ParameterFile &params, double time) {}
