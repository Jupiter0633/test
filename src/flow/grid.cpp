#include "flow/grid.h"

#include <stdexcept>

namespace hcw::flow {

Grid::Grid(Geometry g, int nr_, int nz_, double r0_, double r1_, double z0_, double z1_)
    : geometry(g), nr(nr_), nz(nz_), r0(r0_), r1(r1_), z0(z0_), z1(z1_) {
  if (nr < 2 || nz < 2) throw std::invalid_argument("grid needs at least 2 cells per direction");
  if (!(r1 > r0) || !(z1 > z0)) throw std::invalid_argument("grid extent must be positive");
  if (g == Geometry::Axisymmetric && r0 < 0.0) {
    throw std::invalid_argument("axisymmetric grid must have r0 >= 0");
  }
  dr = (r1 - r0) / nr;
  dz = (z1 - z0) / nz;
}

}  // namespace hcw::flow
