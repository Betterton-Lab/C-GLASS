#ifndef _CGLASS_DEFINITIONS_H_
#define _CGLASS_DEFINITIONS_H_

#include "macros.hpp"
#include "logger.hpp"
#include <math.h>

#define BETTER_ENUMS_DEFAULT_CONSTRUCTOR(Enum) \
 public:                                       \
  Enum() = default;
#include "enum.hpp"

//#if defined(_OPENMP)
//#define ENABLE_OPENMP
//#endif

BETTER_ENUM(species_id, unsigned char, br_bead, filament, rigid_filament,
            spherocylinder, spindle, crosslink, receptor, point_cover, none);
BETTER_ENUM(draw_type, unsigned char, fixed, orientation, bw, none);
BETTER_ENUM(potential_type, unsigned char, none, wca, soft, lj);
BETTER_ENUM(boundary_type, unsigned char, none = 0, box = 1, sphere = 2,
            budding = 3, wall = 4, protrusion = 5);
BETTER_ENUM(poly_state, unsigned char, grow, shrink, pause);
BETTER_ENUM(bind_state, unsigned char, unbound, singly, doubly, free);
BETTER_ENUM(obj_type, unsigned char, generic, bond, site, cortex);
BETTER_ENUM(comp_type, unsigned char, generic, mesh, point_cover);
BETTER_ENUM(density_type, unsigned char, linear, surface, volume);
BETTER_ENUM(shape, unsigned char, sphere, rod, generic)

struct graph_struct {
  double r[3];
  double u[3];
  double length;
  double diameter;
  double color;
  draw_type draw;
};

/* For unit testing */

// namespace unit_test {
// class Tester;
//}

#ifdef TESTS
#define UNIT_TEST       \
  template <typename T> \
  friend class UnitTest;
template <typename T>
class UnitTest {};
#else
#define UNIT_TEST
#endif

#endif
