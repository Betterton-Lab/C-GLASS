#ifndef _CGLASS_RIGID_FILAMENT_SPECIES_H_
#define _CGLASS_RIGID_FILAMENT_SPECIES_H_

#include "minimum_distance.hpp"
#include "rigid_filament.hpp"
#include "species.hpp"

class RigidFilamentSpecies
    : public Species<RigidFilament, species_id::rigid_filament> {
protected:
  double fill_volume_;
  double packing_fraction_;

public:
  RigidFilamentSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();

  void AddMember();

  void Reserve();
  void UpdatePositions();
  void CustomInsert();
  // Redundant for filaments.
  virtual void CenteredOrientedArrangement() {}
};

#endif
