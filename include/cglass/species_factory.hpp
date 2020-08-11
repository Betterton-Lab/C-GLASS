#ifndef _CGLASS_SPECIES_FACTORY_H_
#define _CGLASS_SPECIES_FACTORY_H_

#include "br_bead_species.hpp"
#include "crosslink_species.hpp"
#include "filament_species.hpp"
#include "optical_trap_species.hpp"
#include "rigid_filament_species.hpp"
#include "spherocylinder_species.hpp"
#include "spindle_species.hpp"

class SpeciesFactory {
public:
  SpeciesBase *CreateSpecies(const species_id sid, unsigned long seed) const {
    if (sid == +species_id::filament) {
      return new FilamentSpecies(seed);
    } else if (sid == +species_id::rigid_filament) {
      return new RigidFilamentSpecies(seed);
    } else if (sid == +species_id::br_bead) {
      return new BrBeadSpecies(seed);
    } else if (sid == +species_id::crosslink) {
      return new CrosslinkSpecies(seed);
    } else if (sid == +species_id::spindle) {
      return new SpindleSpecies(seed);
    } else if (sid == +species_id::spherocylinder) {
      return new SpherocylinderSpecies(seed);
    } else if (sid == +species_id::optical_trap) {
      return new OpticalTrapSpecies(seed);
    }
    Logger::Error("Species ID not recognized in SpeciesFactory!");
    return nullptr;
  }
};

#endif
