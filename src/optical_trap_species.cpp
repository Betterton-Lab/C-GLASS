/**
 * @author      : adamlamson (adamlamson@LamsonMacbookPro)
 * @file        : optical_trap_species
 * @created     : Tuesday Jul 07, 2020 09:41:45 MDT
 */

#include "cglass/optical_trap_species.hpp"

OpticalTrapSpecies::OpticalTrapSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::optical_trap);
}

void OpticalTrapSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  attach_species_ = sparams_.attach_species;
}

void OpticalTrapSpecies::InsertOpticalTraps(
    std::vector<SpeciesBase *> *species) {
  for (auto spec : *species) {
    if (attach_species_.compare(spec->GetSpeciesName()) == 0) {
      int num_insert = GetNInsert();
      if (spec->GetNInsert() < num_insert) {
        Logger::Warning("Number of optical traps exceeds %s number. Will go to "
                        "max number of species.",
                        spec->GetSpeciesName().c_str());
        num_insert = spec->GetNInsert();
      }
      std::vector<Object *> attach_objs;
      spec->GetMemberPtrs(attach_objs);

      for (auto &obj : attach_objs) {
        AddMember();
        members_.back().InsertAndAttach(obj);
        //const double *const pos = members_.back().GetScaledPosition();
        //Logger::Info("pos = %f, %f, %f\n", pos[0], pos[1], pos[2]);
      }

      break;
    }
  }
}

/*! \brief Apply all the forces from optical traps to attached objects
 *
 * \return void
 */
void OpticalTrapSpecies::ApplyOpticalTrapForces() {
  if (params_->no_midstep || params_->i_step % 2 == 1) {
    for (auto &otrap : members_) {
      otrap.ApplyOpticalTrapForce();
    }
  }
}

void OpticalTrapSpecies::CleanUp() { members_.clear(); }
