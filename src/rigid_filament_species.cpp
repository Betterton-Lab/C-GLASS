#include "cglass/rigid_filament_species.hpp"

RigidFilamentSpecies::RigidFilamentSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::rigid_filament);
}
void RigidFilamentSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  // fill_volume_ = 0;
  packing_fraction_ = sparams_.packing_fraction;
#ifdef TRACE
  if (packing_fraction_ > 0) {
    Logger::Warning("Simulation run in trace mode with a potentially large "
                    "number of objects in species %s (packing fraction ="
                    " %2.4f)",
                    GetSID()._to_string(), packing_fraction_);
    fprintf(stderr, "Continue anyway? (y/N) ");
    char c;
    if (std::cin.peek() != 'y') {
      Logger::Error("Terminating simulation by user request");
    } else if (!(std::cin >> c)) {
      Logger::Error("Invalid input");
    } else {
      fprintf(stderr, "Resuming simulation\n");
      std::cin.ignore();
    }
  }
#endif

  // double min_length = 2 * sparams_.min_length;
  // if (!sparams_.polydispersity_flag && sparams_.length < min_length) {
  //  Logger::Warning(
  //      "RigidFilament length %2.2f is less than minimum filament length"
  //      " %2.2f for minimum bond length %2.2f. Setting length to "
  //      "minimum length.",
  //      sparams_.length, min_length, min_length);
  //  sparams_.length = min_length;
  //}
}

void RigidFilamentSpecies::CustomInsert() {
  Species::CustomInsert();
  if (sparams_.constrain_motion_flag) {
    if (n_members_ == 2) {

      // Variables to store
      double r_min[3], lambda, mu;
      double dr[3] = {};
      double n_dim = params_->n_dim;
      // Find min distance unit vector between carrier lines to use as
      // constraining vector
      const double *r_1 = members_[0].GetPosition();
      const double *u_1 = members_[0].GetOrientation();
      const double *r_2 = members_[1].GetPosition();
      const double *u_2 = members_[1].GetOrientation();
      for (int i = 0; i < n_dim; ++i) {
        dr[i] = r_2[i] - r_1[i];
      }
      double dr_dot_u_1 = dot_product(n_dim, dr, u_1);
      double dr_dot_u_2 = dot_product(n_dim, dr, u_2);
      double u_1_dot_u_2 = dot_product(n_dim, u_1, u_2);
      double denom = 1.0 - SQR(u_1_dot_u_2);
      if (denom < 1.0e-12) {
        lambda = dr_dot_u_1 / 2.0;
        mu = -dr_dot_u_2 / 2.0;
      } else {
        lambda = (dr_dot_u_1 - u_1_dot_u_2 * dr_dot_u_2) / denom;
        mu = (-dr_dot_u_2 + u_1_dot_u_2 * dr_dot_u_1) / denom;
      }

      /* Calculate minimum distance between two lines. */
      double r_min_mag2 = 0.0;
      for (int i = 0; i < n_dim; ++i) {
        r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
        r_min_mag2 += SQR(r_min[i]);
      }

      normalize_vector(r_min, n_dim);
      Logger::Info("Constraining motion of rigid rods to plane with vector = "
                   "%f, %f, %f \n",
                   r_min[0], r_min[1], r_min[2]);
      //" << std::endl;
      // printf("r_min_vec = %f, %f, %f \n", r_min[0], r_min[1], r_min[2]);

      members_[0].SetConstrainVec(r_min);
      members_[1].SetConstrainVec(r_min);
    } else {
      Logger::Warning("Cannot constrain motion of more than two filaments. "
                      "Constraint not applied!");
      sparams_.constrain_motion_flag = false;
    }
  }
}

void RigidFilamentSpecies::UpdatePositions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  rigid_filament_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  rigid_filament_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    rigid_filament_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i)
      for (auto it = chunks[i].first; it != chunks[i].second; ++it)
        it->UpdatePosition();
  }
#else
  for (rigid_filament_iterator it = members_.begin(); it != members_.end();
       ++it)
    it->UpdatePosition();
#endif
}

void RigidFilamentSpecies::Reserve() {
  int max_insert = GetNInsert();
  if (packing_fraction_ > 0) {
    // double min_length = 2 * sparams_.min_bond_length;
    double min_length = 2 * sparams_.min_length;
    double diameter = sparams_.diameter;
    double min_vol = 0;
    if (params_->n_dim == 2) {
      min_vol = diameter * min_length + 0.25 * M_PI * diameter * diameter;
    }
    if (params_->n_dim == 3) {
      min_vol = 0.25 * M_PI * diameter * diameter * min_length +
                1.0 / 6.0 * M_PI * diameter * diameter * diameter;
    }
    max_insert = (int)ceil(packing_fraction_ * space_->volume / min_vol);
  }
  members_.reserve(max_insert);
  Logger::Debug("Reserving memory for %d members in RigidFilamentSpecies",
                max_insert);
}

void RigidFilamentSpecies::AddMember() {
  Species::AddMember();
  if (packing_fraction_ > 0) {
    double vol = members_.back().GetVolume();
    fill_volume_ += vol;
    /* if we are still short on volume for the target packing fraction, then
       request more members */
    if (fill_volume_ < packing_fraction_ * space_->volume &&
        members_.size() == sparams_.num) {
      sparams_.num++;
    }
  }
}

void RigidFilamentSpecies::PopMember() {
  if (packing_fraction_ > 0) {
    double vol = members_.back().GetVolume();
    fill_volume_ -= vol;
  }
  Species::PopMember();
}

