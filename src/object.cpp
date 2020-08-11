#include "cglass/object.hpp"

Object::Object(unsigned long seed) : rng_(seed) {
  // Initialize object ID, guaranteeing thread safety
  InitOID();
  // Set some defaults
  std::fill(position_, position_ + 3, 0.0);
  std::fill(prev_position_, prev_position_ + 3, 0.0);
  std::fill(prev_orientation_, prev_orientation_ + 3, 0.0);
  std::fill(scaled_position_, scaled_position_ + 3, 0.0);
  std::fill(orientation_, orientation_ + 3, 0.0);
  orientation_[n_dim_ - 1] = 1.0;
  prev_orientation_[n_dim_ - 1] = 1.0;
  std::fill(pos, pos + 3, 0.0);
  std::fill(direction, direction + 3, 0.0);
  std::fill(force_, force_ + 3, 0.0);
  std::fill(torque_, torque_ + 3, 0.0);
  std::fill(dr_zero_, dr_zero_ + 3, 0.0);
  draw_ = draw_type::orientation;
  type_ = obj_type::generic;
  color_ = 0;
  diameter_ = 1;
  length_ = 0;
  p_energy_ = 0;
  interacting_ = true;
  is_mesh_ = false;
  dr_tot_ = 0;
  mesh_id_ = 0;
  polar_order_ = 0;
  contact_number_ = 0;
  n_contact_ = 0;
  has_overlap_ = false;
  gid = oid_;
  interactor_update_ = false;
}

// Set the object OID in a thread-safe way
void Object::InitOID() {
  std::lock_guard<std::mutex> lk(_obj_mtx_);
  SetOID(++_next_oid_);
}

int Object::_next_oid_ = 0;
std::mutex Object::_obj_mtx_;
system_parameters *Object::params_ = nullptr;
space_struct *Object::space_ = nullptr;
int Object::n_dim_ = 0;
double Object::delta_ = 0;

void Object::SetParams(system_parameters *params) { params_ = params; }
void Object::SetSpace(space_struct *space) { space_ = space; }
void Object::SetNDim(int n_dim) { n_dim_ = n_dim; }
void Object::SetDelta(double delta) { delta_ = delta; }
const double Object::GetDelta() { return delta_; }
void Object::SetNextOID(const int next_oid) { _next_oid_ = next_oid; }
const int Object::GetNextOID() { return _next_oid_; }
// Trivial Get/Set functions
int const Object::GetOID() const { return oid_; }
void Object::SetPosition(double const *const new_pos) {
  std::copy(new_pos, new_pos + n_dim_, position_);
}
void Object::SetScaledPosition(double const *const spos) {
  std::copy(spos, spos + n_dim_, scaled_position_);
}
void Object::SetOrientation(double const *const u) {
  std::copy(u, u + n_dim_, orientation_);
}
void Object::SetPrevPosition(double const *const ppos) {
  std::copy(ppos, ppos + n_dim_, prev_position_);
}
void Object::SetPrevOrientation(double const *const pu) {
  std::copy(pu, pu + n_dim_, prev_orientation_);
}
void Object::ResetPreviousPosition() {
  SetPosition(GetPrevPosition());
  SetOrientation(GetPrevOrientation());
  UpdatePeriodic();
}
void Object::SetDiameter(double new_diameter) { diameter_ = new_diameter; }
void Object::SetLength(double new_length) { length_ = new_length; }
void Object::AddForce(double const *const f) {
  for (int i = 0; i < n_dim_; ++i)
    force_[i] += f[i];
}
void Object::SubForce(double const *const f) {
  for (int i = 0; i < n_dim_; ++i)
    force_[i] -= f[i];
}
void Object::SetForce(double const *const f) {
  std::copy(f, f + n_dim_, force_);
}
void Object::AddTorque(double const *const t) {
  for (int i = 0; i < 3; ++i)
    torque_[i] += t[i];
}
void Object::SubTorque(double const *const t) {
  for (int i = 0; i < 3; ++i)
    torque_[i] -= t[i];
}
void Object::SetTorque(double const *const t) { std::copy(t, t + 3, torque_); }
void Object::AddPotential(double const p) { p_energy_ += p; }
void Object::AddPolarOrder(double const po) {
  std::lock_guard<std::mutex> lk(_obj_mtx_);
  polar_order_ += po;
}
void Object::AddContactNumber(double const cn) {
  std::lock_guard<std::mutex> lk(_obj_mtx_);
  contact_number_ += cn;
}
void Object::CalcPolarOrder() {
  polar_order_ = 0;
  contact_number_ = 0;
  for (auto ix = ixs_.begin(); ix != ixs_.end(); ++ix) {
    // Ignore intrafilament bonds
    if (ix->first->obj1->GetMeshID() == ix->first->obj2->GetMeshID())
      continue;
    // if (ix->first->pause_interaction)
    // continue;
    double dr2 = ix->first->dr_mag2;
    if (dr2 < 0) {
      Logger::Error("Object received an interaction that was not parsed");
    }
    double expdr2 = exp(-dr2);
    double const *const u1 = ix->first->obj1->GetInteractorOrientation();
    double const *const u2 = ix->first->obj2->GetInteractorOrientation();
    double u1_dot_u2 = dot_product(n_dim_, u1, u2);
    polar_order_ += u1_dot_u2 * expdr2;
    contact_number_ += expdr2;
  }

  if (contact_number_ > 1e-16) {
    polar_order_ /= contact_number_;
  } else {
    polar_order_ = 0;
  }
}
void Object::ZeroPolarOrder() {
  contact_number_ = 0;
  polar_order_ = 0;
}
void Object::SetInteractor(bool ix) { interacting_ = ix; }
double const *const Object::GetPosition() { return position_; }
double const *const Object::GetPrevPosition() { return prev_position_; }
double const *const Object::GetPrevOrientation() { return prev_orientation_; }
double const *const Object::GetScaledPosition() { return scaled_position_; }
double const *const Object::GetOrientation() { return orientation_; }
double const *const Object::GetForce() { return force_; }
double const *const Object::GetTorque() { return torque_; }
double const Object::GetDiameter() { return diameter_; }
double const Object::GetLength() { return length_; }
double const Object::GetPotentialEnergy() { return p_energy_; }
double const Object::GetPolarOrder() { return polar_order_; }
double const Object::GetContactNumber() { return contact_number_; }
bool const Object::IsInteractor() { return interacting_; }
bool const Object::IsMesh() { return is_mesh_; }
bool const Object::CheckInteractorUpdate() {
  if (interactor_update_) {
    interactor_update_ = false;
    return true;
  }
  return false;
}
void Object::HasOverlap(bool overlap) { has_overlap_ = overlap; }
int const Object::GetMeshID() const { return mesh_id_; }
void Object::SetMeshID(int mid) { mesh_id_ = mid; }
void Object::SetOID(int oid) { oid_ = oid; }
void Object::ToggleIsMesh() { is_mesh_ = !is_mesh_; }
obj_type const Object::GetType() { return type_; }
species_id const Object::GetSID() { return sid_; }
void Object::SetType(obj_type type) { type_ = type; }
void Object::SetSID(species_id sid) { sid_ = sid; }

// Virtual functions
void Object::InsertRandom(double buffer) {
  if (buffer < 0) {
    buffer = diameter_;
  }
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  Logger::Trace("Inserting object %d randomly", GetOID());
  rng_.RandomCoordinate(space_, pos, buffer);
  rng_.RandomUnitVector(n_dim_, u);
  InsertAt(pos, u);
}

void Object::InsertRandomOriented(double const *const u) {
  InsertRandom();
  SetOrientation(u);
  normalize_vector(orientation_, n_dim_);
}

void Object::InsertAt(double const *const new_pos, double const *const u) {
  SetPosition(new_pos);
  SetOrientation(u);
  normalize_vector(orientation_, n_dim_);
  UpdatePeriodic();
  Logger::Trace("Object %d inserted at [%2.2f, %2.2f, %2.2f] with orientation "
                "[%2.2f %2.2f %2.2f]",
                GetOID(), position_[0], position_[1], position_[2],
                orientation_[0], orientation_[1], orientation_[2]);
}

void Object::ZeroForce() {
  std::fill(force_, force_ + 3, 0.0);
  std::fill(torque_, torque_ + 3, 0.0);
  p_energy_ = 0.0;
}

void Object::Draw(std::vector<graph_struct *> &graph_array) {
  std::copy(scaled_position_, scaled_position_ + 3, g_.r);
  for (int i = space_->n_periodic; i < n_dim_; ++i) {
    g_.r[i] = position_[i];
    printf("position_[i] = %f\n", position_[i]);
  }
  std::copy(orientation_, orientation_ + 3, g_.u);
  g_.color = color_;
  g_.diameter = diameter_;
  g_.length = length_;
  g_.draw = draw_;
  graph_array.push_back(&g_);
}

// Updates scaled position, leaving position fixed
void Object::UpdatePeriodic() {
  double s[3], r[3];
  std::copy(scaled_position_, scaled_position_ + 3, s);
  periodic_boundary_conditions(space_->n_dim, space_->n_periodic,
                               space_->unit_cell, space_->unit_cell_inv,
                               position_, s);
  SetScaledPosition(s);
  UpdateKMC();
}

void Object::UpdateKMC() {
  radius = 0.5 * diameter_;
  length = length_;
  for (int i = 0; i < space_->n_periodic; ++i) {
    // pos[i] = space_->unit_cell[n_dim_*i+i]*scaled_position_[i];
    pos[i] = scaled_position_[i];
  }
  for (int i = space_->n_periodic; i < n_dim_; ++i) {
    pos[i] = position_[i];
  }
  for (int i = 0; i < n_dim_; ++i) {
    direction[i] = orientation_[i];
  }
}

void Object::SetColor(double const c, draw_type dtype) {
  color_ = c;
  draw_ = dtype;
}
void Object::ScalePosition() {
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = 0;
    for (int j = 0; j < n_dim_; ++j) {
      position_[i] += space_->unit_cell[n_dim_ * i + j] * scaled_position_[j];
    }
  }
}
int Object::GetCount() { return 1; }
void Object::GetInteractors(std::vector<Object *> &ix) {
  ix.insert(ix.end(), interactors_.begin(), interactors_.end());
}
double const *const Object::GetInteractorPosition() { return GetPosition(); }
double const *const Object::GetInteractorPrevPosition() {
  return GetPrevPosition();
}
double const *const Object::GetInteractorScaledPosition() {
  return GetScaledPosition();
}
double const *const Object::GetInteractorOrientation() {
  return GetOrientation();
}
double const Object::GetInteractorDiameter() { return GetDiameter(); }
double const Object::GetInteractorLength() { return GetLength(); }
double const Object::GetVolume() {
  if (n_dim_ == 2) {
    if (length_ == 0) {
      return 0.25 * M_PI * diameter_ * diameter_;
    } else {
      return length_ * diameter_ + 0.25 * M_PI * diameter_ * diameter_;
    }
  }
  if (n_dim_ == 3) {
    if (length_ == 0) {
      return 1.0 / 6.0 * M_PI * diameter_ * diameter_ * diameter_;
    } else {
      return 1.0 / 6.0 * M_PI * diameter_ * diameter_ * diameter_ +
             0.25 * M_PI * length_ * diameter_ * diameter_;
    }
    return -1;
  }
  return -1;
}
double const Object::GetDrTot() {
  UpdateDrTot();
  return dr_tot_;
}
void Object::ZeroDrTot() {
  std::copy(position_, position_ + 3, dr_zero_);
  dr_tot_ = 0;
}
void Object::UpdateDrTot() {
  for (int i = 0; i < n_dim_; ++i) {
    double dr = position_[i] - dr_zero_[i];
    dr_tot_ += dr * dr;
  }
}
bool Object::HasNeighbor(int other_oid) {
  // Generic objects are not assumed to have neighbors
  return false;
}
void Object::GiveInteraction(object_interaction ix) {
  std::lock_guard<std::mutex> lk(_obj_mtx_);
  ixs_.push_back(ix);
}

void Object::ApplyInteractions() {
  for (auto ix = ixs_.begin(); ix != ixs_.end(); ++ix) {
    if (ix->first->pause_interaction)
      continue;
    if (ix->second) {
      AddForce(ix->first->force);
      AddTorque(ix->first->t1);
      AddPotential(ix->first->pote);
    } else {
      SubForce(ix->first->force);
      SubTorque(ix->first->t2);
      AddPotential(ix->first->pote);
    }
  }
  /* Moved clearing responsibility to global ClearObjectInteractions() function
   in InteractionManager so Objects have an opportunity to inspect their
   interactions (for analysis or other reasons) */
  // ixs_.clear();
}

void Object::FlagDuplicateInteractions() {
  int n_interactions = ixs_.size();
  for (int i = 0; i < n_interactions - 1; ++i) {
    int other_mid = ixs_[i].second ? ixs_[i].first->obj1->GetMeshID()
                                   : ixs_[i].first->obj2->GetMeshID();
    for (int j = i + 1; j < n_interactions; ++j) {
      int other_mid2 = ixs_[j].second ? ixs_[j].first->obj1->GetMeshID()
                                      : ixs_[j].first->obj2->GetMeshID();

      if (other_mid == other_mid2) {
        if (ixs_[i].first->dr_mag2 < ixs_[j].first->dr_mag2) {
          ixs_[j].first->pause_interaction = true;
        } else {
          ixs_[i].first->pause_interaction = true;
        }
      }
    }
  }
}

void Object::GetInteractions(std::vector<object_interaction> &ixs) {
  ixs.insert(ixs.end(), ixs_.begin(), ixs_.end());
}

void Object::ClearInteractions() { ixs_.clear(); }
void Object::Cleanup() {}

// Object I/O functions
void Object::WriteCheckpoint(std::fstream &ocheck) {
  WriteCheckpointHeader(ocheck);
  WriteSpec(ocheck);
}
void Object::WriteCheckpointHeader(std::fstream &ocheck) {
  void *rng_state = rng_.GetState();
  size_t rng_size = rng_.GetSize();
  int oid = GetOID();
  int mid = GetMeshID();
  ocheck.write(reinterpret_cast<char *>(&oid), sizeof(int));
  ocheck.write(reinterpret_cast<char *>(&mid), sizeof(int));
  ocheck.write(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  ocheck.write(reinterpret_cast<char *>(rng_state), rng_size);
}

void Object::ReadCheckpoint(std::fstream &icheck) {
  ReadCheckpointHeader(icheck);
  ReadSpec(icheck);
}

void Object::ReadCheckpointHeader(std::fstream &icheck) {
  if (icheck.eof())
    return;
  void *rng_state = rng_.GetState();
  size_t rng_size;
  int oid;
  int mid;
  icheck.read(reinterpret_cast<char *>(&oid), sizeof(int));
  icheck.read(reinterpret_cast<char *>(&mid), sizeof(int));
  icheck.read(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  icheck.read(reinterpret_cast<char *>(rng_state), rng_size);
  SetOID(oid);
  SetMeshID(mid);
}

void Object::WritePosit(std::fstream &oposit) {
  for (auto &posit : position_)
    oposit.write(reinterpret_cast<char *>(&posit), sizeof(posit));
  for (auto &spos : scaled_position_)
    oposit.write(reinterpret_cast<char *>(&spos), sizeof(spos));
  for (auto &u : orientation_)
    oposit.write(reinterpret_cast<char *>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char *>(&length_), sizeof(length_));
}

void Object::ReadPosit(std::fstream &iposit) {
  if (iposit.eof())
    return;
  for (auto &posit : position_)
    iposit.read(reinterpret_cast<char *>(&posit), sizeof(posit));
  for (auto &spos : scaled_position_)
    iposit.read(reinterpret_cast<char *>(&spos), sizeof(spos));
  for (auto &u : orientation_)
    iposit.read(reinterpret_cast<char *>(&u), sizeof(u));
  iposit.read(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char *>(&length_), sizeof(length_));
}
void Object::WriteSpec(std::fstream &ospec) { WritePosit(ospec); }
void Object::ReadSpec(std::fstream &ispec) { ReadPosit(ispec); }
void Object::ReadPositFromSpec(std::fstream &ispec) { ReadPosit(ispec); }
void Object::GetAvgPosition(double *ap) {
  std::copy(position_, position_ + 3, ap);
}
void Object::GetAvgOrientation(double *au) {
  std::copy(orientation_, orientation_ + 3, au);
}
void Object::SetAvgPosition() {
  // Nothing to be done
}

void Object::Report() {
  fprintf(stderr, "    OID: %d\n", GetOID());
  fprintf(stderr, "      r: {%2.2f %2.2f %2.2f}\n", position_[0], position_[1],
          position_[2]);
  fprintf(stderr, "      u: {%2.2f %2.2f %2.2f}\n", orientation_[0],
          orientation_[1], orientation_[2]);
  fprintf(stderr, "      l: %2.2f\n", length_);
  fprintf(stderr, "      d: %2.2f\n", diameter_);
}
