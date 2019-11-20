#include "anchor.hpp"

Anchor::Anchor(unsigned long seed) : Object(seed) {
  SetSID(species_id::crosslink);
}

void Anchor::Init(crosslink_parameters *sparams) {
  sparams_ = sparams;
  diameter_ = sparams_->diameter;
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  static_flag_ = false;  // Must be explicitly set to true by Crosslink
  Unbind();
  walker_ = sparams_->walker_flag;
  step_direction_ =
      (sparams_->step_direction == 0 ? 0 : SIGNOF(sparams_->step_direction));
  velocity_ = sparams_->velocity;
  max_velocity_ = velocity_;
  k_off_ = sparams_->k_off;
  end_pausing_ = sparams_->end_pausing;
  diffuse_ = sparams_->diffusion_flag;
  f_stall_ = sparams_->f_stall;
  force_dep_vel_flag_ = sparams_->force_dep_vel_flag;
  SetDiffusion();
}

double const Anchor::GetMeshLambda() { return mesh_lambda_; }

double const Anchor::GetBondLambda() { return bond_lambda_; }

void Anchor::SetBondLambda(double l) { bond_lambda_ = l; }
void Anchor::SetMeshLambda(double ml) { mesh_lambda_ = ml; }

void Anchor::SetDiffusion() { diffusion_ = sqrt(24.0 * diameter_ / delta_); }

void Anchor::SetWalker(int dir, double walk_v) {
  if (ABS(dir) != 1) {
    Logger::Error("Walker direction must be set to +/- 1");
  }
  walker_ = true;
  velocity_ = walk_v;
  max_velocity_ = velocity_;
  step_direction_ = dir;
}

void Anchor::UpdateAnchorPositionToMesh() {
  if (!bound_ || static_flag_) return;
  if (!bond_ || !mesh_) {
    Logger::Error("Anchor tried to update position to nullptr bond or mesh");
  }

  /* Use the mesh to determine the bond lengths. The true bond lengths fluctuate
     about this, but should be considered approximations to the ideal mesh. */
  mesh_length_ = mesh_->GetTrueLength();
  /* Use current position along mesh (mesh_lambda) to determine whether the
     anchor fell off the mesh due to dynamic instability */
  if (!CheckMesh()) return;
  // Now figure out which bond we are on in the mesh according to mesh_lambda
  bond_ = mesh_->GetBondAtLambda(mesh_lambda_);

  // Figure out how far we are from the bond tail: bond_lambda
  if (!CalcBondLambda()) {
    return;
  }
  // Update anchor position with respect to bond
  UpdateAnchorPositionToBond();
}

bool Anchor::CalcBondLambda() {
  if (!bond_) {
    Logger::Error(
        "Attempted to calculate bond lambda when not attached to"
        " bond!");
  }
  bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
  bond_length_ = bond_->GetLength();
  if (bond_lambda_ < 0) {
    Bond *bond = bond_->GetNeighborBond(0);
    if (bond) {
      bond_ = bond;
      bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
      bond_length_ = bond_->GetLength();
    } else if (end_pausing_) {
      bond_lambda_ = 0;
    } else {
      Unbind();
      return false;
    }
  } else if (bond_lambda_ > bond_length_) {
    printf("bond_length_ = %f\n", bond_length_);
    Bond *bond = bond_->GetNeighborBond(1);
    if (bond) {
      bond_ = bond;
      bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
      bond_length_ = bond_->GetLength();
    } else if (end_pausing_) {
      bond_lambda_ = bond_length_;
    } else {
      Unbind();
      return false;
    }
  }
  // assert(bond_lambda_ >= 0 && bond_lambda_ <= bond_length_);
  if (bond_lambda_ < -1e-6 || bond_lambda_ > bond_length_ + 1e-6) {
    Logger::Error(
        "Bond lambda out of expected range in UpdateAnchorPositionToMesh, "
        "bond_num: %d, mesh lambda: %2.8f, mesh length: %2.8f, bond lambda: "
        "%2.8f, bond length: %2.8f",
        bond_->GetBondNumber(), mesh_lambda_, mesh_length_, bond_lambda_,
        bond_length_);
  }
  return true;
}
void Anchor::UpdatePosition() {
  // Currently only bound anchors diffuse/walk (no explicit unbound anchors)
  if (!bound_ || static_flag_ || (!diffuse_ && !walker_)) {
    return;
  }
  // Diffuse or walk along the mesh, updating mesh_lambda
  if (diffuse_) {
    Diffuse();
    if (!CheckMesh()) return;
  }
  if (walker_) {
    Walk();
    CheckMesh();
  }
}

void Anchor::ApplyAnchorForces() {
  if (!bound_ || static_flag_) {
    return;
  }
  if (!bond_) {
    Logger::Error("Anchor attempted to apply forces to nullptr bond");
  }
  bond_->AddForce(force_);
  double dlambda[3] = {0};
  for (int i = 0; i < n_dim_; ++i) {
    dlambda[i] = (bond_lambda_ - 0.5 * bond_length_) * orientation_[i];
  }
  cross_product(dlambda, force_, torque_, 3);
  bond_->AddTorque(torque_);
}

void Anchor::Activate() {
  active_ = true;
  step_direction_ = -step_direction_;
}

void Anchor::Deactivate() {
  active_ = false;
  step_direction_ = -step_direction_;
}

void Anchor::Walk() {
  if (force_dep_vel_flag_) {
    double fmag = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      fmag += force_[i] * force_[i];
    }
    fmag = sqrt(fmag);
    // Linear force-velocity relationship
    double fdep = 1 - fmag / f_stall_;
    if (fdep > 1) {
      fdep = 1;
    } else if (fdep < 0) {
      fdep = 0;
    }
    velocity_ = max_velocity_ * fdep;
  }
  double dr = step_direction_ * velocity_ * delta_;
  mesh_lambda_ += dr;
}

// Check that the anchor is still located on the filament mesh
// Returns true if anchor is still on the mesh, false otherwise
bool Anchor::CheckMesh() {
  // Check if we moved off the mesh tail
  if (mesh_lambda_ < 0) {
    // Stick to mesh ends if we have end_pausing
    if (end_pausing_) {
      mesh_lambda_ = 0;
    } else {
      // Otherwise, unbind
      Unbind();
      return false;
    }
  } else if (mesh_lambda_ > mesh_length_) {
    // Same thing for walking off the head
    if (end_pausing_) {
      mesh_lambda_ = mesh_length_;
    } else {
      Unbind();
      return false;
    }
  }
  return true;
}

void Anchor::Unbind() {
  if (static_flag_) {
    Logger::Error("Static anchor attempted to unbind");
  }
  bound_ = false;
  bond_ = nullptr;
  mesh_ = nullptr;
  mesh_n_bonds_ = -1;
  bond_length_ = -1;
  bond_lambda_ = -1;
  mesh_lambda_ = -1;
  active_ = false;
  ClearNeighbors();
  ZeroForce();
  SetMeshID(-1);
  std::fill(position_, position_ + 3, 0.0);
  std::fill(orientation_, orientation_ + 3, 0.0);
}

void Anchor::Diffuse() {
  double kick = rng_.RandomUniform() - 0.5;
  double dr = kick * diffusion_ * delta_ / diameter_;
  mesh_lambda_ += dr;
}

void Anchor::UpdateAnchorPositionToBond() {
  if (!bond_) {
    Logger::Error("Anchor tried to update position relative to nullptr bond");
  }
  double const *const bond_position = bond_->GetPosition();
  double const *const bond_orientation = bond_->GetOrientation();
  for (int i = 0; i < n_dim_; ++i) {
    orientation_[i] = bond_orientation[i];
    position_[i] = bond_position[i] -
                   (0.5 * bond_length_ - bond_lambda_) * bond_orientation[i];
  }
  UpdatePeriodic();
}

void Anchor::Draw(std::vector<graph_struct *> &graph_array) {
  if (!bound_) return;
  std::copy(scaled_position_, scaled_position_ + 3, g_.r);
  for (int i = space_->n_periodic; i < n_dim_; ++i) {
    g_.r[i] = position_[i];
  }
  std::copy(orientation_, orientation_ + 3, g_.u);
  g_.color = color_;
  g_.diameter = diameter_;
  g_.length = length_;
  g_.draw = draw_;
  graph_array.push_back(&g_);
}

void Anchor::AttachObjRandom(Object *o) {
  double length = o->GetLength();
  double lambda = length * rng_.RandomUniform();
  AttachObjLambda(o, lambda);
}

void Anchor::AttachObjLambda(Object *o, double lambda) {
  if (o->GetType() != +obj_type::bond) {
    Logger::Error(
        "Crosslink binding to non-bond objects not yet implemented in "
        "AttachObjLambda.");
  }
  bond_ = dynamic_cast<Bond *>(o);
  if (bond_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a bond!");
  }
  mesh_ = dynamic_cast<Mesh *>(bond_->GetMeshPtr());
  if (mesh_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
  }
  mesh_n_bonds_ = mesh_->GetNBonds();
  bond_length_ = bond_->GetLength();
  mesh_length_ = mesh_->GetTrueLength();
  bond_lambda_ = lambda;

  if (bond_lambda_ < 0 || bond_lambda_ > bond_length_) {
    printf("bond_lambda: %2.2f\n", bond_lambda_);
    Logger::Error(
        "Lambda passed to anchor does not match length of "
        "corresponding bond! lambda: %2.2f, bond_length: %2.2f ",
        bond_lambda_, bond_length_);
  }

  /* Distance anchor is relative to entire mesh length */
  mesh_lambda_ = bond_->GetMeshLambda() + bond_lambda_;
  SetMeshID(bond_->GetMeshID());
  UpdateAnchorPositionToBond();
  ZeroDrTot();
  bound_ = true;
}

void Anchor::AttachObjMeshLambda(Object *o, double mesh_lambda) {
  if (o->GetType() != +obj_type::bond) {
    Logger::Error(
        "Crosslink binding to non-bond objects not yet implemented in "
        "AttachObjMeshLambda.");
  }
  bond_ = dynamic_cast<Bond *>(o);
  if (bond_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a bond!");
  }
  mesh_ = dynamic_cast<Mesh *>(bond_->GetMeshPtr());
  if (mesh_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
  }
  Logger::Trace("Attaching anchor %d to mesh %d", GetOID(), mesh_->GetMeshID());

  bound_ = true;
  mesh_lambda_ = mesh_lambda;
  mesh_n_bonds_ = -1;
  UpdateAnchorPositionToMesh();
  if (!bound_) {
    Logger::Error(
        "Updating anchor to mesh from checkpoint resulted in an unbound "
        "anchor");
  }
  SetMeshID(bond_->GetMeshID());
  ZeroDrTot();
}

void Anchor::BindToPosition(double *bind_pos) {
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = bind_pos[i];
  }
  UpdatePeriodic();
}

bool Anchor::IsBound() { return bound_; }

int const Anchor::GetBoundOID() {
  if (!bond_) {
    return -1;
  }
  return bond_->GetOID();
}

/* Temporary function for setting bound state for singly-bound crosslinks,
   in order to get them to draw while not technically bound to a bond
   ( e.g. bond_ -> null ) */
void Anchor::SetBound() { bound_ = true; }

void Anchor::AddNeighbor(Object *neighbor) { neighbors_.AddNeighbor(neighbor); }

void Anchor::ClearNeighbors() { neighbors_.Clear(); }

const Object *const *Anchor::GetNeighborListMem() {
  return neighbors_.GetNeighborListMem();
}

Object *Anchor::GetNeighbor(int i_neighbor) {
  return neighbors_.GetNeighbor(i_neighbor);
}

const int Anchor::GetNNeighbors() const { return neighbors_.NNeighbors(); }

void Anchor::WriteSpec(std::fstream &ospec) {
  ospec.write(reinterpret_cast<char *>(&bound_), sizeof(bool));
  ospec.write(reinterpret_cast<char *>(&active_), sizeof(bool));
  ospec.write(reinterpret_cast<char *>(&static_flag_), sizeof(bool));
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  ospec.write(reinterpret_cast<char *>(&mesh_lambda_), sizeof(double));
}

void Anchor::ReadSpec(std::fstream &ispec) {
  ispec.read(reinterpret_cast<char *>(&bound_), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&active_), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&static_flag_), sizeof(bool));
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  ispec.read(reinterpret_cast<char *>(&mesh_lambda_), sizeof(double));
  UpdatePeriodic();
  if (active_) step_direction_ = -sparams_->step_direction;
}

void Anchor::SetStatic(bool static_flag) { static_flag_ = static_flag; }
