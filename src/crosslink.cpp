#include "cglass/crosslink.hpp"
#include <iostream>

Crosslink::Crosslink(unsigned long seed) : Object(seed) {
  SetSID(species_id::crosslink);
}

void Crosslink::Init(crosslink_parameters *sparams) {
  sparams_ = sparams;
  length_ = -1;
  diameter_ = sparams_->tether_diameter;
  color_ = sparams_->tether_color;
  draw_ = draw_type::_from_string(sparams_->tether_draw_type.c_str());
  rest_length_ = sparams_->rest_length;
  static_flag_ = sparams_->static_flag;
  k_spring_ = sparams_->k_spring;
  k_spring_compress_ =
      sparams_->k_spring_compress; // spring const for compression if
                                   // asymmetric_spring_flag is true
  // asymmetric_spring_flag_ = sparams_->asymmetric_spring_flag;
  static_flag_ = sparams_->static_flag;
  k_align_ = sparams_->k_align;
  // rcapture_ = sparams_->r_capture;
  bind_site_density_ = sparams_->bind_site_density;
  e_dep_factor_ = sparams_->energy_dep_factor;
  fdep_length_ = sparams_->force_dep_length;
  polar_affinity_ = sparams_->polar_affinity;
  /* TODO generalize crosslinks to more than two anchors */
  Anchor anchor1(rng_.GetSeed());
  Anchor anchor2(rng_.GetSeed());
  anchors_.push_back(anchor1);
  anchors_.push_back(anchor2);
  anchors_[0].Init(sparams_);
  anchors_[1].Init(sparams_);
  SetSingly();
  Logger::Trace("Initializing crosslink %d with anchors %d and %d", GetOID(),
                anchors_[0].GetOID(), anchors_[1].GetOID());
}

void Crosslink::InitInteractionEnvironment(LookupTable *lut) { lut_ = lut; }

/* Function used to set anchor[0] position etc to xlink position etc */
void Crosslink::UpdatePosition() {}

void Crosslink::GetAnchors(std::vector<Object *> &ixors) {
  if (IsUnbound())
    return;
  if (!static_flag_) {
    ixors.push_back(&anchors_[0]);
  }
  if (IsDoubly()) {
    ixors.push_back(&anchors_[1]);
  }
}

/* Perform kinetic monte carlo step of protein with 1 head attached. */
void Crosslink::SinglyKMC() {
  double roll = rng_.RandomUniform();
  int head_bound = 0;
  // Set up KMC objects and calculate probabilities
  double unbind_prob = anchors_[0].GetOffRate() * delta_;
  if (static_flag_) {
    unbind_prob = 0;
  }
  /* Must populate filter with 1 for every neighbor, since KMC
  expects a mask. We already guarantee uniqueness, so we won't overcount. */
  int n_neighbors = anchors_[0].GetNNeighbors();
  std::vector<int> kmc_filter(n_neighbors, 1);
  /* Initialize KMC calculation */
  KMC<Object> kmc_bind(anchors_[0].pos, n_neighbors, delta_, lut_);

  /* Initialize periodic boundary conditions */
  kmc_bind.SetPBCs(n_dim_, space_->n_periodic, space_->unit_cell);

  /* Calculate probability to bind */
  double kmc_bind_prob = 0;
  double k_on_d_prime = anchors_[1].GetOnRate() * bind_site_density_;

  std::vector<double> kmc_bind_factor(n_neighbors, k_on_d_prime);
  if (n_neighbors > 0) {
    if (!static_flag_ && polar_affinity_ < 1) {
      anchors_[0].CalculatePolarAffinity(kmc_bind_factor);
    }
    kmc_bind.LUCalcTotProbsSD(anchors_[0].GetNeighborListMem(), kmc_filter,
                              anchors_[0].GetBoundOID(), kmc_bind_factor);
    kmc_bind_prob = kmc_bind.getTotProb();
  }
  // Find out whether we bind, unbind, or neither.
  int head_activate = choose_kmc_double(unbind_prob, kmc_bind_prob, roll);
  // Change status of activated head
  if (head_activate == 0) {
    // Unbind bound head
    anchors_[0].Unbind();
    SetUnbound();
    Logger::Trace("Crosslink %d came unbound", GetOID());
  } else if (head_activate == 1) {
    // Bind unbound head
    /* Position on rod where protein will bind with respect to center of rod,
     * passed by reference */
    double bind_lambda;
    /* Find out which rod we are binding to */
    int i_bind = kmc_bind.whichRodBindSD(bind_lambda, roll);
    if (i_bind < 0) { // || bind_lambda < 0) {
      printf("i_bind = %d\nbind_lambda = %2.2f\n", i_bind, bind_lambda);
      Logger::Error("kmc_bind.whichRodBindSD in Crosslink::SinglyKMC"
                    " returned an invalid result!");
    }
    Object *bind_obj = anchors_[0].GetNeighbor(i_bind);
    double obj_length = bind_obj->GetLength();
    /* KMC returns bind_lambda to be with respect to center of rod. We want it
       to be specified from the tail of the rod to be consistent */
    bind_lambda += 0.5 * obj_length;
    /* KMC can return values that deviate a very small amount from the true rod
       length. Bind to ends if lambda < 0 or lambda > bond_length. */
    if (bind_lambda > obj_length) {
      bind_lambda = obj_length;
    } else if (bind_lambda < 0) {
      bind_lambda = 0;
    }
    anchors_[1].AttachObjLambda(bind_obj, bind_lambda);
    SetDoubly();
    Logger::Trace("Crosslink %d became doubly bound to obj %d", GetOID(),
                  bind_obj->GetOID());
  }
  kmc_filter.clear();
}

/* Perform kinetic monte carlo step of protein with 2 heads of protein
 * object attached. */
void Crosslink::DoublyKMC() {
  /* Calculate force-dependent unbinding for each head */
  double tether_stretch = length_ - rest_length_;
  // For asymmetric springs apply second spring constant for compression
  double e_dep = e_dep_factor_;
  double f_dep = fdep_length_;
  if (k_spring_compress_ >= 0 && tether_stretch < 0) {
    e_dep *= 0.5 * k_spring_compress_ * SQR(tether_stretch);
    f_dep *= k_spring_compress_ * tether_stretch;
  } else {
    e_dep *= 0.5 * k_spring_ * SQR(tether_stretch);
    f_dep *= k_spring_ * tether_stretch;
  }
  double unbind_prob = anchors_[0].GetOffRate() * delta_ * exp(e_dep + f_dep);
  double roll = rng_.RandomUniform();
  int head_activate = -1;
  if (static_flag_) {
    head_activate = choose_kmc_double(0, unbind_prob, roll);
  } else {
    // Each head has same probability of undbinding.
    // Probability of unbinding follows a poisson process but assume that only
    // one head can unbind during a time step.
    head_activate = choose_kmc_double(unbind_prob, unbind_prob, roll);
  }
  if (head_activate == 0) {
    Logger::Trace("Doubly-bound crosslink %d came unbound from %d", GetOID(),
                  anchors_[0].GetBoundOID());
    anchors_[0] = anchors_[1];
    anchors_[1].Unbind();
    SetSingly();
  } else if (head_activate == 1) {
    Logger::Trace("Doubly-bound crosslink %d came unbound from %d", GetOID(),
                  anchors_[1].GetBoundOID());
    anchors_[1].Unbind();
    SetSingly();
  }
}

void Crosslink::CalculateBinding() {
  if (IsSingly()) {
    SinglyKMC();
  } else if (IsDoubly()) {
    DoublyKMC();
  }
  ClearNeighbors();
}

/* Only singly-bound crosslinks interact */
void Crosslink::GetInteractors(std::vector<Object *> &ixors) {
  ClearNeighbors();
  if (IsSingly()) {
    ixors.push_back(&anchors_[0]);
  }
}

void Crosslink::ClearNeighbors() { anchors_[0].ClearNeighbors(); }

void Crosslink::UpdateAnchorsToMesh() {
  anchors_[0].UpdateAnchorPositionToMesh();
  anchors_[1].UpdateAnchorPositionToMesh();
}

void Crosslink::UpdateAnchorPositions() {
  anchors_[0].UpdatePosition();
  anchors_[1].UpdatePosition();
}

void Crosslink::ApplyTetherForces() {
  if (!IsDoubly())
    return;
  anchors_[0].ApplyAnchorForces();
  anchors_[1].ApplyAnchorForces();
}

void Crosslink::UpdateCrosslinkForces() {
  /* Update anchor positions in space to calculate tether forces */
  UpdateAnchorsToMesh();
  /* Check if an anchor became unbound due to diffusion, etc */
  UpdateXlinkState();
  /* If we are doubly-bound, calculate and apply tether forces */
  CalculateTetherForces();
}

void Crosslink::UpdateCrosslinkPositions() {
  /* Have anchors diffuse/walk along mesh */
  UpdateAnchorPositions();
  /* Check if an anchor became unbound do to diffusion, etc */
  UpdateXlinkState();
  /* Check for binding/unbinding events using KMC */
  CalculateBinding();
}

/* This function ensures that singly-bound crosslinks have anchor[0] bound and
   anchor[1] unbound. */
void Crosslink::UpdateXlinkState() {
  if (!anchors_[0].IsBound() && !anchors_[1].IsBound()) {
    SetUnbound();
    return;
  }
  if (IsDoubly() && !anchors_[1].IsBound()) {
    SetSingly();
  } else if (IsDoubly() && !anchors_[0].IsBound()) {
    anchors_[0] = anchors_[1];
    anchors_[1].Unbind();
    SetSingly();
  }
  if (IsSingly() && anchors_[1].IsBound()) {
    SetDoubly();
  }
}

void Crosslink::ZeroForce() {
  std::fill(force_, force_ + 3, 0.0);
  anchors_[0].ZeroForce();
  anchors_[1].ZeroForce();
}

void Crosslink::CalculateTetherForces() {
  ZeroForce();
  if (!IsDoubly())
    return;
  Interaction ix(&anchors_[0], &anchors_[1]);
  MinimumDistance mindist;
  mindist.ObjectObject(ix);
  length_ = sqrt(ix.dr_mag2);
  double stretch = length_ - rest_length_;
  for (int i = 0; i < params_->n_dim; ++i) {
    orientation_[i] = ix.dr[i] / length_;
    position_[i] = ix.midpoint[i];
  }
  // If compression spring constant is set, use that when spring is compressed.
  if (k_spring_compress_ >= 0 && stretch < 0)
    tether_force_ = k_spring_compress_ * stretch;
  else
    tether_force_ = k_spring_ * stretch;

  for (int i = 0; i < params_->n_dim; ++i) {
    force_[i] = tether_force_ * orientation_[i];
  }
  anchors_[0].AddForce(force_);
  anchors_[1].SubForce(force_);
  // Update xlink's position (for drawing)
  UpdatePeriodic();
}

/* Attach a crosslink anchor to object in a random fashion. Currently only
 * bonds are used, but can be extended to sites (sphere-like objects) */
void Crosslink::AttachObjRandom(Object *obj) {
  /* Attaching to random bond implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if (obj->GetType() == +obj_type::bond) {
    anchors_[0].AttachObjRandom(obj);
    SetMeshID(obj->GetMeshID());
    SetSingly();
  } else {
    /* TODO: add binding to sphere or site-like objects */
    Logger::Error("Crosslink binding to non-bond objects not yet implemented.");
  }
}

void Crosslink::Draw(std::vector<graph_struct *> &graph_array) {
  /* Draw anchors */
  anchors_[0].Draw(graph_array);
  anchors_[1].Draw(graph_array);
  /* Draw tether */
  if (IsDoubly() && length_ > 0) {
    std::copy(scaled_position_, scaled_position_ + 3, g_.r);
    for (int i = space_->n_periodic; i < n_dim_; ++i) {
      g_.r[i] = position_[i];
    }
    // std::copy(position_, position_+3, g_.r);
    std::copy(orientation_, orientation_ + 3, g_.u);
    g_.color = color_;
    if (params_->graph_diameter > 0) {
      g_.diameter = params_->graph_diameter;
    } else {
      g_.diameter = diameter_;
    }
    g_.length = length_;
    g_.draw = draw_;
    graph_array.push_back(&g_);
  }
}

void Crosslink::SetDoubly() {
  state_ = bind_state::doubly;
  SetAnchorStates();
}

void Crosslink::SetSingly() {
  state_ = bind_state::singly;
  SetAnchorStates();
}

void Crosslink::SetUnbound() {
  state_ = bind_state::unbound;
  SetAnchorStates();
}

void Crosslink::SetAnchorStates() {
  anchors_[0].SetState(state_);
  anchors_[1].SetState(state_);
};

const bool Crosslink::IsDoubly() const { return state_ == +bind_state::doubly; }
const bool Crosslink::IsSingly() const { return state_ == +bind_state::singly; }
const bool Crosslink::IsUnbound() const {
  return state_ == +bind_state::unbound;
}

void Crosslink::WriteSpec(std::fstream &ospec) {
  if (IsUnbound()) {
    Logger::Error("Unbound crosslink tried to WriteSpec!");
  }
  bool is_doubly = IsDoubly();
  ospec.write(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ospec.write(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ospec.write(reinterpret_cast<char *>(&length_), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  anchors_[0].WriteSpec(ospec);
  anchors_[1].WriteSpec(ospec);
}

void Crosslink::ReadSpec(std::fstream &ispec) {
  if (ispec.eof())
    return;
  SetSingly();
  bool is_doubly;
  ispec.read(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&length_), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  UpdatePeriodic();
  anchors_[0].ReadSpec(ispec);
  anchors_[1].ReadSpec(ispec);
  if (is_doubly) {
    SetDoubly();
  }
}

void Crosslink::WriteCheckpoint(std::fstream &ocheck) {
  Object::WriteCheckpoint(ocheck);
  anchors_[0].WriteCheckpointHeader(ocheck);
  anchors_[1].WriteCheckpointHeader(ocheck);
}

void Crosslink::ReadCheckpoint(std::fstream &icheck) {
  Object::ReadCheckpoint(icheck);
  anchors_[0].ReadCheckpointHeader(icheck);
  anchors_[1].ReadCheckpointHeader(icheck);
  Logger::Trace("Reloading anchor from checkpoint with mid %d",
                anchors_[0].GetMeshID());
  if (IsDoubly()) {
    Logger::Trace("Reloading anchor from checkpoint with mid %d",
                  anchors_[1].GetMeshID());
  }
}

const double Crosslink::GetDrTot() {
  if (IsSingly()) {
    return anchors_[0].GetDrTot();
  } else if (IsDoubly()) {
    double dr1 = anchors_[0].GetDrTot();
    double dr2 = anchors_[1].GetDrTot();
    if (dr1 > dr2) {
      return dr1;
    } else {
      return dr2;
    }
  } else {
    return 0;
  }
}

void Crosslink::ZeroDrTot() {
  anchors_[0].ZeroDrTot();
  if (IsDoubly()) {
    anchors_[1].ZeroDrTot();
  }
}

void Crosslink::InsertAt(double const *const new_pos, double const *const u) {
  static_flag_ = true;
  anchors_[0].InsertAt(new_pos, u);
  anchors_[0].SetBound();
  anchors_[0].SetStatic(true);
  SetSingly();
}

const int Crosslink::GetNNeighbors() const {
  return anchors_[0].GetNNeighbors();
}

const double *const Crosslink::GetPosition() {
  return anchors_[0].GetPosition();
}

const double *const Crosslink::GetOrientation() {
  return anchors_[0].GetOrientation();
}
