
#include "cglass/interaction_manager.hpp"

void InteractionManager::Init(system_parameters *params,
                              std::vector<SpeciesBase *> *species,
                              space_struct *space, bool processing) {
  // Set up pointer structures
  params_ = params;
  species_ = species;
  space_ = space;
  processing_ = processing;

  // Initialize owned structures
  no_init_ = false;
  static_pnumber_ = params_->static_particle_number;
  n_dim_ = params_->n_dim;
  n_periodic_ = params_->n_periodic;
  n_thermo_ = params_->n_thermo;
  no_interactions_ = !(params_->interaction_flag);
  n_objs_ = -1;
  std::fill(stress_, stress_ + 9, 0);

  // Update dr distance should be half the cell length, and we are comparing the
  // squares of the trajectory distances
  xlink_.Init(params_, space_, &ix_objects_);
  no_boundaries_ = false;
  if (space_->type == +boundary_type::none)
    no_boundaries_ = true;
}

void InteractionManager::InitInteractions() {
  if (processing_ && !run_interaction_analysis_) {
    return;
  }
  CellList::SetMinCellLength(xlink_.GetRCutoff());
  potentials_.InitPotentials(params_);
  CellList::SetMinCellLength(sqrt(potentials_.GetRCut2()));

  CellList::Init(params_->n_dim, params_->n_periodic, params_->system_radius);
  Logger::Info("Constructing cell list data structure");
  clist_.BuildCellList();

  dr_update_ = 0.25 * CellList::GetCellLength() * CellList::GetCellLength();
  MinimumDistance::Init(space_, 2 * dr_update_);
}

void InteractionManager::CalculateInteractions() {
  ClearObjectInteractions();
  if (!no_interactions_) {
    CalculatePairInteractions();
  }
  CalculateBoundaryInteractions();
}

void InteractionManager::ApplyInteractions() {
  if (!no_interactions_) {
    ApplyPairInteractions();
  }
  ApplyBoundaryInteractions();
}

/****************************************
  INTERACT: Loop through interactions and
    apply WCA potential (for now)
*****************************************/

void InteractionManager::Interact() {
  // First check if we need to interact
  if (no_interactions_ && no_boundaries_)
    return;
  // Check if we need to update objects in cell list
  CheckUpdateObjects();
  // Update crosslinks
  xlink_.UpdateCrosslinks();
  // Check if we need to update crosslink interactors
  CheckUpdateXlinks();
  // Check if we need to update pair interactions
  CheckUpdateInteractions();
  // Calculate interactions from potentials and assign interactions to objects
  CalculateInteractions();
  /* Check if forces exceed predetermined limits and if so, abort and reduce
     time step */
  if (params_->dynamic_timestep && potentials_.CheckMaxForce()) {
    decrease_dynamic_timestep_ = true;
    return;
  }
  // Apply forces, torques, and energies on objects
  ApplyInteractions();
  i_update_++;
}

void InteractionManager::CheckUpdateXlinks() {
  /* If we need to update crosslinks */
  if (xlink_.CheckUpdate()) {
    interactors_.clear();
    interactors_.insert(interactors_.end(), ix_objects_.begin(),
                        ix_objects_.end());
    // Add crosslinks as interactors
    std::vector<Object *> xlinks;
    xlink_.GetInteractors(xlinks);
    interactors_.insert(interactors_.end(), xlinks.begin(), xlinks.end());
    Logger::Debug("Updating interactions due to crosslink update. %d steps "
                  "since last update",
                  i_update_);
    UpdateInteractions();
    return;
  }
}

void InteractionManager::ForceUpdate() {
  Logger::Trace("Forcing update in interaction manager");
  UpdateInteractors();
  UpdateInteractions();
}

void InteractionManager::UpdateInteractors() {
  Logger::Trace("Updating interactors");
  ix_objects_.clear();
  interactors_.clear();
  for (auto spec_it = species_->begin(); spec_it != species_->end();
       ++spec_it) {
    (*spec_it)->GetInteractors(ix_objects_);
  }
  // Add crosslinks as interactors
  interactors_.insert(interactors_.end(), ix_objects_.begin(),
                      ix_objects_.end());
  std::vector<Object *> xlinks;
  xlink_.GetInteractors(xlinks);
  interactors_.insert(interactors_.end(), xlinks.begin(), xlinks.end());
  Logger::Trace("Updated interactors: %d objects, %d crosslinks, %d total",
                ix_objects_.size(), xlinks.size(), interactors_.size());
}

/* Checks whether or not the given anchor is supposed to be attached to the
   given bond by checking mesh_id of both the bond and the anchor. If they
   match, the anchor is attached to the mesh using the anchor mesh_lambda */
bool InteractionManager::CheckBondAnchorPair(Object *anchor, Object *bond) {
  // Check that the bond and anchor share a mesh_id
  if (anchor->GetMeshID() == bond->GetMeshID()) {
    Anchor *a = dynamic_cast<Anchor *>(anchor);
    if (a == nullptr) {
      Logger::Error("Object pointer was unsuccessfully dynamically cast to an "
                    "Anchor pointer in CheckBondAnchorPair!");
    }
    if (bond->GetType() != +obj_type::bond) {
      Logger::Error(
          "CheckBondAnchorPair expected Bond object pointer, but object"
          " pointer does not have bond obj_type!");
    }
    // Check that the anchor isn't already attached
    if (a->GetBondLambda() < 0) {
      a->AttachObjMeshLambda(bond, a->GetMeshLambda());
      return true;
    }
  }
  return false;
}

void InteractionManager::PairBondCrosslinks() {
  Logger::Trace("Pairing bound crosslinks and objects");
  ix_objects_.clear();
  interactors_.clear();
  // First get object interactors (non-crosslinks)
  for (auto spec_it = species_->begin(); spec_it != species_->end();
       ++spec_it) {
    (*spec_it)->GetInteractors(ix_objects_);
  }
  interactors_.insert(interactors_.end(), ix_objects_.begin(),
                      ix_objects_.end());

  // Add anchors as interactors
  std::vector<Object *> anchors;
  xlink_.GetAnchorInteractors(anchors);
  interactors_.insert(interactors_.end(), anchors.begin(), anchors.end());
  pair_interactions_.clear();
  clist_.RenewObjectsCells(interactors_);
  clist_.MakePairs(pair_interactions_);
  int n_anchors_attached = 0;
  // for (auto ix = anchors.begin(); ix != anchors.end(); ++ix) {
  // for (auto jx = ix_objects_.begin(); jx != ix_objects_.end(); ++jx) {
  for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
       ++ix) {
    Object *obj1 = ix->obj1;
    Object *obj2 = ix->obj2;
    if (obj1->GetSID() == +species_id::crosslink &&
        obj2->GetType() == +obj_type::bond) {
      if (CheckBondAnchorPair(obj1, obj2))
        // if (CheckBondAnchorPair(*ix, *jx))
        n_anchors_attached++;
    } else if (obj2->GetSID() == +species_id::crosslink &&
               obj1->GetType() == +obj_type::bond) {
      if (CheckBondAnchorPair(obj2, obj1)) {
        n_anchors_attached++;
      }
    }
  }
  /* Check that all anchors found their bond attachments */
  if (n_anchors_attached != anchors.size()) {
    Logger::Error("Not all anchors have found their bond! %d/%lu",
                  n_anchors_attached, anchors.size());
  }
}

void InteractionManager::UpdateInteractions() {
  Logger::Trace("Updating interactions");
  i_update_ = 0;
  UpdatePairInteractions();
  UpdateBoundaryInteractions();
  ZeroDrTot();
}

void InteractionManager::UpdatePairInteractions() {
  if (no_interactions_)
    return;
#ifdef TRACE
  int nix = pair_interactions_.size();
#endif
  pair_interactions_.clear();
  clist_.RenewObjectsCells(interactors_);
  clist_.MakePairs(pair_interactions_);
#ifdef TRACE
  Logger::Trace("Updated interactions: pair interactions: %d -> %d", nix,
                pair_interactions_.size());
#endif
}

void InteractionManager::UpdateBoundaryInteractions() {
  if (no_boundaries_)
    return;
  boundary_interactions_.clear();
  for (auto ixor = interactors_.begin(); ixor != interactors_.end(); ++ixor) {
    Interaction ix(*ixor);
    if (mindist_.CheckBoundaryInteraction(ix)) {
      boundary_interactions_.push_back(ix);
    }
  }
}

int InteractionManager::CountSpecies() {
  int obj_count = 0;
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    obj_count += (*spec)->GetCount();
  }
  return obj_count;
}

const bool InteractionManager::CheckSpeciesInteractorUpdate() const {
  bool result = false;
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->CheckInteractorUpdate()) {
      result = true;
    }
  }
  return result;
}

void InteractionManager::CheckUpdateObjects() {
  /* First check to see if any objects were added to the system;
     If static_pnumber_ is flagged, we know that particle numbers
     never change, so we don't bother counting particles and move on */
  if (processing_)
    return;
  if (static_pnumber_)
    return;
  bool ix_update = CheckSpeciesInteractorUpdate();
  int obj_count = CountSpecies();
  if (obj_count != n_objs_ || ix_update) {
    // reset update count and update number of objects to track
    n_objs_ = obj_count;
    Logger::Debug("Updating interactions due to object update. %d steps since"
                  " last update",
                  i_update_);
    ForceUpdate();
    xlink_.UpdateObjsVolume();
    // XXX This has to be run before the first FlagDuplicateInteractions() call
    // if we check for overlaps at the beginning of the simulation...
    // ClearObjectInteractions();
  }
}

void InteractionManager::ResetCellList() { clist_.ResetNeighbors(); }

void InteractionManager::CheckUpdateInteractions() {
  /* Check if we dynamically changed the timestep last step, if so, always
     update */
  if (decrease_dynamic_timestep_) {
    decrease_dynamic_timestep_ = false;
    UpdateInteractions();
  }
  /* we update nearest neighbors if any particle
     has moved a distance further than dr_update_ */
  double dr_max = GetDrMax();
  if (dr_max > dr_update_) {
    Logger::Debug("Updating interactions due to dr of objects. %d steps since"
                  " last update, dr_max = %2.2f",
                  i_update_, dr_max);
    UpdateInteractions();
  }
}

/* Returns maximum distance traveled by any particle of any species in the
 * system */
double InteractionManager::GetDrMax() {
  double max_dr = 0;
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    double dr = (*spec)->GetDrMax();
    if (dr > max_dr) {
      max_dr = dr;
    }
  }
  double xlink_dr = xlink_.GetDrMax();
  if (xlink_dr > max_dr) {
    max_dr = xlink_dr;
  }
  return max_dr;
}

/* Resets the origin of each particle's tracked trajectory, for tracking
   total distance travelled by particle*/
void InteractionManager::ZeroDrTot() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->ZeroDrTot();
  }
  xlink_.ZeroDrTot();
}

void InteractionManager::ProcessPairInteraction(ix_iterator ix) {
  if (processing_) {
    mindist_.ObjectObject(*ix);
    //Check to see if particles are not close enough to interact
    if (ix->dr_mag2 > potentials_.GetRCut2())
      return;
    //Calculates forces from the potential defined during initialization
    potentials_.CalcPotential(*ix);
    return;
  }
  // Avoid certain types of interactions
  Object *obj1 = ix->obj1;
  Object *obj2 = ix->obj2;
#ifdef TRACE
  Logger::Trace("Processing interaction between %d and %d", obj1->GetOID(),
                obj2->GetOID());
#endif
  // Composite objects do self interact if they want to
  // Check to make sure we aren't self-interacting at the simple object level
  if (obj1->GetOID() == obj2->GetOID()) {
    Logger::Error("Object %d attempted self-interaction!", obj1->GetOID());
  }
  // If we are disallowing like-like interactions, then similar species do not
  // interact...
  // ...so do not interact if same species
  if (!params_->like_like_interactions && obj1->GetSID() == obj2->GetSID()) {
    ix->no_interaction = true;
    return;
  }

  /* If one of the objects are not prime interactors, don't interact */
  if (!obj1->IsInteractor() && !obj2->IsInteractor()) {
    ix->no_interaction = true;
    return;
  }
  // Check that objects are not both motors, which do not interact
  if (obj1->GetSID() == +species_id::crosslink &&
      obj2->GetSID() == +species_id::crosslink) {
    ix->no_interaction = true;
    return;
  }

  // Check that object 1 is part of a mesh, in which case...
  // ...check that object 1 is of the same mesh of object 2, in which case...
  if (obj1->GetMeshID() > 0 && obj1->GetMeshID() == obj2->GetMeshID()) {
    // ..check if object 1 or object 2 are crosslinks, in which case: do not
    // interact
    if (obj1->GetSID() == +species_id::crosslink ||
        obj2->GetSID() == +species_id::crosslink) {
      ix->no_interaction = true;
      return;
    }

    // ...check if object 2 is a neighbor of object 1, in which case: do not
    // interact
    if (obj1->HasNeighbor(obj2->GetOID())) {
      ix->no_interaction = true;
      return;
    }
  }

  /* If one object is a crosslink, add object to crosslink neighbor list */
  if (obj1->GetSID() == +species_id::crosslink) {
    xlink_.AddNeighborToAnchor(obj1, obj2);
    return;
  } else if (obj2->GetSID() == +species_id::crosslink) {
    xlink_.AddNeighborToAnchor(obj2, obj1);
    return;
  }
  mindist_.ObjectObject(*ix);

  // Check for particle overlaps
  if (ix->dr_mag2 < 0.25 * SQR(obj1->GetDiameter() + obj2->GetDiameter())) {
    overlap_ = true;
  }
  /* Check to see if particles are not close enough to interact */
  if (ix->dr_mag2 > potentials_.GetRCut2())
    return;
  /* Calculates forces from the potential defined during initialization */
  potentials_.CalcPotential(*ix);
}

void InteractionManager::ProcessBoundaryInteraction(ix_iterator ix) {
  mindist_.CheckBoundaryInteraction(*ix);
  // Check for particle overlapping with boundary edge
  if (ix->dr_mag2 < 0.25 * SQR(ix->obj1->GetDiameter())) {
    overlap_ = true;
  }
  if (ix->dr_mag2 > potentials_.GetRCut2())
    return;
  potentials_.CalcPotential(*ix);
}

void InteractionManager::CalculateBoundaryInteractions() {
  if (space_->type == +boundary_type::none) {
    return;
  }
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<ix_iterator, ix_iterator>> chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = boundary_interactions_.size() / max_threads;
  auto cur_iter = boundary_interactions_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, boundary_interactions_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto ix = chunks[i].first; ix != chunks[i].second; ++ix) {
        ProcessBoundaryInteraction(ix);
        // Do torque crossproducts
        cross_product(ix->contact1, ix->force, ix->t1, 3);
      }
    }
  }
#else
  for (auto ix = boundary_interactions_.begin();
       ix != boundary_interactions_.end(); ++ix) {
    ProcessBoundaryInteraction(ix);
    // Do torque crossproducts
    cross_product(ix->contact1, ix->force, ix->t1, 3);
  }
#endif
}

void InteractionManager::ClearObjectInteractions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<Object *>::iterator,
                        std::vector<Object *>::iterator>>
      chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = interactors_.size() / max_threads;
  auto cur_iter = interactors_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, interactors_.end()));
#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto it = chunks[i].first; it != chunks[i].second; ++it) {
        (*it)->ClearInteractions();
      }
    }
  }
#else
  for (auto it = interactors_.begin(); it != interactors_.end(); ++it) {
    (*it)->ClearInteractions();
  }
#endif
}

void InteractionManager::CalculatePairInteractions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<ix_iterator, ix_iterator>> chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = pair_interactions_.size() / max_threads;
  auto cur_iter = pair_interactions_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, pair_interactions_.end()));

#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto ix = chunks[i].first; ix != chunks[i].second; ++ix) {
        ix->pause_interaction = false;
        ProcessPairInteraction(ix);
        // Do torque crossproducts
        if (ix->no_interaction)
          continue;
        cross_product(ix->contact1, ix->force, ix->t1, 3);
        cross_product(ix->contact2, ix->force, ix->t2, 3);
        object_interaction oix1 = std::make_pair(&(*ix), true);
        object_interaction oix2 = std::make_pair(&(*ix), false);
        ix->obj1->GiveInteraction(oix1);
        ix->obj2->GiveInteraction(oix2);
      }
    }
  }
#else
  for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
       ++ix) {
    ix->pause_interaction = false;
    ProcessPairInteraction(ix);
    // Do torque crossproducts
    if (ix->no_interaction)
      continue;
    cross_product(ix->contact1, ix->force, ix->t1, 3);
    cross_product(ix->contact2, ix->force, ix->t2, 3);
    object_interaction oix1 = std::make_pair(&(*ix), true);
    object_interaction oix2 = std::make_pair(&(*ix), false);
    ix->obj1->GiveInteraction(oix1);
    ix->obj2->GiveInteraction(oix2);
  }
#endif
  /* After interaction update, remove pairs of interactors who can never
   * interact */
  if (i_update_ == 0) {
#ifdef TRACE
    int nix = pair_interactions_.size();
#endif
    pair_interactions_.erase(
        std::remove_if(pair_interactions_.begin(), pair_interactions_.end(),
                       [](Interaction x) { return x.no_interaction; }),
        pair_interactions_.end());
#ifdef TRACE
    Logger::Trace("Culling pair interactions. Pairs: %d -> %d", nix,
                  pair_interactions_.size());
#endif
  }
}

void InteractionManager::FlagDuplicateInteractions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<Object *>::iterator,
                        std::vector<Object *>::iterator>>
      chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = interactors_.size() / max_threads;
  auto cur_iter = interactors_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, interactors_.end()));
#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto it = chunks[i].first; it != chunks[i].second; ++it) {
        (*it)->FlagDuplicateInteractions();
      }
    }
  }
#else
  for (auto it = interactors_.begin(); it != interactors_.end(); ++it) {
    (*it)->FlagDuplicateInteractions();
  }
#endif
}

void InteractionManager::ApplyPairInteractions() {
  if (params_->remove_duplicate_interactions) {
    FlagDuplicateInteractions();
  }
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  std::vector<std::pair<std::vector<Object *>::iterator,
                        std::vector<Object *>::iterator>>
      chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = interactors_.size() / max_threads;
  auto cur_iter = interactors_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    auto last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, interactors_.end()));
#pragma omp parallel shared(chunks)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto it = chunks[i].first; it != chunks[i].second; ++it) {
        (*it)->ApplyInteractions();
      }
    }
  }
#else
  for (auto it = interactors_.begin(); it != interactors_.end(); ++it) {
    (*it)->ApplyInteractions();
  }
#endif
  if (params_->thermo_flag) {
    for (auto ix = pair_interactions_.begin(); ix != pair_interactions_.end();
         ++ix) {
      if (ix->pause_interaction)
        continue;
      for (int i = 0; i < n_dim_; ++i) {
        for (int j = 0; j < n_dim_; ++j) {
          stress_[n_dim_ * i + j] += ix->stress[n_dim_ * i + j];
        }
      }
    }
  }
}

void InteractionManager::ApplyBoundaryInteractions() {
  for (auto ix = boundary_interactions_.begin();
       ix != boundary_interactions_.end(); ++ix) {
    Object *obj1 = ix->obj1;
    obj1->AddForce(ix->force);
    obj1->AddTorque(ix->t1);
    obj1->AddPotential(ix->pote);
    for (int i = 0; i < n_dim_; ++i) {
      for (int j = 0; j < n_dim_; ++j) {
        stress_[n_dim_ * i + j] += ix->stress[n_dim_ * i + j];
      }
    }
  }
}

// Compute pressure tensor after n_thermo_ steps
void InteractionManager::CalculatePressure() {
  double inv_V = 1.0 / space_->volume;
  std::fill(space_->pressure_tensor, space_->pressure_tensor + 9, 0);
  // Calculate pressure tensor from stress tensor (only physical for periodic
  // subspace)
  for (int i = 0; i < n_dim_; ++i) {
    for (int j = 0; j < n_dim_; ++j) {
      // Add particle density along principle axes
      if (i == j) {
        space_->pressure_tensor[n_dim_ * i + j] += n_objs_ * inv_V;
      }
      // Add time-averaged virial component
      space_->pressure_tensor[n_dim_ * i + j] +=
          stress_[n_dim_ * i + j] * inv_V / n_thermo_;
    }
  }
  // Calculate isometric pressure
  space_->pressure = 0;
  for (int i = 0; i < n_dim_; ++i) {
    space_->pressure += space_->pressure_tensor[n_dim_ * i + i];
  }
  space_->pressure /= n_dim_;
  // Reset local stress tensor
  std::fill(stress_, stress_ + 9, 0);
}

bool InteractionManager::CheckOverlap(std::vector<Object *> &ixors) {
  overlap_ = false;
  /* Only consider objects (not crosslinks) for overlaps */
  for (auto ixor = ixors.begin(); ixor != ixors.end(); ++ixor) {
    pair_interactions_.clear();
    clist_.PairSingleObject(**ixor, pair_interactions_);
    CalculatePairInteractions();
    if (overlap_)
      break;
  }
  ClearObjectInteractions();
  if (potentials_.CheckMaxForce()) {
    return true;
  }
  return overlap_;
}

bool InteractionManager::CheckBoundaryConditions(std::vector<Object *> &ixors) {
  bool outside_boundary = false;
  for (auto ixor = ixors.begin(); ixor != ixors.end(); ++ixor) {
    if (mindist_.CheckOutsideBoundary(**ixor)) {
      outside_boundary = true;
    }
  }
  return outside_boundary;
}

/* Here, I want to calculate P(r, phi), which tells me the probability
   of finding an object at a position (r, phi) in its reference frame. */
void InteractionManager::InteractionAnalysis() {
  if (run_interaction_analysis_) {
    ForceUpdate();
    CalculateInteractions();
  }
}

void InteractionManager::Clear() {
  if (no_init_)
    return;
  clist_.Clear();
  xlink_.Clear();
}

/* Only used during species insertion */
void InteractionManager::Reset() {
  pair_interactions_.clear();
  ix_objects_.clear();
  interactors_.clear();
  clist_.ClearCellObjects();
}

/* Only used during species insertion */
void InteractionManager::AddInteractors(std::vector<Object *> &ixs) {
  clist_.AssignObjectsCells(ixs);
  ix_objects_.insert(ix_objects_.end(), ixs.begin(), ixs.end());
}

void InteractionManager::DrawInteractions(
    std::vector<graph_struct *> &graph_array) {
  xlink_.Draw(graph_array);
}

void InteractionManager::WriteOutputs() { xlink_.WriteOutputs(); }

void InteractionManager::InitOutputs(bool reading_inputs,
                                     run_options *run_opts) {
  xlink_.InitOutputs(reading_inputs, run_opts);
}

void InteractionManager::ReadInputs() { xlink_.ReadInputs(); }

void InteractionManager::InitCrosslinkSpecies(sid_label &slab,
                                              ParamsParser &parser,
                                              unsigned long seed) {
  xlink_.InitSpecies(slab, parser, seed);
}

void InteractionManager::LoadCrosslinksFromCheckpoints(
    std::string run_name, std::string checkpoint_run_name) {
  xlink_.LoadCrosslinksFromCheckpoints(run_name, checkpoint_run_name);
  PairBondCrosslinks();
  ForceUpdate();
}

void InteractionManager::InsertCrosslinks() { xlink_.InsertCrosslinks(); }

bool InteractionManager::CheckDynamicTimestep() {
  if (decrease_dynamic_timestep_) {
    for (auto spec_it = species_->begin(); spec_it != species_->end();
         ++spec_it) {
      (*spec_it)->ResetPreviousPositions();
    }
    return true;
  }
  return false;
}
