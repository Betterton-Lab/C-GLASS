#include "cglass/crosslink_manager.hpp"

void CrosslinkManager::Init(system_parameters *params, space_struct *space,
                            std::vector<Object *> *objs) {
  objs_ = objs;
  update_ = false;
  obj_volume_ = 0;
  params_ = params;
  space_ = space;
}

void CrosslinkManager::InitSpecies(sid_label &slab, ParamsParser &parser,
                                   unsigned long seed) {
  if (xlink_species_.size() == 0) {
    xlink_species_.reserve(parser.GetNCrosslinkSpecies());
  }
  xlink_species_.push_back(new CrosslinkSpecies(seed));
  xlink_species_.back()->Init(slab.second, parser);
  if (xlink_species_.back()->GetNInsert() <= 0) {
    delete xlink_species_.back();
    xlink_species_.pop_back();
  } else {
    xlink_species_.back()->InitInteractionEnvironment(objs_, &obj_volume_,
                                                      &update_);
    rcutoff_ = xlink_species_.back()->GetRCutoff();
  }
}

/* Keep track of volume of objects in the system. Affects the
 * probability of a free crosslink binding to an object. */
void CrosslinkManager::UpdateObjsVolume() {
  obj_volume_ = 0;
  for (auto it = objs_->begin(); it != objs_->end(); ++it) {
    // obj_volume_ += (*it)->GetVolume(); //XXX: Binding based off length not
    // volume
    obj_volume_ += (*it)->GetLength();
  }
}

/* Whether to reinsert anchors into the interactors list */
bool CrosslinkManager::CheckUpdate() {
  if (update_) {
    update_ = false;
    return true;
  }
  return false;
}

/* Return singly-bound anchors, for finding neighbors to bind to */
void CrosslinkManager::GetInteractors(std::vector<Object *> &ixors) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->GetInteractors(ixors);
  }
}

/* Returns all anchors, not just singly-bound anchors. Used for reassigning
   bound anchors to bonds upon a checkpoint reload */
void CrosslinkManager::GetAnchorInteractors(std::vector<Object *> &ixors) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->GetAnchorInteractors(ixors);
  }
}

void CrosslinkManager::UpdateCrosslinks() {
  update_ = false;
  UpdateObjsVolume();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->UpdatePositions();
  }
}

void CrosslinkManager::InsertCrosslinks() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->InsertCrosslinks();
  }
}

void CrosslinkManager::Clear() {
  output_mgr_.Close();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->CleanUp();
    delete (*it);
  }
  xlink_species_.clear();
}

void CrosslinkManager::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->Draw(graph_array);
  }
}

void CrosslinkManager::AddNeighborToAnchor(Object *anchor, Object *neighbor) {
  if (anchor->GetSID() != +species_id::crosslink) {
    Logger::Error(
        "AddNeighborToAnchor expected crosslink object, got generic object.");
  }
  Anchor *a = dynamic_cast<Anchor *>(anchor);
  a->AddNeighbor(neighbor);
}

void CrosslinkManager::InitOutputs(bool reading_inputs, run_options *run_opts) {
  output_mgr_.Init(params_, &xlink_species_, space_, reading_inputs, run_opts);
}

void CrosslinkManager::WriteOutputs() { output_mgr_.WriteOutputs(); }

void CrosslinkManager::ZeroDrTot() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->ZeroDrTot();
  }
}

const int CrosslinkManager::GetDoublyBoundCrosslinkNumber() const {
  int num = 0;
  for (const auto &xl_spec : xlink_species_) {
    num += xl_spec->GetDoublyBoundCrosslinkNumber();
  }
  return num;
}

const double CrosslinkManager::GetDrMax() {
  double dr_max = 0;
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    double dr = (*it)->GetDrMax();
    if (dr > dr_max) {
      dr_max = dr;
    }
  }
  return dr_max;
}
void CrosslinkManager::LoadCrosslinksFromCheckpoints(
    std::string run_name, std::string checkpoint_run_name) {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->LoadFromCheckpoints(run_name, checkpoint_run_name);
  }
}

void CrosslinkManager::ReadInputs() { output_mgr_.ReadInputs(); }
