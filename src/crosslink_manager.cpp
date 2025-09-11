#include "cglass/crosslink_manager.hpp"

void CrosslinkManager::Init(system_parameters *params, SpaceBase *space, 
                            Tracker *tracker, std::vector<Object *> *objs, unsigned long seed) {
  objs_ = objs;
  update_ = false;
  obj_size_ = 0.0;
  params_ = params;
  space_ = space;
  tracker_ = tracker;
  rng_ = new RNG(seed);
}

void CrosslinkManager::InitSpecies(sid_label &slab, ParamsParser &parser,
                                   unsigned long seed) {
  if (xlink_species_.size() == 0) {
    xlink_species_.reserve(parser.GetNCrosslinkSpecies());
  }
  xlink_species_.push_back(new CrosslinkSpecies(seed));
  xlink_species_.back()->Init(slab.second, parser);
  Logger::Trace("Crosslinker species %s was added", (xlink_species_.back()->GetSpeciesName()).c_str());
  species_map_[(xlink_species_.back()->GetSpeciesName()).c_str()] = xlink_species_.back();
  if (xlink_species_.back()->GetNInsert() <= 0) {
    delete xlink_species_.back();
    xlink_species_.pop_back();
  } else {
    xlink_species_.back()->InitInteractionEnvironment(objs_, &obj_size_, tracker_,
                                                      &update_, &bound_curr_);
    rcutoff_ = xlink_species_.back()->GetRCutoff();
  }
}

/* Keep track of volume of objects in the system. Affects the
 * probability of a free crosslink binding to an object. */
void CrosslinkManager::UpdateObjsSize() {
  obj_size_ = 0;
  for (auto it = objs_->begin(); it != objs_->end(); ++it) {
    switch ((*it)->GetShape()) {
      case shape::rod:
        obj_size_ += (*it)->GetLength();
        break;
      case shape::sphere:
        if ((*it)->GetNAnchored() == 0) obj_size_ += (*it)->GetArea();
        break;
      default:
        break;
    }
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
  UpdateObjsSize();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->UpdatePositions();
    (*it)->UpdateBindRate();
  }
  if (!bound_curr_.empty()){
    Knockout();
  }
  //Use global check for cross so we don't need to check individual crosslinkers every time step
  if (global_check_for_cross == true) {
    CheckForCross();
    global_check_for_cross = false;
  }
}

//Check if any crosslinkers are crossing
void CrosslinkManager::CheckForCross() {
  Logger::Trace("Checking if crosslinkers are crossing");
  //Sorry this is a mess to look at, this is going over over crosslink in every species 
  //and comparing it to every other crosslink in every species to see if the crosslinks
  //are crossing
  for (auto spec_one = xlink_species_.begin(); spec_one != xlink_species_.end(); ++spec_one) {
  int member_num_ = (*spec_one) -> GetNMembers(); 
    for (int i = 0; i <= (member_num_-1); i++) {
      Crosslink* link_one = (*spec_one) -> GetCrosslink(i);
next_cl:
      // Check is we need to check for crossing for this crosslinker 
      if (link_one -> IsDoubly() && link_one -> ReturnCheckForCross() == true) {
        for (auto spec_two = xlink_species_.begin(); spec_two != xlink_species_.end(); ++spec_two) {
          int member_num_two_ = (*spec_two) -> GetNMembers();  
          for (int j = 0; j <= (member_num_two_-1); j++) {
            Crosslink* link_two = (*spec_two) -> GetCrosslink(j);
            //If link two is double bound and the crosslinkers aren't the same crosslinker
            if(link_two -> IsDoubly() && link_one -> GetOID() != link_two -> GetOID()) {
              //Check if anchors still double bound 
              bool link_one_db = link_one -> StillDBound();
              bool link_two_db = link_two -> StillDBound();
              if (!(link_one_db && link_two_db)) {
                continue;
              }                         
              //Get how far the crosslinker anchors are along the filament
              std::vector<double> link_one_s = link_one -> GetAnchorS();
              std::vector<double> link_two_s = link_two -> GetAnchorS();
              //Get the IDs of the PCs crosslink one is attatched to
              std::vector<int> rec_ids_one = link_one -> GetReceptorPCIDs();
              //Get the IDs of the PCs crosslink two is attatched to
              std::vector<int> rec_ids_two = link_two -> GetReceptorPCIDs();
              //If same index doesn't refer to sites on the same filament switch indexs so 
              //so the same index refers to anchors on the same microtubule
              if (rec_ids_one[0] != rec_ids_two[0]) {
                std::swap(link_two_s[0], link_two_s[1]);
                std::swap(rec_ids_two[0], rec_ids_two[1]);
              } 
              //If the sites are still on different filaments then crosslinkers aren't binding the same microtubules
              if (rec_ids_one[0] == rec_ids_two[0]) {
                same_microtubules = true;
              }
              else {
                same_microtubules = false;
              }
              //If crosslinkers are crossing and connected same microtubules, then unbind most rect crosslinker
              if (same_microtubules && ((link_one_s[0] > link_two_s[0] && link_one_s[1] < link_two_s[1]) 
                     || (link_one_s[0] < link_two_s[0] && link_one_s[1] > link_two_s[1]))) { 
                link_one -> UnbindCrossing();
                goto next_cl;
              } else {
                link_one -> SetCheckForCross();
              } 
            }
          }
        }    
      }
    }
  }
}

// Loop over all sites crosslinkers tried to bind to in the last dt. If multiple xlinks want to bind
// to a site, roll and select based on their relative probabilities, and unbind
void CrosslinkManager::Knockout() {

  //If trace is on print everything in bound_curr
  #ifdef TRACE
	  Logger::Warning("\n");
    Logger::Warning("Printing bound_curr, size is %i", bound_curr_.size());
    for (auto& kv : bound_curr_) {
      Sphere* receptor = kv.first;
      std::pair<std::vector<double>, std::vector<std::pair<Anchor*, std::string> > >&  receptor_info = kv.second;
      Logger::Trace("Receptor id is, %d", receptor->GetOID());
      int events = receptor_info.first.size();
      while (events>0) {
        events-=1;
        if (receptor_info.second[events].second == "forward step" || receptor_info.second[events].second == "back step" || receptor_info.second[events].second == "single to double") {
         Logger::Trace("Prop is %f, Anchor is %i, Bind type is %s", receptor_info.first[events], receptor_info.second[events].first->GetOID(), receptor_info.second[events].second.c_str());
        }
      }
    }
    Logger::Trace("\n");
  #endif

  for (auto& kv : bound_curr_) {
    Sphere* receptor = kv.first;
    std::pair<std::vector<double>, std::vector<std::pair<Anchor*, std::string> > >&  receptor_info = kv.second;
    // If one anchor wants to bind to a receptor bind it
    if (receptor_info.first.size() == 1) {
      KnockoutBind(receptor, 0);
    }

    //If two or three anchors want to bind to a receptor choose which one binds 
    else if (receptor_info.first.size() >= 3 ) {
      double sum_of_probs = std::accumulate(receptor_info.first.begin(), receptor_info.first.end(), 0.0);
      Logger::Warning("%i anchors want to bind to site %i with total probability %f. If this warning is common, or if probability is greater than 1 reduce time step.",
                     receptor_info.first.size(), receptor->GetOID(), sum_of_probs);
      double roll = sum_of_probs*rng_->RandomUniform();
      int count = 0;
      std::vector<double> prob_list = receptor_info.first;
      for (auto prob = prob_list.begin(); prob != prob_list.end(); ++prob) {
        if (*prob>roll) {
          KnockoutBind(receptor, count);
          Logger::Trace("chose %s, from roll, %f", receptor_info.second[count].second.c_str(), roll); 
          break;
        }
        else {
          roll -= *prob;
          count+= 1;
        }
      }
    }

    // If more than 3 anchors want to bind to a receptor throw error
    else if (receptor_info.first.size()>3) {
      Logger::Error("More than 3 anchors trying to bind to receptor, time step too large");
    }
    for (auto it = receptor_info.second.begin(); it != receptor_info.second.end(); ++it){
      if (it->second == "forward step" || it->second == "back step") {
        it->first->ResetChangedThisStep();
      }
    }
  }
  AddKnockoutCrosslinks();
  bound_curr_.clear();
  Logger::Trace("bound_curr cleared");
}

void CrosslinkManager::KnockoutBind(Sphere* receptor, int winner) {

  //Forward Step
  if (bound_curr_[receptor].second[winner].second == "forward step") {
    bound_curr_[receptor].second[winner].first->StepForward();
  }
  //Back Step
  else if (bound_curr_[receptor].second[winner].second == "back step") {
    bound_curr_[receptor].second[winner].first->StepBack();
  }

  //Single to double bind
  else if (bound_curr_[receptor].second[winner].second == "single to double") {
    Logger::Trace("single to double %i", bound_curr_[receptor].second[winner].first->GetOID());
    Object* ob_pointer = bound_curr_[receptor].second[winner].first->GetCrosslinkPointer();
    Crosslink* cl_pointer = dynamic_cast<Crosslink*>(ob_pointer);
    if (!cl_pointer) {
      Logger::Warning("Dynamic cast failed, cl_pointer doesn't exist");
    }    
    if (cl_pointer){
    cl_pointer->SetDoubly();
    bound_curr_[receptor].second[winner].first->AttachObjCenter(receptor);
    } 
 }

  //Free to single bind
  else if (bound_curr_[receptor].second[winner].second == "free to single") {
  Logger::Info("Free to single bind");
    Object* ob_pointer = bound_curr_[receptor].second[winner].first->GetCrosslinkPointer();
    Crosslink* cl_pointer = dynamic_cast<Crosslink*>(ob_pointer);
    if (!cl_pointer) {
      Logger::Error("Dynamic cast failed, cl_ponter doesn't exist");
    }
    cl_pointer->SetSingly(0);
    bound_curr_[receptor].second[winner].first->AttachObjCenter(receptor);
  }

  //Bind from soultion, bind type is set to species name if binding from solution 
  else {
   CrosslinkSpecies* cl_species_ = species_map_[bound_curr_[receptor].second[winner].second.c_str()];
   std::pair<CrosslinkSpecies*, Sphere*> cl_rec_pair_;
   cl_rec_pair_.first = cl_species_;
   cl_rec_pair_.second = receptor;
   cl_to_add_this_step_.push_back(cl_rec_pair_);
  }
}

void CrosslinkManager::AddKnockoutCrosslinks() {
  for (auto it = cl_to_add_this_step_.begin(); it != cl_to_add_this_step_.end(); ++it) {
    Logger::Trace("Added crosslink from solution during knockout");
    (it->first)->KnockoutBind(it->second);
  }
  cl_to_add_this_step_.clear();
}

void CrosslinkManager::InsertCrosslinks() {
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->InsertCrosslinks();
    (*it)->SetGlobalCheckForCrossPointer(&global_check_for_cross);
  } 
}

void CrosslinkManager::InsertAttachedCrosslinks(std::vector<std::vector<Object *>> receptor_list) {
  // Need to do this for GetRandomObject to work with spheres
  UpdateObjsSize();
  for (auto it = xlink_species_.begin(); it != xlink_species_.end(); ++it) {
    (*it)->UpdateBindRate();
    (*it)->InsertAttachedCrosslinksSpecies(receptor_list);
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
void CrosslinkManager::Convert() { output_mgr_.Convert(); }
