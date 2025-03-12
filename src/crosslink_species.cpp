#include "cglass/crosslink_species.hpp"

CrosslinkSpecies::CrosslinkSpecies(unsigned long seed) : Species(seed), bind_param_map_(2),
                                                         default_bind_params_(2) {
  SetSID(species_id::crosslink);
}
void CrosslinkSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  k_on_ = sparams_.anchors[0].k_on_s + sparams_.anchors[1].k_on_s;
  k_off_ = sparams_.anchors[0].k_off_s + sparams_.anchors[1].k_off_s;
  bind_site_density_ = sparams_.bind_site_density;
  begin_with_bound_crosslinks_ = sparams_.begin_with_bound_crosslinks;
  xlink_concentration_ = sparams_.concentration;
  infinite_reservoir_flag_ = sparams_.infinite_reservoir_flag;
  step_count_=0;
  if (sparams_.use_number==false) {
    sparams_.num = (int)round(sparams_.concentration * space_->volume);
  }
  printf("number is %d \n", sparams_.num);
  std::vector<std::string> bind_file = {sparams_.anchors[0].bind_file, sparams_.anchors[1].bind_file};
  
  // Create a default set of specific binding parameters
  InitializeBindParams();

  if (bool(bind_file[0].compare("none"))!=bool(bind_file[1].compare("none"))) {
    Logger::Error("Cannot have bind file=none for one anchor and a bind file for the other.");
  }

  // If bind file is not "none", load bind file parameters and calulate the number of 
  // expected binding events per file.
  if (bind_file[0].compare("none")) {
    LoadBindingSpecies();
    use_bind_file_ = true;
  }
}

// Set default params to whatever this crosslink species has listed (if not listed,
// these will be the default params from the config file).
void CrosslinkSpecies::InitializeBindParams() {
  for (int i = 0; i < 2; ++i) {
    default_bind_params_[i].use_partner = sparams_.anchors[i].use_partner;
    default_bind_params_[i].k_on_s = sparams_.anchors[i].k_on_s;
    default_bind_params_[i].partner_on_s = sparams_.anchors[i].partner_on_s;
    default_bind_params_[i].k_on_d = sparams_.anchors[i].k_on_d;
    default_bind_params_[i].partner_on_d = sparams_.anchors[i].partner_on_d;
    default_bind_params_[i].k_off_s = sparams_.anchors[i].k_off_s;
    default_bind_params_[i].k_off_d = sparams_.anchors[i].k_off_d;
    default_bind_params_[i].dens_type = density_type::linear;
    default_bind_params_[i].bind_site_density = sparams_.bind_site_density;
    default_bind_params_[i].single_occupancy = false;
  }
}

void CrosslinkSpecies::LoadBindingSpecies() {
  YAML::Node bnode;
  for (int anchor_index = 0; anchor_index < 2; ++anchor_index) {
    try {
      bnode = YAML::LoadFile(sparams_.anchors[anchor_index].bind_file);
    } catch (...) {
      Logger::Error("Failed to load binding species file in crosslink_species.cpp");
    }
    for (auto it = bnode.begin(); it != bnode.end(); ++it) {
      std::string param_name = it->first.as<std::string>();
      // Error if there is already a map entry with given string name
      if (bind_param_map_[anchor_index].find(param_name) != bind_param_map_[anchor_index].end()) {
        Logger::Error("Duplicate species names found in Crosslink bind file.");
      }
      // Initialize to default
      bind_param_map_[anchor_index][param_name] = default_bind_params_[anchor_index];
      if (it->second["use_partner"]) {
        bind_param_map_[anchor_index][param_name].use_partner = it->second["use_partner"].as<bool>();
      }if (it->second["k_on_s"]) {
        bind_param_map_[anchor_index][param_name].k_on_s = it->second["k_on_s"].as<double>();
      }
      if (it->second["k_on_d"]) {
        bind_param_map_[anchor_index][param_name].k_on_d = it->second["k_on_d"].as<double>();
      }
      if (it->second["k_off_s"]) {
        bind_param_map_[anchor_index][param_name].k_off_s = it->second["k_off_s"].as<double>();
      }
      if (it->second["k_off_d"]) {
        bind_param_map_[anchor_index][param_name].k_off_d = it->second["k_off_d"].as<double>();
      }
      if (it->second["partner_on_s"]) {
        bind_param_map_[anchor_index][param_name].partner_on_s = it->second["partner_on_s"].as<double>();
      }
      if (it->second["partner_on_d"]) {
        bind_param_map_[anchor_index][param_name].partner_on_d = it->second["partner_on_d"].as<double>();
      }
      if (it->second["density_type"]) {
        bind_param_map_[anchor_index][param_name].dens_type = 
                        density_type::_from_string(it->second["density_type"].Scalar().c_str());
      } else {
        Logger::Warning("No density type found for species name %s, defaulting to linear.", param_name.c_str());
      }
      if (it->second["bind_site_density"]) {
        bind_param_map_[anchor_index][param_name].bind_site_density = it->second["bind_site_density"].as<double>();
      }
      if (it->second["single_occupancy"]) {
        bind_param_map_[anchor_index][param_name].single_occupancy = it->second["single_occupancy"].as<bool>();
      }
    }
  }
}

void CrosslinkSpecies::AddMember() {
  Species::AddMember();
  members_.back().InitInteractionEnvironment(&lut_, tracker_, bound_curr_);
  members_.back().SetObjSize(obj_size_);
  members_.back().SetBindRate(&bind_rate_);
  members_.back().SetBindParamMap(&bind_param_map_);
  *update_ = true;
}

void CrosslinkSpecies::InitInteractionEnvironment(std::vector<Object *> *objs,
                                                  double *obj_size,
                                                  Tracker *tracker, bool *update,
                       std::map<Sphere *, std::pair<std::vector<double>, std::vector<std::pair<Anchor*, std::string> > > > *bound_curr) {
  objs_ = objs;
  obj_size_ = obj_size;
  bound_curr_ = bound_curr;
  update_ = update;
  tracker_ = tracker;
  LUTFiller *lut_filler_ptr = MakeLUTFiller();
  lut_ = LookupTable(lut_filler_ptr);
  if (sparams_.use_binding_volume) {
    lut_.setBindVol(lut_filler_ptr->getBindingVolume());
  }
  /* TODO: Add time testing right here <24-06-20, ARL> */
  //TestKMCStepSize();
  delete lut_filler_ptr;
}

/*! \brief Test simulation step size to make sure KMC will work properly.
 *
 * \return void, will end program if dt is too high.
 */
void CrosslinkSpecies::TestKMCStepSize() {
  /* TODO: Make this acceptable for C-GLASS <24-06-20, ARL> */
  KMC<int, char> kmc_diag(params_->delta, &lut_);

  const double prob_thresh = 1e-4;

  double probs[4];

  // Constant rate factors for (un)binding. Change if bind model changes.
  double u_s_fact =
      xlink_concentration_ * k_on_ * bind_site_density_; //per length basis
  double s_u_fact = k_off_;
  double s_d_fact = k_on_ * bind_site_density_;
  if (sparams_.use_binding_volume)
    s_d_fact /= lut_.getBindVolume();
  double d_s_fact = k_off_;

  if (sparams_.energy_dep_factor == 0 && sparams_.force_dep_length == 0) {
    kmc_diag.Diagnostic(u_s_fact, s_u_fact, s_d_fact, d_s_fact, probs);
  } else {
    kmc_diag.DiagnosticUnBindDep(u_s_fact, s_u_fact, s_d_fact, probs);
  }

  // If any probabilities of a double binding event occuring is too high,
  // throw a warning.
  bool throw_error = false;
  if (probs[0] > prob_thresh || std::isnan(probs[0])) {
    Logger::Warning(
        "Probability of double event (U->S->U = %f) is too high. Try "
        "decreasing delta, concentration, or singly (un)binding parameters.",
        probs[0]);
    throw_error = true;
  }
  if (probs[1] > prob_thresh || std::isnan(probs[1])) {
    Logger::Warning(
        "Probability of double event (U->S->D = %f) is too high. Try "
        "decreasing delta, concentration, or binding parameters.",
        probs[1]);
    throw_error = true;
  }
  if (probs[2] > prob_thresh || std::isnan(probs[2])) {
    Logger::Warning("Probability of double event (S->D->S = %f) is too high. "
                    "Try decreasing delta or doubly (un)binding parameters.",
                    probs[2]);
    throw_error = true;
  }
  if (probs[3] > prob_thresh || std::isnan(probs[3])) {
    Logger::Warning(" !!!Probability of double event (D->S->U = %f) is too "
                    "high. Try decreasing delta or unbinding parameters.",
                    probs[3]);
    throw_error = true;
  }

  if (throw_error) {
    Logger::Error("The likelyhood of a double KMC event is too high. Fix "
                  "before continuing.");
  }
}

LUTFiller *CrosslinkSpecies::MakeLUTFiller() {
  int grid_num = sparams_.lut_grid_num;
  if (sparams_.k_spring_compress >= 0) {
    Logger::Warning("!!!Asymmetric springs are being used. This has not been "
                    "fully tested. Use at your own risk!");
    LUTFillerAsym *lut_filler_ptr = new LUTFillerAsym(grid_num, grid_num);
    lut_filler_ptr->Init(sparams_.k_spring_compress, sparams_.k_spring,
                         sparams_.energy_dep_factor, sparams_.force_dep_length,
                         sparams_.rest_length, 1);
    return lut_filler_ptr;

  } else if (sparams_.force_dep_length == 0) {
    LUTFillerEdep *lut_filler_ptr = new LUTFillerEdep(grid_num, grid_num);
    lut_filler_ptr->Init(sparams_.k_spring * .5 *
                             (1. - sparams_.energy_dep_factor),
                         sparams_.rest_length, 1);
    return lut_filler_ptr;
  } else {
    LUTFillerFdep *lut_filler_ptr = new LUTFillerFdep(grid_num, grid_num);
    lut_filler_ptr->Init(sparams_.k_spring, sparams_.energy_dep_factor,
                         sparams_.force_dep_length, sparams_.rest_length, 1);
    return lut_filler_ptr;
  }
}

void CrosslinkSpecies::InsertCrosslinks() {
  // Random insertion (default) implies crosslinks in solution
  if (sparams_.insertion_type.compare("random") == 0) {
    sparams_.static_flag = false;
    return;
  } else if (sparams_.insertion_type.compare("free") == 0) {
    Logger::Info("Inserting %i unbound crosslinks", sparams_.num);
if (space_->type != +boundary_type::sphere  && space_->type != +boundary_type::protrusion) {
    Logger::Error("Tracking unbound crosslinkers not set up for current boundry type");
    } else if (space_->type == +boundary_type::sphere || space_->type == +boundary_type::protrusion) {
      for (int i = 0; i < sparams_.num; ++i) {
        AddMember();
        members_.back().InsertRandom(); 
        if (params_->n_dim != 3) {
          Logger::Error("Tracking unbound crosslinkers only set up for 3 dimensions");
        } else {
          double projected_pos[3] = {0};
          double same_u[3] = {0};
          const double *const pos = members_.back().GetPosition();
          const double *const u = members_.back().GetOrientation();
          for (int j = 0; j < 3; ++j) {
            projected_pos[j] = pos[j];
            same_u[j] = u[j];
          }
        members_.back().InsertFree(projected_pos, same_u);
        }
      }
    } 
  } else if (sparams_.insertion_type.compare("centered") == 0) {
    sparams_.num = 1;
    sparams_.static_flag = true;
    sparams_.infinite_reservoir_flag = false;
    AddMember();
    double pos[3] = {0, 0, 0};
    double u[3] = {0, 0, 0};
    u[params_->n_dim - 1] = 1.0;
    members_.back().InsertAt(pos, u);
  } else if (sparams_.insertion_type.compare("random_grid") == 0) {
    sparams_.static_flag = true;
    sparams_.infinite_reservoir_flag = false;
    if (params_->n_dim == 3) {
      if (space_->type == +boundary_type::none ||
          space_->type == +boundary_type::box) {
        sparams_.num = (int)round(4 * space_->radius * space_->radius *
                                  xlink_concentration_);
      } else if (space_->type == +boundary_type::sphere || space_->type == +boundary_type::protrusion) {
        sparams_.num = (int)round(M_PI * space_->radius * space_->radius *
                                  xlink_concentration_);
      } else if (space_->type == +boundary_type::budding) {
        double R = space_->radius;
        double r = space_->bud_radius;
        double d = space_->bud_height;
        sparams_.num = (int)round(
            M_PI * SQR(R) + M_PI * SQR(r) -
            SQR(r) * acos((SQR(d) + SQR(r) - SQR(R)) / (2 * d * r)) -
            SQR(R) * acos((SQR(d) + SQR(R) - SQR(r)) / (2 * d * R)) +
            0.5 * sqrt((R + r - d) * (r + d - R) * (R + d - r) * (R + r + d)));
      } else {
        Logger::Error("Boundary type not recognized in CrosslinkSpecies");
      }
    }
    for (int i = 0; i < sparams_.num; ++i) {
      AddMember();
      members_.back().InsertRandom();
      /* If in 3D, zero out the third dimension so the anchor is on a plane */
      if (params_->n_dim == 3) {
        double projected_pos[3] = {0};
        double same_u[3] = {0};
        const double *const pos = members_.back().GetPosition();
        const double *const u = members_.back().GetOrientation();
        for (int j = 0; j < 3; ++j) {
          projected_pos[j] = pos[j];
          same_u[j] = u[j];
        }
        projected_pos[2] = 0;
        members_.back().InsertAt(projected_pos, same_u);
      }
    }
  } else if (sparams_.insertion_type.compare("random_boundary") == 0) {
    sparams_.infinite_reservoir_flag = false;
    sparams_.num = (int)round(space_->BoundaryArea() * xlink_concentration_);        
    double pos[3] = {0};
    double u[3] = {1, 0, 0};
    for (int i = 0; i < sparams_.num; ++i) {
      rng_.RandomBoundaryCoordinate(space_, pos);
      AddMember();
      members_.back().InsertAt(pos, u);
    } 
  } else {
    Logger::Error("Insertion type %s not implemented yet for crosslinks",
                  sparams_.insertion_type.c_str());
  }
}

//Adds in Crosslinkers for begin_with_bound_crosslinks flag
void CrosslinkSpecies::InsertAttachedCrosslinksSpecies(std::vector<std::vector<Object *>> receptor_list) {
  if (begin_with_bound_crosslinks_<=0) {
    return;
  }
  Crosslink xlink(rng_.GetSeed());
  xlink.Init(&sparams_);
  xlink.InitInteractionEnvironment(&lut_, tracker_, bound_curr_);
  xlink.SetObjSize(obj_size_);
  xlink.SetBindRate(&bind_rate_);
  xlink.SetBindParamMap(&bind_param_map_);
  xlink.SetSID(GetSID());
  members_.resize(begin_with_bound_crosslinks_, xlink);
  UpdateBoundCrosslinks();
  //If crosslinkers are starting singly bound
  if (sparams_.begin_double_bound == false) {
    for (int i=0; i < begin_with_bound_crosslinks_; ++i) {
      BindCrosslink();
      members_.back().SetGlobalCheckForCross(global_check_for_cross_);
    }
  }
  //If crosslinkers are starting doubly bound 
  else {
    for (int i=0; i < begin_with_bound_crosslinks_; ++i) {
      BindDoubly(receptor_list[0][i], receptor_list[1][i]);
      members_.back().SetGlobalCheckForCross(global_check_for_cross_);
    }
  }
     
}

// Calculate and bind crosslinkers from solution implicitly
void CrosslinkSpecies::CalculateBindingFree() {
  /* Static crosslinks are never free */
  if (sparams_.static_flag) {
    return;
  }
  /* Check crosslink binding */
  double free_concentration;
  if (infinite_reservoir_flag_) { // Have a constant concentration of
                                  // crosslinkers binding from solution
    free_concentration = xlink_concentration_;
  } else { // Have a constant number of crosslinkers in a space
    //printf("Number is %i \n",sparams_.num - n_members_);
    free_concentration = (sparams_.num - n_members_) / space_->volume;
  }
  if (use_bind_file_) {
    int num_to_bind = rng_.RandomPoisson(free_concentration * params_->delta * bind_rate_);
    for (int i = 0; i < num_to_bind; ++i) {
      BindCrosslink();
    }
  } else {
    double expected_bind_n = bind_site_density_ * free_concentration *
                       (*obj_size_) * k_on_ * params_->delta;
    // Use a Poisson distribution to calculate the number of particles
    int bind_num = rng_.RandomPoisson(expected_bind_n);
    // Track US probabilities
    tracker_->TrackUS(expected_bind_n);
    tracker_->BindUS(bind_num);
    for (int i = 0; i < bind_num; ++i) {
      /* Create a new crosslink and bind an anchor to a random rod
       * in the system */
      BindCrosslink();
    }
  }
}

/* Returns a random object with selection probability proportional to object
   bind rate. */
std::pair <Object*, int> CrosslinkSpecies::GetRandomObject() {
  if (use_bind_file_) {
    double bind_rate_sum = 0;
    std::string name = "";
    // Roll to choose object
    double roll = bind_rate_*rng_.RandomUniform();
    // Count up bind rates to pick correct object
    for (int anchor_index = 0; anchor_index < 2; ++anchor_index) {
      for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
        name = (*obj)->GetName();
        auto bind_param_it = bind_param_map_[anchor_index].find(name);
        // Do not count objects that are not in bind_param_map
        if (bind_param_it == bind_param_map_[anchor_index].end()) continue;
        // obj_amount is area if surface density used, length if linear density used
        double obj_amount = (bind_param_it->second.dens_type == +density_type::linear) 
                            ? (*obj)->GetLength() : (*obj)->GetArea();
        // Do not contribute area/length if object is already occupied and single occupancy is
        // selected
        if (bind_param_it->second.single_occupancy && ((*obj)->GetNAnchored() > 0)) {
          obj_amount = 0;
        }
        bind_rate_sum += bind_param_it->second.k_on_s * 
                         bind_param_it->second.bind_site_density * obj_amount;
        if (bind_rate_sum > roll) {
          Logger::Trace("Binding free crosslink to random object: xl %d -> obj %d",
                        members_.back().GetOID(), (*obj)->GetOID());
          random_obj_probability_ = bind_rate_sum;
          return std::make_pair(*obj, anchor_index);
        }
      }
    }
    Logger::Error("Crosslinks tried to bind to more sites than available- lower timestep.");
    return std::make_pair(nullptr, -1);
  } else {
    // 2x factor comes from 2 anchors
    double roll = 2*rng_.RandomUniform()*(*obj_size_);
    double vol = 0;

    // Search through interactors to find an object of shape sh
    for (int anchor_index = 0; anchor_index < 2; ++anchor_index) {
      for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
        switch((*obj)->GetShape()) {
          case shape::rod:
            vol += (*obj)->GetLength();
            break;
          case shape::sphere:{
            std::string name = (*obj)->GetName();
            auto bind_param_it = bind_param_map_[anchor_index].find(name);
            if (bind_param_it->second.single_occupancy == 0 || ((*obj)->GetNAnchored()==0)) {
              vol += (*obj)->GetArea();
            }
          }
          break;
          default:
          break;
        }
        if (vol > roll) {
          Logger::Trace("Binding free crosslink to random object: xl %d -> obj %d",
                        members_.back().GetOID(), (*obj)->GetOID());
          return std::make_pair(*obj, anchor_index);
        }
      }
    }
    Logger::Error("Crosslinks tried to bind to more sites than available- lower timestep.");
    return std::make_pair(nullptr, -1);
  }
}

//Same function as GetRandomObj().
//Get GetRandomObj() expects that a new crosslink has been added before getting an object to bind it to, this doesn't.
std::pair <Object*, int> CrosslinkSpecies::GetRandomObjectKnockout() {
  if (use_bind_file_) {
    double bind_rate_sum = 0;
    // Roll to choose object
    double roll = bind_rate_*rng_.RandomUniform();
    // Count up bind rates to pick correct object
    for (int anchor_index = 0; anchor_index < 2; ++anchor_index) {
      for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
        std::string name = (*obj)->GetName();
        auto bind_param_it = bind_param_map_[anchor_index].find(name);
        // Do not count objects that are not in bind_param_map
        if (bind_param_it == bind_param_map_[anchor_index].end()) continue;
        // obj_amount is area if surface density used, length if linear density used
        double obj_amount = (bind_param_it->second.dens_type == +density_type::linear)
                            ? (*obj)->GetLength() : (*obj)->GetArea();
        // Do not contribute area/length if object is already occupied and single occupancy is
        // selected
        if (bind_param_it->second.single_occupancy && ((*obj)->GetNAnchored() > 0)) {
          obj_amount = 0;
        }
        bind_rate_sum += bind_param_it->second.k_on_s *
                         bind_param_it->second.bind_site_density * obj_amount;
        if (bind_rate_sum > roll) {
          random_obj_probability_ = bind_rate_sum;
          return std::make_pair(*obj, anchor_index);
        }
      }
    }
    Logger::Error("Crosslinks tried to bind to more sites than available- lower timestep.");
    return std::make_pair(nullptr, -1);
  } else {
    // 2x factor comes from 2 anchors
    double roll = 2*rng_.RandomUniform()*(*obj_size_);
    double vol = 0;

    // Search through interactors to find an object of shape sh
    for (int anchor_index = 0; anchor_index < 2; ++anchor_index) {
      for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
        switch((*obj)->GetShape()) {
          case shape::rod:
            vol += (*obj)->GetLength();
            break;
          case shape::sphere:{
            std::string name = (*obj)->GetName();
            auto bind_param_it = bind_param_map_[anchor_index].find(name); 
            if (bind_param_it->second.single_occupancy == 0 || (*obj)->GetNAnchored()==0) {
              vol += (*obj)->GetArea();
            }
          }
          break;
          default:
          break;
        }
        if (vol > roll) {
         return std::make_pair(*obj, anchor_index);
        }
      }
    }
    Logger::Error("Crosslinks tried to bind to more sites than available- lower timestep.");
    return std::make_pair(nullptr, -1);
  }
}

/* A crosslink binds to an object from solution */
void CrosslinkSpecies::BindCrosslink() {  
  std::pair <Object*, int>  RandomObj = GetRandomObjectKnockout();
  //If binding to sphere, at to knockout
  if (RandomObj.first->GetShape() == +shape::sphere) {
    Sphere* RandomSphere = dynamic_cast<Sphere*>(RandomObj.first);
    (*bound_curr_)[RandomSphere].first.push_back(random_obj_probability_);
    std::pair<Anchor*, std::string> anchor_and_bind_type;
    anchor_and_bind_type.first = nullptr;
    anchor_and_bind_type.second = this->GetSpeciesName();
    (*bound_curr_)[RandomSphere].second.push_back(anchor_and_bind_type);

  //If binding to rod, bind now
  } else {
    AddMember();
    members_.back().AttachObjRandom(RandomObj);
    members_.back().SetGlobalCheckForCross(global_check_for_cross_);
  }
}

//Bind during knockout
void CrosslinkSpecies::KnockoutBind(Sphere* receptor) {
    AddMember();
    std::pair <Sphere*, int>  Receptor_pair;
    //receptor to bind to
    Receptor_pair.first = receptor;
    //anchor 0
    Receptor_pair.second = 0;
    members_.back().AttachSphere(Receptor_pair);
    members_.back().SetGlobalCheckForCross(global_check_for_cross_);
}

//Insert crosslink attatched to two receptors
void CrosslinkSpecies::BindDoubly(Object* receptor_one, Object* receptor_two) {
  AddMember();
  members_.back().DoublyCenter(receptor_one, receptor_two);
}

/* Return singly-bound anchors, for finding neighbors to bind to */
void CrosslinkSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->GetInteractors(ixors);
  }
}

/* Returns all anchors, not just singly-bound anchors. Used for reassigning
   bound anchors to objects upon a checkpoint reload */
void CrosslinkSpecies::GetAnchorInteractors(std::vector<Object *> &ixors) {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->GetAnchors(ixors);
  }
}

Crosslink* CrosslinkSpecies::GetCrosslink(int i) {
  return &members_[i];
}

void CrosslinkSpecies::UpdatePositions() {
  /* Only do this every other step (assuming flexible filaments with midstep)
   */
  if (!params_->on_midstep) {
    /* First update bound crosslinks state and positions */
    UpdateBoundCrosslinks();
    if (sparams_.no_binding == false && sparams_.no_solution_binding == false && sparams_.exist_while_unbound == false && step_count_>sparams_.Add_at){
      /* Calculate implicit binding of crosslinks from solution */
      CalculateBindingFree();
    }
    step_count_+=1;
  } else {
    /* Apply tether forces from doubly-bound crosslinks onto anchored objects.
       We do this every half step only, because the fullstep tether forces are
       handled by the crosslink update on the full step */
    ApplyCrosslinkTetherForces();
    // Add a step to clear neighbors on the midstep- otherwise they only get
    // cleared on the full step and are double-counted
    ClearNeighbors();
  }
  // for (auto it=members_.begin(); it!=members_.end(); ++it) {
  // it->SanityCheck();
  //}
}

void CrosslinkSpecies::UpdateBindRate() {
  if (!use_bind_file_) return;  
  bind_rate_ = 0;
  for (int anchor_index = 0; anchor_index < 2; ++anchor_index) {
    for (auto obj = objs_->begin(); obj != objs_->end(); ++obj) {
      std::string name = (*obj)->GetName();
      // Find if name of object is in the map
      auto bind_param_it = bind_param_map_[anchor_index].find(name);
      if (bind_param_it != bind_param_map_[anchor_index].end()) {
        // obj_amount is area if surface density used, length if linear density used
        double obj_amount = (bind_param_it->second.dens_type == +density_type::linear) 
                            ? (*obj)->GetLength() : (*obj)->GetArea();
        // Do not contribute area/length if object is already occupied and single occupancy is
        // selected
        if (bind_param_it->second.single_occupancy && ((*obj)->GetNAnchored() > 0)) {
          obj_amount = 0;
        }
        bind_rate_ += bind_param_it->second.k_on_s * 
                      bind_param_it->second.bind_site_density * obj_amount;
      }
    }
  }
}

void CrosslinkSpecies::UpdateBoundCrosslinks() {
  n_members_ = 0;
  /* Update anchor positions to their attached meshes and calculate anchor
     forces */
  UpdateBoundCrosslinkForces();
  /* Apply anchor forces on bound objects sequentially */
  ApplyCrosslinkTetherForces();
  /* Update anchor positions from diffusion, walking */
  UpdateBoundCrosslinkPositions();
  /* Remove crosslinks that came unbound */
  if (!sparams_.static_flag) {
    members_.erase(std::remove_if(members_.begin(), members_.end(),
                                  [](Crosslink x) { return x.IsUnbound(); }),
                   members_.end());
  }
  /* Get the number of bound crosslinks so we know what the current
     concentration of free crosslinks is */
  n_members_ = members_.size();
}

/* This must be done sequentially to avoid racy conditions when accessing
   bound object's forces */
void CrosslinkSpecies::ApplyCrosslinkTetherForces() {
  for (auto xlink = members_.begin(); xlink != members_.end(); ++xlink) {
    xlink->ApplyTetherForces();
  }
}

void CrosslinkSpecies::UpdateBoundCrosslinkForces() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  xlink_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  xlink_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    xlink_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks, update_)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto xlink = chunks[i].first; xlink != chunks[i].second; ++xlink) {
        bool init_state = xlink->IsSingly();
        if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
          continue;
        }
        xlink->UpdateCrosslinkForces();
        if (xlink->IsSingly() != init_state) {
          *update_ = true;
        }
      }
    }
  }
#else
  for (xlink_iterator xlink = members_.begin(); xlink != members_.end();
       ++xlink) {
    bool init_state = xlink->IsSingly();
    if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
      continue;
    }
    xlink->UpdateCrosslinkForces();
    if (xlink->IsSingly() != init_state) {
      *update_ = true;
    }
  }
#endif
}

void CrosslinkSpecies::UpdateBoundCrosslinkPositions() {
#ifdef ENABLE_OPENMP
  int max_threads = omp_get_max_threads();
  xlink_chunk_vector chunks;
  chunks.reserve(max_threads);
  size_t chunk_size = members_.size() / max_threads;
  xlink_iterator cur_iter = members_.begin();
  for (int i = 0; i < max_threads - 1; ++i) {
    xlink_iterator last_iter = cur_iter;
    std::advance(cur_iter, chunk_size);
    chunks.push_back(std::make_pair(last_iter, cur_iter));
  }
  chunks.push_back(std::make_pair(cur_iter, members_.end()));

#pragma omp parallel shared(chunks, update_)
  {
#pragma omp for
    for (int i = 0; i < max_threads; ++i) {
      for (auto xlink = chunks[i].first; xlink != chunks[i].second; ++xlink) {
        bool init_state = xlink->IsSingly();
        if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
          continue;
        }
        xlink->UpdateCrosslinkPositions();
        /* Xlink is no longer bound, return to solution */
        if (xlink->IsUnbound()) {
          if (sparams_.static_flag) {
            Logger::Error("Static crosslinks became unbound");
          }
          *update_ = true;
          /* If a crosslink enters or leaves the singly state, we need to
           * update xlink interactors */
        } else if (xlink->IsSingly() != init_state) {
          *update_ = true;
        }
      }
    }
  }
#else
  for (xlink_iterator xlink = members_.begin(); xlink != members_.end();
       ++xlink) {
    bool init_state = xlink->IsSingly();
    if (sparams_.static_flag && init_state && xlink->GetNNeighbors() == 0) {
      continue;
    }
    xlink->UpdateCrosslinkPositions();
    /* Xlink is no longer bound, return to solution */
    if (xlink->IsUnbound()) {
      if (sparams_.static_flag) {
        Logger::Error("Static crosslinks became unbound");
      }
      *update_ = true;
      /* If a crosslink enters or leaves the singly state, we need to update
       * xlink interactors */
    } else if (xlink->IsSingly() != init_state) {
      *update_ = true;
    }
  }
#endif
}

void CrosslinkSpecies::CleanUp() { members_.clear(); }

void CrosslinkSpecies::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->Draw(graph_array);
  }
}

void CrosslinkSpecies::ReadSpecs() {
  if (ispec_file_.eof()) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (!ispec_file_.is_open()) {
    Logger::Warning("ERROR. Spec file unexpectedly not open! Exiting early.");
    early_exit = true;
    return;
  }
  n_members_ = -1;
  ispec_file_.read(reinterpret_cast<char *>(&n_members_), sizeof(int));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (n_members_ == -1) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (n_members_ == 0) {
    members_.clear();
  } else if (n_members_ != members_.size()) {
    Crosslink xlink(rng_.GetSeed());
    xlink.Init(&sparams_);
    xlink.InitInteractionEnvironment(&lut_, tracker_, bound_curr_);
    xlink.SetObjSize(obj_size_);
    xlink.SetBindRate(&bind_rate_);
    xlink.SetBindParamMap(&bind_param_map_);
    xlink.SetSID(GetSID());
    members_.resize(n_members_, xlink);
  }
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ReadSpec(ispec_file_);
  }
}

void CrosslinkSpecies::ClearNeighbors() {
  // Clear all anchor neighborlists.
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ClearNeighbors();
  }
}

const int CrosslinkSpecies::GetDoublyBoundCrosslinkNumber() const {
  int num = 0;
  for (const auto &xl : members_) {
    if (xl.IsDoubly())
      ++num;
  }
  return num;
}

const double CrosslinkSpecies::GetConcentration() const {
  return sparams_.concentration;
}

//Set pointer to global check for cross flag so individual crosslinks can change it
void CrosslinkSpecies::SetGlobalCheckForCrossPointer(bool* global_check_for_cross) {
  global_check_for_cross_ = global_check_for_cross;
}

const double CrosslinkSpecies::GetRCutoff() const { return lut_.getLUCutoff(); }
