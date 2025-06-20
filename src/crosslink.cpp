#include "cglass/crosslink.hpp"
#include <iostream>

Crosslink::Crosslink(unsigned long seed) : Object(seed) {
  SetSID(species_id::crosslink);
}

void Crosslink::Init(crosslink_parameters *sparams) {
  sparams_ = sparams;
  name_ = sparams_->name;
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
  use_bind_file_ = sparams_->anchors[0].bind_file.compare("none");
  Anchor anchor1(rng_.GetSeed());
  Anchor anchor2(rng_.GetSeed());
  anchors_.push_back(anchor1);
  anchors_.push_back(anchor2);
  anchors_[0].Init(sparams_, 0);
  anchors_[1].Init(sparams_, 1);
  SetSingly(bound_anchor_);
  free_diffuse_kick_ = sqrt(2*sparams_->diffusion_free*params_->delta);
  Logger::Trace("Initializing crosslink %d with anchors %d and %d", GetOID(),
                anchors_[0].GetOID(), anchors_[1].GetOID());
}

void Crosslink::InitInteractionEnvironment(LookupTable *lut, Tracker *tracker, 
                                           std::map<Sphere *, std::pair<std::vector<double>, 
                                           std::vector<std::pair<Anchor*, std::string> > > > *bound_curr) { 
  lut_ = lut;
  tracker_ = tracker;
  bound_curr_ = bound_curr; 
  anchors_[0].SetBoundCurr(bound_curr);
  anchors_[1].SetBoundCurr(bound_curr);
}

/* Function used to set anchor[0] position etc to xlink position etc */
void Crosslink::UpdatePosition() {}

void Crosslink::GetAnchors(std::vector<Object *> &ixors) {
  if (IsUnbound())
    return;
  if (!static_flag_) {
    ixors.push_back(&anchors_[bound_anchor_]);
  }
  if (IsDoubly()) {
    ixors.push_back(&anchors_[(int)!bound_anchor_]);
  }
}

//Perform kinetic monte carlo step for an unbound crosslinker, only used if exist_while_unbound flag is on
void Crosslink::FreeKMC () {
  double roll = rng_.RandomUniform();
  std::vector<double> prob_list;
  double total_bind_prop = 0;
  double bind_factor = 0;
  const std::vector<const Sphere*>& sphere_nbr_list = anchors_[0].GetNeighborListMemSpheres();
  for (int i = 0; i < sphere_nbr_list.size(); ++i) {
    double const *const  anchor_pos = anchors_[0].GetPosition();
    double const *const sphere_pos = sphere_nbr_list[i]->GetPosition();
    double squared_distance = 0;
    for (int i = 0; i < params_->n_dim; ++i) {
      squared_distance += SQR(anchor_pos[i]-sphere_pos[i]);
    }
    double distance = sqrt(squared_distance);
    double prob_factor = 0;
    if (distance <= sparams_->f_to_s_radius) {
      prob_factor = sparams_->f_to_s_rate*delta_;
    }
    total_bind_prop+=prob_factor;
    prob_list.push_back(prob_factor);
  }
  const std::vector<const Rod*>& rod_nbr_list = anchors_[0].GetNeighborListMemRods();
  for (int i = 0; i < rod_nbr_list.size(); ++i) {
    double const *const anchor_pos = anchors_[0].GetPosition();
    double const *const rod_pos = rod_nbr_list[i]->GetPosition();
    Rod *bind_obj = anchors_[0].GetRodNeighbor(i);
    double obj_length = bind_obj->GetLength();
    double perp_distance_sq = 0;
    for (int i = 1; i < params_->n_dim; ++i) {
      perp_distance_sq += SQR(anchor_pos[i]-rod_pos[i]);
    }
    double perp_distance = sqrt(perp_distance_sq);
    double prob_factor = 0;
    double bind_radius=sparams_->f_to_s_radius;
    //Calculate chance to bind if rod is within bind radius of crosslink, currently only set up for rods parallel to the xaxis like in the protrusion simulations
    if (perp_distance <= bind_radius) {
      //Don't bind if past the edges of the rod
      if (anchor_pos[0]<=(rod_pos[0]-.5*obj_length-bind_radius) || anchor_pos[0]>=(rod_pos[0]+.5*obj_length+bind_radius) ) {
        prob_factor=0;
      }
      //Binding to interior section of rod
      else if (anchor_pos[0]>=(rod_pos[0]-.5*obj_length+bind_radius) || anchor_pos[0]<=(rod_pos[0]+.5*obj_length+bind_radius) ) {
        double length_in_sphere=2*sqrt(SQR(bind_radius)-SQR(perp_distance));
        prob_factor=length_in_sphere*sparams_->f_to_s_rate*bind_site_density_*delta_;
      }
      //Binding where the crosslink is at the edge of the rod
      else if (anchor_pos[0]<=(rod_pos[0]-.5*obj_length+bind_radius) || anchor_pos[0]>=(rod_pos[0]-.5*obj_length-bind_radius) ) {
        double sphere_rod_int=anchor_pos[0]+sqrt(SQR(bind_radius)-SQR(perp_distance));
        double end_of_rod=rod_pos[0]-bind_radius;
        double length_in_sphere=sphere_rod_int-end_of_rod;
        prob_factor=length_in_sphere*sparams_->f_to_s_rate*bind_site_density_*delta_;
      }
      else if (anchor_pos[0]>=(rod_pos[0]+.5*obj_length-bind_radius) || anchor_pos[0]<=(rod_pos[0]+.5*obj_length+bind_radius) ) {
        double sphere_rod_int=anchor_pos[0]-sqrt(SQR(bind_radius)-SQR(perp_distance));
        double end_of_rod=rod_pos[0]+bind_radius;
        double length_in_sphere=end_of_rod-sphere_rod_int;
        prob_factor=length_in_sphere*sparams_->f_to_s_rate*bind_site_density_*delta_;
      }
    }
    total_bind_prop+=prob_factor;
    prob_list.push_back(prob_factor);
    }
  if (total_bind_prop>=roll) {
    for (int i = 0; i < sphere_nbr_list.size(); ++i) {
      if(prob_list[i]>=roll) {
        Sphere *bind_obj = anchors_[0].GetSphereNeighbor(i);
        (*bound_curr_)[bind_obj].first.push_back(prob_list[i]);
        std::pair<Anchor*, std::string> anchor_and_bind_type;
        anchor_and_bind_type.first = &anchors_[0];
        anchor_and_bind_type.second = "free to single";
        (*bound_curr_)[bind_obj].second.push_back(anchor_and_bind_type);
        anchors_[0].SetCrosslinkPointer(this);
        Logger::Trace("Unbound to singly bind added for object %i", this->GetOID());
        return;
        }
      else{
         roll -= prob_list[i];
      }
    }
    for (int i = 0; i < rod_nbr_list.size(); ++i) {
      if(prob_list[i]>=roll) {
        Rod *bind_obj = anchors_[0].GetRodNeighbor(i);
        double obj_length = bind_obj->GetLength();
        double const *const anchor_pos = anchors_[0].GetPosition();
        double const *const rod_pos = rod_nbr_list[i]->GetPosition();
        double bind_lambda = anchor_pos[0]-rod_pos[0];
        bind_lambda+=.5*obj_length;
        if (bind_lambda > obj_length) {
          bind_lambda = obj_length;
        } else if (bind_lambda < 0) {
          bind_lambda = 0;
        }
        anchors_[0].AttachObjLambda(bind_obj, bind_lambda);
        SetSingly(0);
        Logger::Trace("Unbound to singly bind added for object %i", anchors_[0].GetOID());
        return;
      }
      else{
        roll -= prob_list[i];
      }
    }
  }
}

/* Perform kinetic monte carlo step of protein with 1 head attached. */
void Crosslink::SinglyKMC() {
  double roll = rng_.RandomUniform();
  int head_bound = 0;
  // Set up KMC objects and calculate probabilities
  double unbind_prob = anchors_[bound_anchor_].GetOffRate() * delta_;
  if (static_flag_ || sparams_->no_solution_binding) {
    unbind_prob = 0;
  }
  tracker_->TrackSU(unbind_prob);
  int n_neighbors_rod = anchors_[bound_anchor_].GetNNeighborsRod();
  int n_neighbors_sphere = anchors_[bound_anchor_].GetNNeighborsSphere();
  int n_neighbors = n_neighbors_rod + n_neighbors_sphere;
 
  /* Initialize KMC calculation */
  KMC<Rod, Sphere> kmc_bind(anchors_[bound_anchor_].pos, n_neighbors_rod, n_neighbors_sphere, delta_, lut_);

  /* Initialize periodic boundary conditions */
  kmc_bind.SetPBCs(n_dim_, space_->n_periodic, space_->unit_cell);

  /* Calculate probability to bind */
  double kmc_bind_prob = 0;
  std::vector<double> bind_factors(n_neighbors);
  const std::vector<const Rod*>& rod_nbr_list = anchors_[bound_anchor_].GetNeighborListMemRods();
  const std::vector<const Sphere*>& sphere_nbr_list = anchors_[bound_anchor_].GetNeighborListMemSpheres();
  if (use_bind_file_) {
    for (int i = 0; i < rod_nbr_list.size(); ++i) {
      std::string name = rod_nbr_list[i]->GetName();
      bind_factors[i] = bind_param_map_->at((int)!bound_anchor_)[name].k_on_d
                        * bind_param_map_->at((int)!bound_anchor_)[name].bind_site_density;
    }
    for (int i = 0; i < sphere_nbr_list.size(); ++i) {
      std::string name = sphere_nbr_list[i]->GetName();
      bind_factors[rod_nbr_list.size() + i] =
               bind_param_map_->at((int)!bound_anchor_)[name].k_on_d
               * bind_param_map_->at((int)!bound_anchor_)[name].bind_site_density;
    }
  } else {
  double bind_factor_rod = anchors_[(int)!bound_anchor_].GetOnRate() * bind_site_density_;
  double bind_factor_sphere = anchors_[(int)!bound_anchor_].GetOnRate() * bind_site_density_;

  /* Fill vector of bind factors with rod factors, then sphere factors */
    std::fill(bind_factors.begin(), bind_factors.begin() + n_neighbors_rod, bind_factor_rod);
    std::fill(bind_factors.begin() + n_neighbors_rod, bind_factors.end(), bind_factor_sphere);
  }
 
  if (n_neighbors > 0) {
    if (!static_flag_ && polar_affinity_ != 1.0) {
      anchors_[bound_anchor_].CalculatePolarAffinity(bind_factors);
    }
    /* Use auto-filter populated with 1's for every neighbor.
    We already guarantee uniqueness, so we won't overcount. */
    kmc_bind.LUCalcTotProbsSD(anchors_[bound_anchor_].GetNeighborListMemRods(), 
                              anchors_[bound_anchor_].GetNeighborListMemSpheres(), 
                              anchors_[bound_anchor_].GetBoundOID(), bind_factors); 

    kmc_bind_prob = kmc_bind.getTotProb();
    tracker_->TrackSD(kmc_bind_prob);
  } 
  // Find out whether we bind, unbind, or neither.
  int head_activate = choose_kmc_double(unbind_prob, kmc_bind_prob, roll);
  // Change status of activated head
  if (head_activate == 0 && anchors_[bound_anchor_].GetChangedThisStep()==false) {
    // Unbind bound head
    // Track unbinding
    tracker_->UnbindSU();
    if (anchors_[bound_anchor_].GetReachedPlusEnd()) {
      anchors_[bound_anchor_].SubtractFilEndProteins(true);
      anchors_[bound_anchor_].SetReachedPlusEnd(false);
    }
    //If unbound crosslinkers are tracked, release crosslink instead of deleting it
    if (free_flag_) {
      SetFree(bound_anchor_);
      anchors_[bound_anchor_].Unbind();
      anchors_[0].SetBound();
      anchors_[1].SetUnbound();
      return;
    }
    else {
      anchors_[bound_anchor_].Unbind();
      SetUnbound();
      Logger::Trace("Crosslink %i with anchor %i came unbound", GetOID(), anchors_[bound_anchor_].GetOID()); 
    }
  } else if (head_activate == 1) {
    //Logger::Info("Single Binding");
    // Bind unbound head
    // Track binding
    tracker_->BindSD();
    /* Position on rod where protein will bind with respect to center of rod,
     * passed by reference */
    double bind_lambda;
    /* Find out which rod we are binding to */
    int i_bind = kmc_bind.whichObjBindSD(bind_lambda, roll);
    if (i_bind < 0) {
      printf("i_bind = %d\nbind_lambda = %2.2f\n", i_bind, bind_lambda);
      Logger::Error("kmc_bind.whichRodBindSD in Crosslink::SinglyKMC"
                    " returned an invalid result!");
    }
    if (anchors_[bound_anchor_].GetReachedPlusEnd()) {
      anchors_[bound_anchor_].SubtractFilEndProteins(false);
      anchors_[bound_anchor_].SetReachedPlusEnd(false);
    }
    if (i_bind < n_neighbors_rod) {
      Rod *bind_obj = anchors_[bound_anchor_].GetRodNeighbor(i_bind);
      double obj_length = bind_obj->GetLength();
      /* KMC returns bind_lambda to be with respect to center of rod. We want 
      it to be specified from the tail of the rod to be consistent */
      bind_lambda += 0.5 * obj_length;
      /* KMC can return values that deviate a very small amount from the true 
      rod length. Bind to ends if lambda < 0 or lambda > bond_length. */
      if (bind_lambda > obj_length) {
        bind_lambda = obj_length;
      } else if (bind_lambda < 0) {
        bind_lambda = 0;
      }
      anchors_[(int)!bound_anchor_].AttachObjLambda(bind_obj, bind_lambda);
      SetDoubly();
      Logger::Trace("Crosslink %d became doubly bound to obj %d", GetOID(),
                  bind_obj->GetOID());
    } else {
      Sphere *bind_obj = anchors_[bound_anchor_].GetSphereNeighbor(i_bind - n_neighbors_rod);
      //if (bind_obj -> GetNAnchored() == 0) {
      (*bound_curr_)[bind_obj].first.push_back(kmc_bind.getProb(i_bind));
      std::pair<Anchor*, std::string> anchor_and_bind_type;
      anchors_[(int)!bound_anchor_].SetCrosslinkPointer(this);
      anchor_and_bind_type.first = &anchors_[(int)!bound_anchor_];
      anchor_and_bind_type.second = "single to double";
      //Add anchor to bound_curr, will be deciding if it binds during knockout
      (*bound_curr_)[bind_obj].second.push_back(anchor_and_bind_type);
      Logger::Trace("Crosslink %d, with anchor %d, became doubly bound to obj %d", GetOID(), anchors_[(int)!bound_anchor_].GetOID(),
                  bind_obj->GetOID());
      //If crosslinkers can't cross check if newly bound crosslinker is crossing
      if (sparams_ -> cant_cross == true) {
        //Logger::Error("cant cross flag needs repairing");
        check_for_cross = true;
        last_bound_ = (int)!bound_anchor_;
        if (*global_check_for_cross_ == true) {
          //Logger::Warning("Two crosslinks bound during same time step");
        } else {
          *global_check_for_cross_ = true;
        }
      }
      //bind_obj -> IncrementNAnchored();
     //}
    }
  }
}

//If a crosslinker was bound during this step it is flaged to check to see if 
//it is crossing another crosslinker
bool Crosslink::ReturnCheckForCross() {
  return check_for_cross;
}

//After a crosslinker is checked if its crossing set check_for_cross back to false
void Crosslink::SetCheckForCross() {
  check_for_cross = false;
  *global_check_for_cross_ = false;
}

void Crosslink::SetGlobalCheckForCross(bool* global_check_for_cross){
  global_check_for_cross_ = global_check_for_cross;
}

//Get index of most recently bound anchor (if crosslinker unbinds do to crosslinking 
//most recently bound anchor needs to unbind)
int Crosslink::GetLastBound() {
  return last_bound_;
}

//Unbind anchor if unbinding is due to crosslinkers crossing
void Crosslink::UnbindCrossing() {
  Logger::Trace("Crosslinker %f came unbound because it was crossing another crosslinker", GetOID());
  int head_activate = last_bound_;
  ClearNeighbors();
  UpdateXlinkState();
  tracker_ -> UnbindDS();
  anchors_[head_activate].Unbind();
  SetSingly((int)!head_activate);
  SetCheckForCross();
}
  

/* Perform kinetic monte carlo step of protein with 2 heads of protein
 * object attached. */
void Crosslink::DoublyKMC() {
  /* Calculate force-dependent unbinding for each head */
  double tether_stretch = length_ - rest_length_;
  // For asymmetric springs apply second spring constant for compression
  double e_dep = e_dep_factor_;
  double f_dep = fdep_length_;
  if (sparams_->motor_off==true){
    e_dep=0;
    f_dep=1;
  }
  //double factor = 1;

  //if (force_[1]<0){   
        //factor=(1-force_[0]/30.38*.375);
  //   factor=(1-force_[1]*0.066);
  //   } else {
  //   factor=(1+force_[1]*0.066);
        //factor=(1+force_[0]/30.38*.8);
  //    }
  //if (force_[0]<0){   
        //factor=(1-force_[0]/30.38*.375);
  //   factor+=0;
  //   } else {
  //   factor+=force_[0]*0.033;
        //factor=(1+force_[0]/30.38*.8);
  //    }
  //f_dep *= factor;
  double s_fac =0;
  if (k_spring_compress_ >= 0 && tether_stretch < 0) {
    e_dep *= 0.5 * k_spring_compress_ * SQR(tether_stretch);
  } else {
    e_dep *= 0.5 * k_spring_ * SQR(tether_stretch);
  }
  double separation = anchors_[0].GetPosition()[1]-anchors_[1].GetPosition()[1];
  if (separation<0) {
    separation*=-1;
  }
  if (tether_stretch < 0 and separation > 0.6) {
    s_fac=1;
  }
  if (tether_stretch > 0 and separation < 0.6) {
    s_fac=1;
  }
  //printf("[%f,%f]", f_dep);

  std::vector<double> unbind_prob;
  double k1_0 = sparams_->k_0; //0.005; //0.05;
  double x1 = 0.04; //0.0587;
  double y1 =0.018; //0.027; //0.018;
  double k2_0 = 0.356;
  double x2=0.0072;
  double y2=0;
  for (int i = 0; i < 2; i++) {
    //Logger::Info("Off is %f", anchors_[i].GetOffRate());
    if (sparams_->motor_off==false){
      unbind_prob.push_back(anchors_[i].GetOffRate() * delta_ * exp(e_dep + f_dep));
    } else{
      if (off_calc_fx_<0) {
        off_calc_fx_=off_calc_fx_*-1;
      }
      double k1=k1_0*exp(off_calc_fx_*x1+off_calc_fy_*y1); 
      double k2=k2_0*exp(off_calc_fx_*x2+off_calc_fy_*y2);
      //double off_total=k1*k2/(k1+k2)*s_fac;
      double off_total=k1*s_fac;
      //printf("Force x = %f, force y = %f, off = %f \n", off_calc_fx_, off_calc_fy_, off_total);

      unbind_prob.push_back(off_total* delta_);
    }
  }
 /*
  std::vector<double> unbind_prob;
  double off_rate = sparams_->k_0;
  double base_off = 0.0249;
  double off_total = 0;
  for (int i = 0; i < 2; i++) {
    //Logger::Info("Off is %f", anchors_[i].GetOffRate());
    if (sparams_->motor_off==false){
      unbind_prob.push_back(anchors_[i].GetOffRate() * delta_ * exp(e_dep + f_dep));
    } else{
      if (off_calc_fx_<0) {
        off_total=base_off-off_calc_fx_*0.000469;
      } else{
        off_total=base_off+off_calc_fx_*0.00257;
      }
      off_total*=off_rate/base_off;
      unbind_prob.push_back(off_total* delta_*s_fac);
    }
  }
  */


  tracker_->TrackDS(unbind_prob[0]); // Richelle modify to track full prob
  double roll = rng_.RandomUniform();
  int head_activate = -1;
  if (static_flag_) {
    head_activate = choose_kmc_double(0, unbind_prob[1], roll);
  } else {
    // Probability of unbinding follows a poisson process but assume that only
    // one head can unbind during a time step.
    head_activate = choose_kmc_double(unbind_prob[0], unbind_prob[1], roll);
  }
  if (head_activate > -1 && anchors_[head_activate].GetChangedThisStep() == false) {
    tracker_->UnbindDS();
    Logger::Trace("Doubly-bound crosslink %d came unbound from %d", GetOID(),
                  anchors_[head_activate].GetBoundOID());
    anchors_[head_activate].Unbind();
    SetSingly((int)!head_activate);
    //SetCheckForCross();
  }
}

void Crosslink::CalculateBinding() {
  if (IsSingly()) {
    SinglyKMC();
  } else if (IsDoubly()) {
    DoublyKMC();
  } else if (IsFree()) {
    FreeKMC();
  }
  ClearNeighbors();
}

/* Only singly-bound crosslinks interact */
void Crosslink::GetInteractors(std::vector<Object *> &ixors) {
  ClearNeighbors();
  if (IsSingly() || IsFree()) {
    ixors.push_back(&anchors_[bound_anchor_]);
  }
}

//Get how far anchors are along the filaments
std::vector<double> Crosslink::GetAnchorS() {
  std::vector<double> s_values_;
  s_values_.push_back(anchors_[0].GetRecS());
  s_values_.push_back(anchors_[1].GetRecS());
  return s_values_;
}

//Get how far anchors are along the filaments
bool Crosslink::StillDBound() {
  bool double_bound;
  double_bound = (anchors_[0].StillBound() && anchors_[1].StillBound());
  return double_bound;
}

//Get the IDs of the filaments that the receptors the anchors are connected to are connected to
std::vector<int> Crosslink::GetReceptorPCIDs() {
  std::vector<int> rec_ids_;
  rec_ids_.push_back(anchors_[0].GetPCID());  
  rec_ids_.push_back(anchors_[1].GetPCID()); 
  return rec_ids_;
} 

void Crosslink::ClearNeighbors() { anchors_[bound_anchor_].ClearNeighbors(); }

void Crosslink::UpdateAnchorsToMesh() {
  if (!IsFree()) {
  anchors_[0].UpdateAnchorPositionToMesh();
  anchors_[1].UpdateAnchorPositionToMesh();
  }
}

void Crosslink::UpdateAnchorPositions() {
  if (!sparams_->stationary_flag and !IsFree()) {
    anchors_[0].UpdatePosition();   
    anchors_[1].UpdatePosition();
  }
}

void Crosslink::ApplyTetherForces() {
  if (!IsDoubly())
    return;
  anchors_[0].ApplyAnchorForces();
  anchors_[1].ApplyAnchorForces();
}

void Crosslink::UpdateCrosslinkForces() {
  if (!IsFree()) {
    /* Update anchor positions in space to calculate tether forces */
    UpdateAnchorsToMesh();
  }
  /* Check if an anchor became unbound due to diffusion, etc */
  UpdateXlinkState();
  if (!IsFree()) {
    /* If we are doubly-bound, calculate and apply tether forces */
    CalculateTetherForces();
  }
}

void Crosslink::UpdateCrosslinkPositions() {
  if (!IsFree()) {
    /* Have anchors diffuse/walk along mesh */
    UpdateAnchorPositions();
  } else {
    DiffuseFree();
    UpdatePeriodic();
  }
  /* Check if an anchor became unbound do to diffusion, etc */
  UpdateXlinkState();
  if (sparams_->no_binding == false){
    /* Check for binding/unbinding events using KMC */
    CalculateBinding();
  }
}

//Diffuse for unbound, tracked crosslinkers
void Crosslink::DiffuseFree() {;
  for (int i = 0; i < params_->n_dim; ++i) {
    double dr = free_diffuse_kick_ * rng_.RandomNormal(1);
    position_[i]+=dr;
  }
  if (params_->reflect_at_boundary == true) {
    ReflectAtBoundary();
  }
  //ZeroForce();
  //Move anchors to new position
  anchors_[0].BindToPosition(position_);
  anchors_[1].BindToPosition(position_);
}

void Crosslink::ReflectAtBoundary() {
  if (space_->type != +boundary_type::sphere && space_->type != +boundary_type::protrusion){
    Logger::Error("Reflecting boundary only set up sphere and protrusion boundaries");
  }
  //Check is crosslinker is outside of the simulation boundry
  MinimumDistance mindist;
  bool outside = mindist.CheckOutsideBoundary(*this);
  if (outside == true) {
    double new_radius = mindist.GetNewRadius();
    if (space_->type == +boundary_type::sphere) {
      //Scale position to new position inside boundry
      double new_radius = mindist.GetNewRadius();
      double radius_ratio = new_radius/(sqrt( SQR(position_[1]) + SQR(position_[2])+ SQR(position_[0])) );
      position_[1] = position_[1]*radius_ratio;
      position_[2] = position_[2]*radius_ratio;
      position_[0] = position_[0]*radius_ratio;
    } else if (space_->type == +boundary_type::protrusion){
      double new_radius = mindist.GetNewRadius();
      double pro_start = space_->pro_start; //x coordinate where the protrusion starts
      //If past the end of the cylinder
      if (position_[0] < (pro_start-space_->pro_length)) {
        double new_x = mindist.GetNewXValue();
        position_[0] = new_x;
      }
      //If radially outside the protrusion section
      else if (position_[0] < pro_start) {
        double radius_ratio = new_radius/(sqrt( SQR(position_[1]) + SQR(position_[2])) );
        position_[1] = position_[1]*radius_ratio;
        position_[2] = position_[2]*radius_ratio;
      }
      //If inside the sphere section
      else {
        double radius_ratio = new_radius/(sqrt( SQR(position_[1]) + SQR(position_[2])+ SQR(position_[0])) );
        position_[1] = position_[1]*radius_ratio;
        position_[2] = position_[2]*radius_ratio;
        position_[0] = position_[0]*radius_ratio;
      }
    }
  }
}

/* This function ensures that singly-bound crosslinks have anchor[0] bound and
   anchor[1] unbound. */
void Crosslink::UpdateXlinkState() {
  if (!anchors_[0].IsBound() && !anchors_[1].IsBound()) {
    SetUnbound();
    return;
  }
  if (IsDoubly() && !anchors_[1].IsBound()) {
    SetSingly(0);
  } else if (IsDoubly() && !anchors_[0].IsBound()) {
    SetSingly(1);
  }
  if (IsSingly() && (anchors_[1].IsBound() && anchors_[0].IsBound())) {
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
  //printf("orientation = %f, %f, force = mag %f, %f, %f, anchor locations are %f, %f\n", orientation_[0], orientation_[1], tether_force_, force_[0], force_[1], anchors_[0].GetPosition()[1], anchors_[1].GetPosition()[1]);
  anchors_[0].AddForce(force_);
  anchors_[1].SubForce(force_);
 
  if (anchors_[0].GetPosition()[1]>anchors_[1].GetPosition()[1]){
    off_calc_fx_=-force_[0];
    off_calc_fy_=force_[1];
  } else {
    off_calc_fx_=force_[0];
    off_calc_fy_=-force_[1];
  }

  //If anchors are double bound and hopping between receptors then we need to calculate 
  //the energy to the neighboring sites to determine hopping rates
  if ((anchors_[0].IsBoundToSphere() == true) && (anchors_[1].IsBoundToSphere() == true)) {
    Sphere* sphere_zero_ = anchors_[0].GetBoundPointer();
    Sphere* sphere_one_ = anchors_[1].GetBoundPointer();
 
    //Calculate current length of crosslinker
    Interaction ix_zo(sphere_zero_, sphere_one_);
    MinimumDistance mindist_zo;
    mindist_zo.ObjectObject(ix_zo);
    double length_zo_ = sqrt(ix_zo.dr_mag2);
    if (anchors_[0].IsWalker() == false) {
      anchors_[0].SetCrosslinkLength(length_zo_);
    }
    if (anchors_[1].IsWalker() == false) {
      anchors_[1].SetCrosslinkLength(length_zo_);
    }
    
    if(anchors_[1].IsWalker() == false) {
      //receptor in plus direction of anchor zero
      Sphere* one_plus_n_ = sphere_one_ -> GetPlusNeighbor();
      Sphere* one_minus_n_ = sphere_one_ -> GetMinusNeighbor();

      //If anchor is at plus end of microtubule it won't have a plus neighbor
      if (one_plus_n_ == nullptr) {
        anchors_[1].SetLengthAtPlus(-1);
      }
      else {
        //Get length between anchor zero and plus neighbor of anchor one
        Interaction ix_zp(sphere_zero_, one_plus_n_);
        MinimumDistance mindist_zp;
        mindist_zp.ObjectObject(ix_zp);
        double length_zp_ = sqrt(ix_zp.dr_mag2);
        anchors_[1].SetLengthAtPlus(length_zp_);  
      }

      //If anchor is at minus end of microtubule it won't have a plus neighbor
      if (one_minus_n_ == nullptr) {
        anchors_[1].SetLengthAtMinus(-1);
      }
      else { 
        //Get length between anchor zero and minus neighbor of anchor one
        Interaction ix_zm(sphere_zero_, one_minus_n_);
        MinimumDistance mindist_zm;
        mindist_zm.ObjectObject(ix_zm);
        double length_zm_ = sqrt(ix_zm.dr_mag2);
        anchors_[1].SetLengthAtMinus(length_zm_);
      }
   }

   if (anchors_[0].IsWalker() == false) {
     Sphere* zero_plus_n_ = sphere_zero_ -> GetPlusNeighbor();
     Sphere* zero_minus_n_ = sphere_zero_ -> GetMinusNeighbor();

     //If anchor is at plus end of microtubule it won't have a plus neighbor
     if (zero_plus_n_ == nullptr) {
       anchors_[0].SetLengthAtPlus(-1);
     }

     else { 

       //Get length between anchor one and plus neighbor of anchor zero
       Interaction ix_op(sphere_one_, zero_plus_n_);
       MinimumDistance mindist_op;
       mindist_op.ObjectObject(ix_op);
       double length_op_ = sqrt(ix_op.dr_mag2); 
       anchors_[0].SetLengthAtPlus(length_op_);
      }
      //If anchor is at plus end of microtubule it won't have a plus neighbor
      if (zero_minus_n_ == nullptr) {
        anchors_[0].SetLengthAtMinus(-1);
      }
      else {  
        //Get length between anchor one and minus neighbor of anchor zero
        Interaction ix_om(sphere_one_, zero_minus_n_);
        MinimumDistance mindist_om;
        mindist_om.ObjectObject(ix_om);
        double length_om_ = sqrt(ix_om.dr_mag2); 
        anchors_[0].SetLengthAtMinus(length_om_);
      }
    }
  }
  // If one anchor induces catastrophe and the other is attached to a filament, depolymerize
  // attached filament.
  if (anchors_[0].InducesCatastrophe() && anchors_[1].AttachedToFilamentLastBond()) {
    anchors_[1].InduceCatastrophe();
  } else if (anchors_[1].InducesCatastrophe() && anchors_[0].AttachedToFilamentLastBond()) {
    anchors_[0].InduceCatastrophe();
  }
  
  // Update xlink's position (for drawing)
  UpdatePeriodic();
}

/* Attach a crosslink anchor to object in a random fashion */
void Crosslink::AttachObjRandom(std::pair<Object*, int> obj_index) {
  /* Attaching to random obj implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if ((obj_index.first->GetShape() == +shape::rod) || (obj_index.first->GetShape() == +shape::sphere)) {
    bound_anchor_ = obj_index.second;
    anchors_[obj_index.second].AttachObjRandom(obj_index.first);
    SetCompID(obj_index.first->GetCompID());
  } else {
    Logger::Error("Crosslink binding to %s shaped objects not yet implemented.", obj_index.first->GetShape()._to_string());
  }
}


/* Attach a crosslink anchor to object in a random fashion */
void Crosslink::AttachSphere(std::pair<Sphere*, int> obj_index) {
  /* Attaching to random obj implies first anchor binding from solution, so
   * this crosslink should be new and should not be singly or doubly bound */
  if (obj_index.first->GetShape() == +shape::sphere) {
    bound_anchor_ = obj_index.second;
    anchors_[obj_index.second].AttachObjRandom(obj_index.first);
    SetCompID(obj_index.first->GetCompID());
  } else {
    Logger::Error("Crosslink binding to %s shaped objects not yet implemented.", obj_index.first->GetShape()._to_string());
  }
}

//Attatch crosslinker connected to two receptors, only used when crosslinkers
//start doubly bound
void Crosslink::DoublyCenter(Object* receptor_one, Object* receptor_two) {
  anchors_[0].AttachObjCenter(receptor_one);
  SetCompID(receptor_one->GetCompID());
  anchors_[1].AttachObjCenter(receptor_two);
  SetCompID(receptor_two->GetCompID());
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

void Crosslink::SetSingly(int bound_anchor) {
  state_ = bind_state::singly;
  bound_anchor_ = bound_anchor;
  SetAnchorStates();
}

//SetFree is for when crosslinkers unbind and they are being tracked while unbound
void Crosslink::SetFree(bool a) {
  state_ = bind_state::free;
  const double* anchor_position_= anchors_[a].GetPosition();
  for (int i = 0; i < params_->n_dim; ++i) {
    position_[i] = anchor_position_[i];
  }
  anchors_[0].SetState(bind_state::unbound);
  anchors_[1].SetState(bind_state::unbound);
  bound_anchor_ = 0;
}

//SetUnbound is for crosslinkers unbind and aren't being tracked while unbound
void Crosslink::SetUnbound() {
  state_ = bind_state::unbound;
  SetAnchorStates();
}

void Crosslink::SetAnchorStates() {
  anchors_[0].SetState(state_);
  anchors_[1].SetState(state_);
};

const bool Crosslink::IsFree() const {  return state_ == +bind_state::free; }
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
  bool is_free = IsFree();
  ospec.write(reinterpret_cast<char *>(&is_free), sizeof(bool));
  ospec.write(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ospec.write(reinterpret_cast<char *>(&length_), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ospec.write(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  ospec.write(reinterpret_cast<char *>(&oid_), sizeof(int));
  anchors_[0].WriteSpec(ospec);
  anchors_[1].WriteSpec(ospec);
}

void Crosslink::WriteSpecTextHeader(std::fstream &otext) {
  otext << "is_doubly diameter length position[0] position[1] position[2] "
        << "orientation[0] orientation[1] orientation[2] oid_" << std::endl;
}

void Crosslink::ConvertSpec(std::fstream &ispec, std::fstream &otext) {
  if (ispec.eof())
    return;
  bool is_doubly;
  bool is_free;
  double diameter, length;
  double position[3], orientation[3];
  int oid;
  // Read in all data from spec file ispec
  ispec.read(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  ispec.write(reinterpret_cast<char *>(&is_free), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&diameter), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&length), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&position[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&orientation[i]), sizeof(double));
  }
  ispec.read(reinterpret_cast<char *>(&oid), sizeof(int));
  // Write out data to SpecText file otext
  otext << is_doubly << " " << is_free << " " << diameter << " " << length << " " << position[0] << " "
        << position[1] << " " << position[2] << " " << orientation[0] 
        << " " << orientation[1] << " " << orientation[2] << " " << oid << std::endl;
  // Convert anchor data
  Anchor::WriteSpecTextHeader(otext);
  for (int i = 0; i < 2; ++i) {
    Anchor::ConvertSpec(ispec, otext);
  }
}

void Crosslink::ReadSpec(std::fstream &ispec) {
  if (ispec.eof())
    return;
  SetSingly(bound_anchor_);
  bool is_doubly;
  ispec.read(reinterpret_cast<char *>(&is_doubly), sizeof(bool));
  bool is_free = IsFree();
  ispec.write(reinterpret_cast<char *>(&is_free), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&diameter_), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&length_), sizeof(double));
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  }
  ispec.read(reinterpret_cast<char *>(&oid_), sizeof(int));
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
  Logger::Trace("Reloading anchor from checkpoint with cid %d",
                anchors_[bound_anchor_].GetCompID());
  if (IsDoubly()) {
    Logger::Trace("Reloading anchor from checkpoint with cid %d",
                  anchors_[(int)!bound_anchor_].GetCompID());
  }
}

const double Crosslink::GetDrTot() {
  if (IsSingly() || IsFree()) {
    return anchors_[bound_anchor_].GetDrTot();
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
  anchors_[bound_anchor_].ZeroDrTot();
  if (IsDoubly()) {
    anchors_[(int)!bound_anchor_].ZeroDrTot();
  }
}

void Crosslink::InsertAt(double const *const new_pos, double const *const u) {
  static_flag_ = true;
  anchors_[bound_anchor_].InsertAt(new_pos, u);
  anchors_[bound_anchor_].SetBound();
  anchors_[bound_anchor_].SetStatic(true);
  SetSingly(bound_anchor_);
}

//Insert crosslink unbound
void Crosslink::InsertFree(double const *const new_pos, double const *const u) {
  free_flag_ = true;
  static_flag_ = false;
  anchors_[0].SetStatic(false);
  anchors_[1].SetStatic(false);
  SetFree(0);
  for (int i = 0; i < params_->n_dim; ++i) {
    position_[i] = new_pos[i];
  }
}

void Crosslink::SetBindParamMap(std::vector<std::map<std::string, bind_params> > *bind_param_map) {
  bind_param_map_ = bind_param_map;
  anchors_[0].SetBindParamMap(bind_param_map_);
  anchors_[1].SetBindParamMap(bind_param_map_);
}

void Crosslink::SetObjSize(double *obj_size) {
  if (!obj_size) Logger::Error("Crosslink received nullptr obj_size");
  obj_size_ = obj_size;
  anchors_[0].SetObjSize(obj_size);
  anchors_[1].SetObjSize(obj_size);
}

void Crosslink::SetBindRate(double *bind_rate) {
  if (!bind_rate) Logger::Warning("Crosslink received nullptr bind_rate");
  bind_rate_ = bind_rate;
  anchors_[0].SetBindRate(bind_rate);
  anchors_[1].SetBindRate(bind_rate);
}

const double* const Crosslink::GetObjSize() {
  if (!obj_size_) Logger::Warning("Crosslink sent nullptr obj_size");
  return obj_size_;
}

const int Crosslink::GetNNeighbors() const {
  return anchors_[bound_anchor_].GetNNeighbors();
}

const double *const Crosslink::GetPosition() {
  return anchors_[bound_anchor_].GetPosition();
}

const double *const Crosslink::GetOrientation() {
  return anchors_[bound_anchor_].GetOrientation();
}
