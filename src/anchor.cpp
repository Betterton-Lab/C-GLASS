#include "cglass/anchor.hpp"

Anchor::Anchor(unsigned long seed) : Object(seed) {
  SetSID(species_id::crosslink);
}

void Anchor::Init(crosslink_parameters *sparams, int index) {
  index_ = index;
  sparams_ = sparams;
  name_ = sparams_->name;
  diameter_ = sparams_->diameter;
  color_ = sparams_->anchors[index_].color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  static_flag_ = false; // Must be explicitly set to true by Crosslink
  Unbind();
  step_direction_ =
      (sparams_->step_direction == 0 ? 0 : SIGNOF(sparams_->step_direction));
  max_velocity_s_ = sparams->anchors[index_].velocity_s;
  max_velocity_d_ = sparams->anchors[index_].velocity_d;
  //If anchor parameters are not set and at their default values use the crosslinker diffusion.
  //If anchor parameters are set use the anchor parameters. 
  if (sparams_->anchors[index_].diffusion_s == -1) {
    diffusion_s_ = sparams->diffusion_s;
  } else {
    diffusion_s_ = sparams_->anchors[index_].diffusion_s;
  }
  if (sparams_->anchors[index_].diffusion_d == -1) {
    diffusion_d_ = sparams->diffusion_d;
  } else {
    diffusion_d_ = sparams_->anchors[index_].diffusion_d;
  }
  use_partner_ = sparams_->anchors[index_].use_partner;
  k_on_s_ = sparams_->anchors[index_].k_on_s;
  partner_on_s_ = sparams_->anchors[index_].partner_on_s;
  k_on_d_ = sparams_->anchors[index_].k_on_d;
  partner_on_d_ = sparams_->anchors[index_].partner_on_d;
  k_off_s_ = sparams_->anchors[index_].k_off_s;
  k_off_d_ = sparams_->anchors[index_].k_off_d;
  plus_end_pausing_ = sparams_->plus_end_pausing;
  minus_end_pausing_ = sparams_->minus_end_pausing;
  f_stall_ = sparams_->f_stall;
  force_dep_vel_flag_ = sparams_->force_dep_vel_flag;
  polar_affinity_ = sparams_->polar_affinity;
  use_bind_file_ = sparams_->anchors[index_].bind_file.compare("none");
  assert(polar_affinity_ >= 0 && polar_affinity_ <= 1);
  SetDiffusion();
  changed_this_step_ = false;
}

void Anchor::SetBindParamMap(std::vector<std::map<std::string, bind_params> > *bind_param_map) {
  bind_param_map_ = bind_param_map;
}

double const Anchor::GetMeshLambda() { return mesh_lambda_; }

double const Anchor::GetBondLambda() { return bond_lambda_; }

void Anchor::SetRodLambda(double l) { bond_lambda_ = l; }
void Anchor::SetMeshLambda(double ml) { mesh_lambda_ = ml; }

void Anchor::SetDiffusion() {
  // Solve them now so you do not have to keep taking sqrts
  kick_amp_s_ = sqrt(2. * diffusion_s_ * delta_);
  kick_amp_d_ = sqrt(2. * diffusion_d_ * delta_);
}

void Anchor::SetReachedPlusEnd(bool plus_end) {
  reached_plus_end_ = plus_end;
}

void Anchor::UpdateAnchorPositionToMesh() {
  if (!bound_ || static_flag_)
    return;
  if (!mesh_) {
    UpdateAnchorPositionToObj();
    return;
  }

  /* Use the mesh to determine the rod lengths. The true rod lengths fluctuate
     about this, but should be considered approximations to the ideal mesh. */
  mesh_length_ = mesh_->GetTrueLength();
  /* Use current position along mesh (mesh_lambda) to determine whether the
     anchor fell off the mesh due to dynamic instability */
  if (!CheckMesh())
    return;
  // Now figure out which rod we are on in the mesh according to mesh_lambda
  rod_ = mesh_->GetBondAtLambda(mesh_lambda_);
  if (rod_->GetType() == +obj_type::bond) bond_ = dynamic_cast<Bond*>(rod_);

  // Figure out how far we are from the rod tail: bond_lambda
  if (!CalcRodLambda()) {
    return;
  }
  // Update anchor position with respect to bond
  UpdateAnchorPositionToObj();
}

bool Anchor::CalcRodLambda() {
  if (!rod_) {
    Logger::Error("Attempted to calculate rod lambda when not attached to"
                  " rod!");
  }
  bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
  rod_length_ = bond_->GetLength();
  if (bond_lambda_ < 0) {
    Bond *bond = bond_->GetNeighborBond(0);
    if (bond) {
      rod_ = bond;
      bond_ = bond;
      bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
      rod_length_ = bond_->GetLength();
    } else if (minus_end_pausing_) {
      bond_lambda_ = 0;
    } else {
      Unbind();
      return false;
    }
  } else if (bond_lambda_ > rod_length_) {
    Bond *bond = bond_->GetNeighborBond(1);
    if (bond) {
      rod_ = bond;
      bond_ = bond;
      bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
      rod_length_ = rod_->GetLength();
    } else if (plus_end_pausing_) {
      bond_lambda_ = rod_length_;
      if (!reached_plus_end_) {
        AddFilEndProteins();
        reached_plus_end_ = true;
      }
    } else {
      Unbind();
      return false;
    }
  }
  // assert(bond_lambda_ >= 0 && bond_lambda_ <= rod_length_);
  if (bond_lambda_ < -1e-6 || bond_lambda_ > rod_length_ + 1e-6) {
    Logger::Error(
        "Bond lambda out of expected range in UpdateAnchorPositionToMesh, "
        "rod_num: %d, mesh lambda: %2.8f, mesh length: %2.8f, bond lambda: "
        "%2.8f, bond length: %2.8f",
        bond_->GetBondNumber(), mesh_lambda_, mesh_length_, bond_lambda_,
        rod_length_);
  }
  return true;
}
void Anchor::UpdatePosition() {
  // Currently only bound anchors diffuse/walk (no explicit unbound anchors)
  bool diffuse = GetDiffusionConst() > 0 ? true : false;
  bool walker = abs(GetMaxVelocity()) > input_tol ? true : false;
  if (!bound_ || static_flag_ || !(rod_ || sphere_) || (!diffuse && !walker)) {
    return;
  }
  // Diffuse or walk along the mesh, updating mesh_lambda
  // If anchors are walking directy along a rod
  if (rod_) {
    if (diffuse) {
      Diffuse();
     if (!CheckMesh())
       return;
    }
    if (walker) {
      Walk();
      CheckMesh();
    }
  }
  //If anchors are hopping between receptors
  if (sphere_) {
    double discrete_diffusion_ = 0;
    double discrete_velocity_ = 0;
    if (diffuse) {
      discrete_diffusion_ = DiscreteDiffuse();
    }
    if (walker) {
      discrete_velocity_ = DiscreteWalk();
      DecideToStepMotor(discrete_diffusion_, discrete_velocity_);
    }
    if (!walker) {
       switch (state_){
       case +bind_state::singly:
         DecideToStepCrosslinkSingly(discrete_diffusion_);
       break;
       case +bind_state::doubly:
         DecideToStepCrosslink(discrete_diffusion_);
       break; 
       }
    }  
  }
}

void Anchor::DecideToStepMotor(double discrete_diffusion_, double discrete_velocity_) {
  double roll = rng_.RandomUniform();
  double vel_ = discrete_velocity_;
  double step_size_ = sphere_ -> GetStepSize();
  double chance_forward_ = 0;
  double chance_back_ = 0;
  double D = discrete_diffusion_;
  chance_forward_ = vel_/step_size_*delta_;
  chance_back_ = 0;
  if (chance_forward_>roll) {
    PrepareToStepForward(chance_forward_);
  }
  else if (chance_back_>(1-roll)) {
    PrepareToStepBack(chance_back_);
  }
  if ( (chance_back_+chance_forward_) > 1) {
    //Logger::Warning("Chance of anchor, %i, hopping sites greater than one chance back %f, chance forward %f", this->GetOID(), chance_back_, chance_forward_); 
  }
}

void Anchor::DecideToStepCrosslinkSingly(double discrete_diffusion_) {
  double roll = rng_.RandomUniform();  
  double D = discrete_diffusion_;
  double chance_forward_ = 0;
  double chance_back_ = 0;
  double step_size_ = sphere_ -> GetStepSize();

  //Calculating the chance the crosslinker will diffuse toward plus end
  //If distance has been set to -1 this means the anchor is already at
  //the plus end of the microtubule and can't diffuse towards the plus end
  chance_forward_ = (D/pow(step_size_,2))*delta_;
  chance_back_ = (D/pow(step_size_,2))*delta_;
  if (chance_forward_>roll) {
    PrepareToStepForward(chance_forward_);
  }
  else if (chance_back_>(1-roll)) {
    PrepareToStepBack(chance_back_);
  }
  if ( (chance_back_+chance_forward_) > 1) {
    //Logger::Warning("Chance of anchor ,%i,hopping sites greater than one (chance back %f, chance forward %f)", this->GetOID(),chance_back_, chance_forward_);
  }
}

/* Crosslinker decides if it's going to step forward, backwards, or not step*/
void Anchor::DecideToStepCrosslink(double discrete_diffusion_) {
  double roll = rng_.RandomUniform();  
  double D = discrete_diffusion_;
  double chance_forward_ = 0;
  double chance_back_ = 0;
  double k = 0;
  double k_stretch = sparams_ -> k_spring;
  double k_comp = sparams_ -> k_spring_compress;
  double r_l = sparams_ -> rest_length;
  double r_a = sparams_ -> rest_angle;
  double k_align = sparams_-> k_align;
  double step_size_ = sphere_ -> GetStepSize();
  double del_l = cl_length_ - r_l;
  double del_a = cl_angle_ - r_a; 
  if (del_l<0){
    k=k_comp;
  } else {
    k=k_stretch;
  }

  //Calculate the current energy
  double energy = 0.5 * k * pow((del_l), 2);
  double energy_angle = 0.5 * k_align * pow((del_a), 2);
  energy += energy_angle;


  //Calculating the chance the crosslinker will diffuse toward plus end
  //If distance has been set to -1 this means the anchor is already at
  //the plus end of the microtubule and can't diffuse towards the plus end
  if (distance_to_plus_ == -1) {
    chance_forward_ = 0;
  }
  else {
    //Calculate the energy change between current length/angle and length/angle at plus
    del_l = distance_to_plus_ - r_l;
    del_a = angle_at_plus_ - r_a;
    if (del_l<0){
      k=k_comp;
    } else {
      k=k_stretch;
    }
    double energy_to_plus = 0.5 * k * pow(del_l, 2); 
    double energy_angle_p = 0.5 * k_align * pow((del_a), 2);
    //printf("angular energy %f, angle p energy %f \n", energy_angle, energy_angle_p);
    energy_to_plus += energy_angle_p;  
    double e_change_to_p = energy_to_plus - energy;
    //Calculate Boltz factor assuming lambda = 1/2
    double boltz_factor_p = exp(-0.5 * e_change_to_p);
    //modify diffusion rate using Boltzmann factor
    double plus_diffusion = boltz_factor_p * D;
    //See  equation 7.30 and 7.31 from "Molecular motors: thermodynamics and
    //the random walk" (Thomas et al. 2001). Equation rearranged to solve for
    //k+ and k-
    chance_forward_ = (plus_diffusion/pow(step_size_,2))*delta_;
  } 

  //Calculate the chance the crosslinker will diffuse toward minus end
  //If at minus end chance to diffuse further is zero
  if (distance_to_minus_ == -1) {
    chance_back_ = 0;
  }
  else {
    del_l = distance_to_minus_ - r_l;
    del_a = angle_at_minus_ - r_a;
    if (del_l<0){
      k=k_comp;
    } else {
      k=k_stretch;
    }
    //Calculate the energy change between current length and length at minus
    double energy_to_minus = 0.5 * k * pow((del_l), 2); 
    double energy_angle_m = 0.5 * k_align * pow((del_a), 2);
    //printf("angle m energy %f, angles, at base %f, plue %f, at minus %f \n", energy_angle_m, cl_angle_, angle_at_plus_, angle_at_minus_);
    energy_to_minus+= energy_angle_m;   //Calculate Boltz factor assuming lambda = 1/2
    double e_change_to_m = energy_to_minus - energy;
    double boltz_factor_m = exp(-0.5 * e_change_to_m);
    //modify diffusion rate using Boltzmann factor
    double minus_diffusion = boltz_factor_m * D;
    //See  equation 7.30 and 7.31 from "Molecular motors: thermodynamics and
    //the random walk" (Thomas et al. 2001). Equation rearranged to solve for
    //k+ and k-
    chance_back_ = (minus_diffusion/pow(step_size_,2))*delta_;
  } 
  //printf("Chance of anchor  %i hopping sites (chance back %f, chance forward %f) \n", this->GetOID(),chance_back_, chance_forward_);
 //printf("cl angle is %f \n", cl_angle_);

 //Calculate if crosslinker steps forward or backwords
 if (chance_forward_>roll) {
    PrepareToStepForward(chance_forward_);
  }
  else if (chance_back_>(1-roll)) {
    PrepareToStepBack(chance_back_);
  }
  if ( (chance_back_+chance_forward_) > 1) {
    //Logger::Warning("Chance of anchor ,%i,hopping sites greater than one (chance back %f, chance forward %f)", this->GetOID(),chance_back_, chance_forward_);
  }
}

//Add Anchor that has decided to step forward to bound_curr
void Anchor::PrepareToStepForward(double prob) {
  Sphere* next_receptor_ = nullptr;
  next_receptor_ = sphere_->GetPlusNeighbor();
  bool single_occupancy = true;
  if (next_receptor_ != NULL) {
    std::string name = next_receptor_->GetName();
    single_occupancy = bind_param_map_->at(index_)[name].single_occupancy;
  } else if (plus_end_pausing_== false){
    Unbind(); 
  }
  //If receptor is trying to move to an open receptor move.
  //If null that means the receptor is on edge of filament already.
  //If NAnchored is 0, the next receptor isn't occupied.
  if((next_receptor_ != NULL) && (next_receptor_ -> GetNAnchored() == 0 || single_occupancy == false)){
    (*bound_curr_)[next_receptor_].first.push_back(prob);
    std::pair<Anchor*, std::string> anchor_and_bind_type;
    anchor_and_bind_type.first = this;
    anchor_and_bind_type.second = "forward step";
    (*bound_curr_)[next_receptor_].second.push_back(anchor_and_bind_type);
    SetChangedThisStep();
    Logger::Trace("Anchor %i, added to bound curr (Stepping Forward)", this->GetOID());
  }
}

//Add anchor that has decided to step back to bound_curr
void Anchor::PrepareToStepBack(double prob) {
 Sphere* last_receptor_ = nullptr;
 last_receptor_ = sphere_->GetMinusNeighbor();
 bool single_occupancy = true;
  if (last_receptor_ != NULL) {
   std::string name = last_receptor_->GetName();
    single_occupancy = bind_param_map_->at(index_)[name].single_occupancy;
  } else if (minus_end_pausing_== false){
    Unbind(); 
  }
  //If receptor is trying to move to an open receptor move
  //If null that means the receptor is on edge of tube already
  if((last_receptor_ != NULL) && (last_receptor_ -> GetNAnchored() == 0 || single_occupancy == false)){
   (*bound_curr_)[last_receptor_].first.push_back(prob);
   std::pair<Anchor*, std::string> anchor_and_bind_type;
   anchor_and_bind_type.first = this;
   anchor_and_bind_type.second = "back step";
   (*bound_curr_)[last_receptor_].second.push_back(anchor_and_bind_type);
   SetChangedThisStep();
   Logger::Trace("Anchor %i, added to bound_curr (Stepping Back)", this->GetOID());
  }
}

//Set the distance to the plus neighbor of the other head of
//the crosslinker
void Anchor::SetLengthAtPlus(double distance) {
  distance_to_plus_ = distance;
}

//Set the distance to the minus neighbor of the other head of
//the crosslinker
void Anchor::SetLengthAtMinus(double distance) {
  distance_to_minus_ = distance;
}

//Set current length of crosslink anchor is a part off
void Anchor::SetCrosslinkLength(double cl_length) {
  cl_length_ = cl_length;
}


//Set the angle to the plus neighbor of the other head of
//the crosslinker
void Anchor::SetAngleAtPlus(double angle) {
  angle_at_plus_ = angle;
}

//Set the angle to the minus neighbor of the other head of
//the crosslinker
void Anchor::SetAngleAtMinus(double angle) {
  angle_at_minus_ = angle;
}

//Set current length of crosslink anchor is a part off
void Anchor::SetCrosslinkAngle(double cl_angle) {
  cl_angle_ = cl_angle;
}

//Set pointer to the crosslink anchor is a part of
void Anchor::SetCrosslinkPointer(Object* cl_pointer) {
  cl_pointer_ = cl_pointer;
  Logger::Trace("crosslink  %d pointer set for anchor %d", cl_pointer_->GetOID(), this->GetOID());
}

//Anchor steps in the minus direction
void Anchor::StepBack() {
  Sphere* last_receptor_ = nullptr;
  last_receptor_ = sphere_->GetMinusNeighbor();
  bool single_occupancy = true;
  if (last_receptor_ != NULL) { 
   std::string name = last_receptor_->GetName();
    single_occupancy = bind_param_map_->at(index_)[name].single_occupancy; 
  }
  //If receptor is trying to move to an open receptor move
  //If null that means the receptor is on edge of filament already
  if((last_receptor_ != NULL) && (last_receptor_ -> GetNAnchored() == 0 || single_occupancy == false)){
    Unbind();
    AttachObjCenter(last_receptor_);
  }
}

//Anchor steps in the plus direction
void Anchor::StepForward() {
  Sphere* next_receptor_ = nullptr;
  next_receptor_ = sphere_->GetPlusNeighbor();
  bool single_occupancy = true;
  if (next_receptor_ != NULL) {
    std::string name = next_receptor_->GetName();
    single_occupancy = bind_param_map_->at(index_)[name].single_occupancy;
  }
  //If receptor is trying to move to an open receptor move
  //If null that means the receptor is on edge of tube already
  if((next_receptor_ != NULL) && (next_receptor_ -> GetNAnchored() == 0 || single_occupancy == false)){
    Unbind();
    AttachObjCenter(next_receptor_);
  }
}

void Anchor::ApplyAnchorForces() {
  if (!bound_ || static_flag_) {
    return;
  }
  // Spheres don't calculate torque/receptors calculate torque themselves
  if (sphere_) {
    if (sphere_->IsFixed()) return; 
    sphere_->AddForce(force_);
    sphere_->AddTorque(torque_);
  } else if (rod_) {
    double dlambda[3] = {0};
    for (int i = 0; i < n_dim_; ++i) {
      dlambda[i] = (bond_lambda_ - 0.5 * rod_length_) * orientation_[i];
    }
    cross_product(dlambda, force_, torque_, 3);
    rod_->AddForce(force_);
    rod_->AddTorque(torque_);
  } else Logger::Error("Anchor attempted to apply forces to nullptr");
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
  double vel = GetMaxVelocity();
  if (force_dep_vel_flag_) {
    // Only consider projected force in direction of stepping
    double const force_proj =
        step_direction_ * dot_product(n_dim_, force_, orientation_);
    // Linear force-velocity relationship
    double fdep = 1. + (force_proj / f_stall_);
    if (fdep > 1) {
      fdep = 1;
    } else if (fdep < 0) {
      fdep = 0.;
    }
    vel *= fdep;
  }
  double dr = step_direction_ * vel * delta_;
  mesh_lambda_ += dr;
  // Should this also add to bond lambda?
}

double Anchor::DiscreteWalk() {
  double vel = GetMaxVelocity();
   if (force_dep_vel_flag_) {
    //Want oriontation of filament sphere is on, not orientation of sphere
    double const *const rod_orientation_ = (sphere_ ->GetPCObjectForSphere()) -> GetOrientation();
    double const force_proj =
      step_direction_ * dot_product(n_dim_, force_, rod_orientation_);
    double fdep = 1. + (force_proj / f_stall_);
    //double const *const position = this->GetPosition();
    //Logger::Warning("rod orientation x is %f, force x is %f, for anchor %i at x location %f", rod_orientation_[0], force_proj, this->GetOID(), position[0]);
    if (fdep >1) {
      fdep = 1;
    } else if (fdep<0) {
      fdep = 0.;
    }
    vel *= fdep;
  }
  return vel;
}

double Anchor::GetRodOrientation() {
  double const *const rod_orientation_ = (sphere_ ->GetPCObjectForSphere()) -> GetOrientation();
  return rod_orientation_[0];
}

// Check that the anchor is still located on the filament mesh
// Returns true if anchor is still on the mesh, false otherwise
bool Anchor::CheckMesh() {
  // Check if we moved off the mesh tail
  if (mesh_lambda_ < 0) {
    // Stick to mesh ends if we have end_pausing
    if (minus_end_pausing_) {
      mesh_lambda_ = 0;
    } else {
      // Otherwise, unbind
      Unbind();
      return false;
    }
  } else if (mesh_lambda_ > mesh_length_) {
    // Same thing for walking off the head
    if (plus_end_pausing_) {
      mesh_lambda_ = mesh_length_;
    } else {
      Unbind();
      return false;
    }
  }
  return true;
}

void Anchor::Unbind() {
  Logger::Trace("Anchor %i unbound", this->GetOID());
  Logger::Trace("DeAttathed to %i", sphere_);
  if (static_flag_) {
    Logger::Error("Static anchor attempted to unbind");
  }
  if (sphere_) sphere_->DecrementNAnchored();
  AddBackBindRate();
  bound_ = false;
  rod_ = nullptr;
  sphere_ = nullptr;
  comp_ = nullptr;
  mesh_ = nullptr;
  bond_ = nullptr;
  mesh_n_bonds_ = -1;
  rod_length_ = -1;
  bond_lambda_ = -1;
  mesh_lambda_ = -1;
  active_ = false;
  reached_plus_end_ = false;
  if (sparams_->exist_while_unbound == false) {
    ClearNeighbors();
  }
  ZeroForce();
  SetCompID(-1);
  std::fill(position_, position_ + 3, 0.0);
  std::fill(orientation_, orientation_ + 3, 0.0);
}

// Get the bind rate just for the object attached to this anchor
double Anchor::CalcSingleBindRate() {
  double single_bind_rate = 0.0;
  Object* o;
  if (sphere_) o = sphere_;
  else if (rod_) o = rod_;
  // Will sometimes be called w/ unbound anchor during initialization/clearing steps
  else return 0.0;
  std::string name = o->GetName();
  if (bind_param_map_->at(index_)[name].single_occupancy) {
    double obj_amount = (bind_param_map_->at(index_)[name].dens_type == +density_type::linear) 
                        ? o->GetLength() : o->GetArea();
    // Sum both anchor rates because during binding and unbinding the bind rate will 
    // increased for both floating anchors.
    for (int i = 0; i < 2; i++) {
      single_bind_rate += bind_param_map_->at(index_)[name].k_on_s 
                          * bind_param_map_->at(index_)[name].bind_site_density * obj_amount;
    }
  }
  return single_bind_rate;
}

// Increase the bind rate if an occupied object became unbound
void Anchor::AddBackBindRate() {
  if (use_bind_file_ && bind_rate_) {
    *bind_rate_ += CalcSingleBindRate();
  } else {
    if (sphere_) *obj_size_ += sphere_->GetArea();
  }
}

void Anchor::Diffuse() {
  // Motion from thermal kicks
  double dr = GetKickAmplitude() * rng_.RandomNormal(1);
  bool walker = abs(GetMaxVelocity()) > input_tol ? true: false;

  // Force dependence diffusion depends on the mobility (D/kBT) and the force
  // applied along the direction of the filament. Force velocity relation for
  // motors is already taken into account for walkers in Walk().
  if (force_dep_vel_flag_ && !walker) {
    double force_proj = dot_product(n_dim_, force_, orientation_);
    dr += GetDiffusionConst() * force_proj * delta_;
  }
  mesh_lambda_ += dr;
  // Should this also add to bond lambda?
}

double Anchor::DiscreteDiffuse() {
  return GetDiffusionConst();
} 

void Anchor::UpdateAnchorPositionToObj() {
  if (rod_) {
    if (rod_->IsFixed()) return;
    double const *const rod_position = rod_->GetPosition();
    double const *const rod_orientation = rod_->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      orientation_[i] = rod_orientation[i];
      position_[i] = rod_position[i] -
                     (0.5 * rod_length_ - bond_lambda_) * rod_orientation[i];
    }
  } else if (sphere_) {
    if (sphere_->IsFixed()) return;
    /* When attached to sphere, position is just the sphere position */
    std::copy(sphere_->GetPosition(), sphere_->GetPosition() + n_dim_, position_);
  } else {
    Logger::Error("Anchor tried to change position but wasn't bound to sphere_ or rod_.");
  }
  UpdatePeriodic();
}

/*Creates Vector that has different binding rates for parallel and anti-parallel
 * bonds*/
void Anchor::CalculatePolarAffinity(std::vector<double> &doubly_binding_rates) {
  double orientation[3]; 
  if (rod_) std::copy(rod_->GetOrientation(), rod_->GetOrientation() + 3, orientation);
  else if (sphere_) std::copy(sphere_->GetOrientation(), sphere_->GetOrientation() + 3, orientation);
  else
    Logger::Error("Rod and sphere are nullptr in Anchor::CalculatePolarAffinity"); 
  for (int i = 0; i < doubly_binding_rates.size(); ++i) {
    Object *obj = neighbors_.GetNeighbor(i);
    if (obj->GetShape() == +shape::sphere) continue;
    double const *const i_orientation = obj->GetOrientation();
    double alignment = dot_product(n_dim_, orientation, i_orientation);
    if (alignment < 0) {
      doubly_binding_rates[i] *= polar_affinity_;
    }
  }
}

void Anchor::Draw(std::vector<graph_struct *> &graph_array) {
  if ((!bound_) && (state_ != +bind_state::free))
    return;
  std::copy(scaled_position_, scaled_position_ + 3, g_.r);
  for (int i = space_->n_periodic; i < n_dim_; ++i) {
    g_.r[i] = position_[i];
  }
  std::copy(orientation_, orientation_ + 3, g_.u);
  g_.color = color_;
  g_.diameter = .32; //diameter_;
  g_.length = length_;
  g_.draw = draw_;
  graph_array.push_back(&g_);
}

void Anchor::AttachObjRandom(Object *o) {
  switch (o->GetShape()) {
    case shape::rod: {
      double length = o->GetLength();
      double lambda = length * rng_.RandomUniform();
      AttachObjLambda(o, lambda);
      break;
    }
    case shape::sphere: {
      // Automatically let spheres have occupancy 1 if no bind file is used.
      if (!use_bind_file_) {
        if (o->GetNAnchored() > 0) {
          Logger::Error("Xlink tried to bind to already occupied sphere!");
        }
        *obj_size_ -= o->GetArea();
      }
      AttachObjCenter(o); 
      break;
    }
    default: {
      Logger::Error("Crosslink binding to %s type objects not yet implemented in" 
                    " Anchor::AttachObjRandom", o->GetType()._to_string());
    }
  }
  if (use_bind_file_ && bind_rate_) {
    *bind_rate_ -= CalcSingleBindRate();
  }
}
    
void Anchor::AttachObjLambda(Object *o, double lambda) {
  if (o->GetShape() != +shape::rod) {
    Logger::Error(
        "Crosslink binding to %s objects not implemented in "
        "AttachObjLambda.", o->GetShape()._to_string());
  }
  if (use_bind_file_) SetRatesFromBindFile(o->GetName());
  rod_ = dynamic_cast<Rod *>(o);
  bond_ = dynamic_cast<Bond *>(o);
  if (bond_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a bond!");
  }
  comp_ = dynamic_cast<Composite *>(rod_->GetCompPtr());
  if (comp_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a composite!");
  }
  if (comp_->GetCompType() == +comp_type::mesh) {
    mesh_ = dynamic_cast<Mesh *>(rod_->GetCompPtr());
    if (mesh_ == nullptr) {
      Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
    }
    mesh_n_bonds_ = mesh_->GetNBonds();
    mesh_length_ = mesh_->GetTrueLength();
  }
  rod_length_ = rod_->GetLength();
  bond_lambda_ = lambda;

  if (bond_lambda_ < 0 || bond_lambda_ > rod_length_) {
    printf("bond_lambda: %2.2f\n", bond_lambda_);
    Logger::Error("Lambda passed to anchor does not match length of "
                  "corresponding bond! lambda: %2.2f, rod_length: %2.2f ",
                  bond_lambda_, rod_length_);
  }

  /* Distance anchor is relative to entire mesh length */
  mesh_lambda_ = bond_->GetMeshLambda() + bond_lambda_;
  SetCompID(rod_->GetCompID());
  UpdateAnchorPositionToObj();
  ZeroDrTot();
  bound_ = true;
}

/* Attach object in center of site. Site binding likelihood weighted by 
 * surface area, but binding places crosslinks in center regardless. */
void Anchor::AttachObjCenter(Object *o) {
  //Logger::Info("Attathed to %i", o);
  //Logger::Info("Anchor %i unbound", this->GetOID());
  o->IncrementNAnchored();
  if (use_bind_file_) SetRatesFromBindFile(o->GetName());
  if (o->GetShape() != +shape::sphere) {
    Logger::Error(
        "Crosslink binding to non-sphere objects not implemented in "
        "AttachObjCenter.");
  }
  sphere_ = dynamic_cast<Sphere *>(o);
  if (sphere_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a sphere!");
  }
  comp_ = dynamic_cast<Composite *>(sphere_->GetCompPtr());
  if (comp_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a composite!");
  }
  if (comp_->GetCompType() == +comp_type::mesh) {
    mesh_ = dynamic_cast<Mesh *>(sphere_->GetCompPtr());
    if (mesh_ == nullptr) {
      Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
    }
  }

  mesh_lambda_ = -1; // not used for sites
  SetCompID(sphere_->GetCompID());
  std::copy(sphere_->GetPosition(), sphere_->GetPosition() + 3, position_); 
  std::copy(sphere_->GetOrientation(), sphere_->GetOrientation() + 3, orientation_);
  UpdatePeriodic();
  ZeroDrTot();
  bound_ = true;
}

void Anchor::AttachObjMeshLambda(Object *o, double mesh_lambda) {
  if (o->GetShape() != +shape::rod) {
    Logger::Error(
        "Crosslink binding to non-rod objects not allowed in "
        "AttachObjMeshLambda.");
  }
  rod_ = dynamic_cast<Rod *>(o);
  bond_ = dynamic_cast<Bond *>(o);
  if (bond_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a bond!");
  }
  if (rod_->GetType() == +obj_type::bond) bond_ = dynamic_cast<Bond *>(o);
  comp_ = dynamic_cast<Composite *>(rod_->GetCompPtr());
  if (comp_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a composite!");
  }
  if (comp_->GetCompType() == +comp_type::mesh) {
    mesh_ = dynamic_cast<Mesh *>(rod_->GetCompPtr());
    if (mesh_ == nullptr) {
      Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
    }
  }
  Logger::Trace("Attaching anchor %d to comp %d", GetOID(), comp_->GetCompID());

  bound_ = true;
  mesh_lambda_ = mesh_lambda;
  mesh_n_bonds_ = -1;
  UpdateAnchorPositionToMesh();
  if (!bound_) {
    Logger::Error(
        "Updating anchor to mesh from checkpoint resulted in an unbound "
        "anchor");
  }
  SetCompID(rod_->GetCompID());
  ZeroDrTot();
}

void Anchor::AttachObjMeshCenter(Object *o) {
  o->IncrementNAnchored();
  if (o->GetType() != +obj_type::site) {
    Logger::Error(
        "Crosslink binding to non-bond objects not allowed in "
        "AttachObjMeshCenter.");
  }
  sphere_ = dynamic_cast<Sphere *>(o);
  if (sphere_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a sphere!");
  }
  comp_ = dynamic_cast<Composite *>(sphere_->GetCompPtr());
  if (comp_ == nullptr) {
    Logger::Error("Object ptr passed to anchor was not referencing a composite!");
  }
  if (comp_->GetCompType() == +comp_type::mesh) {
    mesh_ = dynamic_cast<Mesh *>(sphere_->GetCompPtr());
    if (mesh_ == nullptr) {
      Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
    }
  }
  
  Logger::Trace("Attaching anchor %d to comp %d", GetOID(), comp_->GetCompID());

  bound_ = true;
  bond_lambda_ = 0;
  mesh_n_bonds_ = -1;
  SetCompID(sphere_->GetCompID());
  ZeroDrTot();
}

// Save binding parameters based on species name
void Anchor::SetRatesFromBindFile(const std::string &name) {
  k_on_d_ = bind_param_map_->at(index_)[name].k_on_d;
  k_on_s_ = bind_param_map_->at(index_)[name].k_on_s;
  k_off_d_ = bind_param_map_->at(index_)[name].k_off_d;
  k_off_s_ = bind_param_map_->at(index_)[name].k_off_s;
}

void Anchor::BindToPosition(double *bind_pos) {
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] = bind_pos[i];
  }
  UpdatePeriodic();
}

Sphere* Anchor::GetBoundPointer() {
  if (sphere_) {
    return sphere_;
  }
  else {
    Logger::Error("GetBoundPointer only set up for spheres/receptors");
    return nullptr;
  }
}

bool Anchor::IsBound() { return bound_; }
 
bool Anchor::IsBoundToSphere() {
  if (sphere_) {
    return true;
  }
  else {
    return false;
  } 
}
 
int const Anchor::GetBoundOID() {
  if (sphere_) return sphere_->GetOID();
  else if (rod_) return rod_->GetOID();
  else return -1;
}

//Get ID for point cover anchor is on
int Anchor::GetPCID() {
  if (sphere_) {
    return (sphere_ ->GetPCObjectForSphere()) -> GetOID();
  } else {
    Logger::Error("Tried to get PC ID for a non-receptor");
    return 0; 
  }
}

/* Temporary function for setting bound state for singly-bound crosslinks,
   in order to get them to draw while not technically bound to a bond
   ( e.g. rod_ -> null ) */
void Anchor::SetBound() { bound_ = true; }

//Set pointer to map that contains the binding events from this step
//Pointer is set so it can be added to within anchor.cpp 
void Anchor::SetBoundCurr(std::map<Sphere *, std::pair<std::vector<double>, std::vector<std::pair<Anchor*, std::string> > > > *bound_curr) {
  bound_curr_ = bound_curr;
}

//Set if this anchor bound this step, unbinding events won't happen on the same turn
void Anchor::SetChangedThisStep() {
  changed_this_step_ = true;
}

//Set changed this step back to false
void Anchor::ResetChangedThisStep() {
  changed_this_step_ = false;
}

bool Anchor::GetChangedThisStep() {
  return changed_this_step_;
}

void Anchor::AddNeighbor(Object *neighbor) { neighbors_.AddNeighbor(neighbor); }

void Anchor::ClearNeighbors() { neighbors_.Clear(); }

const Object *const *Anchor::GetNeighborListMem() {
  return neighbors_.GetNeighborListMem();
}

const std::vector<const Sphere*>& Anchor::GetNeighborListMemSpheres() {
  return neighbors_.GetNeighborListMemSpheres();
}

const std::vector<const Rod*> &Anchor::GetNeighborListMemRods() {
  return neighbors_.GetNeighborListMemRods();
}

Object *Anchor::GetNeighbor(int i_neighbor) {
  return neighbors_.GetNeighbor(i_neighbor);
}

Sphere *Anchor::GetSphereNeighbor(int i_neighbor) {
  return neighbors_.GetSphereNeighbor(i_neighbor);
}

//Get how far the anchor's receptor is on the filament (currently only set up for parallel and antiparallel) 
double Anchor::GetRecS() {
  if (sphere_) {    
    double const *const rod_orientation_ = (sphere_ ->GetPCObjectForSphere()) -> GetOrientation();
    if (rod_orientation_[0]>0) {
      return sphere_ -> GetSphereS();
    }
    else {
      double s = sphere_-> GetSphereS();
      return -s;
    }
  } else {
    Logger::Error("Anchor is not attatched to Sphere");
  return 0;
  }
}

//Check if anchor is still bound to receptor 
bool Anchor::StillBound() {
  if (sphere_) {    
    return true;
  } else {
    return false;
  }
}

Rod *Anchor::GetRodNeighbor(int i_neighbor) {
  return neighbors_.GetRodNeighbor(i_neighbor);
}
const int Anchor::GetNNeighbors() const { return neighbors_.NNeighbors(); }
const int Anchor::GetNNeighborsRod() const { return neighbors_.NNeighborsRod(); }
const int Anchor::GetNNeighborsSphere() const { return neighbors_.NNeighborsSphere(); }



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
  int attached_comp_id = comp_ != nullptr ? comp_->GetCompID() : -1;
  ospec.write(reinterpret_cast<char *>(&attached_comp_id), sizeof(int));
}

void Anchor::WriteSpecTextHeader(std::fstream &otext) {
  otext << "bound active static_flag position[0] position[1] position[2] "
        << "orientation[0] orientation[1] orientation[2] mesh_lambda "
        << "comp_id" << std::endl;
}

void Anchor::ConvertSpec(std::fstream &ispec, std::fstream &otext) {
  bool bound, active, static_flag;
  double position[3], orientation[3];
  double mesh_lambda;
  int comp_id;
  if (ispec.eof())
    return;
  ispec.read(reinterpret_cast<char *>(&bound), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&active), sizeof(bool));
  ispec.read(reinterpret_cast<char *>(&static_flag), sizeof(bool));
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&position[i]), sizeof(double));
  }
  for (int i = 0; i < 3; ++i) {
    ispec.read(reinterpret_cast<char *>(&orientation[i]), sizeof(double));
  }
  ispec.read(reinterpret_cast<char *>(&mesh_lambda), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&comp_id), sizeof(int));
  otext << bound << " " << active << " " << static_flag << " " << position[0] << " " 
        << position[1] << " " << position[2] << " " << orientation[0] << " " << orientation[1] 
        << " " << orientation[2] << " " << mesh_lambda << " " << comp_id << std::endl;
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
  int attached_comp_id = -1; // Just a place holder at the moment
  ispec.read(reinterpret_cast<char *>(&attached_comp_id), sizeof(int));
  UpdatePeriodic();
  if (active_)
    step_direction_ = -sparams_->step_direction;
}

void Anchor::SetStatic(bool static_flag) { static_flag_ = static_flag; }
void Anchor::SetState(bind_state state) { state_ = state; }
void Anchor::SetUnbound() {bound_ = false;}

const double Anchor::GetOnRate() const {
  switch (state_) {
  case +bind_state::unbound:
    return k_on_s_;
    break;
  case +bind_state::singly:
    return k_on_d_;
    break;
  case +bind_state::doubly:
    Logger::Error(
        "Crosslinker is already doubly bound. No on rate exists because both "
        "anchors are bound");
    return 0;
    break;
  default:
    Logger::Error("State of anchor is not a bind_state enum.");
    return 0;
  }
}

const double Anchor::GetOffRate() const {
  switch (state_) {
  case +bind_state::singly:
    return k_off_s_;
    break;
  case +bind_state::doubly:
    return k_off_d_;
    break;
  case +bind_state::unbound:
    Logger::Error("Crosslinker is already unbound. No off rate exists because "
                  "both anchors are off filament");
    return 0;
    break;
  default:
    Logger::Error("State of anchor is not a bind_state enum.");
    return 0;
  }
}

const double Anchor::GetMaxVelocity() const {
  switch (state_) {
  case +bind_state::singly:
    return max_velocity_s_;
    break;
  case +bind_state::doubly:
    return max_velocity_d_;
    break;
  case +bind_state::unbound:
    // Logger::Error(
    //    "Crosslinker is unbound. Anchors cannot walk if not attached.");
    return 0;
    break;
  default:
    Logger::Error("State of anchor is not a bind_state enum");
    return 0;
  }
}

//Returns whether anchor is walker or not
bool Anchor::IsWalker () { 
  bool walker = abs(GetMaxVelocity()) > input_tol ? true : false;
  return walker;
}

const double Anchor::GetDiffusionConst() const {
  switch (state_) {
  case +bind_state::singly:
    return diffusion_s_;
    break;
  case +bind_state::doubly:
    return diffusion_d_;
    break;
  case +bind_state::unbound:
    // Logger::Error(
    //    "Crosslinker is unbound. Anchors cannot diffuse on objects if not "
    //    "attached.");
    return 0;
    break;
  default:
    Logger::Error("State of anchor is not a bind_state enum");
    return 0;
  }
}

const double Anchor::GetKickAmplitude() const {
  switch (state_) {
  case +bind_state::singly:
    return kick_amp_s_;
    break;
  case +bind_state::doubly:
    return kick_amp_d_;
    break;
  case +bind_state::unbound:
    // Logger::Error(
    //    "Crosslinker is unbound. Anchors cannot diffuse on objects if not "
    //    "attached.");
    return 0;
    break;
  default:
    Logger::Error("State of anchor is not a bind_state enum.");
    return 0;
  }
}

const double* const Anchor::GetObjSize() {
  if (!obj_size_) Logger::Warning("Anchor passed nullptr obj_size");
  return obj_size_;
}

void Anchor::SetObjSize(double* obj_size) {
  if (!obj_size) Logger::Error("Anchor received nullptr obj_size");
  obj_size_ = obj_size;
}

const double* const Anchor::GetBindRate() {
  if (!bind_rate_) Logger::Warning("Anchor passed nullptr bind_rate");
  return bind_rate_;
}

void Anchor::SetBindRate(double* bind_rate) {
  if (!bind_rate) Logger::Warning("Anchor received nullptr bind_rate");
  bind_rate_ = bind_rate;
}

const double Anchor::GetKonS() const {
  return k_on_s_;
}

// Returns true if the anchor is a catastrophe-inducer (it is bound
// to a receptor that had the induce_catastrophe flag checked)
bool Anchor::InducesCatastrophe() {
  if (!sphere_ || !(sphere_->InducesCatastrophe())) return false;
  return true;
}

// Returns true if anchor is attached to a RigidFilament or Filament last bond.
bool Anchor::AttachedToFilamentLastBond() {

  // Check if attached to bond of a filament
  if (!rod_ || !bond_ || !mesh_ || (mesh_->GetSID() != +species_id::filament)) {
    return false;
  }
  // Check if attached to last bond of filament
  Bond *bond = bond_->GetNeighborBond(1);
  if (bond) return false;
  return true;
}

// Returns true if anchor is attached to a RigidFilament or Filament last bond.
bool Anchor::GetReachedPlusEnd() {
  return reached_plus_end_;
}

/* Decrement xlink and neighbor end count if anchor dettached from singly bound filament end */
void Anchor::SubtractFilEndProteins(bool singly) {
  if (use_bind_file_) {
    if (bind_param_map_->at(index_)[bond_->GetName()].use_partner) {
      bond_->DecrementNEndXlinks();
      if (singly) {
        bond_->SubNPartners(bind_param_map_->at(index_)[bond_->GetName()].partner_on_s);
      } else {
        bond_->SubNPartners(bind_param_map_->at(index_)[bond_->GetName()].partner_on_d);
      }
    }
  } else {
    if (use_partner_) {
      bond_->DecrementNEndXlinks();
      bond_->SubNPartners(partner_on_s_);
    }
  }
}

/* Decrement xlink and neighbor end count if anchor dettached from doubly bound filament end */
void Anchor::AddFilEndProteins() {
  if (use_bind_file_) {
    if (bind_param_map_->at(index_)[bond_->GetName()].use_partner) {
      bond_->IncrementNEndXlinks();
      bond_->AddNPartners(bind_param_map_->at(index_)[bond_->GetName()].partner_on_s);
    }
  } else {
    if (use_partner_) {
      bond_->IncrementNEndXlinks();
      bond_->AddNPartners(partner_on_s_);
    }
  }
}

// Depolymerize attached anchor
void Anchor::InduceCatastrophe() {
  Filament* fil = dynamic_cast<Filament*>(mesh_);
  fil->Depolymerize();
}
