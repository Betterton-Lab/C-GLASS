#include "cglass/filament.hpp"

Filament::Filament(unsigned long seed) : Mesh(seed) {
  SetSID(species_id::filament);
}
void Filament::SetParameters() {
  /* Read parameters from filament parameters */
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams_->length;
  /* Bending_stiffness = persistence_length*kT. Bending_stiffness is used in
     equations below for clarity. */
  bending_stiffness_ = sparams_->persistence_length;
  diameter_ = sparams_->diameter;
  max_length_ = sparams_->max_length;
  min_length_ = sparams_->min_length;
  min_bond_length_ = sparams_->min_bond_length;
  if (min_length_ < 2 * min_bond_length_)
    min_length_ = 2 * min_bond_length_;
  dynamic_instability_flag_ = sparams_->dynamic_instability_flag;
  polydispersity_flag_ = sparams_->polydispersity_flag;
  spiral_init_flag_ = sparams_->spiral_init_flag;
  force_induced_catastrophe_flag_ = sparams_->force_induced_catastrophe_flag;
  p_g2s_ = sparams_->f_grow_to_shrink * delta_;
  p_g2p_ = sparams_->f_grow_to_pause * delta_;
  p_s2p_ = sparams_->f_shrink_to_pause * delta_;
  p_s2g_ = sparams_->f_shrink_to_grow * delta_;
  p_p2s_ = sparams_->f_pause_to_shrink * delta_;
  p_p2g_ = sparams_->f_pause_to_grow * delta_;
  v_depoly_ = sparams_->v_depoly;
  v_poly_ = sparams_->v_poly;
  driving_factor_ = sparams_->driving_factor;
  nematic_driving_ = sparams_->nematic_driving;
  p_driving_switch_ = sparams_->nematic_driving_freq * delta_;
  peclet_number_ = sparams_->peclet_number;
  flexure_number_ = sparams_->flexure_number;
  if (dynamic_instability_flag_ &&
      (peclet_number_ >= 0 || flexure_number_ >= 0)) {
    /* Dynamic instability and Peclet number do not mix because a Peclet number
       assumes a fixed length, and the driving factor should not change as a
       filament grows or shrinks */
    Logger::Error("Dynamic instability and Peclet number are not compatible");
  }
  driving_factor_ = sparams_->driving_factor;
  friction_ratio_ = sparams_->friction_ratio;
  zero_temperature_ = params_->zero_temperature; // don't include thermal forces
  eq_steps_ = sparams_->n_equil;
  // Don't bother with equilibriating if this is a reloaded simulation
  if (params_->n_load > 0) {
    eq_steps_ = 0;
  }
  optical_trap_spring_ = sparams_->optical_trap_spring;
  optical_trap_flag_ = sparams_->optical_trap_flag;
  optical_trap_fixed_ = sparams_->optical_trap_fixed;
  cilia_trap_flag_ = sparams_->cilia_trap_flag;
  fic_factor_ = sparams_->fic_factor;
  tip_force_ = 0.0;
  /* Intrinsic curvature is given in the format of d_theta/d_s where s is the
     arc length, then the angle between each bond must be d_theta/d_s *
     bond_length_. The additional factor of 1/2 is due to the fact that
     curvature is the adjusted angle for each bond when calculating the bending
     forces. If using radius of curvature, it should have already been converted
     by filament species initialization*/
  curvature_ = 0.5 * (sparams_->intrinsic_curvature +
                      rng_.RandomNormal(sparams_->intrinsic_curvature_sig));
  if (sparams_->randomize_intrinsic_curvature_handedness) {
    int sign = (rng_.RandomUniform() < 0.5 ? -1 : 1);
    curvature_ = sign * curvature_;
  }
  if (ABS(curvature_) < 0.5 * sparams_->intrinsic_curvature_min) {
    curvature_ = SIGNOF(curvature_) * 0.5 * sparams_->intrinsic_curvature_min;
  }
  flagella_flag_ = sparams_->flagella_flag;
  flagella_freq_ = sparams_->flagella_freq;
  flagella_period_ = sparams_->flagella_period;
  flagella_amplitude_ = sparams_->flagella_amplitude;
  /* Default site in optical trap is the tail */
  trapped_site_ = 0;
  custom_set_tail_ = sparams_->custom_set_tail;
  error_analysis_ = sparams_->error_analysis;
  /* Refine parameters */
  if (dynamic_instability_flag_) {
    /* Since dynamic instability requires a bond number that is a power of two,
       the maximum bond length must be twice the min bond length */
    max_bond_length_ = 2 * min_bond_length_;
  } else {
    /* For a static filament length, any real valued length greater than
     * 2 x min_bond_length (the minimum filament length) can be represented
     * with an integer number of bonds with bond lengths between 1 and 1.5 times
     * the min bond length */
    max_bond_length_ = 1.5 * min_bond_length_;
  }
  if (sparams_->highlight_handedness) {
    draw_ = draw_type::fixed;
    color_ = (curvature_ > 0 ? sparams_->color : sparams_->color + M_PI);
  }
}

void Filament::Init(filament_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  InitFilamentLength();
  Reserve();
  InsertFilament();
  SetDiffusion();

  if (optical_trap_flag_) {
    double const *const r0 = sites_[0].GetPosition();
    std::copy(r0, r0 + 3, optical_trap_pos_);
    double const *const r1 = sites_[1].GetPosition();
    std::copy(r1, r1 + 3, optical_trap_pos2_);
  }
  poly_ = poly_state::grow;
}

/* Returns number of bonds to initialize */
void Filament::InitFilamentLength() {
  if (max_length_ < min_length_) {
    Logger::Warning("Minimum filament length larger than max length -- setting "
                    "max_length_ = min_length_");
    max_length_ = min_length_;
  }

  if (polydispersity_flag_) {
    ExponentialDist expon;
    expon.Init(length_, min_length_, max_length_);
    double roll = rng_.RandomUniform();
    length_ = expon.Rand(roll);
    if (length_ > max_length_ + 1e-6 || length_ < min_length_ - 1e-6) {
      Logger::Error(
          "Filament polydispersity distribution returned a value out of "
          "range:\nAttempted length: %2.2f, Min length: %2.2f, Max length: "
          "%2.2f",
          length_, min_length_, max_length_);
    }
  }
  if (length_ > max_length_) {
    Logger::Warning(
        "Filament length larger than max length -- setting length = "
        "max_length");
    length_ = max_length_;
  } else if (length_ < min_length_) {
    Logger::Warning("Filament length less than min length -- setting length = "
                    "min_length");
    length_ = min_length_;
  }
  /* Calculate number of bonds in filament at initialization */
  if (dynamic_instability_flag_) {
    n_bonds_ = 2;
    bond_length_ = length_ / n_bonds_;
    while (dynamic_instability_flag_ && bond_length_ > max_bond_length_) {
      n_bonds_ *= 2;
      bond_length_ = length_ / n_bonds_;
    }
  } else {
    n_bonds_ = (int)floor(length_ / min_bond_length_);
    /* Minimum number of bonds must be 2 */
    n_bonds_ = (n_bonds_ == 1 ? 2 : n_bonds_);
    bond_length_ = length_ / n_bonds_;
  }
  /* Determine maximum number of bonds we may have */
  if (dynamic_instability_flag_) {
    int max_bonds = (int)ceil(max_length_ / min_bond_length_);
    n_bonds_max_ = 2;
    while (n_bonds_max_ < max_bonds) {
      n_bonds_max_ *= 2;
    }
  } else {
    /* Static bond number */
    n_bonds_max_ = n_bonds_;
  }

  Logger::Trace("Filament initialized with length %2.2f with %d bonds, mesh_id:"
                " %d",
                length_, n_bonds_, GetMeshID());
  true_length_ = length_;
  if (flexure_number_ >= 0) {
    peclet_number_ = flexure_number_ * bending_stiffness_ / length_;
  }
  if (peclet_number_ >= 0) {
    driving_factor_ = peclet_number_ / SQR(length_);
  }
}

void Filament::Reserve() {
  // Initialize mesh
  Mesh::Reserve();
  // Allocate control structures
  int n_sites_max = n_bonds_max_ + 1;
  tensions_.resize(n_sites_max - 1);                    // max_sites -1
  g_mat_lower_.resize(n_sites_max - 2);                 // max_sites-2
  g_mat_upper_.resize(n_sites_max - 2);                 // max_sites-2
  g_mat_diag_.resize(n_sites_max - 1);                  // max_sites-1
  det_t_mat_.resize(n_sites_max + 1);                   // max_sites+1
  det_b_mat_.resize(n_sites_max + 1);                   // max_sites+1
  g_mat_inverse_.resize(n_sites_max - 2);               // max_sites-2
  k_eff_.resize(n_sites_max - 2);                       // max_sites-2
  h_mat_diag_.resize(n_sites_max - 1);                  // max_sites-1
  h_mat_upper_.resize(n_sites_max - 2);                 // max_sites-2
  h_mat_lower_.resize(n_sites_max - 2);                 // max_sites-2
  gamma_inverse_.resize(n_sites_max * n_dim_ * n_dim_); // max_sites*ndim*ndim
  cos_thetas_.resize(n_sites_max - 2);                  // max_sites-2
}

void Filament::InsertFirstBond() {
  if (sparams_->insertion_type.compare("random") == 0) {
    InitRandomBond(diameter_);
  } else if (sparams_->insertion_type.compare("random_nematic") == 0) {
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = (rng_.RandomUniform() > 0.5 ? 1.0 : -1.0);
    InitRandomBondOriented(orientation_, diameter_);
  } else if (sparams_->insertion_type.compare("random_polar") == 0) {
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
    InitRandomBondOriented(orientation_, diameter_);
  } else if (sparams_->insertion_type.compare("centered_oriented") == 0) {
    std::fill(orientation_, orientation_ + 3, 0.0);
    orientation_[n_dim_ - 1] = 1.0;
    for (int i = 0; i < n_dim_; ++i) {
      position_[i] = -0.5 * length_ * orientation_[i];
    }
    InitSiteAt(position_, diameter_);
    AddBondToTip(orientation_, bond_length_);
  } else if (sparams_->insertion_type.compare("centered_random") == 0) {
    rng_.RandomUnitVector(n_dim_, orientation_);
    for (int i = 0; i < n_dim_; ++i) {
      position_[i] = -0.5 * length_ * orientation_[i];
    }
    InitSiteAt(position_, diameter_);
    AddBondToTip(orientation_, bond_length_);
  } else {
    // Assume custom arrangement for now
    InitSiteAt(position_, diameter_);
    AddBondToTip(orientation_, bond_length_);
  }
}

void Filament::InsertFilament() {
  int n_bonds_init = n_bonds_;
  Clear();
  InsertFirstBond();
  SetOrientation(bonds_.back().GetOrientation());
  bool probable_orientation =
      (sparams_->insertion_type.compare("simple_crystal") != 0 &&
       sparams_->insertion_type.compare("random_polar") != 0 &&
       sparams_->insertion_type.compare("random_nematic") != 0);
  for (int i = 0; i < n_bonds_init - 1; ++i) {
    if (probable_orientation) {
      GenerateProbableOrientation();
    }
    AddBondToTip(orientation_, bond_length_);
  }
  for (bond_iterator bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
    bond->SetColor(color_, draw_);
  }
  if (spiral_init_flag_) {
    InitSpiral2D();
  }
  UpdateBondPositions();
  UpdatePrevPositions();
  CalculateAngles();
}

void Filament::InsertAt(const double *const new_pos, const double *const u) {
  Logger::Trace("Inserting filament at [%2.1f, %2.1f, %2.1f] with orientation"
                "[%2.1f, %2.1f, %2.1f]",
                new_pos[0], new_pos[1], new_pos[2], u[0], u[1], u[2]);
  if (custom_set_tail_) {
    std::copy(new_pos, new_pos + n_dim_, position_);
    for (int i = 0; i < n_dim_; ++i) {
      position_[i] += u[i] * length_ * .5;
    }
    Logger::Trace("Insert file locations are the locations of tails");
    RelocateMesh(position_, u);
  } else {
    RelocateMesh(new_pos, u);
  }
  UpdatePrevPositions();
  CalculateAngles();
  SetDiffusion();
  if (optical_trap_flag_) {
    trapped_site_ = 0;
    /* For cilia flag, if filament is oriented along a negative
     * dimension, assume the filament needs to be fixed at the plus
     * end */
    if (cilia_trap_flag_) {
      for (int i = 0; i < n_dim_; ++i) {
        if (u[i] < 0) {
          trapped_site_ = n_sites_ - 1;
        }
      }
    }
    int trapped_2 = (trapped_site_ == 0 ? 1 : n_sites_ - 2);
    const double *const r0 = sites_[trapped_site_].GetPosition();
    const double *const r1 = sites_[trapped_2].GetPosition();
    std::copy(r0, r0 + 3, optical_trap_pos_);
    std::copy(r1, r1 + 3, optical_trap_pos2_);
  }
  poly_ = poly_state::grow;
}

// Place a spool centered at the origin
void Filament::InitSpiral2D() {
  if (n_dim_ > 2)
    Logger::Error("3D Spirals not coded yet.");
  double prev_pos[3] = {0, 0, 0};
  std::fill(position_, position_ + 3, 0);
  sites_[n_sites_ - 1].SetPosition(prev_pos);
  for (auto site = sites_.begin(); site != sites_.end(); ++site) {
    site->SetDiameter(diameter_);
    site->SetLength(bond_length_);
  }
  double step = diameter_ / M_PI;
  double theta = bond_length_ / step;
  for (int i = 2; i < n_sites_ + 1; ++i) {
    double move = step * theta;
    double angle = theta + 2 * M_PI;
    position_[0] = move * cos(angle);
    position_[1] = move * sin(angle);
    theta += bond_length_ / move;
    // Set current site position and orientation
    sites_[n_sites_ - i].SetPosition(position_);
    for (int j = 0; j < 3; ++j) {
      orientation_[j] = (prev_pos[j] - position_[j]) / bond_length_;
      prev_pos[j] = position_[j];
    }
    sites_[n_sites_ - i].SetOrientation(orientation_);
  }
  // Set last site orientation
  sites_[n_sites_ - 1].SetOrientation(sites_[n_sites_ - 2].GetOrientation());
}

void Filament::SetDiffusion() {
  /* Using the friction parameters of Montesi et al. for slender filaments.
     Gives good results for MSD, but not for VCF. */
  // double eps = 1.0 / log(2.0 * length_ / diameter_);
  // double gamma_0 = 4.0 / 3.0 * eps *
  //                 ((1 + 0.64 * eps) / (1 - 1.15 * eps) + 1.659 * SQR(eps));
  // friction_perp_ = bond_length_ * gamma_0;
  // friction_par_ = friction_perp_ / friction_ratio_;
  /* Using the friction coefficients given by Doi & Edwards, 1987. Derived from
     the Smoluchowski eqn by the Kirkwood theory. Gives expected results for MSD
     and VCF. */
  double logLD = log(length_ / diameter_);
  friction_perp_ = 4.0 * length_ / (3.0 * n_sites_ * logLD);
  friction_par_ = friction_perp_ / friction_ratio_;
  /* Using the friction parameters of HSK et al. The MSD agrees with theory but
     the VCF does not. The rotational diffusion of the filament appears to be
     artificially inflated by a factor of 10 or so. The friction ratio here is
     fixed to be 2. */
  // double a = length_ / diameter_ + 1;// + 1;
  // double lna = log(a);
  // friction_par_ =
  // 2.0 / 3.0 * (length_ + 1) / (lna - 0.207 + 0.980 / a - 0.133 / (a * a));
  // friction_perp_ =
  // 4.0 / 3.0 * (length_ + 1) / (lna + 0.839 + 0.185 / a + 0.233 / (a * a));
  // friction_par_ /= n_sites_;
  // friction_perp_ /= n_sites_;
  rand_sigma_perp_ = sqrt(24.0 * friction_perp_ / delta_);
  rand_sigma_par_ = sqrt(24.0 * friction_par_ / delta_);
}

void Filament::GenerateProbableOrientation() {
  /* This updates the current orientation with a generated probable
  orientation where we generate random theta pulled from probability
  distribution P(th) = exp(k cos(th)) where k is the bending stiffness of the
  filament. If k is too large, there is enormous imprecision in this
  calculation since sinh(k) is very large so to fix this I introduce an
  approximate distribution that is valid for large k */
  double theta;
  if (bending_stiffness_ == 0) {
    theta = rng_.RandomUniform() * M_PI;
  } else if (bending_stiffness_ < 100) {
    theta = acos(log(exp(-bending_stiffness_ / bond_length_) +
                     2.0 * rng_.RandomUniform() *
                         sinh(bending_stiffness_ / bond_length_)) /
                 (bending_stiffness_ / bond_length_));
  } else {
    theta = acos((log(2.0 * rng_.RandomUniform()) - log(2.0) +
                  bending_stiffness_ / bond_length_) /
                 (bending_stiffness_ / bond_length_));
  }
  double new_orientation[3] = {0, 0, 0};
  if (n_dim_ == 2) {
    theta = (rng_.RandomInt(2) == 0 ? -1 : 1) * theta;
    new_orientation[0] = cos(theta);
    new_orientation[1] = sin(theta);
  } else {
    double phi = rng_.RandomUniform() * 2.0 * M_PI;
    new_orientation[0] = sin(theta) * cos(phi);
    new_orientation[1] = sin(theta) * sin(phi);
    new_orientation[2] = cos(theta);
  }
  rotate_vector_relative(n_dim_, new_orientation, orientation_);
  std::copy(new_orientation, new_orientation + 3, orientation_);
}

double const Filament::GetVolume() {
  if (n_dim_ == 2) {
    return diameter_ * length_ + 0.25 * M_PI * diameter_ * diameter_;
  } else {
    return 0.25 * M_PI * diameter_ * diameter_ * length_ +
           1.0 / 6.0 * M_PI * diameter_ * diameter_ * diameter_;
  }
}

void Filament::UpdatePosition(bool midstep) {
  midstep_ = midstep;
  ApplyForcesTorques();
  Integrate();
  UpdateAvgPosition();
  DynamicInstability();
}

/*******************************************************************************
  BD algorithm for inextensible wormlike chains with anisotropic friction
  Montesi, Morse, Pasquali. J Chem Phys 122, 084903 (2005).
********************************************************************************/
void Filament::Integrate() {
  CalculateAngles();
  CalculateTangents();
  if (midstep_) {
    ConstructUnprojectedRandomForces();
    GeometricallyProjectRandomForces();
    UpdatePrevPositions();
  }
  AddRandomForces();
  CalculateBendingForces();
  CalculateTensions();
  UpdateSitePositions();
  UpdateBondPositions();
  if (sparams_->reference_frame_flag) {
    /* Rotate/translate filament into COM reference frame coordinates */
    RotateToReferenceFrame();
  }
}

void Filament::CalculateAngles() {
  for (int i_site = 0; i_site < n_sites_ - 2; ++i_site) {
    double const *const u1 = sites_[i_site].GetOrientation();
    double const *const u2 = sites_[i_site + 1].GetOrientation();
    double cos_angle = dot_product(n_dim_, u1, u2);
    cos_thetas_[i_site] = cos_angle;
  }
  if (spiral_init_flag_) {
    CalculateSpiralNumber();
  }
}

void Filament::CalculateTangents() {
  for (auto it = sites_.begin(); it != sites_.end(); ++it) {
    it->CalcTangent();
  }
}

void Filament::CalculateSpiralNumber() {
  // This will calculate the angle
  double u_bond[3] = {0, 0, 0};
  double const *const u_head = bonds_[n_bonds_ - 1].GetOrientation();
  double const *const p_head = sites_[n_sites_ - 1].GetPosition();
  std::copy(u_head, u_head + 3, u_bond);
  for (int i = 0; i < 3; ++i) {
    u_bond[i] = -u_bond[i];
  }
  spiral_number_ = 0;
  for (int i = n_sites_ - 3; i >= 0; --i) {
    // for (int i=n_sites_-2; i>=0; --i) {
    double u[3] = {0, 0, 0};
    double const *const p = sites_[i].GetPosition();
    for (int i = 0; i < n_dim_; ++i) {
      u[i] = p[i] - p_head[i];
    }
    normalize_vector(u, n_dim_);
    double dp = dot_product(n_dim_, u_bond, u);
    /* Note that floating point errors occasionally lead to
       NANs in the acos function, these fix that */
    if (dp > 1)
      dp = 1;
    else if (dp < -1)
      dp = -1;
    double angle = acos(dp);
    double temp[3];
    cross_product(u_bond, u, temp, 3);
    int sign = SIGNOF(temp[2]);
    for (int i = 0; i < n_dim_; ++i) {
      u_bond[i] = u[i];
    }
    spiral_number_ += sign * angle;
  }
}

double Filament::GetSpiralNumber() { return 0.5 * spiral_number_ / M_PI; }

void Filament::GetNematicOrder(double *nematic_order_tensor) {
  for (auto it = bonds_.begin(); it != bonds_.end(); ++it) {
    double const *const u_bond = it->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      for (int j = 0; j < n_dim_; ++j) {
        if (i == j) {
          nematic_order_tensor[3 * i + j] +=
              (2 * u_bond[i] * u_bond[j] - 1) / n_bonds_;
        } else {
          nematic_order_tensor[3 * i + j] +=
              2 * u_bond[i] * u_bond[j] / n_bonds_;
        }
      }
    }
  }
}

void Filament::GetPolarOrder(double *polar_order_vector) {
  double p_i[3] = {0, 0, 0};
  for (auto it = bonds_.begin(); it != bonds_.end(); ++it) {
    double const *const u_bond = it->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      p_i[i] += u_bond[i];
    }
  }
  for (int i = 0; i < n_dim_; ++i) {
    polar_order_vector[i] += p_i[i] / n_bonds_;
  }
}

void Filament::ConstructUnprojectedRandomForces() {
  // Create unprojected forces, see J. Chem. Phys. 122, 084903 (2005),
  // eqn. 40. xi is the random force vector with elements that are
  // uncorrelated and randomly distributed uniformly between -0.5 and 0.5,
  // xi_term is the outer product of the tangent vector u_tan_i u_tan_i acting
  // on the vector xi
  if (zero_temperature_)
    return;
  double xi[3], xi_term[3], f_rand[3];
  for (int i_site = 0; i_site < n_sites_; ++i_site) {
    double const *const utan = sites_[i_site].GetTangent();
    for (int i = 0; i < n_dim_; ++i)
      xi[i] = rng_.RandomUniform() - 0.5;
    if (n_dim_ == 2) {
      xi_term[0] = SQR(utan[0]) * xi[0] + utan[0] * utan[1] * xi[1];
      xi_term[1] = SQR(utan[1]) * xi[1] + utan[0] * utan[1] * xi[0];
    } else if (n_dim_ == 3) {
      xi_term[0] = SQR(utan[0]) * xi[0] + utan[0] * utan[1] * xi[1] +
                   utan[0] * utan[2] * xi[2];
      xi_term[1] = SQR(utan[1]) * xi[1] + utan[0] * utan[1] * xi[0] +
                   utan[1] * utan[2] * xi[2];
      xi_term[2] = SQR(utan[2]) * xi[2] + utan[0] * utan[2] * xi[0] +
                   utan[1] * utan[2] * xi[1];
    }
    for (int i = 0; i < n_dim_; ++i) {
      f_rand[i] = rand_sigma_perp_ * xi[i] +
                  (rand_sigma_par_ - rand_sigma_perp_) * xi_term[i];
    }
    sites_[i_site].SetRandomForce(f_rand);
  }
}

void Filament::GeometricallyProjectRandomForces() {
  if (zero_temperature_)
    return;
  double f_rand_temp[3];
  for (int i_site = 0; i_site < n_sites_ - 1; ++i_site) {
    // Use the tensions vector to calculate the hard components of the random
    // forces These are not the same as the tensions, they will be calculated
    // later
    double const *const f_rand1 = sites_[i_site].GetRandomForce();
    double const *const f_rand2 = sites_[i_site + 1].GetRandomForce();
    double const *const u_site = sites_[i_site].GetOrientation();
    for (int i = 0; i < n_dim_; ++i)
      f_rand_temp[i] = f_rand2[i] - f_rand1[i];
    tensions_[i_site] = dot_product(n_dim_, f_rand_temp, u_site);
    // Then get the G arrays (for inertialess case where m=1, see
    // ref. 15 of above paper)
    g_mat_diag_[i_site] = 2;
    if (i_site > 0) {
      g_mat_upper_[i_site - 1] = -cos_thetas_[i_site - 1];
      g_mat_lower_[i_site - 1] = -cos_thetas_[i_site - 1];
    }
  }
  // Now solve using tridiagonal solver
  tridiagonal_solver(&g_mat_lower_, &g_mat_diag_, &g_mat_upper_, &tensions_,
                     n_sites_ - 1);
  // Update to the projected brownian forces
  // First the end sites:
  double f_proj[3];
  for (int i = 0; i < n_dim_; ++i) {
    f_proj[i] = sites_[0].GetRandomForce()[i] +
                tensions_[0] * sites_[0].GetOrientation()[i];
  }
  sites_[0].SetRandomForce(f_proj);
  for (int i = 0; i < n_dim_; ++i) {
    f_proj[i] =
        sites_[n_sites_ - 1].GetRandomForce()[i] -
        tensions_[n_sites_ - 2] * sites_[n_sites_ - 2].GetOrientation()[i];
  }
  sites_[n_sites_ - 1].SetRandomForce(f_proj);
  // Then the rest
  for (int i_site = 1; i_site < n_sites_ - 1; ++i_site) {
    double const *const u1 = sites_[i_site - 1].GetOrientation();
    double const *const u2 = sites_[i_site].GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      f_proj[i] = sites_[i_site].GetRandomForce()[i] +
                  tensions_[i_site] * u2[i] - tensions_[i_site - 1] * u1[i];
    }
    sites_[i_site].SetRandomForce(f_proj);
  }
}

void Filament::AddRandomForces() {
  if (zero_temperature_)
    return;
  for (auto site = sites_.begin(); site != sites_.end(); ++site)
    site->AddRandomForce();
}

void Filament::CalculateBendingForces() {
  /* Metric forces give the appropriate equilibrium behavior at zero
   * persistence
   * length: all angles have an equal probability of being sampled */
  { // Apply metric forces
    det_t_mat_[0] = 1;
    det_t_mat_[1] = 2;
    det_b_mat_[n_sites_] = 1;
    det_b_mat_[n_sites_ - 1] = 2;
    for (int i = 2; i < n_sites_; ++i) {
      det_t_mat_[i] =
          2 * det_t_mat_[i - 1] - SQR(-cos_thetas_[i - 2]) * det_t_mat_[i - 2];
      det_b_mat_[n_sites_ - i] =
          2 * det_b_mat_[n_sites_ - i + 1] -
          SQR(-cos_thetas_[n_sites_ - i - 1]) * det_b_mat_[n_sites_ - i + 2];
    }
    double det_g = det_t_mat_[n_sites_ - 1];
    for (int i = 0; i < n_sites_ - 2; ++i) {
      g_mat_inverse_[i] =
          cos_thetas_[i] * det_t_mat_[i] * det_b_mat_[i + 3] / det_g;
    }
  }
  // In the case of no metric forces
  // for (int i = 0; i < n_sites_ - 2; ++i) {
  // g_mat_inverse_[i] = 0;
  //}

  // Now calculate the effective rigidities
  for (int i = 0; i < n_sites_ - 2; ++i) {
    k_eff_[i] = (bending_stiffness_ + bond_length_ * g_mat_inverse_[i]) /
                SQR(bond_length_);
  }
  /* The following algorithm calculates the bending forces on each of
   * the sites.
   *
   * These calculations were done by hand in order to be expressed in
   * a form that is efficient to compute, but is not human readable.
   * The correct equilibrium behavior has been verified for low
   * persistence lengths by sampling of bond angle probabilities and
   * at high persistence lengths by looking at the filament
   * mean-square end-to-end distances (before the implicit curvature
   * code was added).
   *
   * WARNING: Be very careful when changing the following code. If
   * this code ever breaks, you need to either check through the
   * indices very carefully or completely redo the calculation by
   * hand! See Pasquali and Morse, J. Chem. Phys. Vol 116, No 5
   * (2002) */
  double f_site[3] = {0, 0, 0};
  if (n_dim_ == 2) {
    double zvec[3] = {0, 0, 1};
    double u1[3] = {0, 0, 0};
    double u2[3] = {0, 0, 0};
    double curve_mag = 0.5 * flagella_amplitude_ * M_PI / length_;
    double theta_t = 0;
    if (flagella_freq_ != 0) {
      theta_t = 2 * M_PI * n_step_ * delta_ / flagella_freq_;
    }
    double theta_x = 2 * M_PI * flagella_period_ * bond_length_ / length_;
    double curve = curvature_;
    for (int k_site = 0; k_site < n_sites_; ++k_site) {
      std::fill(f_site, f_site + 3, 0.0);
      if (k_site > 1) {
        if (flagella_flag_) {
          curve = curve_mag * sin(theta_x * (k_site - 1) - theta_t);
        }
        double const *const u1_temp = sites_[k_site - 2].GetOrientation();
        double const *const u2_temp = sites_[k_site - 1].GetOrientation();
        std::copy(u1_temp, u1_temp + 3, u1);
        std::copy(u2_temp, u2_temp + 3, u2);
        if (curvature_ != 0 || flagella_flag_) {
          rotate_vector(u1, zvec, curve * bond_length_, n_dim_);
          rotate_vector(u2, zvec, -curve * bond_length_, n_dim_);
        }
        f_site[0] += k_eff_[k_site - 2] *
                     ((1 - SQR(u2[0])) * u1[0] - u2[0] * u2[1] * u1[1]);
        f_site[1] += k_eff_[k_site - 2] *
                     ((1 - SQR(u2[1])) * u1[1] - u2[0] * u2[1] * u1[0]);
      }
      if (k_site > 0 && k_site < n_sites_ - 1) {
        if (flagella_flag_) {
          curve = curve_mag * sin(theta_x * (k_site)-theta_t);
        }
        double const *const u1_temp = sites_[k_site - 1].GetOrientation();
        double const *const u2_temp = sites_[k_site].GetOrientation();
        std::copy(u1_temp, u1_temp + 3, u1);
        std::copy(u2_temp, u2_temp + 3, u2);
        if (curvature_ != 0 || flagella_flag_) {
          rotate_vector(u1, zvec, curve * bond_length_, n_dim_);
          rotate_vector(u2, zvec, -curve * bond_length_, n_dim_);
        }
        f_site[0] += k_eff_[k_site - 1] *
                     ((1 - SQR(u1[0])) * u2[0] - u1[0] * u1[1] * u2[1] -
                      ((1 - SQR(u2[0])) * u1[0] - u2[0] * u2[1] * u1[1]));
        f_site[1] += k_eff_[k_site - 1] *
                     ((1 - SQR(u1[1])) * u2[1] - u1[0] * u1[1] * u2[0] -
                      ((1 - SQR(u2[1])) * u1[1] - u2[0] * u2[1] * u1[0]));
      }
      if (k_site < n_sites_ - 2) {
        if (flagella_flag_) {
          curve = curve_mag * sin(theta_x * (k_site + 1) - theta_t);
        }
        double const *const u1_temp = sites_[k_site].GetOrientation();
        double const *const u2_temp = sites_[k_site + 1].GetOrientation();
        std::copy(u1_temp, u1_temp + 3, u1);
        std::copy(u2_temp, u2_temp + 3, u2);
        if (curvature_ != 0 || flagella_flag_) {
          rotate_vector(u1, zvec, curve * bond_length_, n_dim_);
          rotate_vector(u2, zvec, -curve * bond_length_, n_dim_);
        }
        f_site[0] -=
            k_eff_[k_site] * ((1 - SQR(u1[0])) * u2[0] - u1[0] * u1[1] * u2[1]);
        f_site[1] -=
            k_eff_[k_site] * ((1 - SQR(u1[1])) * u2[1] - u1[0] * u1[1] * u2[0]);
      }
      sites_[k_site].AddForce(f_site);
    }
  } else if (n_dim_ == 3) {
    for (int k_site = 0; k_site < n_sites_; ++k_site) {
      std::fill(f_site, f_site + 3, 0.0);
      if (k_site > 1) {
        double const *const u1 = sites_[k_site - 2].GetOrientation();
        double const *const u2 = sites_[k_site - 1].GetOrientation();
        f_site[0] += k_eff_[k_site - 2] *
                     ((1 - SQR(u2[0])) * u1[0] - u2[0] * u2[1] * u1[1] -
                      u2[0] * u2[2] * u1[2]);
        f_site[1] += k_eff_[k_site - 2] *
                     ((1 - SQR(u2[1])) * u1[1] - u2[1] * u2[0] * u1[0] -
                      u2[1] * u2[2] * u1[2]);
        f_site[2] += k_eff_[k_site - 2] *
                     ((1 - SQR(u2[2])) * u1[2] - u2[2] * u2[0] * u1[0] -
                      u2[2] * u2[1] * u1[1]);
      }
      if (k_site > 0 && k_site < n_sites_ - 1) {
        double const *const u1 = sites_[k_site - 1].GetOrientation();
        double const *const u2 = sites_[k_site].GetOrientation();
        f_site[0] += k_eff_[k_site - 1] *
                     ((1 - SQR(u1[0])) * u2[0] - u1[0] * u1[1] * u2[1] -
                      u1[0] * u1[2] * u2[2] -
                      ((1 - SQR(u2[0])) * u1[0] - u2[0] * u2[1] * u1[1] -
                       u2[0] * u2[2] * u1[2]));
        f_site[1] += k_eff_[k_site - 1] *
                     ((1 - SQR(u1[1])) * u2[1] - u1[1] * u1[0] * u2[0] -
                      u1[1] * u1[2] * u2[2] -
                      ((1 - SQR(u2[1])) * u1[1] - u2[1] * u2[0] * u1[0] -
                       u2[1] * u2[2] * u1[2]));
        f_site[2] += k_eff_[k_site - 1] *
                     ((1 - SQR(u1[2])) * u2[2] - u1[2] * u1[0] * u2[0] -
                      u1[2] * u1[1] * u2[1] -
                      ((1 - SQR(u2[2])) * u1[2] - u2[2] * u2[0] * u1[0] -
                       u2[1] * u2[2] * u1[1]));
      }
      if (k_site < n_sites_ - 2) {
        double const *const u1 = sites_[k_site].GetOrientation();
        double const *const u2 = sites_[k_site + 1].GetOrientation();
        f_site[0] -=
            k_eff_[k_site] * ((1 - SQR(u1[0])) * u2[0] - u1[0] * u1[1] * u2[1] -
                              u1[0] * u1[2] * u2[2]);
        f_site[1] -=
            k_eff_[k_site] * ((1 - SQR(u1[1])) * u2[1] - u1[1] * u1[0] * u2[0] -
                              u1[1] * u1[2] * u2[2]);
        f_site[2] -=
            k_eff_[k_site] * ((1 - SQR(u1[2])) * u2[2] - u1[2] * u1[0] * u2[0] -
                              u1[2] * u1[1] * u2[1]);
      }
      sites_[k_site].AddForce(f_site);
    }
  }
}

void Filament::CalculateTensions() {
  // Calculate friction_inverse matrix
  int site_index = 0;
  int next_site = n_dim_ * n_dim_;
  for (int i_site = 0; i_site < n_sites_; ++i_site) {
    int gamma_index = 0;
    double const *const utan = sites_[i_site].GetTangent();
    for (int i = 0; i < n_dim_; ++i) {
      for (int j = 0; j < n_dim_; ++j) {
        gamma_inverse_[site_index + gamma_index] =
            1.0 / friction_par_ * (utan[i] * utan[j]) +
            1.0 / friction_perp_ * ((i == j ? 1 : 0) - utan[i] * utan[j]);
        gamma_index++;
      }
    }
    site_index += next_site;
  }
  // Populate the H matrix and Q vector using tensions (p_vec) array
  double temp_a, temp_b;
  double f_diff[3];
  double utan1_dot_u2, utan2_dot_u2;
  site_index = 0;
  for (int i_site = 0; i_site < n_sites_ - 1; ++i_site) {
    // f_diff is the term in par_entheses in equation 29 of J. Chem. Phys.
    // 122, 084903 (2005)
    double const *const f1 = sites_[i_site].GetForce();
    double const *const f2 = sites_[i_site + 1].GetForce();
    double const *const u2 = sites_[i_site].GetOrientation();
    double const *const utan1 = sites_[i_site].GetTangent();
    double const *const utan2 = sites_[i_site + 1].GetTangent();
    for (int i = 0; i < n_dim_; ++i) {
      temp_a = gamma_inverse_[site_index + n_dim_ * i] * f1[0] +
               gamma_inverse_[site_index + n_dim_ * i + 1] * f1[1];
      if (n_dim_ == 3)
        temp_a += gamma_inverse_[site_index + n_dim_ * i + 2] * f1[2];
      temp_b = gamma_inverse_[site_index + next_site + n_dim_ * i] * f2[0] +
               gamma_inverse_[site_index + next_site + n_dim_ * i + 1] * f2[1];
      if (n_dim_ == 3)
        temp_b +=
            gamma_inverse_[site_index + next_site + n_dim_ * i + 2] * f2[2];
      f_diff[i] = temp_b - temp_a;
    }
    tensions_[i_site] = dot_product(n_dim_, u2, f_diff);
    utan1_dot_u2 = dot_product(n_dim_, utan1, u2);
    utan2_dot_u2 = dot_product(n_dim_, utan2, u2);
    h_mat_diag_[i_site] =
        2.0 / friction_perp_ + (1.0 / friction_par_ - 1.0 / friction_perp_) *
                                   (SQR(utan1_dot_u2) + SQR(utan2_dot_u2));
    if (i_site > 0) {
      double const *const u1 = sites_[i_site - 1].GetOrientation();
      h_mat_upper_[i_site - 1] =
          -1.0 / friction_perp_ * dot_product(n_dim_, u2, u1) -
          (1.0 / friction_par_ - 1.0 / friction_perp_) *
              (dot_product(n_dim_, utan1, u1) * dot_product(n_dim_, utan1, u2));
      h_mat_lower_[i_site - 1] = h_mat_upper_[i_site - 1];
    }
    site_index += next_site;
  }
  tridiagonal_solver(&h_mat_lower_, &h_mat_diag_, &h_mat_upper_, &tensions_,
                     n_sites_ - 1);
}

void Filament::UpdateSitePositions() {
  double delta = (midstep_ ? 0.5 * delta_ : delta_);
  double f_site[3];
  // First get total forces
  // Handle end sites first
  for (int i = 0; i < n_dim_; ++i)
    f_site[i] = tensions_[0] * sites_[0].GetOrientation()[i];
  sites_[0].AddForce(f_site);
  for (int i = 0; i < n_dim_; ++i)
    f_site[i] =
        -tensions_[n_sites_ - 2] * sites_[n_sites_ - 2].GetOrientation()[i];
  sites_[n_sites_ - 1].AddForce(f_site);
  // and then the rest
  for (int i_site = 1; i_site < n_sites_ - 1; ++i_site) {
    double const *const u_site1 = sites_[i_site - 1].GetOrientation();
    double const *const u_site2 = sites_[i_site].GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      f_site[i] =
          tensions_[i_site] * u_site2[i] - tensions_[i_site - 1] * u_site1[i];
    }
    sites_[i_site].AddForce(f_site);
  }
  // Now update positions
  double f_term[3], r_new[3];
  int site_index = 0;
  int next_site = n_dim_ * n_dim_;
  for (int i_site = 0; i_site < n_sites_; ++i_site) {
    double const *const f_site1 = sites_[i_site].GetForce();
    double const *const r_prev = sites_[i_site].GetPrevPosition();
    for (int i = 0; i < n_dim_; ++i) {
      f_term[i] = gamma_inverse_[site_index + n_dim_ * i] * f_site1[0] +
                  gamma_inverse_[site_index + n_dim_ * i + 1] * f_site1[1];
      if (n_dim_ == 3)
        f_term[i] += gamma_inverse_[site_index + n_dim_ * i + 2] * f_site1[2];
      r_new[i] = r_prev[i] + delta * f_term[i];
    }
    sites_[i_site].SetPosition(r_new);
    site_index += next_site;
  }
  // Next, update orientation vectors
  double u_mag, r_diff[3];
  for (int i_site = 0; i_site < n_sites_ - 1; ++i_site) {
    double const *const r_site1 = sites_[i_site].GetPosition();
    double const *const r_site2 = sites_[i_site + 1].GetPosition();
    u_mag = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_diff[i] = r_site2[i] - r_site1[i];
      u_mag += SQR(r_diff[i]);
    }
    u_mag = sqrt(u_mag);
    for (int i = 0; i < n_dim_; ++i)
      r_diff[i] /= u_mag;
    sites_[i_site].SetOrientation(r_diff);
  }
  sites_[n_sites_ - 1].SetOrientation(sites_[n_sites_ - 2].GetOrientation());
  // Finally, normalize site positions, making sure the sites are still
  // rod-length apart
  if (CheckBondLengths()) {
    Logger::Debug(
        "Renormalizing bond lengths: %d steps since last renormalization",
        n_normalize_);
    if (error_analysis_) {
      error_rates_.push_back(n_normalize_);
    }
    normalize_switch_ = !normalize_switch_;
    n_normalize_ = 0;
    if (normalize_switch_) {
      // Normalize from tail to head
      for (int i_site = 1; i_site < n_sites_; ++i_site) {
        double const *const r_site1 = sites_[i_site - 1].GetPosition();
        double const *const u_site1 = sites_[i_site - 1].GetOrientation();
        for (int i = 0; i < n_dim_; ++i)
          r_diff[i] = r_site1[i] + bond_length_ * u_site1[i];
        sites_[i_site].SetPosition(r_diff);
      }
    } else {
      // Normalize from head to tail
      for (int i_site = n_sites_ - 1; i_site > 0; --i_site) {
        double const *const r_site1 = sites_[i_site].GetPosition();
        double const *const u_site1 = sites_[i_site - 1].GetOrientation();
        for (int i = 0; i < n_dim_; ++i)
          r_diff[i] = r_site1[i] - bond_length_ * u_site1[i];
        sites_[i_site - 1].SetPosition(r_diff);
      }
    }
  }
  n_normalize_++;
}

void Filament::GetErrorRates(std::vector<int> &rates) {
  rates.insert(rates.end(), error_rates_.begin(), error_rates_.end());
  error_rates_.clear();
}

bool Filament::CheckBondLengths() {
  bool renormalize = false;
  for (int i_site = 1; i_site < n_sites_; ++i_site) {
    double const *const r_site1 = sites_[i_site - 1].GetPosition();
    double const *const r_site2 = sites_[i_site].GetPosition();
    double a = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      double temp = r_site2[i] - r_site1[i];
      a += temp * temp;
    }
    a = sqrt(a);
    double err = ABS(bond_length_ - a) / bond_length_;
    if (err > 1e-3) {
      renormalize = true;
    }
  }
  return renormalize;
}

void Filament::UpdateAvgPosition() {
  std::fill(position_, position_ + 3, 0.0);
  std::fill(orientation_, orientation_ + 3, 0.0);
  for (auto site_it = sites_.begin(); site_it != sites_.end(); ++site_it) {
    double const *const site_pos = site_it->GetPosition();
    double const *const site_u = site_it->GetOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      position_[i] += site_pos[i];
      orientation_[i] += site_u[i];
    }
  }
  normalize_vector(orientation_, n_dim_);
  for (int i = 0; i < n_dim_; ++i) {
    position_[i] /= n_sites_;
  }
}

void Filament::ApplyForcesTorques() {
  ApplyInteractionForces();
  if (optical_trap_flag_) {
    double f_trap1[3] = {0};
    double f_trap2[3] = {0};
    double const *const r0 = sites_[trapped_site_].GetPosition();
    int trap2 = (trapped_site_ == 0 ? 1 : n_sites_ - 2);
    double const *const r1 = sites_[trap2].GetPosition();
    for (int i = 0; i < n_dim_; ++i) {
      f_trap1[i] = optical_trap_spring_ * (optical_trap_pos_[i] - r0[i]);
      if (optical_trap_fixed_) {
        f_trap2[i] = optical_trap_spring_ * (optical_trap_pos2_[i] - r1[i]);
      }
    }
    sites_[trapped_site_].AddForce(f_trap1);
    if (optical_trap_fixed_) {
      sites_[trap2].AddForce(f_trap2);
    }
  }
  // if (anchored_) ApplyAnchorForces();
}

void Filament::ApplyInteractionForces() {
  double pure_torque[3] = {0, 0, 0};
  double site_force[3] = {0, 0, 0};
  double linv = 1.0 / bond_length_;
  if (!sparams_->drive_from_bond_center) {
    // Driving originating from the site tangents
    CalculateTangents();
  }
  if (nematic_driving_ && rng_.RandomUniform() < p_driving_switch_) {
    driving_factor_ = -driving_factor_;
  }
  for (int i = 0; i < n_bonds_; ++i) {
    double const *const f = bonds_[i].GetForce();
    double const *const t = bonds_[i].GetTorque();
    double const *const u = sites_[i].GetOrientation();
    if (i == n_bonds_ - 1) {
      tip_force_ = -dot_product(n_dim_, u, f);
    }
    AddPotential(bonds_[i].GetPotentialEnergy());
    // Convert torques into forces at bond ends
    // u x t / bond_length = pure torque force at tail of bond
    cross_product(u, t, pure_torque, 3);
    for (int i = 0; i < n_dim_; ++i) {
      pure_torque[i] *= linv;
      site_force[i] = 0.5 * f[i];
    }
    // Add translational forces and pure torque forces at bond ends
    sites_[i].AddForce(site_force);
    sites_[i].AddForce(pure_torque);
    for (int j = 0; j < n_dim_; ++j)
      pure_torque[j] *= -1;
    sites_[i + 1].AddForce(site_force);
    sites_[i + 1].AddForce(pure_torque);
    // The driving factor is a force per unit length,
    // so need to multiply by bond length to get f_dr on bond
    if (params_->i_step > eq_steps_) {
      double f_dr[3] = {};
      double mag = 0.5 * driving_factor_ * bond_length_;
      if (sparams_->drive_from_bond_center) {
        // Add driving (originating from the com of the bond)
        // Driving toward filament center (for fun)
        // if (i >= n_bonds_/2) {
        // mag *= -1;
        //} else if (n_bonds_%2!=0 && i == floor(n_bonds_/2)) {
        // mag = 0;
        //}
        for (int j = 0; j < n_dim_; ++j)
          f_dr[j] = mag * u[j];
        sites_[i].AddForce(f_dr);
        sites_[i + 1].AddForce(f_dr);
      } else {
        // Driving from sites
        double const *const u_tan1 = sites_[i].GetTangent();
        double const *const u_tan2 = sites_[i + 1].GetTangent();
        for (int j = 0; j < n_dim_; ++j) {
          f_dr[j] = mag * u_tan1[j];
        }
        sites_[i].AddForce(f_dr);
        for (int j = 0; j < n_dim_; ++j) {
          f_dr[j] = mag * u_tan2[j];
        }
        sites_[i + 1].AddForce(f_dr);
      }
    }
  }
}

void Filament::DynamicInstability() {
  if (midstep_ || !dynamic_instability_flag_)
    return;
  UpdatePolyState();
  GrowFilament();
  SetDiffusion();
}

void Filament::GrowFilament() {
  // If the filament is paused, do nothing
  if (poly_ == +poly_state::pause)
    return;
  // Otherwise, adjust filament length due to polymerization
  double delta_length = 0;
  if (poly_ == +poly_state::grow) {
    delta_length = v_poly_ * delta_;
  } else if (poly_ == +poly_state::shrink) {
    if (n_bonds_ == 2 && bond_length_ <= min_bond_length_) {
      return;
    }
    delta_length = -v_depoly_ * delta_;
  }
  length_ += delta_length;
  RescaleBonds();
  if (bond_length_ > max_bond_length_) {
    DoubleGranularityLinear();
    // RebindMotors();
  } else if (bond_length_ < min_bond_length_ && n_bonds_ > 2) {
    HalfGranularityLinear();
    // RebindMotors();
  }
}

void Filament::RescaleBonds() {
  UpdatePrevPositions();
  double old_bond_length = bond_length_;
  bond_length_ = length_ / n_bonds_;
  double k, dl;
  double r2[3] = {0, 0, 0};
  if (poly_ == +poly_state::shrink) {
    // Old code
    double const *const r0 = GetTailPosition();
    double const *const u0 = GetTailOrientation();
    for (int i = 0; i < n_dim_; ++i) {
      r2[i] = r0[i] + u0[i] * bond_length_;
    }
    sites_[1].SetPosition(r2);
    dl = old_bond_length - bond_length_;
    for (int i_site = 2; i_site < n_sites_; ++i_site) {
      double const *const r1 = sites_[i_site - 1].GetPrevPosition();
      double const *const u1 = sites_[i_site - 1].GetOrientation();
      k = SQR(dl * cos_thetas_[i_site - 2]) - SQR(dl) + SQR(bond_length_);
      k = (k > 0 ? k : 0);
      k = -cos_thetas_[i_site - 2] * dl + sqrt(k);
      for (int i = 0; i < n_dim_; ++i)
        r2[i] = r1[i] + k * u1[i];
      sites_[i_site].SetPosition(r2);
      dl = old_bond_length - k;
    }
  } else if (poly_ == +poly_state::grow) {
    dl = old_bond_length;
    for (int i_site = 1; i_site < n_sites_ - 1; ++i_site) {
      double const *const r_old = sites_[i_site].GetPrevPosition();
      double const *const u = sites_[i_site].GetOrientation();
      k = SQR(dl * cos_thetas_[i_site - 1]) - SQR(dl) + SQR(bond_length_);
      k = (k > 0 ? k : 0);
      k = -cos_thetas_[i_site - 1] * dl + sqrt(k);
      for (int i = 0; i < n_dim_; ++i) {
        r2[i] = r_old[i] + k * u[i];
      }
      sites_[i_site].SetPosition(r2);
      dl = old_bond_length - k;
    }
    double const *const u = sites_[n_sites_ - 2].GetOrientation();
    double const *const r_old = sites_[n_sites_ - 2].GetPosition();
    for (int i = 0; i < n_dim_; ++i) {
      r2[i] = r_old[i] + bond_length_ * u[i];
    }
    sites_[n_sites_ - 1].SetPosition(r2);
  }
  UpdateBondPositions();
  CalculateAngles();
}

void Filament::UpdatePolyState() {
  double p_g2s = p_g2s_;
  double p_p2s = p_p2s_;
  double roll = rng_.RandomUniform();
  double p_norm;
  // Modify catastrophe probabilities if the end of the filament is under a
  // load
  if (force_induced_catastrophe_flag_ && tip_force_ > 0) {
    double p_factor = exp(fic_factor_ * tip_force_);
    p_g2s = (p_g2s + p_g2p_) * p_factor;
    p_p2s = p_p2s * p_factor;
  }
  // Filament shrinking
  if (poly_ == +poly_state::shrink) {
    p_norm = p_s2g_ + p_s2p_;
    if (p_norm > 1.0)
      poly_ = (roll < p_s2g_ / p_norm ? poly_state::grow : poly_state::pause);
    else {
      if (roll < p_s2g_)
        poly_ = poly_state::grow;
      else if (roll < (p_s2g_ + p_s2p_))
        poly_ = poly_state::pause;
    }
  }
  // Filament growing
  else if (poly_ == +poly_state::grow) {
    p_norm = p_g2s + p_g2p_;
    if (p_norm > 1.0)
      poly_ = (roll < p_g2s / p_norm ? poly_state::shrink : poly_state::pause);
    else {
      if (roll < p_g2s)
        poly_ = poly_state::shrink;
      else if (roll < (p_g2s + p_g2p_))
        poly_ = poly_state::pause;
    }
  }
  // Filament paused
  else if (poly_ == +poly_state::pause) {
    p_norm = p_p2g_ + p_p2s;
    if (p_norm > 1)
      poly_ = (roll < p_p2g_ / p_norm ? poly_state::grow : poly_state::shrink);
    else {
      if (roll < p_p2g_)
        poly_ = poly_state::grow;
      else if (roll < (p_p2g_ + p_p2s))
        poly_ = poly_state::shrink;
    }
  }
  // Check to make sure the filament lengths stay in the correct ranges
  if (length_ < min_length_)
    poly_ = poly_state::grow;
  else if (length_ > max_length_)
    poly_ = poly_state::shrink;
}

void Filament::CheckFlocking() {
  if (!sparams_->polar_order_analysis) {
    /* Only need to call if we aren't calculating polar order in
       PolarOrderAnalysis*/
    CalcPolarOrder();
  }
  double avg_polar_order = 0;
  double avg_contact_number = 0;
  for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
    /* Check if filament satisfies the flock condition by checking whether
       the average bond polar order is above cutoff. Also see whether the
       average contact number exceeds the cutoff for interior/exterior
       filament. */
    avg_polar_order += bond->GetPolarOrder();
    avg_contact_number += bond->GetContactNumber();
  }
  avg_polar_order /= n_bonds_;
  avg_contact_number /= n_bonds_;

  int in_flock_prev = in_flock_;
  in_flock_ = 0;
  flock_change_state_ = 0;
  if (avg_polar_order >= sparams_->flock_polar_min) {
    // Filament is in a flock
    if (in_flock_prev == 0) {
      // Filament joined flock this timestep
      flock_change_state_ = 1;
    }
    if (avg_contact_number >= sparams_->flock_contact_min) {
      // Filament is in flock interior
      in_flock_ = 2;
    } else {
      // Filament is in flock exterior
      in_flock_ = 1;
    }
  } else if (in_flock_prev > 0) {
    // Filament left flock this timestep
    flock_change_state_ = 2;
  }
}

void Filament::Draw(std::vector<graph_struct *> &graph_array) {
  if (sparams_->curvature_cluster_analysis && cluster_ > 0) {
    double color = sparams_->color + cluster_ * 0.1 * M_PI * M_PI;
    for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
      bond->SetColor(color, draw_type::fixed);
    }
  } else if (sparams_->highlight_curvature) {
    draw_ = draw_type::fixed;
    color_ = sparams_->color + M_PI;
    double avg_cos_theta = 0;
    // Average over bond angles
    for (int i = 0; i < n_bonds_ - 1; ++i) {
      avg_cos_theta += cos_thetas_[i];
    }
    avg_cos_theta /= n_bonds_ - 1;
    double scale = 4 * M_PI;
    double nominal_cos_theta = cos(2 * curvature_ * bond_length_);
    double color_diff = scale * (avg_cos_theta - nominal_cos_theta);
    for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
      bond->SetColor(color_ + color_diff, draw_);
    }
  } else if (sparams_->highlight_flock) {
    if (!sparams_->flocking_analysis) {
      // We already call this during FlockingAnalysis
      CheckFlocking();
    }
    if (in_flock_ == 1) {
      // Part of flock exterior
      for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
        bond->SetColor(sparams_->flock_color_ext, draw_type::fixed);
      }
    } else if (in_flock_ == 2) {
      // Part of flock interior
      for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
        bond->SetColor(sparams_->flock_color_int, draw_type::fixed);
      }
    } else {
      for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
        bond->SetColor(color_, draw_);
      }
    }
  } else {
    for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
      bond->SetColor(color_, draw_);
    }
  }
  for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
    bond->Draw(graph_array);
  }
  if (sparams_->draw_center_of_curvature) {
    double pos[3] = {0, 0, 0};
    double roc = GetCenterOfCurvature(pos);
    std::copy(pos, pos + 3, g_.r);
    std::copy(orientation_, orientation_ + 3, g_.u);
    g_.color = sparams_->color;
    g_.diameter = 0.5;
    g_.length = 0;
    g_.draw = draw_type::bw;
    graph_array.push_back(&g_);
  }
}

// Scale bond and site positions from new unit cell
void Filament::ScalePosition() {
  // scale first bond position using new unit cell
  bonds_[0].ScalePosition();
  // then reposition sites based on first bond position
  // handle first site
  double r[3];
  double const *const bond_r = bonds_[0].GetPosition();
  double const *const bond_u = bonds_[0].GetOrientation();
  for (int i = 0; i < n_dim_; ++i)
    r[i] = bond_r[i] - 0.5 * bond_length_ * bond_u[i];
  sites_[0].SetPosition(r);
  // then handle remaining sites
  for (int i_site = 1; i_site < n_sites_; ++i_site) {
    double const *const prev_r = sites_[i_site - 1].GetPosition();
    double const *const prev_u = sites_[i_site - 1].GetOrientation();
    for (int i = 0; i < n_dim_; ++i)
      r[i] = prev_r[i] + bond_length_ * prev_u[i];
    sites_[i_site].SetPosition(r);
  }
  // update remaining bond positions
  UpdateBondPositions();
}

/* Returns average radius of curvature and position of center of curvature */
const double Filament::GetCenterOfCurvature(double *center) {
  std::fill(center, center + 3, 0.0);
  double avg_theta = 0;
  for (int i = 0; i < n_bonds_ - 1; ++i) {
    avg_theta += acos(cos_thetas_[i]);
  }
  avg_theta /= (n_bonds_ - 1);
  double avg_roc = bond_length_ / avg_theta;
  if (avg_roc > params_->system_radius || avg_roc != avg_roc) {
    avg_roc = 0;
  }
  // avg_roc = sparams_->radius_of_curvature;
  double z_vec[3] = {0, 0, 1};
  double mid_pos[3] = {0, 0, 0};
  double mid_tan[3] = {0, 0, 0};
  /* If odd number of sites, use middle site for center of filament, otherwise,
     use center of middle bond */
  if (n_sites_ % 2 != 0) {
    int mid_site = (int)floor(n_sites_ / 2);
    const double *const pos = sites_[mid_site].GetPosition();
    sites_[mid_site].CalcTangent();
    const double *const tan = sites_[mid_site].GetTangent();
    std::copy(tan, tan + 3, mid_tan);
    std::copy(pos, pos + 3, mid_pos);
  } else {
    int mid_bond = (int)floor(n_bonds_ / 2);
    const double *const pos = bonds_[mid_bond].GetPosition();
    const double *const tan = bonds_[mid_bond].GetOrientation();
    std::copy(tan, tan + 3, mid_tan);
    std::copy(pos, pos + 3, mid_pos);
  }
  /* Rotate vector 90 degrees toward the direction of curvature. Use the spiral
     number to determine handedness of curvature */
  CalculateSpiralNumber();
  rotate_vector(mid_tan, z_vec, -0.5 * SIGNOF(GetSpiralNumber()) * M_PI,
                n_dim_);
  for (int i = 0; i < n_dim_; ++i) {
    center[i] = mid_pos[i] + avg_roc * mid_tan[i];
  }
  std::copy(center, center + 3, position_);
  UpdatePeriodic();
  std::copy(scaled_position_, scaled_position_ + 3, center);
  return avg_roc;
}

void Filament::ReportAll() {
  printf("tensions:\n  {");
  for (int i = 0; i < n_sites_ - 1; ++i)
    printf(" %5.5f ", tensions_[i]);
  printf("}\n");
  printf("cos_thetas:\n  {");
  for (int i = 0; i < n_sites_ - 2; ++i)
    printf(" %5.5f ", cos_thetas_[i]);
  printf("}\n");
  printf("g_mat_lower:\n  {");
  for (int i = 0; i < n_sites_ - 2; ++i)
    printf(" %5.5f ", g_mat_lower_[i]);
  printf("}\n");
  printf("g_mat_upper:\n  {");
  for (int i = 0; i < n_sites_ - 2; ++i)
    printf(" %5.5f ", g_mat_upper_[i]);
  printf("}\n");
  printf("g_mat_diag:\n  {");
  for (int i = 0; i < n_sites_ - 1; ++i)
    printf(" %5.5f ", g_mat_diag_[i]);
  printf("}\n");
  printf("det_t_mat:\n  {");
  for (int i = 0; i < n_sites_ + 1; ++i)
    printf(" %5.5f ", det_t_mat_[i]);
  printf("}\n");
  printf("det_b_mat:\n  {");
  for (int i = 0; i < n_sites_ + 1; ++i)
    printf(" %5.5f ", det_b_mat_[i]);
  printf("}\n");
  printf("h_mat_diag:\n  {");
  for (int i = 0; i < n_sites_ - 1; ++i)
    printf(" %5.5f ", h_mat_diag_[i]);
  printf("}\n");
  printf("h_mat_upper:\n  {");
  for (int i = 0; i < n_sites_ - 2; ++i)
    printf(" %5.5f ", h_mat_upper_[i]);
  printf("}\n");
  printf("h_mat_lower:\n  {");
  for (int i = 0; i < n_sites_ - 2; ++i)
    printf(" %5.5f ", h_mat_lower_[i]);
  printf("}\n");
  printf("k_eff:\n  {");
  for (int i = 0; i < n_sites_ - 2; ++i)
    printf(" %5.5f ", k_eff_[i]);
  printf("}\n\n\n");
}

/* The spec output for one filament is:
   (from Mesh::WriteSpec)
   double diameter
   double length
   double bond_length
   int n_sites
   double[3*n_sites] site_positions
   (from Filament::WriteSpec)
   double bending_stiffness
   double curvature (half the intrinsic curvature)
   uchar polymerization_state
    */
void Filament::WriteSpec(std::fstream &ospec) {
  Logger::Trace("Writing filament specs, object id: %d", GetOID());
  Mesh::WriteSpec(ospec);
  ospec.write(reinterpret_cast<char *>(&bending_stiffness_), sizeof(double));
  ospec.write(reinterpret_cast<char *>(&curvature_), sizeof(double));
  ospec.write(reinterpret_cast<char *>(&poly_), sizeof(unsigned char));
}

void Filament::ReadSpec(std::fstream &ispec) {
  if (ispec.eof())
    return;
  Mesh::ReadSpec(ispec);
  ispec.read(reinterpret_cast<char *>(&bending_stiffness_), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&curvature_), sizeof(double));
  ispec.read(reinterpret_cast<char *>(&poly_), sizeof(unsigned char));
  CalculateAngles();
  if (sparams_->highlight_handedness) {
    draw_ = draw_type::fixed;
    color_ = (curvature_ > 0 ? sparams_->color : sparams_->color + M_PI);
  }
}

/* double[3] avg_pos
   double[3] avg_scaled_pos
   double[3] avg_orientation
   double diameter
   double length
*/
void Filament::WritePosit(std::fstream &oposit) {
  double avg_pos[3], avg_u[3];
  GetAvgPosition(avg_pos);
  GetAvgOrientation(avg_u);
  std::copy(avg_pos, avg_pos + 3, position_);
  UpdatePeriodic();
  for (auto &pos : position_)
    oposit.write(reinterpret_cast<char *>(&pos), sizeof(pos));
  for (auto &spos : scaled_position_)
    oposit.write(reinterpret_cast<char *>(&spos), sizeof(spos));
  for (auto &u : avg_u)
    oposit.write(reinterpret_cast<char *>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char *>(&length_), sizeof(length_));
}

/* double[3] avg_pos
   double[3] avg_scaled_pos
   double[3] avg_orientation
   double diameter
   double length
*/
void Filament::ReadPosit(std::fstream &iposit) {
  if (iposit.eof())
    return;
  posits_only_ = true;
  double avg_pos[3], avg_u[3], s_pos[3];
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&avg_pos[i]), sizeof(double));
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&s_pos[i]), sizeof(double));
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&avg_u[i]), sizeof(double));
  iposit.read(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char *>(&length_), sizeof(length_));
  // Initialize first bond position
  std::copy(avg_pos, avg_pos + 3, position_);
  for (int i = 0; i < n_dim_; ++i)
    avg_pos[i] = avg_pos[i] - 0.5 * length_ * avg_u[i];
  for (int i_bond = 0; i_bond < n_bonds_; ++i_bond) {
    sites_[i_bond].SetPosition(avg_pos);
    sites_[i_bond].SetOrientation(avg_u);
    for (int i = 0; i < n_dim_; ++i) {
      avg_pos[i] += 0.5 * bond_length_ * avg_u[i];
    }
    bonds_[i_bond].SetPosition(avg_pos);
    bonds_[i_bond].SetOrientation(avg_u);
    bonds_[i_bond].SetDiameter(diameter_);
    bonds_[i_bond].UpdatePeriodic();
    // Set next bond position
    for (int i = 0; i < n_dim_; ++i) {
      avg_pos[i] += 0.5 * bond_length_ * avg_u[i];
    }
  }
  sites_[n_bonds_].SetPosition(avg_pos);
  sites_[n_bonds_].SetOrientation(avg_u);
  SetOrientation(avg_u);
  UpdatePeriodic();
  CalculateAngles();
}

void Filament::WriteCheckpoint(std::fstream &ocheck) {
  Mesh::WriteCheckpoint(ocheck);
}

void Filament::ReadCheckpoint(std::fstream &icheck) {
  Mesh::ReadCheckpoint(icheck);
}

void Filament::RotateToReferenceFrame() {
  /* Rotates filament site positions relative filament COM and translates
     filament to COM coordinates*/
  if (n_dim_ != 2)
    // TODO: Add 3D reference frame
    Logger::Error("Reference frame does not work for 3D"
                  " filaments yet!");
  /* If we have an even number of sites, then find the center bond and rotate
     filament relative to bond orientation and center filament on bond */
  double u0[3] = {0};
  double r0[3] = {0};
  if (n_sites_ % 2 == 0) {
    double const *const u = bonds_[n_sites_ / 2 - 1].GetOrientation();
    double const *const r = bonds_[n_sites_ / 2 - 1].GetPosition();
    std::copy(u, u + 3, u0);
    std::copy(r, r + 3, r0);
  } else {
    /* If we have an odd number of sites, center filament with respect to middle
       site and rotate filament relative to vector tangent to site */
    sites_[n_sites_ / 2].CalcTangent();
    double const *const u = sites_[n_sites_ / 2].GetTangent();
    double const *const r = sites_[n_sites_ / 2].GetPosition();
    std::copy(u, u + 3, u0);
    std::copy(r, r + 3, r0);
  }
  double new_pos[3] = {0};
  /* Only 2D rotation for now */
  for (auto site = sites_.begin(); site != sites_.end(); ++site) {
    const double *const pos = site->GetPosition();
    new_pos[0] = u0[1] * (pos[0] - r0[0]) - u0[0] * (pos[1] - r0[1]);
    new_pos[1] = u0[0] * (pos[0] - r0[0]) + u0[1] * (pos[1] - r0[1]);
    new_pos[2] = 0;
    site->SetPosition(new_pos);
  }

  /* Now update site orientation vectors */
  double u_mag, r_diff[3];
  for (int i_site = 0; i_site < n_sites_ - 1; ++i_site) {
    double const *const r_site1 = sites_[i_site].GetPosition();
    double const *const r_site2 = sites_[i_site + 1].GetPosition();
    u_mag = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_diff[i] = r_site2[i] - r_site1[i];
      u_mag += SQR(r_diff[i]);
    }
    u_mag = sqrt(u_mag);
    for (int i = 0; i < n_dim_; ++i)
      r_diff[i] /= u_mag;
    sites_[i_site].SetOrientation(r_diff);
  }
  sites_[n_sites_ - 1].SetOrientation(sites_[n_sites_ - 2].GetOrientation());

  /* Finally, update bond positions for visualization */
  UpdateBondPositions();
}
