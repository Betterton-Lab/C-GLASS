#include "cglass/rigid_filament.hpp"

RigidFilament::RigidFilament(unsigned long seed) : Mesh(seed) {
  std::fill(body_frame_, body_frame_ + 6, 0.0);
  SetSID(species_id::rigid_filament);
}

void RigidFilament::SetParameters() {
  /* Read parameters from filament parameters */
  color_ = sparams_->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams_->length;
  diameter_ = sparams_->diameter;
  max_length_ = sparams_->max_length;
  min_length_ = sparams_->min_length;
  zero_temperature_ = params_->zero_temperature; // include thermal forces
  eq_steps_count_ = 0;
  pseudo1D_ = sparams_->pseudo1D;

  /* Refine parameters */
}

void RigidFilament::Init(rigid_filament_parameters *sparams) {
  sparams_ = sparams;
  SetParameters();
  InitRigidFilamentLength();
  Reserve();
  InsertRigidFilament(sparams_->insertion_type, -1);
  SetDiffusion();
}

/* Returns number of bonds to initialize */
void RigidFilament::InitRigidFilamentLength() {
  if (max_length_ < min_length_) {
    Logger::Warning(
        "Minimum rigid filament length larger than max length -- setting "
        "max_length_ = min_length_");
    max_length_ = min_length_;
  }

  if (length_ > max_length_) {
    Logger::Warning(
        "RigidFilament length larger than max length -- setting length = "
        "max_length");
    length_ = max_length_;
  } else if (length_ < min_length_) {
    Logger::Warning(
        "RigidFilament length less than min length -- setting length = "
        "min_length");
    length_ = min_length_;
  }
  Logger::Trace(
      "RigidFilament initialized with length %2.2f with %d bonds, mesh_id:"
      " %d",
      length_, n_bonds_, GetMeshID());
  true_length_ = length_;
  bond_length_ = length_;
  n_bonds_max_ = 1;
}

void RigidFilament::InsertRigidFilament(std::string insertion_type,
                                        double buffer) {
  if (buffer < 0) {
    buffer = length_;
  }
  if (insertion_type.compare("random") == 0 ||
      insertion_type.compare("custom") == 0) {
    AddRandomBondAnywhere(length_, diameter_);
    SetPosition(bonds_.back().GetPosition());
    SetOrientation(bonds_.back().GetOrientation());
    UpdateBondPositions();
  } else if (insertion_type.compare("random_oriented") == 0) {
    AddRandomBondAnywhere(length_, diameter_);
    double orient[3] = {0};
    orient[n_dim_ - 1] = 1.0;
    bonds_.back().SetOrientation(orient);
    SetPosition(bonds_.back().GetPosition());
    SetOrientation(bonds_.back().GetOrientation());
    UpdateBondPositions();
  } else {
    Logger::Error("Rigid Filament insertion type not recognized!");
  }
}

/* Integration scheme taken from Yu-Guo Tao,
   J Chem Phys 122 244903 (2005)
   Explicit calculation of the friction tensor acting on force vector,

   r(t+dt) = r(t) + (Xi^-1 . F_s(t)) * dt + dr(t),

   where friction tensor Xi = gamma_par * |u><u| + gamma_perp * (I - |u><u|),
   u is the orientation, F_s is the force from interactions, and dr(t) is the
   random displacement due to random force,

   dr(t) = Xi^-1 . F_r * dt,

   which is treated separately as a random displacement with std dev
   sqrt(2*kT*dt/gamma_(par/perp)) along par/perp unit vectors
   relative to rod. */
void RigidFilament::Integrate() {
  // Explicit calculation of Xi.F_s
  double fric_mat[9] = {}; // Friction matrix for filament
  double mob_mat[9] = {};  // Mobility matrix for filament
  double force_eff[3];
  double torque_eff[3];
  std::copy(force_, force_ + 3, force_eff);
  std::copy(torque_, torque_ + 3, torque_eff);

  // Constrain the motion of filaments to a plane defined by constrain_vec_
  // by subtracting force components out of the plane and keeping torque
  // components only along constrain_vec_ direction.
  if (sparams_->constrain_motion_flag) {
    double f_dot_u = dot_product(3, force_, constrain_vec_);
    double t_dot_u = dot_product(3, torque_, constrain_vec_);
    for (int i = 0; i < 3; ++i) {
      force_eff[i] -= f_dot_u * constrain_vec_[i];
      torque_eff[i] = t_dot_u * constrain_vec_[i];
    }
  }
  // Construct mobility matrix
  for (int i = 0; i < n_dim_; ++i) {
    for (int j = i + 1; j < n_dim_; ++j) {
      fric_mat[i * n_dim_ + j] = fric_mat[j * n_dim_ + i] =
          (gamma_par_ - gamma_perp_) * orientation_[i] * orientation_[j];
    }
    fric_mat[i * n_dim_ + i] += gamma_perp_;
  }
  if (n_dim_ == 2)
    invert_sym_2d_matrix(fric_mat, mob_mat);
  else if (n_dim_ == 3)
    invert_sym_3d_matrix(fric_mat, mob_mat);

  for (int i = 0; i < n_dim_; ++i) {
    for (int j = 0; j < n_dim_; ++j) {
      position_[i] += mob_mat[i * n_dim_ + j] * force_eff[j] * delta_;
    }
  }
  // Reorientation due to external torques
  double du[3];
  cross_product(torque_eff, orientation_, du, 3); // ndim=3 since torques
  for (int i = 0; i < n_dim_; ++i) {
    orientation_[i] += du[i] * delta_ / gamma_rot_;
  }
  normalize_vector(orientation_, n_dim_);
  if (!zero_temperature_) {
    // Add the random displacement dr(t)
    AddRandomDisplacement();
    // Update the orientation due to torques and random rotation
    if (!pseudo1D_) {
      AddRandomReorientation();
    }
  }
  // double f_mag = sqrt(dot_product(n_dim_, force_, force_));
  // printf("f_mag = %f\n", f_mag);

  UpdatePeriodic();
  UpdateSitePositions();
  UpdateBondPositions();
}

/* Calculates body frame, which returns the vector(s) orthogonal
   to u(t), then applies random displacements along each
   orthogonal vector and along u(t) pulled from a distribution
   with std dev sqrt(2*kT*dt/gamma) where gamma is the friction
   coefficient along that direction */
void RigidFilament::AddRandomDisplacement() {
  double mag = 0;
  // Get vector(s) orthogonal to orientation
  GetBodyFrame();
  // First handle the parallel component
  if (!pseudo1D_) {
    mag = rng_.RandomNormal(diffusion_par_);
    for (int i = 0; i < n_dim_; ++i) position_[i] += mag * orientation_[i];
  }
  // Then the perpendicular component(s)  
  for (int j = 0; j < n_dim_ - 1; ++j) {
    mag = rng_.RandomNormal(diffusion_perp_);
    for (int i = 0; i < n_dim_; ++i)
      position_[i] += mag * body_frame_[n_dim_ * j + i];
  }
  // Handle the random orientation update after updating orientation from
  // interaction torques
}

/* The orientation update is also from Yu-Guo Tao,

   u(t+dt) = u(t) + gamma_rot^-1 * T_s(t) x u(t) * dt + du(t)

   where similar to above, du(t) is the reorientation due to
   random forces, and is treated as random displacement vector(s)
   orthogonal to u(t) with std dev sqrt(2*kT*dt/gamma_rot) */
void RigidFilament::AddRandomReorientation() {
  // Now handle the random orientation update
  for (int j = 0; j < n_dim_ - 1; ++j) {
    double mag = rng_.RandomNormal(diffusion_rot_);
    for (int i = 0; i < n_dim_; ++i) {
      orientation_[i] += mag * body_frame_[n_dim_ * j + i];
    }
  }
  normalize_vector(orientation_, n_dim_);
}

void RigidFilament::GetBodyFrame() {
  if (n_dim_ == 2) {
    body_frame_[0] = orientation_[1];
    body_frame_[1] = -orientation_[0];
  } else {
    double vect1[3] = {1.0, 0.0, 0.0};
    double vect2[3] = {0.0, 1.0, 0.0};
    if (1.0 - ABS(orientation_[0]) > 1e-2)
      cross_product(orientation_, vect1, &(body_frame_[0]), n_dim_);
    else
      cross_product(orientation_, vect2, &(body_frame_[0]), n_dim_);
    normalize_vector(&(body_frame_[0]), n_dim_);
    cross_product(orientation_, &(body_frame_[0]), &(body_frame_[3]), n_dim_);
  }
}

void RigidFilament::InsertAt(const double *const new_pos,
                             const double *const u) {
  Logger::Trace(
      "Inserting rigid filament at [%2.1f, %2.1f, %2.1f] with orientation"
      "[%2.1f, %2.1f, %2.1f]",
      new_pos[0], new_pos[1], new_pos[2], u[0], u[1], u[2]);
  RelocateMesh(new_pos, u);
  std::copy(new_pos, new_pos + 3, position_);
  SetDiffusion();
}

void RigidFilament::SetDiffusion() {
  // Sets the diffusion in accordance to Lowen, Phys. Rev. E, 1994.
  double L = length_ + diameter_;
  double p = L / diameter_;
  double log_p = log(p);
  gamma_par_ = 2.0 * L / 3.0 / (log_p - 0.207 + 0.980 / p - 0.133 / SQR(p));
  gamma_perp_ = 4.0 * L / 3.0 / (log_p + 0.839 + 0.185 / p + 0.233 / SQR(p));
  gamma_rot_ = CUBE(L) / 9.0 / (log_p - 0.662 + 0.917 / p - 0.050 / SQR(p));
  diffusion_par_ = sqrt(2 * delta_ / gamma_par_);
  diffusion_perp_ = sqrt(2 * delta_ / gamma_perp_);
  diffusion_rot_ = sqrt(2 * delta_ / gamma_rot_);
}

double const RigidFilament::GetVolume() {
  if (n_dim_ == 2) {
    return diameter_ * length_ + 0.25 * M_PI * diameter_ * diameter_;
  } else {
    return 0.25 * M_PI * diameter_ * diameter_ * length_ +
           1.0 / 6.0 * M_PI * diameter_ * diameter_ * diameter_;
  }
}

void RigidFilament::UpdatePosition() {
  ApplyForcesTorques();
  if (!sparams_->stationary_flag)
    Integrate();
  eq_steps_count_++;
}

void RigidFilament::GetNematicOrder(double *nematic_order_tensor) {
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

void RigidFilament::GetPolarOrder(double *polar_order_vector) {
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

void RigidFilament::UpdateSitePositions() {
  double s0_pos[3] = {0};
  double s1_pos[3] = {0};
  for (int i = 0; i < n_dim_; ++i) {
    s0_pos[i] = position_[i] - .5 * length_ * orientation_[i];
    s1_pos[i] = position_[i] + .5 * length_ * orientation_[i];
  }
  sites_[0].SetPosition(s0_pos);
  sites_[1].SetPosition(s1_pos);
}

void RigidFilament::ApplyForcesTorques() {
  const double *force = bonds_.back().GetForce();
  const double *torque = bonds_.back().GetTorque();
  if (!pseudo1D_) {
    for (int i = 0; i < 3; ++i) {
      force_[i] = force[i];
      torque_[i] = torque[i];
    }
  } else {
      force_[0] = 0;
      force_[1] = force[1];
      force_[2] = 0;
      for (int i =0; i < 3; ++i) {
        torque_[i]=0;
      }
    } 
}

void RigidFilament::ApplyInteractionForces() {
  double pure_torque[3] = {0, 0, 0};
  double site_force[3] = {0, 0, 0};
  double linv = 1.0 / bond_length_;
  for (int i = 0; i < n_bonds_; ++i) {
    double const *const f = bonds_[i].GetForce();
    double const *const t = bonds_[i].GetTorque();
    double const *const u = sites_[i].GetOrientation();
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
  }
}

void RigidFilament::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto bond = bonds_.begin(); bond != bonds_.end(); ++bond) {
    bond->Draw(graph_array);
  }
}

// Scale bond and site positions from new unit cell
void RigidFilament::ScalePosition() {
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

void RigidFilament::ReportAll() {}

void RigidFilament::WriteSpec(std::fstream &ospec) {
  Logger::Trace("Writing rigid filament specs, object id: %d", GetOID());
  Mesh::WriteSpec(ospec);
}

void RigidFilament::ReadSpec(std::fstream &ispec) {
  if (ispec.eof())
    return;
  Mesh::ReadSpec(ispec);
}

/* double[3] avg_pos
   double[3] avg_scaled_pos
   double[3] avg_orientation
   double diameter
   double length
*/
void RigidFilament::WritePosit(std::fstream &oposit) {
  UpdatePeriodic();
  for (auto &pos : position_)
    oposit.write(reinterpret_cast<char *>(&pos), sizeof(pos));
  for (auto &spos : scaled_position_)
    oposit.write(reinterpret_cast<char *>(&spos), sizeof(spos));
  for (auto &u : orientation_)
    oposit.write(reinterpret_cast<char *>(&u), sizeof(u));
  oposit.write(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  oposit.write(reinterpret_cast<char *>(&length_), sizeof(length_));
  auto mesh_id = GetMeshID();
  oposit.write(reinterpret_cast<char *>(&mesh_id), sizeof(mesh_id));
}

/* double[3] avg_pos
   double[3] avg_scaled_pos
   double[3] avg_orientation
   double diameter
   double length
*/
void RigidFilament::ReadPosit(std::fstream &iposit) {
  if (iposit.eof())
    return;
  posits_only_ = true;
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&position_[i]), sizeof(double));
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&scaled_position_[i]), sizeof(double));
  for (int i = 0; i < 3; ++i)
    iposit.read(reinterpret_cast<char *>(&orientation_[i]), sizeof(double));
  iposit.read(reinterpret_cast<char *>(&diameter_), sizeof(diameter_));
  iposit.read(reinterpret_cast<char *>(&length_), sizeof(length_));
  int mesh_id = 0;
  iposit.read(reinterpret_cast<char *>(mesh_id), sizeof(int));
  SetMeshID(mesh_id);
  UpdateBondPositions();
}

void RigidFilament::WriteCheckpoint(std::fstream &ocheck) {
  Mesh::WriteCheckpoint(ocheck);
}

void RigidFilament::ReadCheckpoint(std::fstream &icheck) {
  Mesh::ReadCheckpoint(icheck);
}
