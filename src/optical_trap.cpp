/**
 * @author      : Adam Lamson (adam.lamson@colorado.edu)
 * @file        : optical_trap
 * @created     : Friday Jul 03, 2020 14:49:10 MDT
 */

#include "cglass/optical_trap.hpp"

OpticalTrap::OpticalTrap(unsigned long seed) : Object(seed), bead_(seed) {
  SetSID(species_id::optical_trap);
}

/*! \brief Initialize optical trap object
 *
 *  Detailed description
 *
 * \param *sparams Parameter description
 * \return Return parameter description
 */
void OpticalTrap::Init(optical_trap_parameters *sparams) {
  sparams_ = sparams;
  trap_spring_ = sparams_->trap_spring;
  diameter_ = sparams_->trap_diameter;
  color_ = sparams_->trap_color;

  bead_.SetSpace(space_);
  bead_.SetNDim(n_dim_);
  bead_.SetDiameter(sparams_->bead_diameter);
  bead_.SetColor(sparams_->bead_color, draw_type::fixed);
  //bead_=Object(0)
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
}

/*! \brief Set position of optical trap relative to obj minus end
 *
 *  Detailed description
 *
 * \param Object *obj Parameter description
 * \return Return parameter description
 */
void OpticalTrap::InsertAndAttach(Object *obj) {
  attach_obj_ = obj;
  SetMeshID(-1);

  AttachObjRelLambda(0); //TODO attach to minus end for now
  UpdateBeadPosition();
  // Bead and trap are initially set to same position
  SetPosition(bead_.GetPosition());
  UpdatePeriodic();

  return;
}

void OpticalTrap::CalculateOpticalTrapForce() {
  ZeroForce();
  UpdateBeadPosition();
  double const *const bead_pos = bead_.GetPosition();
  // Calculate the force on bead from trap
  for (int i = 0; i < n_dim_; ++i) {
    force_[i] = trap_spring_ * (position_[i] - bead_pos[i]);
  }
}

/*! \brief Apply force from optical trap to bond
 *
 * \return void
 */
void OpticalTrap::ApplyOpticalTrapForce() {
  if (!bond_) {
    Logger::Error("Optical trap attempted to apply forces to nullptr bond");
  }
  CalculateOpticalTrapForce();
  bond_->AddForce(force_);
  double dlambda[3] = {0};
  double const *const bond_orientation = bond_->GetOrientation();
  for (int i = 0; i < n_dim_; ++i) {
    dlambda[i] = (bond_lambda_ - 0.5 * bond_length_) * bond_orientation[i];
  }
  //printf("dlambda = (%f, %f, %f)\n", dlambda[0], dlambda[1], dlambda[2]);
  cross_product(dlambda, force_, torque_, 3);
  //printf("torque = (%f, %f, %f)\n", torque_[0], torque_[1], torque_[2]);
  bond_->AddTorque(torque_);
}

void OpticalTrap::AttachObjRelLambda(double lambda) {
  //if (attach_obj_->GetType() != +obj_type::bond) {
  //  Logger::Error(
  //      "Optical traps bound to non-bond objects not yet implemented in "
  //      "AttachObjLambda.");
  //}

  if (attach_obj_->GetType() == +obj_type::bond) {
    bond_ = dynamic_cast<Bond *>(attach_obj_);
    bond_length_ = bond_->GetLength();
    if (bond_ == nullptr) {
      Logger::Error(
          "Object ptr passed to optical trap was not referencing a bond!");
    }
    mesh_ = dynamic_cast<Mesh *>(bond_->GetMeshPtr());
    if (mesh_ != nullptr) {
      mesh_n_bonds_ = mesh_->GetNBonds();
      mesh_length_ = mesh_->GetTrueLength();
      mesh_lambda_ = mesh_length_ * lambda;
      bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
      SetMeshID(bond_->GetMeshID());
      //  Logger::Error("Object ptr passed to anchor was not referencing a mesh!");
    } else {
      bond_lambda_ = lambda * bond_length_;
    }
  } else if (attach_obj_->IsMesh()) {
    mesh_ = dynamic_cast<Mesh *>(attach_obj_);
    mesh_n_bonds_ = mesh_->GetNBonds();
    mesh_length_ = mesh_->GetTrueLength();
    mesh_lambda_ = mesh_length_ * lambda;
    bond_ = mesh_->GetBondAtLambda(mesh_lambda_);
    bond_length_ = bond_->GetLength();
    bond_lambda_ = mesh_lambda_ - bond_->GetMeshLambda();
    SetMeshID(mesh_->GetMeshID());
  } else {
    Logger::Error(
        "Optical traps for non-bond or non-mesh objects not yet implemented in "
        "AttachObjLambda.");
  }

  if (bond_lambda_ < 0 || bond_lambda_ > bond_length_) {
    printf("bond_lambda: %2.2f\n", bond_lambda_);
    Logger::Error("Lambda passed to anchor does not match length of "
                  "corresponding bond! lambda: %2.2f, bond_length: %2.2f ",
                  bond_lambda_, bond_length_);
  }
}

void OpticalTrap::UpdateBeadPosition() {
  double bead_pos[3] = {};
  double const *const bond_position = bond_->GetPosition();
  double const *const bond_orientation = bond_->GetOrientation();
  for (int i = 0; i < n_dim_; ++i) {
    bead_pos[i] = bond_position[i] +
                  (bond_lambda_ - 0.5 * bond_length_) * bond_orientation[i];
  }
  bead_.SetPosition(bead_pos);
  bead_.UpdatePeriodic();
}

void OpticalTrap::Draw(std::vector<graph_struct *> &graph_array) {
  Object::Draw(graph_array);
  UpdateBeadPosition();
  bead_.Draw(graph_array);
}

void OpticalTrap::WriteSpec(std::fstream &ospec) {
  Logger::Trace("Writing optical trap specs, mesh id: %d", GetMeshID());
  Object::WriteSpec(ospec);
  UpdateBeadPosition();
  double const *const bead_pos = bead_.GetPosition();
  double const *const bead_spos = bead_.GetScaledPosition();
  for (int i = 0; i < 3; ++i) {
    double bpos = bead_pos[i];
    ospec.write(reinterpret_cast<char *>(&bpos), sizeof(bpos));
  }
  for (int i = 0; i < 3; ++i) {
    double bspos = bead_spos[i];
    ospec.write(reinterpret_cast<char *>(&bspos), sizeof(bspos));
  }
  /* TODO: Make this an object ID one day <10-08-20, ARL> */
  int attach_id = attach_obj_->GetMeshID();
  ospec.write(reinterpret_cast<char *>(&attach_id), sizeof(attach_id));
}

void OpticalTrap::ReadSpec(std::fstream &ispec) {
  if (ispec.eof())
    return;
  Object::ReadSpec(ispec);
  double bead_pos[3], bead_spos[3];
  for (auto &bpos : bead_pos)
    ispec.read(reinterpret_cast<char *>(&bpos), sizeof(bpos));
  for (auto &bspos : bead_spos)
    ispec.read(reinterpret_cast<char *>(&bspos), sizeof(bspos));
  int attach_id;
  ispec.read(reinterpret_cast<char *>(&attach_id), sizeof(int));
  UpdatePeriodic();
};
