/**
 * @author      : Adam Lamson (adam.lamson@colorado.edu)
 * @file        : optical_trap
 * @created     : Friday Jul 03, 2020 14:49:10 MDT
 */

#ifndef OPTICAL_TRAP_HPP

#define OPTICAL_TRAP_HPP

#include "mesh.hpp"

class OpticalTrap : public Object {
private:
  optical_trap_parameters *sparams_;

  Object *attach_obj_;

  Bond *bond_ = nullptr;
  double bond_length_ = -1;
  double bond_lambda_ = -1;

  Mesh *mesh_ = nullptr;
  int mesh_n_bonds_ = -1;
  double mesh_length_ = -1;
  double mesh_lambda_ = -1;

  double trap_spring_;

  Object bead_; //Purely for graphing right now

public:
  OpticalTrap(unsigned long seed);
  void Init(optical_trap_parameters *sparams);
  void InsertAndAttach(Object *obj);
  void AttachObjRelLambda(double lambda);
  void ApplyOpticalTrapForce();
  void CalculateOpticalTrapForce();
  void UpdateBeadPosition();
  double const GetMeshLambda() { return mesh_lambda_; };
  double const GetBondLambda() { return bond_lambda_; };
  void Draw(std::vector<graph_struct *> &graph_array);
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);

  //virtual ~optical_trap();
};

#endif /* end of include guard OPTICAL_TRAP_HPP */

