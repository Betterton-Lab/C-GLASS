#ifndef _CGLASS_ANCHOR_H_
#define _CGLASS_ANCHOR_H_

#include "mesh.hpp"
#include "neighbor_list.hpp"

/* Class for bound crosslink heads (called anchors). Tracks and updates its
   absolute position in space and relative position to its bound object. */
class Anchor : public Object {
 private:
  bool bound_;
  bool static_flag_;
  bool active_;
  bool end_pausing_;
  crosslink_parameters *sparams_;
  int step_direction_;

  double bond_length_;
  double bond_lambda_;
  double mesh_length_;
  double mesh_lambda_;
  double max_velocity_s_;
  double max_velocity_d_;
  double diffusion_s_;
  double diffusion_d_;
  double kick_amp_s_;
  double kick_amp_d_;
  double k_on_s_;
  double k_on_d_;
  double k_off_s_;
  double k_off_d_;
  double polar_affinity_;
  double f_stall_;
  double force_dep_vel_flag_;

  bind_state state_;

  NeighborList neighbors_;

  Bond *bond_;
  Site *site_ = nullptr;
  Mesh *mesh_;

  double *obj_area_ = nullptr;

  int mesh_n_bonds_;

  void UpdateAnchorPositionToBond();
  void Diffuse();
  void Walk();
  bool CheckMesh();
  bool CalcBondLambda();

 public:
  Anchor(unsigned long seed);
  void Init(crosslink_parameters *sparams);
  bool IsBound();
  void UpdatePosition();
  void Activate();
  void Deactivate();
  void ApplyAnchorForces();
  void UpdateAnchorPositionToMesh();
  void SetDiffusion();
  void AttachObjRandom(Object *o);
  void AttachObjLambda(Object *o, double lambda);
  void AttachObjCenter(Object *o);
  void AttachObjMeshLambda(Object *o, double mesh_lambda);
  void AttachObjMeshCenter(Object *o);
  void CalculatePolarAffinity(std::vector<double> &doubly_binding_rates);
  void SetBondLambda(double l);
  void SetMeshLambda(double ml);
  void SetBound();
  void Unbind();
  int const GetBoundOID();
  void Draw(std::vector<graph_struct *> &graph_array);
  void AddNeighbor(Object *neighbor);
  void ClearNeighbors();
  const Object *const *GetNeighborListMem();
  const std::vector<Bond*>& GetNeighborListMemBonds();
  const std::vector<Site*>& GetNeighborListMemSites();
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);
  void BindToPosition(double *bind_pos);
  void SetStatic(bool static_flag);
  void SetState(bind_state state);

  double const GetMeshLambda();
  double const GetBondLambda();
  Object *GetNeighbor(int i_neighbor);
  Site *GetSiteNeighbor(int i_neighbor);
  Bond *GetBondNeighbor(int i_neighbor);
  const int GetNNeighbors() const;
  const int GetNNeighborsSite() const;
  const int GetNNeighborsBond() const;
  const double GetOnRate() const;
  const double GetOffRate() const;
  const double GetMaxVelocity() const;
  const double GetDiffusionConst() const;
  const double GetKickAmplitude() const;
  const double* const GetObjArea();
  void SetObjArea(double* obj_area);
};

#endif
