#ifndef _CGLASS_RIGID_FILAMENT_H_
#define _CGLASS_RIGID_FILAMENT_H_

#include "exponential_dist.hpp"
#include "flory_schulz.hpp"
#include "mesh.hpp"

class RigidFilament : public Mesh {
private:
  rigid_filament_parameters *sparams_;
  double gamma_par_ = 0;
  double gamma_perp_ = 0;
  double gamma_rot_ = 0;
  double diffusion_par_ = 0;
  double diffusion_perp_ = 0;
  double diffusion_rot_ = 0;
  double body_frame_[6];
  double min_length_ = 0;
  double max_length_ = 1000;
  double constrain_vec_[3] = {}; // Rigid filaments move only in plane defined
                                 // by this unit vector.

  bool zero_temperature_ = false;
  bool pseudo1D_ = false;
  int n_step_ = 0;
  int eq_steps_ = 0;
  int eq_steps_count_ = 0;

  void UpdateSiteBondPositions();
  void SetDiffusion();
  void GenerateProbableOrientation();
  void AddRandomForces();

  void UpdateSitePositions();

  void ApplyForcesTorques();
  void ApplyInteractionForces();

  void SetParameters();
  void InitRigidFilamentLength();
  void UpdateAvgPosition();
  void DynamicInstability();
  void UpdatePolyState();
  void GrowFilament();
  void ReportAll();

  void CalculateBinding();

protected:
  void InsertRigidFilament(std::string insertion_type, double buffer = -1);
  void GetBodyFrame();
  void AddRandomDisplacement();
  void AddRandomReorientation();

public:
  RigidFilament(unsigned long seed);
  virtual void Init(rigid_filament_parameters *sparams);
  virtual void InsertAt(const double *const new_pos, const double *const u);
  virtual void Integrate();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void UpdatePosition();
  double const GetLength() { return length_; }
  double const GetTrueLength() const { return length_; }
  void CheckFlocking();

  void GetNematicOrder(double *nematic_order_tensor);
  void GetPolarOrder(double *polar_order_vector);
  double GetTipZ() { return sites_[n_sites_ - 1].GetOrientation()[n_dim_ - 1]; }
  double const *const GetHeadPosition() {
    return sites_[n_sites_ - 1].GetPosition();
  }
  void SetConstrainVec(const double *constrain_vec) {
    std::copy(constrain_vec, constrain_vec + 3, constrain_vec_);
  }
  void WritePosit(std::fstream &oposit);
  void ReadPosit(std::fstream &iposit);
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);
  void WriteCheckpoint(std::fstream &ocheck);
  void ReadCheckpoint(std::fstream &icheck);
  void ScalePosition();
  double const GetVolume();
};

typedef std::vector<RigidFilament>::iterator rigid_filament_iterator;
typedef std::vector<std::pair<std::vector<RigidFilament>::iterator,
                              std::vector<RigidFilament>::iterator>>
    rigid_filament_chunk_vector;

#endif // _CGLASS_RIGID_FILAMENT_H_
