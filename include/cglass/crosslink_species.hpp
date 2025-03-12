#ifndef _CGLASS_CROSSLINK_SPECIES_H_
#define _CGLASS_CROSSLINK_SPECIES_H_

#include "crosslink.hpp"
#include "species.hpp"
#include <KMC/kmc.hpp>

typedef std::vector<std::pair<std::vector<Crosslink>::iterator,
                              std::vector<Crosslink>::iterator>>
    xlink_chunk_vector;
typedef std::vector<Crosslink>::iterator xlink_iterator;

class CrosslinkSpecies : public Species<Crosslink, species_id::crosslink> {
private:
  bool *update_;
  std::string checkpoint_file_;
  double *obj_size_; // Total area/length of all the objects in the system
  double xlink_concentration_;
  int begin_with_bound_crosslinks_;
  double bind_site_density_;
  bool infinite_reservoir_flag_;
  double k_on_;
  double k_off_;
  double random_obj_probability_;//Made global so probability can be added to bound_curr
  bool static_flag_;
  bool use_bind_file_ = false; // Whether to use a file with species names and parameters associated
  std::vector<std::string> bind_species_;
  Tracker *tracker_ = nullptr;
  LookupTable lut_;
  std::vector<Object *> *objs_;
  std::map<Sphere *, std::pair<std::vector<double>, std::vector<std::pair<Anchor*, std::string> > > > *bound_curr_;
  
  // use a map of species names to binding parameters to
  // store binding parameters for specific species for each anchor
  std::vector<std::map<std::string, bind_params> > bind_param_map_;

  // A default binding parameter list to initialize bind_param_map to
  std::vector<bind_params> default_bind_params_;

  // Binding factor- k_on_s * bind_site_density * object_amount
  double bind_rate_;
  int step_count_;
  void InitializeBindParams();
  LUTFiller *MakeLUTFiller();
  void CalculateBindingFree();
  void BindCrosslink();
  void BindDoubly(Object* receptor_one, Object* receptor_two);
  void UpdateBoundCrosslinks();
  void UpdateBoundCrosslinkForces();
  void UpdateBoundCrosslinkPositions();
  void ApplyCrosslinkTetherForces();
  std::pair<Object*, int> GetRandomObject();
  std::pair<Object*, int> GetRandomObjectKnockout();
  bool* global_check_for_cross_;

public:
  CrosslinkSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void LoadBindingSpecies();
  void InitInteractionEnvironment(std::vector<Object *> *objs, double *obj_size, Tracker *tracker, bool *update,
                                  std::map<Sphere *, std::pair<std::vector<double>, 
                                  std::vector<std::pair<Anchor*, std::string> > > > *bound_curr);
  void TestKMCStepSize();
  void SetGlobalCheckForCrossPointer(bool* check);
  void GetInteractors(std::vector<Object *> &ixors);
  void UpdatePositions();
  void UpdateBindRate();
  void CleanUp();
  void ClearNeighbors();
  void Draw(std::vector<graph_struct *> &graph_array);
  void BindCrosslinkObj(Object *obj);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void AddMember();
  void KnockoutBind(Sphere* receptor);
  void InsertAttachedCrosslinksSpecies(std::vector<std::vector<Object *>> receptor_list);
  void GetAnchorInteractors(std::vector<Object *> &ixors);
  void ReadSpecs();
  void InsertCrosslinks();
  Crosslink* GetCrosslink(int i); 
  const int GetDoublyBoundCrosslinkNumber() const;
  const double GetConcentration() const;
  const double GetRCutoff() const;
};

#endif
