#ifndef _CGLASS_CROSSLINK_MANAGER_H_
#define _CGLASS_CROSSLINK_MANAGER_H_

#include "output_manager.hpp"
#include "species_factory.hpp"

class CrosslinkOutputManager : public OutputManagerBase<CrosslinkSpecies> {
  /* Do not write thermo - base output manager will handle that */
  void WriteThermo() {}
  void ReadThermo() {}
  void InitThermo(std::string fname) {}
};

class CrosslinkManager {
 private:
  system_parameters *params_;
  CrosslinkOutputManager output_mgr_;
  double obj_volume_;
  double rcutoff_ = 0;  // Cutoff for binding any crosslink and bond
  bool update_;
  std::vector<CrosslinkSpecies *> xlink_species_;
  std::vector<Object *> *objs_;
  space_struct *space_;

 public:
  void Init(system_parameters *params, space_struct *space,
            std::vector<Object *> *objs);
  void GetInteractors(std::vector<Object *> &ixors);
  void UpdateCrosslinks();
  void UpdateObjsVolume();
  bool CheckUpdate();
  void Clear();
  void Draw(std::vector<graph_struct *> &graph_array);
  void AddNeighborToAnchor(Object *anchor, Object *neighbor);
  void WriteOutputs();
  void InitOutputs(bool reading_inputs = false,
                   run_options *run_opts = nullptr);
  void GetAnchorInteractors(std::vector<Object *> &ixors);
  void InitSpecies(sid_label &slab, ParamsParser &parser, unsigned long seed);
  void LoadCrosslinksFromCheckpoints(std::string run_name,
                                     std::string checkpoint_run_name);
  const int GetDoublyBoundCrosslinkNumber() const;
  void ZeroDrTot();
  const double GetDrMax();
  void ReadInputs();
  void InsertCrosslinks();
  const double GetRCutoff() const {
    Logger::Trace("Crosslink rcutoff is %2.2f", rcutoff_);
    return rcutoff_;
  }
};

#endif
