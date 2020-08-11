#ifndef _CGLASS_SPECIES_BASE_H_
#define _CGLASS_SPECIES_BASE_H_

#include "object.hpp"
#include "params_parser.hpp"

class SpeciesBase {
private:
  species_id sid_;

protected:
  int n_members_ = 0;
  int spec_file_iterator_ = -1;
  static const system_parameters *params_;
  static const space_struct *space_;
  RNG rng_;
  std::fstream oposit_file_;
  std::fstream iposit_file_;
  std::fstream ospec_file_;
  std::fstream ispec_file_;
  std::string checkpoint_file_;
  std::vector<std::string> spec_file_names_;

public:
  SpeciesBase(unsigned long seed);
  static void SetParams(system_parameters *params);
  static void SetSpace(space_struct *space);
  void SetSID(species_id sid);
  virtual void UpdatePositions() {}
  virtual void Draw(std::vector<graph_struct *> &graph_array) {}
  virtual void Init(std::string spec_name, ParamsParser &parser) {}
  virtual void InitMembers() {}
  virtual void ZeroForces() {}
  virtual void GetInteractors(std::vector<Object *> &ix) {}
  virtual void GetLastInteractors(std::vector<Object *> &ix) {}
  virtual double GetPotentialEnergy() { return 0; }
  virtual void ScalePositions() {}
  virtual void AddMember() {}
  virtual void GetMemberPtrs(std::vector<Object *> &objs) {}
  virtual void SetLastMemberPosition(double const *const pos) {}
  virtual void ResetPreviousPositions() {}
  virtual void PopMember() {}
  virtual void PopAll() {}
  virtual const double GetSpecDiameter() const { return -1; }
  virtual const double GetSpecLength() const { return -1; }
  virtual void ArrangeMembers() {}
  virtual const bool CanOverlap() const { return -1; }
  species_id const GetSID() const { return sid_; }
  virtual void Report() {}
  virtual int const GetNMembers() const { return n_members_; }
  virtual int const GetNInsert() const { return -1; }
  virtual int const GetNPosit() const { return -1; }
  virtual int const GetNSpec() const { return -1; }
  virtual int const GetNCheckpoint() const { return params_->n_checkpoint; }
  virtual bool const GetPositFlag() const { return false; }
  virtual bool const GetSpecFlag() const { return false; }
  virtual bool const GetCheckpointFlag() const {
    return params_->checkpoint_flag;
  }
  virtual std::string GetInsertionType() const { return ""; }
  virtual bool CheckInteractionAnalysis() { return false; }
  virtual int GetCount() const { return 0; }
  virtual void WriteOutputs(std::string run_name) {}
  virtual void WritePosits() {}
  virtual void WriteSpecs() {}
  virtual void WriteCheckpoints() {}
  virtual void ReadSpecs() {}
  virtual void ReadCheckpoints() {}
  virtual void ReadPosits() {}
  virtual void ReadPositsFromSpecs() {}
  virtual void InitAnalysis() {}
  virtual void RunAnalysis() {}
  virtual void FinalizeAnalysis() {}
  virtual void InitOutputFiles(std::string run_name);
  virtual void InitPositFile(std::string run_name);
  virtual void InitSpecFile(std::string run_name);
  virtual void InitPositFileInput(std::string run_name);
  virtual void InitSpecFileInput(std::string run_name);
  virtual bool InitSpecFileInputFromFile(std::string run_name);
  virtual bool HandleEOF();
  virtual void InitInputFiles(std::string run_name, bool posits_only,
                              bool with_reloads);
  virtual void InitCheckpoints(std::string run_name);
  virtual void LoadFromCheckpoints(std::string run_name,
                                   std::string checkpoint_run_name);
  virtual int OutputIsOpen() { return oposit_file_.is_open(); }
  virtual int InputIsOpen() { return iposit_file_.is_open(); }
  virtual void CloseFiles();
  virtual void CleanUp() {}
  virtual void Reserve() {}
  virtual double const GetVolume() { return 0; }
  virtual double const GetDrMax() { return 0; }
  virtual void ZeroDrTot() {}
  virtual void CustomInsert() {}
  virtual const bool CheckInteractorUpdate() { return false; }
  virtual const std::string GetSpeciesName() const { return "base"; }
};

#endif
