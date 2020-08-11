#ifndef _CGLASS_SPECIES_H_
#define _CGLASS_SPECIES_H_

#include "species_base.hpp"
#include <stdexcept>
#ifdef ENABLE_OPENMP
#include "omp.h"
#endif
#include "analysis.hpp"

template <typename T, unsigned char S>
class Species : public SpeciesBase {
protected:
  std::vector<T> members_;
  species_parameters<S> sparams_;
  std::vector<Analysis<T, S> *> analysis_;
  virtual void LoadAnalysis() {
    /* Check parameters for analyses and load them into analysis_ here */
  }

public:
  Species(unsigned long seed) : SpeciesBase(seed) {}
  // Initialize function for setting it up on the first pass
  virtual void Init(std::string spec_name, ParamsParser &parser) {
    SetSID(species_id::_from_integral(S));
    species_base_parameters *sparams =
        parser.GetNewSpeciesParameters(GetSID(), spec_name);
    sparams_ = *dynamic_cast<species_parameters<S> *>(sparams);
    delete sparams;
    Logger::Debug("Initializing %s %s", GetSID()._to_string(),
                  GetSpeciesName().c_str());
  }
  // Virtual functions
  virtual const double GetSpecDiameter() const { return sparams_.diameter; }
  virtual const double GetSpecLength() const { return sparams_.length; }
  virtual const bool CanOverlap() const { return sparams_.overlap; }
  virtual const int GetNInsert() const { return sparams_.num; }
  virtual const int GetNPosit() const { return sparams_.n_posit; }
  virtual const int GetNSpec() const { return sparams_.n_spec; }
  virtual const bool GetPositFlag() const { return sparams_.posit_flag; }
  virtual const bool GetSpecFlag() const { return sparams_.spec_flag; }
  virtual const std::string GetSpeciesName() const { return sparams_.name; }
  std::string GetInsertionType() const { return sparams_.insertion_type; }

  virtual void AddMember();
  virtual void AddMember(T newmem);
  virtual void PopMember();
  virtual void PopAll();
  // For any post-species insertion operations
  virtual void InitMembers() {}
  virtual void ArrangeMembers();
  virtual void CrystalArrangement();
  virtual void CenteredOrientedArrangement();
  void SetLastMemberPosition(double const *const pos);
  virtual void ResetPreviousPositions();
  virtual void Draw(std::vector<graph_struct *> &graph_array);
  virtual void UpdatePositions();
  virtual void GetInteractors(std::vector<Object *> &ix);
  virtual void GetLastInteractors(std::vector<Object *> &ix);
  virtual double GetPotentialEnergy();
  virtual void ZeroForces();
  virtual void Report();
  virtual int GetCount();
  virtual void WritePosits();
  virtual void WriteSpecs();
  virtual void WriteCheckpoints();
  virtual void ReadPosits();
  virtual void ReadPositsFromSpecs();
  virtual void ReadSpecs();
  virtual void ReadCheckpoints();
  virtual void ScalePositions();
  virtual const std::vector<T> &GetMembers() { return members_; }
  /* XXX: Dangerous. This need to grab members from species_base ptr. 
   *    If there is a better way, please let me know. <14-07-20, ARL> */
  virtual void GetMemberPtrs(std::vector<Object *> &objs) {
    for (auto &obj : members_) {
      objs.push_back(&obj);
    }
  }
  virtual void CleanUp();
  virtual void Reserve();
  virtual double const GetVolume();
  virtual double const GetDrMax();
  virtual void ZeroDrTot();
  virtual void CustomInsert();
  virtual const bool CheckInteractorUpdate();
  virtual void InitAnalysis() {
    LoadAnalysis();
    for (auto it = analysis_.begin(); it != analysis_.end(); ++it) {
      (*it)->Init(members_, sparams_);
    }
  }
  virtual void RunAnalysis() {
    for (auto it = analysis_.begin(); it != analysis_.end(); ++it) {
      (*it)->Run();
    }
  }
  virtual void FinalizeAnalysis() {
    for (auto it = analysis_.begin(); it != analysis_.end(); ++it) {
      (*it)->End();
      delete (*it);
    }
  }
  virtual bool CheckInteractionAnalysis() {
    for (auto it = analysis_.begin(); it != analysis_.end(); ++it) {
      if ((*it)->CheckInteractionAnalysis()) {
        return true;
      }
    }
    return false;
  }
};

template <typename T, unsigned char S>
double const Species<T, S>::GetVolume() {
  double vol = 0.0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    vol += it->GetVolume();
  }
  return vol;
}

template <typename T, unsigned char S>
void Species<T, S>::Reserve() {
  members_.reserve(GetNInsert());
  Logger::Debug("Reserving memory for %d members in %s %s", GetNInsert(),
                GetSID()._to_string(), GetSpeciesName().c_str());
}

template <typename T, unsigned char S>
void Species<T, S>::AddMember() {
  T newmember(rng_.GetSeed());
  Logger::Trace("Adding member to %s %s, member number %d, member id %d",
                GetSID()._to_string(), GetSpeciesName().c_str(), n_members_ + 1,
                newmember.GetOID());
  members_.push_back(newmember);
  members_.back().SetSID(GetSID());
  members_.back().Init(&sparams_);
  n_members_++;
}

template <typename T, unsigned char S>
double const Species<T, S>::GetDrMax() {
  double max_dr = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    double dr = it->GetDrTot();
    if (dr > max_dr) {
      max_dr = dr;
    }
  }
  // if (max_dr > params_->system_radius * params_->system_radius) {
  if (max_dr > 10000) {
    Logger::Warning("Oddly large dr (%2.2f) in %s %s", max_dr,
                    GetSID()._to_string(), GetSpeciesName().c_str());
    throw std::runtime_error("Unphysical result");
  }
  return max_dr;
}

template <typename T, unsigned char S>
void Species<T, S>::ZeroDrTot() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ZeroDrTot();
  }
}

template <typename T, unsigned char S>
void Species<T, S>::ResetPreviousPositions() {
  Logger::Trace("Resetting previous positions for species %s %s",
                GetSID()._to_string(), GetSpeciesName().c_str());
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ResetPreviousPosition();
  }
}

template <typename T, unsigned char S>
void Species<T, S>::AddMember(T newmem) {
  Logger::Trace("Adding preexisting member to %s %s", GetSID()._to_string(),
                GetSpeciesName().c_str());
  members_.push_back(newmem);
  newmem.SetSID(GetSID());
  n_members_++;
}

template <typename T, unsigned char S>
void Species<T, S>::PopMember() {
  Logger::Trace("Removing last member of %s %s", GetSID()._to_string(),
                GetSpeciesName().c_str());
  members_.back().Cleanup();
  members_.pop_back();
  n_members_--;
}

template <typename T, unsigned char S>
void Species<T, S>::PopAll() {
  while (n_members_ > 0) {
    PopMember();
  }
}

template <typename T, unsigned char S>
void Species<T, S>::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->Draw(graph_array);
  }
  for (auto it = analysis_.begin(); it != analysis_.end(); ++it) {
    (*it)->Draw(graph_array);
  }
}

template <typename T, unsigned char S>
void Species<T, S>::UpdatePositions() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->UpdatePosition();
  }
}

template <typename T, unsigned char S>
void Species<T, S>::GetInteractors(std::vector<Object *> &ixors) {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->GetInteractors(ixors);
  }
}

template <typename T, unsigned char S>
void Species<T, S>::GetLastInteractors(std::vector<Object *> &ix) {
  if (members_.size() == 0) {
    Logger::Error("Called for last interactors of %s %s, but species has zero "
                  "members",
                  GetSID()._to_string(), GetSpeciesName().c_str());
  }
  members_.back().GetInteractors(ix);
}

template <typename T, unsigned char S>
double Species<T, S>::GetPotentialEnergy() {
  double pe = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    pe += it->GetPotentialEnergy();
  }
  /* The total potential energy is going to be half of the
   potential energy felt by each particle. Potential energy
   is shared, so I need to avoid double counting. */
  return 0.5 * pe;
}

template <typename T, unsigned char S>
void Species<T, S>::ZeroForces() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->ZeroForce();
  }
}

template <typename T, unsigned char S>
void Species<T, S>::Report() {
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->Report();
  }
}

template <typename T, unsigned char S>
int Species<T, S>::GetCount() {
  int count = 0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    count += it->GetCount();
  }
  return count;
}

template <typename T, unsigned char S>
void Species<T, S>::WritePosits() {
  int size = members_.size();
  oposit_file_.write(reinterpret_cast<char *>(&size), sizeof(int));
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->WritePosit(oposit_file_);
}

template <typename T, unsigned char S>
void Species<T, S>::WriteSpecs() {
  int size = members_.size();
  ospec_file_.write(reinterpret_cast<char *>(&size), sizeof(int));
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->WriteSpec(ospec_file_);
}

template <typename T, unsigned char S>
void Species<T, S>::WriteCheckpoints() {
  Logger::Trace("Writing checkpoints for %s %s", GetSID()._to_string(),
                GetSpeciesName().c_str());
  std::fstream ocheck_file(checkpoint_file_, std::ios::out | std::ios::binary);
  if (!ocheck_file.is_open()) {
    Logger::Error("Output %s file did not open", checkpoint_file_.c_str());
  }
  void *rng_state = rng_.GetState();
  size_t rng_size = rng_.GetSize();
  int next_oid = Object::GetNextOID();
  int size = members_.size();
  ocheck_file.write(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  ocheck_file.write(reinterpret_cast<char *>(rng_state), rng_size);
  ocheck_file.write(reinterpret_cast<char *>(&next_oid), sizeof(int));
  ocheck_file.write(reinterpret_cast<char *>(&size), sizeof(int));
  if (size > 0) {
    for (auto it = members_.begin(); it != members_.end(); ++it)
      it->WriteCheckpoint(ocheck_file);
  }
  ocheck_file.close();
}

template <typename T, unsigned char S>
void Species<T, S>::ReadPosits() {
  if (iposit_file_.eof()) {
    Logger::Info("EOF reached while reading posits");
    early_exit = true;
    return;
  }
  if (!iposit_file_.is_open()) {
    Logger::Warning("ERROR Posit file unexpectedly not open! Exiting early.");
    early_exit = true;
    return;
  }
  int size = -1;
  iposit_file_.read(reinterpret_cast<char *>(&size), sizeof(int));
  // Hacky workaround FIXME
  // This prevents strange errors that occasionally crop up when reading inputs
  if (size == -1) {
    early_exit = true;
    return;
  }
  if (size != n_members_) {
    T member(rng_.GetSeed());
    members_.resize(size, member);
    n_members_ = size;
  }
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->ReadPosit(iposit_file_);
}

template <typename T, unsigned char S>
void Species<T, S>::ReadPositsFromSpecs() {
  ReadSpecs();
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->SetAvgPosition();
  }
}

template <typename T, unsigned char S>
void Species<T, S>::ReadCheckpoints() {
  Logger::Trace("Reading checkpoints for %s %s", GetSID()._to_string(),
                GetSpeciesName().c_str());
  std::fstream icheck_file(checkpoint_file_, std::ios::in | std::ios::binary);
  if (!icheck_file.is_open()) {
    Logger::Error("Output %s file did not open", checkpoint_file_.c_str());
  }
  void *rng_state = rng_.GetState();
  size_t rng_size;
  int next_oid = -1;
  int size = -1;
  icheck_file.read(reinterpret_cast<char *>(&rng_size), sizeof(size_t));
  icheck_file.read(reinterpret_cast<char *>(rng_state), rng_size);
  icheck_file.read(reinterpret_cast<char *>(&next_oid), sizeof(int));
  icheck_file.read(reinterpret_cast<char *>(&size), sizeof(int));
  if (size == 0) {
    members_.clear();
  } else {
    if (members_.size() == 0) {
      AddMember();
    }
    members_.resize(size, members_[0]);
    n_members_ = size;
    for (auto it = members_.begin(); it != members_.end(); ++it)
      it->ReadCheckpoint(icheck_file);
  }
  icheck_file.close();
  Object::SetNextOID(next_oid);
}

template <typename T, unsigned char S>
void Species<T, S>::ReadSpecs() {
  if (ispec_file_.eof()) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (!ispec_file_.is_open()) {
    Logger::Warning("ERROR: Spec file unexpectedly not open! Exiting early.");
    early_exit = true;
    return;
  }
  n_members_ = -1;
  ispec_file_.read(reinterpret_cast<char *>(&n_members_), sizeof(int));
  /* For some reason, we can't catch the EOF above. If size == -1 still, then
     we caught a EOF here */
  if (n_members_ == -1) {
    if (HandleEOF()) {
      return;
    } else {
      Logger::Info("EOF reached in spec file for %s %s", GetSID()._to_string(),
                   GetSpeciesName().c_str());
      early_exit = true;
      return;
    }
  }
  if (n_members_ != members_.size()) {
    T member(rng_.GetSeed());
    member.Init(&sparams_);
    member.SetSID(GetSID());
    members_.resize(n_members_, member);
  }
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->ReadSpec(ispec_file_);
}

template <typename T, unsigned char S>
void Species<T, S>::ScalePositions() {
  for (auto it = members_.begin(); it != members_.end(); ++it)
    it->ScalePosition();
}

template <typename T, unsigned char S>
void Species<T, S>::CleanUp() {
  members_.clear();
}

template <typename T, unsigned char S>
void Species<T, S>::ArrangeMembers() {
  if (GetInsertionType().compare("custom") == 0)
    CustomInsert();
  else if (GetInsertionType().compare("simple_crystal") == 0)
    CrystalArrangement();
  else if (GetInsertionType().compare("centered_oriented") == 0)
    CenteredOrientedArrangement();
  else
    Logger::Error(
        "Arrangement not recognized and ArrangeMembers not overwritten by "
        "species!");
}

template <typename T, unsigned char S>
void Species<T, S>::SetLastMemberPosition(double const *const pos) {
  members_.back().SetPosition(pos);
  members_.back().UpdatePeriodic();
}

template <typename T, unsigned char S>
void Species<T, S>::CenteredOrientedArrangement() {
  // This is redundant for filaments, since they have already inserted
  // themselves properly.
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  u[params_->n_dim - 1] = 1.0;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    it->SetPosition(pos);
    it->SetOrientation(u);
  }
}

template <typename T, unsigned char S>
void Species<T, S>::CrystalArrangement() {
  double d = members_[0].GetDiameter();
  double l = members_[0].GetLength();
  int n_dim = params_->n_dim;
  double R = space_->radius;
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  u[n_dim - 1] = 1;
  double nX_max = floor((2.0 * R / d) - 1);
  double nY_max = floor((2.0 * R - d) / (l + d));
  double nZ_max = (n_dim == 3 ? nX_max : 1);
  double N = nX_max * nY_max * nZ_max;
  if (n_members_ > N) {
    Logger::Error("Number of members in species exceeds maximum possible for "
                  "crystal in "
                  "system radius! Max possible: %d",
                  (int)N);
  }
  double fraction = N / n_members_;
  if (fraction < 1)
    fraction = 1;
  if (l == 0) {
    if (n_dim == 2)
      nY_max = floor(pow(n_members_, 1.0 / n_dim));
    else if (n_dim == 3)
      nY_max = ceil(pow(n_members_, 1.0 / n_dim));
  }
  if (n_dim == 2) {
    nX_max = ceil(n_members_ / nY_max);
  } else if (n_dim == 3) {
    nX_max = floor(sqrt(n_members_ / nY_max));
    if (nX_max == 0)
      nX_max = 1;
    nZ_max = ceil(n_members_ / (nX_max * nY_max));
  }
  for (int i = 0; i < n_dim; ++i) {
    pos[i] = -R + 0.5 * d;
  }
  pos[n_dim - 1] += 0.5 * l;
  int inserted = 0;
  int insert_x = 0;
  int insert_y = 0;
  int insert_z = 0;
  bool shift = false;
  double diff_y = (2 * R - nY_max * (l + d)) / nY_max;
  double diff_x = (2 * R - nX_max * d) / nX_max;
  double diff_z = (2 * R - nZ_max * d) / nZ_max;
  int u0 = 1;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    // Check crystal orientation type
    if (params_->uniform_crystal == 0) {
      // Random orientations
      u[n_dim - 1] = (rng_.RandomInt(2) == 0 ? 1 : -1);
    } else if (params_->uniform_crystal == 2) {
      // Stagger orientations
      u[n_dim - 1] = -u[n_dim - 1];
    }
    // Insert object
    it->InsertAt(pos, u);
    inserted++;
    // Update next position
    pos[n_dim - 1] += l + d + diff_y;
    if (++insert_y == nY_max) {
      if (params_->uniform_crystal == 2 && shift) {
        // Stagger orientations row-wise
        u[n_dim - 1] = u0;
        u0 = -u0;
      }
      // Stagger positions column-wise
      shift = !shift;
      insert_y = 0;
      pos[n_dim - 1] =
          -R + 0.5 * (l + d) + (shift ? 0.5 * (l + d + diff_y) : 0);
      pos[0] += d + diff_x;
      if (++insert_x == nX_max) {
        if (inserted < n_members_ && n_dim == 2) {
          Logger::Error(
              "Ran out of room while arranging crystal in 2D! Arranged %d/%d",
              inserted, n_members_);
        } else if (n_dim == 2)
          continue;
        insert_x = 0;
        pos[0] = -R + 0.5 * d;
        pos[1] += d + diff_z;
        if (++insert_z == nZ_max && inserted < n_members_) {
          Logger::Error(
              "Ran out of room while arranging crystal in 3D! Arranged %d/%d",
              inserted, n_members_);
        }
      }
    }
  }
}

template <typename T, unsigned char S>
void Species<T, S>::CustomInsert() {
  YAML::Node inode;
  std::string spec_name = GetSID()._to_string();
  try {
    inode = YAML::LoadFile(sparams_.insert_file);
  } catch (...) {
    std::cout << "Failed to load custom insert file " << sparams_.insert_file
              << " for species " << spec_name << "\n";
    Logger::Error("");
  }
  if (!inode[spec_name] || inode[spec_name].size() != n_members_) {
    std::cout << "Custom insert file for species " << spec_name
              << " was invalid: \n";
    if (!inode[spec_name]) {
      std::cout << "  Species ID header was missing from file\n";
    } else if (inode[spec_name].size() != n_members_) {
      std::cout << "  There were " << inode[spec_name].size()
                << " positions specified for " << n_members_
                << " species members\n";
    }
    Logger::Error("");
  }
  for (YAML::const_iterator it = inode.begin(); it != inode.end(); ++it) {
    if (it->first.as<std::string>().compare(spec_name) != 0)
      continue;
    if (!it->second.IsSequence()) {
      Logger::Error("Custom insert file positions not specified as sequence");
    }
    int i_member = 0;
    for (YAML::const_iterator jt = it->second.begin(); jt != it->second.end();
         ++jt) {
      double pos[3], u[3];
      if (!jt->IsSequence()) {
        Logger::Error("Custom insert position not specified as sequence");
      }
      if (jt->size() != 2 || (*jt)[0].size() != 3 || (*jt)[1].size() != 3) {
        Logger::Error("Custom insert position not in format "
                      "[[pos_x,pos_y,pos_z],[u_x,u_y,u_z]]");
      }
      for (int i = 0; i < 3; ++i) {
        pos[i] = (*jt)[0][i].as<double>();
        u[i] = (*jt)[1][i].as<double>();
      }
      members_[i_member++].InsertAt(pos, u);
    }
  }
}

template <typename T, unsigned char S>
const bool Species<T, S>::CheckInteractorUpdate() {
  bool result = false;
  for (auto it = members_.begin(); it != members_.end(); ++it) {
    if (it->CheckInteractorUpdate()) {
      result = true;
    }
  }
  return result;
}

#endif // _CGLASS_SPECIES_H_
