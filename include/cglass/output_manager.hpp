#ifndef _CGLASS_OUTPUT_MANAGER_H_
#define _CGLASS_OUTPUT_MANAGER_H_

#include "parse_flags.hpp"
#include "species.hpp"

template <class T>
class OutputManagerBase {
private:
  bool posit_flag_ = false;
  bool spec_flag_ = false;
  bool checkpoint_flag_ = false;
  bool thermo_flag_ = false;
  bool posits_only_ = false;
  bool reduce_flag_ = false;
  bool thermo_analysis_ = false;
  bool with_reloads_ = false;
  int n_posit_;
  int n_spec_;
  int n_checkpoint_;
  int n_thermo_;
  int reduce_factor_ = 1;
  std::string run_name_;
  system_parameters *params_;
  std::vector<T *> *species_;
  space_struct *space_;
  void WritePosits();
  void WriteSpecs();
  void WriteCheckpoints();
  virtual void WriteThermo();
  virtual void ReadThermo();
  virtual void InitThermo(std::string fname);
  void WriteReduce();
  void WriteTime();
  std::fstream othermo_file_;
  std::fstream ithermo_file_;
  std::fstream time_file_;

public:
  OutputManagerBase() {}
  void Init(system_parameters *params, std::vector<T *> *species,
            space_struct *space, bool reading_inputs = false,
            run_options *run_opts = nullptr);
  int GetNPosit() { return n_posit_; }
  int GetNSpec() { return n_spec_; }
  int GetNCheckpoint() { return n_checkpoint_; }
  void WriteOutputs();
  void InitInputs();
  void ReadInputs();
  void Close();
};

template <class T>
void OutputManagerBase<T>::Init(system_parameters *params,
                                std::vector<T *> *species, space_struct *space,
                                bool reading_inputs, run_options *run_opts) {
  params_ = params;
  run_name_ = params_->run_name;
  species_ = species;
  space_ = space;
  n_posit_ = n_spec_ = n_checkpoint_ = ABS((int)params_->n_steps);
  n_thermo_ = params_->n_thermo;
  thermo_flag_ = params_->thermo_flag;

  if (run_opts) {
    reduce_flag_ = run_opts->reduce_flag;
    reduce_factor_ = run_opts->reduce_factor;
    with_reloads_ = run_opts->with_reloads;
    posits_only_ = run_opts->use_posits;
  }
  if (thermo_flag_ && !reading_inputs && !reduce_flag_) {
    InitThermo(run_name_);
  } else if (thermo_flag_ && !reading_inputs && reduce_flag_) {
    std::string red_file_name =
        run_name_ + "_reduced" + std::to_string(reduce_factor_);
    InitThermo(red_file_name);
  }

  for (auto it = species_->begin(); it != species_->end(); ++it) {
    if (!params_->load_checkpoint && reading_inputs) {
      (*it)->InitInputFiles(run_name_, posits_only_, with_reloads_);
      if (reduce_flag_) {
        std::string red_file_name =
            run_name_ + "_reduced" + std::to_string(reduce_factor_);
        (*it)->InitOutputFiles(red_file_name);
      } else if (params->checkpoint_from_spec) {
        (*it)->InitCheckpoints(run_name_);
      }
    } else {
      (*it)->InitOutputFiles(run_name_);
    }
    if ((*it)->GetPositFlag()) {
      posit_flag_ = true;
      if ((*it)->GetNPosit() < n_posit_) {
        n_posit_ = (*it)->GetNPosit();
      }
    }
    if ((*it)->GetSpecFlag()) {
      spec_flag_ = true;
      if ((*it)->GetNSpec() < n_spec_) {
        n_spec_ = (*it)->GetNSpec();
      }
    }
    if ((*it)->GetCheckpointFlag()) {
      checkpoint_flag_ = true;
      if ((*it)->GetNCheckpoint() < n_checkpoint_) {
        n_checkpoint_ = (*it)->GetNCheckpoint();
      }
    }
  }
}

template <class T>
void OutputManagerBase<T>::WriteOutputs() {
  if (posit_flag_ && (params_->i_step % n_posit_ == 0)) {
    WritePosits();
  }
  if (spec_flag_ && params_->i_step != params_->prev_step &&
      (params_->i_step % n_spec_ == 0)) {
    Logger::Debug("Writing spec files");
    WriteSpecs();
  }
  if (checkpoint_flag_ && params_->i_step != params_->prev_step &&
      (params_->i_step % n_checkpoint_ == 0)) {
    Logger::Debug("Writing checkpoint files");
    WriteCheckpoints();
  }
  if (thermo_flag_ && params_->i_step != params_->prev_step &&
      (params_->i_step % n_thermo_ == 0)) {
    Logger::Debug("Writing thermo file");
    WriteThermo();
  }
}

template <class T>
void OutputManagerBase<T>::WritePosits() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetPositFlag() && params_->i_step != params_->prev_step &&
        params_->i_step % (*spec)->GetNPosit() == 0) {
      (*spec)->WritePosits();
    }
  }
}

template <class T>
void OutputManagerBase<T>::WriteSpecs() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetSpecFlag() && params_->i_step != params_->prev_step &&
        params_->i_step % (*spec)->GetNSpec() == 0) {
      (*spec)->WriteSpecs();
    }
  }
}

template <class T>
void OutputManagerBase<T>::WriteCheckpoints() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetCheckpointFlag() && params_->i_step != params_->prev_step &&
        params_->i_step % (*spec)->GetNCheckpoint() == 0) {
      (*spec)->WriteCheckpoints();
    }
  }
}

template <class T>
void OutputManagerBase<T>::InitThermo(std::string fname) {
  fname.append(".thermo");
  othermo_file_.open(fname, std::ios::out | std::ios::binary);
  if (!othermo_file_.is_open()) {
    Logger::Error("Thermo file failed to open!");
  }
  othermo_file_.write(reinterpret_cast<char *>(&(params_->n_steps)),
                      sizeof(int));
  othermo_file_.write(reinterpret_cast<char *>(&(n_thermo_)), sizeof(int));
  othermo_file_.write(reinterpret_cast<char *>(&(params_->delta)),
                      sizeof(double));
  othermo_file_.write(reinterpret_cast<char *>(&(params_->n_dim)), sizeof(int));
}

template <class T>
void OutputManagerBase<T>::WriteThermo() {
  for (int i = 0; i < 9; ++i) {
    othermo_file_.write(reinterpret_cast<char *>(&(space_->unit_cell[i])),
                        sizeof(double));
  }
  for (int i = 0; i < 9; ++i) {
    othermo_file_.write(reinterpret_cast<char *>(&(space_->pressure_tensor[i])),
                        sizeof(double));
  }
  othermo_file_.write(reinterpret_cast<char *>(&(space_->pressure)),
                      sizeof(double));
  othermo_file_.write(reinterpret_cast<char *>(&(space_->volume)),
                      sizeof(double));
}

template <class T>
void OutputManagerBase<T>::ReadThermo() {
  for (int i = 0; i < 9; ++i) {
    ithermo_file_.read(reinterpret_cast<char *>(&(space_->unit_cell[i])),
                       sizeof(double));
  }
  for (int i = 0; i < 9; ++i) {
    ithermo_file_.read(reinterpret_cast<char *>(&(space_->pressure_tensor[i])),
                       sizeof(double));
  }
  ithermo_file_.read(reinterpret_cast<char *>(&(space_->pressure)),
                     sizeof(double));
  ithermo_file_.read(reinterpret_cast<char *>(&(space_->volume)),
                     sizeof(double));
}

template <class T>
void OutputManagerBase<T>::ReadInputs() {
  // if ( params_->i_step % n_posit_ != 0 && params_->i_step % n_spec_ != 0) {
  // return;
  //}
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if (posits_only_ && (*spec)->GetPositFlag() &&
        params_->i_step != params_->prev_step &&
        params_->i_step % (*spec)->GetNPosit() == 0) {
      (*spec)->ReadPosits();
    }
    // In the case that we only want to consider posits, but we only have spec
    // files, use specs to read average positions
    else if (posits_only_ && !(*spec)->GetPositFlag() &&
             (*spec)->GetSpecFlag() && params_->i_step != params_->prev_step &&
             params_->i_step % (*spec)->GetNSpec() == 0) {
      (*spec)->ReadPositsFromSpecs();
    } else if (!posits_only_ && (*spec)->GetSpecFlag() &&
               params_->i_step != params_->prev_step &&
               params_->i_step % (*spec)->GetNSpec() == 0) {
      (*spec)->ReadSpecs();
      if (params_->checkpoint_from_spec &&
          params_->i_step % (*spec)->GetNCheckpoint() == 0 &&
          params_->i_step != params_->prev_step) {
        (*spec)->WriteCheckpoints();
      }
    }
  }
  if (thermo_flag_) {
    ReadThermo();
  }
  if (reduce_flag_) {
    WriteReduce();
  }
}

template <class T>
void OutputManagerBase<T>::WriteReduce() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    if ((*spec)->GetPositFlag() && params_->i_step != params_->prev_step &&
        params_->i_step % (reduce_factor_ * (*spec)->GetNPosit()) == 0) {
      (*spec)->WritePosits();
    }
    if (!posits_only_ && (*spec)->GetSpecFlag() &&
        params_->i_step != params_->prev_step &&
        params_->i_step % (reduce_factor_ * (*spec)->GetNSpec()) == 0) {
      (*spec)->WriteSpecs();
    }
  }
  if (thermo_flag_ && params_->i_step != params_->prev_step &&
      params_->i_step % (reduce_factor_ * n_thermo_) == 0) {
    WriteThermo();
  }
}

template <class T>
void OutputManagerBase<T>::Close() {
  for (auto spec = species_->begin(); spec != species_->end(); ++spec) {
    (*spec)->CloseFiles();
  }
  if (ithermo_file_.is_open()) {
    ithermo_file_.close();
  }
  if (othermo_file_.is_open()) {
    othermo_file_.close();
  }
  if (params_->time_analysis) {
    WriteTime();
  }
}

template <class T>
void OutputManagerBase<T>::WriteTime() {
  std::string fname = run_name_;
  fname.append(".time");
  time_file_.open(fname, std::ios::out);
  time_file_ << "i_step n_steps delta\n";
  time_file_ << params_->i_step << " " << params_->n_steps << " "
             << params_->delta << "\n";
  time_file_.close();
}

class OutputManager : public OutputManagerBase<SpeciesBase> {};

#endif //_CGLASS_OUTPUT_MANAGER_H_
