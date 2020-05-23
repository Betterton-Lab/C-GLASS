#ifndef _CGLASS_SIMULATION_MANAGER_H_
#define _CGLASS_SIMULATION_MANAGER_H_

#include "simulation.hpp"

class SimulationManager {
 private:
  UNIT_TEST
  bool make_movie_ = false;
  bool run_analysis_ = false;
  int n_runs_ = 1;               // Number of runs per parameters set
  int n_var_ = 1;                // Number of parameter variations
  int n_random_ = 1;             // Number of random params
  std::string run_name_ = "sc";  // simulation batch name
  std::vector<std::string> pfiles_;
  Simulation *sim_;   // New sim created and destroyed for every set of
                      // parameters
  YAML::Node pnode_;  // Main node to initialize pvector
  std::vector<YAML::Node> pvector_;
  system_parameters params_;
  run_options run_opts_;
  RNG *rng_ = nullptr;  // ptr so we can instantiate it after SetSeed
  void CheckAppendParams();
  void AppendParams(YAML::Node app_node);
  void LoadDefaultParams();
  void CountVariations();
  void SetRandomParams();
  double GetRandomParam(std::string rtype, double min, double max);
  void CheckRandomParams();
  void GenerateParameters();
  void WriteParams();
  void RunSimulations();
  // void ParseParams(std::string file_name);
  void ProcessOutputs();
  void InitLogger();

 public:
  SimulationManager() { early_exit = false; }
  ~SimulationManager() {
    /* rng_ is a ptr here because I have to instantiate the RNG after
       initializing the static seed */
    if (rng_ != nullptr) {
      delete rng_;
    }
  }
  void InitManager(run_options run_opts);
  void RunManager();
};

#endif  // _CGLASS_SIMULATION_MANAGER_H_
