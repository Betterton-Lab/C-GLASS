#ifndef _CGLASS_PARAMS_PARSER_H_
#define _CGLASS_PARAMS_PARSER_H_

#include "auxiliary.hpp"
#include "yaml-cpp/yaml.h"

typedef std::pair<std::string, std::string> sid_label;

class ParamsParser {
private:
  YAML::Node sim_node_;
  YAML::Node species_node_;
  std::vector<sid_label> labels_;
  void CheckDuplicateLabels();
  void ParseSpeciesParameters();
  std::vector<species_base_parameters *> species_params_;
  int n_species_ = 0;
  int n_crosslinks_ = 0;
  int n_optical_traps_ = 0;

public:
  void Init(YAML::Node sim_params);
  system_parameters ParseSystemParameters();
  std::vector<sid_label> &GetSpeciesLabels() { return labels_; }
  species_base_parameters *GetNewSpeciesParameters(species_id sid,
                                                   std::string spec_name);
  const int GetNCrosslinkSpecies() { return n_crosslinks_; }
  const int GetNOpticalTrapSpecies() { return n_optical_traps_; }
  const int GetNSpecies() { return n_species_; }
};

#endif
