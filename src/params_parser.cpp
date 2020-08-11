#include "cglass/params_parser.hpp"
#include "cglass/parse_params.hpp"

void ParamsParser::Init(YAML::Node sim_params) {
  Logger::Trace("Initializing yaml parameter parser");
  sim_node_ = sim_params;
  ParseSpeciesParameters();
}

void ParamsParser::CheckDuplicateLabels() {
  Logger::Trace("Checking for duplicate species labels");
  int n_labels = labels_.size();
  for (int i = 0; i < n_labels - 1; ++i) {
    for (int j = i + 1; j < n_labels; ++j) {
      if (labels_[i] == labels_[j]) {
        Logger::Error("Duplicate species label found! %s %s",
                      labels_[i].first.c_str(), labels_[i].second.c_str());
      }
    }
  }
  for (int i = 0; i < n_labels; ++i) {
    if (species_id::_from_string(labels_[i].first.c_str()) ==
        +species_id::crosslink) {
      n_crosslinks_++;
    } else if (species_id::_from_string(labels_[i].first.c_str()) ==
               +species_id::optical_trap) {
      n_optical_traps_++;
    } else {
      n_species_++;
    }
  }
}

void ParamsParser::ParseSpeciesParameters() {
  Logger::Trace("Parsing species parameters");
  for (auto it = sim_node_.begin(); it != sim_node_.end(); ++it) {
    if (it->second.IsMap()) {
      if (it->first.as<std::string>() == "species") {
        continue;
      }
      species_node_[it->first.as<std::string>()] = it->second;
      labels_.push_back(std::make_pair(it->first.as<std::string>(), "species"));
    } else if (it->second.IsSequence()) {
      species_node_[it->first.as<std::string>()] = it->second;
      for (int i = 0; i < it->second.size(); ++i) {
        YAML::Node subnode = it->second[i];
        for (YAML::const_iterator jt = subnode.begin(); jt != subnode.end();
             ++jt) {
          if (jt->first.as<std::string>().compare("name") == 0) {
            labels_.push_back(std::make_pair(it->first.as<std::string>(),
                                             jt->second.as<std::string>()));
          }
        }
      }
    }
  }
  CheckDuplicateLabels();
}

system_parameters ParamsParser::ParseSystemParameters() {
  Logger::Trace("Building system parameters");
  return parse_system_params(sim_node_);
}

species_base_parameters *
ParamsParser::GetNewSpeciesParameters(species_id sid, std::string spec_name) {
  std::string sid_str(sid._to_string());
  Logger::Trace("Building parameters for %s %s", sid_str.c_str(),
                spec_name.c_str());
  for (auto it = species_node_.begin(); it != species_node_.end(); ++it) {
    if (it->first.as<std::string>().compare(sid_str) == 0) {
      if (it->second.IsMap() && spec_name.compare("species") == 0) {
        return parse_species_params(sid_str, it->second, sim_node_);
      } else if (it->second.IsSequence()) {
        for (int i = 0; i < it->second.size(); ++i) {
          YAML::Node subnode = it->second[i];
          if (subnode["name"] &&
              subnode["name"].as<std::string>().compare(spec_name) == 0) {
            return parse_species_params(sid_str, subnode, sim_node_);
          }
        }
      }
    }
  }
  Logger::Error("ParamsParser did not find species id/label combination %s/%s"
                " in species params YAML node!",
                sid_str.c_str(), spec_name.c_str());
  return nullptr;
}
