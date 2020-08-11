#ifndef _CGLASS_PARSE_PARAMS_H_
#define _CGLASS_PARSE_PARAMS_H_

#include "yaml-cpp/yaml.h"
#include "auxiliary.hpp"

system_parameters parse_system_params(YAML::Node &node) {
  system_parameters params;
  for (auto it=node.begin(); it!= node.end(); ++it) {
    if (!it->second.IsScalar()) {
     continue;
    }
    std::string param_name = it->first.as<std::string>();
    if (false) {
    } else if (param_name.compare("seed")==0) {
    params.seed = it->second.as<long>();
    } else if (param_name.compare("n_runs")==0) {
    params.n_runs = it->second.as<int>();
    } else if (param_name.compare("n_random")==0) {
    params.n_random = it->second.as<int>();
    } else if (param_name.compare("run_name")==0) {
    params.run_name = it->second.as<std::string>();
    } else if (param_name.compare("n_dim")==0) {
    params.n_dim = it->second.as<int>();
    } else if (param_name.compare("n_periodic")==0) {
    params.n_periodic = it->second.as<int>();
    } else if (param_name.compare("boundary")==0) {
    params.boundary = it->second.as<int>();
    } else if (param_name.compare("system_radius")==0) {
    params.system_radius = it->second.as<double>();
    } else if (param_name.compare("n_steps")==0) {
    params.n_steps = it->second.as<int>();
    } else if (param_name.compare("i_step")==0) {
    params.i_step = it->second.as<int>();
    } else if (param_name.compare("prev_step")==0) {
    params.prev_step = it->second.as<int>();
    } else if (param_name.compare("delta")==0) {
    params.delta = it->second.as<double>();
    } else if (param_name.compare("dynamic_timestep")==0) {
    params.dynamic_timestep = it->second.as<bool>();
    } else if (param_name.compare("dynamic_timestep_ramp")==0) {
    params.dynamic_timestep_ramp = it->second.as<double>();
    } else if (param_name.compare("graph_flag")==0) {
    params.graph_flag = it->second.as<bool>();
    } else if (param_name.compare("n_graph")==0) {
    params.n_graph = it->second.as<int>();
    } else if (param_name.compare("graph_diameter")==0) {
    params.graph_diameter = it->second.as<double>();
    } else if (param_name.compare("object_opacity")==0) {
    params.object_opacity = it->second.as<double>();
    } else if (param_name.compare("invert_background")==0) {
    params.invert_background = it->second.as<bool>();
    } else if (param_name.compare("draw_boundary")==0) {
    params.draw_boundary = it->second.as<bool>();
    } else if (param_name.compare("load_checkpoint")==0) {
    params.load_checkpoint = it->second.as<bool>();
    } else if (param_name.compare("checkpoint_run_name")==0) {
    params.checkpoint_run_name = it->second.as<std::string>();
    } else if (param_name.compare("n_load")==0) {
    params.n_load = it->second.as<int>();
    } else if (param_name.compare("movie_flag")==0) {
    params.movie_flag = it->second.as<bool>();
    } else if (param_name.compare("movie_directory")==0) {
    params.movie_directory = it->second.as<std::string>();
    } else if (param_name.compare("time_analysis")==0) {
    params.time_analysis = it->second.as<bool>();
    } else if (param_name.compare("bud_height")==0) {
    params.bud_height = it->second.as<double>();
    } else if (param_name.compare("bud_radius")==0) {
    params.bud_radius = it->second.as<double>();
    } else if (param_name.compare("lj_epsilon")==0) {
    params.lj_epsilon = it->second.as<double>();
    } else if (param_name.compare("wca_eps")==0) {
    params.wca_eps = it->second.as<double>();
    } else if (param_name.compare("wca_sig")==0) {
    params.wca_sig = it->second.as<double>();
    } else if (param_name.compare("ss_a")==0) {
    params.ss_a = it->second.as<double>();
    } else if (param_name.compare("ss_rs")==0) {
    params.ss_rs = it->second.as<double>();
    } else if (param_name.compare("ss_eps")==0) {
    params.ss_eps = it->second.as<double>();
    } else if (param_name.compare("f_cutoff")==0) {
    params.f_cutoff = it->second.as<double>();
    } else if (param_name.compare("constant_pressure")==0) {
    params.constant_pressure = it->second.as<bool>();
    } else if (param_name.compare("constant_volume")==0) {
    params.constant_volume = it->second.as<bool>();
    } else if (param_name.compare("target_pressure")==0) {
    params.target_pressure = it->second.as<double>();
    } else if (param_name.compare("target_radius")==0) {
    params.target_radius = it->second.as<double>();
    } else if (param_name.compare("pressure_time")==0) {
    params.pressure_time = it->second.as<int>();
    } else if (param_name.compare("compressibility")==0) {
    params.compressibility = it->second.as<double>();
    } else if (param_name.compare("zero_temperature")==0) {
    params.zero_temperature = it->second.as<bool>();
    } else if (param_name.compare("thermo_flag")==0) {
    params.thermo_flag = it->second.as<bool>();
    } else if (param_name.compare("n_thermo")==0) {
    params.n_thermo = it->second.as<int>();
    } else if (param_name.compare("insert_radius")==0) {
    params.insert_radius = it->second.as<double>();
    } else if (param_name.compare("interaction_flag")==0) {
    params.interaction_flag = it->second.as<bool>();
    } else if (param_name.compare("remove_duplicate_interactions")==0) {
    params.remove_duplicate_interactions = it->second.as<bool>();
    } else if (param_name.compare("coarse_grained_mesh_interactions")==0) {
    params.coarse_grained_mesh_interactions = it->second.as<bool>();
    } else if (param_name.compare("mesh_coarsening")==0) {
    params.mesh_coarsening = it->second.as<int>();
    } else if (param_name.compare("species_insertion_failure_threshold")==0) {
    params.species_insertion_failure_threshold = it->second.as<int>();
    } else if (param_name.compare("species_insertion_reattempt_threshold")==0) {
    params.species_insertion_reattempt_threshold = it->second.as<int>();
    } else if (param_name.compare("uniform_crystal")==0) {
    params.uniform_crystal = it->second.as<bool>();
    } else if (param_name.compare("n_steps_equil")==0) {
    params.n_steps_equil = it->second.as<int>();
    } else if (param_name.compare("n_steps_target")==0) {
    params.n_steps_target = it->second.as<int>();
    } else if (param_name.compare("static_particle_number")==0) {
    params.static_particle_number = it->second.as<bool>();
    } else if (param_name.compare("checkpoint_from_spec")==0) {
    params.checkpoint_from_spec = it->second.as<bool>();
    } else if (param_name.compare("potential")==0) {
    params.potential = it->second.as<std::string>();
    } else if (param_name.compare("soft_potential_mag")==0) {
    params.soft_potential_mag = it->second.as<double>();
    } else if (param_name.compare("soft_potential_mag_target")==0) {
    params.soft_potential_mag_target = it->second.as<double>();
    } else if (param_name.compare("like_like_interactions")==0) {
    params.like_like_interactions = it->second.as<bool>();
    } else if (param_name.compare("auto_graph")==0) {
    params.auto_graph = it->second.as<bool>();
    } else if (param_name.compare("local_order_analysis")==0) {
    params.local_order_analysis = it->second.as<bool>();
    } else if (param_name.compare("local_order_width")==0) {
    params.local_order_width = it->second.as<double>();
    } else if (param_name.compare("local_order_bin_width")==0) {
    params.local_order_bin_width = it->second.as<double>();
    } else if (param_name.compare("local_order_n_analysis")==0) {
    params.local_order_n_analysis = it->second.as<int>();
    } else if (param_name.compare("density_analysis")==0) {
    params.density_analysis = it->second.as<int>();
    } else if (param_name.compare("density_bin_width")==0) {
    params.density_bin_width = it->second.as<double>();
    } else if (param_name.compare("density_com_only")==0) {
    params.density_com_only = it->second.as<bool>();
    } else if (param_name.compare("overlap_analysis")==0) {
    params.overlap_analysis = it->second.as<bool>();
    } else if (param_name.compare("highlight_overlaps")==0) {
    params.highlight_overlaps = it->second.as<bool>();
    } else if (param_name.compare("reduced")==0) {
    params.reduced = it->second.as<bool>();
    } else if (param_name.compare("reload_reduce_switch")==0) {
    params.reload_reduce_switch = it->second.as<bool>();
    } else if (param_name.compare("checkpoint_flag")==0) {
    params.checkpoint_flag = it->second.as<bool>();
    } else if (param_name.compare("n_checkpoint")==0) {
    params.n_checkpoint = it->second.as<int>();
    } else if (param_name.compare("no_midstep")==0) {
    params.no_midstep = it->second.as<bool>();
    } else {
      Logger::Warning("Unrecognized parameter '%s'", param_name.c_str());
    }
  }
  return params;
}

void parse_species_base_params(species_base_parameters &params,
                               YAML::Node &node) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    if (it->first.as<std::string>().compare("species") == 0) {
      if (!it->second.IsMap()) {
        Logger::Error("Species base params yaml node should always be a map!");
      }
      for (auto jt = it->second.begin(); jt != it->second.end(); ++jt) {
        std::string param_name = jt->first.as<std::string>();
        if (false) {
        } else if (param_name.compare("name")==0) {
        params.name = jt->second.as<std::string>();
        } else if (param_name.compare("num")==0) {
        params.num = jt->second.as<int>();
        } else if (param_name.compare("diameter")==0) {
        params.diameter = jt->second.as<double>();
        } else if (param_name.compare("length")==0) {
        params.length = jt->second.as<double>();
        } else if (param_name.compare("insertion_type")==0) {
        params.insertion_type = jt->second.as<std::string>();
        } else if (param_name.compare("insert_file")==0) {
        params.insert_file = jt->second.as<std::string>();
        } else if (param_name.compare("overlap")==0) {
        params.overlap = jt->second.as<bool>();
        } else if (param_name.compare("draw_type")==0) {
        params.draw_type = jt->second.as<std::string>();
        } else if (param_name.compare("color")==0) {
        params.color = jt->second.as<double>();
        } else if (param_name.compare("posit_flag")==0) {
        params.posit_flag = jt->second.as<bool>();
        } else if (param_name.compare("spec_flag")==0) {
        params.spec_flag = jt->second.as<bool>();
        } else if (param_name.compare("n_posit")==0) {
        params.n_posit = jt->second.as<int>();
        } else if (param_name.compare("n_spec")==0) {
        params.n_spec = jt->second.as<int>();
        } else {
          Logger::Warning("Species base parameter %s not recognized!", param_name.c_str());
        }
      }
      return;
    }
  }
}

species_base_parameters *parse_species_params(std::string sid,
                                              YAML::Node &subnode,
                                              YAML::Node &node) {
  if (false) {
  } else if (sid.compare("rigid_filament") == 0) {
    rigid_filament_parameters params;
    parse_species_base_params(params, node);
    for (auto jt = subnode.begin(); jt != subnode.end(); ++jt) {
      std::string param_name = jt->first.as<std::string>();
      if (false) {
      } else if (param_name.compare("name")==0) {
      params.name = jt->second.as<std::string>();
      } else if (param_name.compare("num")==0) {
      params.num = jt->second.as<int>();
      } else if (param_name.compare("diameter")==0) {
      params.diameter = jt->second.as<double>();
      } else if (param_name.compare("length")==0) {
      params.length = jt->second.as<double>();
      } else if (param_name.compare("insertion_type")==0) {
      params.insertion_type = jt->second.as<std::string>();
      } else if (param_name.compare("insert_file")==0) {
      params.insert_file = jt->second.as<std::string>();
      } else if (param_name.compare("overlap")==0) {
      params.overlap = jt->second.as<bool>();
      } else if (param_name.compare("draw_type")==0) {
      params.draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("color")==0) {
      params.color = jt->second.as<double>();
      } else if (param_name.compare("posit_flag")==0) {
      params.posit_flag = jt->second.as<bool>();
      } else if (param_name.compare("spec_flag")==0) {
      params.spec_flag = jt->second.as<bool>();
      } else if (param_name.compare("n_posit")==0) {
      params.n_posit = jt->second.as<int>();
      } else if (param_name.compare("n_spec")==0) {
      params.n_spec = jt->second.as<int>();
      } else if (param_name.compare("max_length")==0) {
      params.max_length = jt->second.as<double>();
      } else if (param_name.compare("min_length")==0) {
      params.min_length = jt->second.as<double>();
      } else if (param_name.compare("stationary_flag")==0) {
      params.stationary_flag = jt->second.as<bool>();
      } else if (param_name.compare("constrain_motion_flag")==0) {
      params.constrain_motion_flag = jt->second.as<bool>();
      } else if (param_name.compare("packing_fraction")==0) {
      params.packing_fraction = jt->second.as<double>();
      } else if (param_name.compare("n_equil")==0) {
      params.n_equil = jt->second.as<int>();
      } else {
        Logger::Warning("Unrecognized %s parameter: '%s'", sid.c_str(), param_name.c_str());
      }
    }
    return new rigid_filament_parameters(params);
  } else if (sid.compare("filament") == 0) {
    filament_parameters params;
    parse_species_base_params(params, node);
    for (auto jt = subnode.begin(); jt != subnode.end(); ++jt) {
      std::string param_name = jt->first.as<std::string>();
      if (false) {
      } else if (param_name.compare("name")==0) {
      params.name = jt->second.as<std::string>();
      } else if (param_name.compare("num")==0) {
      params.num = jt->second.as<int>();
      } else if (param_name.compare("diameter")==0) {
      params.diameter = jt->second.as<double>();
      } else if (param_name.compare("length")==0) {
      params.length = jt->second.as<double>();
      } else if (param_name.compare("insertion_type")==0) {
      params.insertion_type = jt->second.as<std::string>();
      } else if (param_name.compare("insert_file")==0) {
      params.insert_file = jt->second.as<std::string>();
      } else if (param_name.compare("overlap")==0) {
      params.overlap = jt->second.as<bool>();
      } else if (param_name.compare("draw_type")==0) {
      params.draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("color")==0) {
      params.color = jt->second.as<double>();
      } else if (param_name.compare("posit_flag")==0) {
      params.posit_flag = jt->second.as<bool>();
      } else if (param_name.compare("spec_flag")==0) {
      params.spec_flag = jt->second.as<bool>();
      } else if (param_name.compare("n_posit")==0) {
      params.n_posit = jt->second.as<int>();
      } else if (param_name.compare("n_spec")==0) {
      params.n_spec = jt->second.as<int>();
      } else if (param_name.compare("packing_fraction")==0) {
      params.packing_fraction = jt->second.as<double>();
      } else if (param_name.compare("persistence_length")==0) {
      params.persistence_length = jt->second.as<double>();
      } else if (param_name.compare("perlen_ratio")==0) {
      params.perlen_ratio = jt->second.as<double>();
      } else if (param_name.compare("polydispersity_flag")==0) {
      params.polydispersity_flag = jt->second.as<bool>();
      } else if (param_name.compare("max_length")==0) {
      params.max_length = jt->second.as<double>();
      } else if (param_name.compare("min_length")==0) {
      params.min_length = jt->second.as<double>();
      } else if (param_name.compare("min_bond_length")==0) {
      params.min_bond_length = jt->second.as<double>();
      } else if (param_name.compare("driving_factor")==0) {
      params.driving_factor = jt->second.as<double>();
      } else if (param_name.compare("n_equil")==0) {
      params.n_equil = jt->second.as<int>();
      } else if (param_name.compare("nematic_driving")==0) {
      params.nematic_driving = jt->second.as<bool>();
      } else if (param_name.compare("nematic_driving_freq")==0) {
      params.nematic_driving_freq = jt->second.as<double>();
      } else if (param_name.compare("peclet_number")==0) {
      params.peclet_number = jt->second.as<double>();
      } else if (param_name.compare("flexure_number")==0) {
      params.flexure_number = jt->second.as<double>();
      } else if (param_name.compare("radius_of_curvature")==0) {
      params.radius_of_curvature = jt->second.as<double>();
      } else if (param_name.compare("intrinsic_curvature")==0) {
      params.intrinsic_curvature = jt->second.as<double>();
      } else if (param_name.compare("intrinsic_curvature_sig")==0) {
      params.intrinsic_curvature_sig = jt->second.as<double>();
      } else if (param_name.compare("randomize_intrinsic_curvature_handedness")==0) {
      params.randomize_intrinsic_curvature_handedness = jt->second.as<bool>();
      } else if (param_name.compare("intrinsic_curvature_min")==0) {
      params.intrinsic_curvature_min = jt->second.as<double>();
      } else if (param_name.compare("highlight_handedness")==0) {
      params.highlight_handedness = jt->second.as<bool>();
      } else if (param_name.compare("highlight_curvature")==0) {
      params.highlight_curvature = jt->second.as<bool>();
      } else if (param_name.compare("draw_center_of_curvature")==0) {
      params.draw_center_of_curvature = jt->second.as<bool>();
      } else if (param_name.compare("error_analysis")==0) {
      params.error_analysis = jt->second.as<bool>();
      } else if (param_name.compare("theta_analysis")==0) {
      params.theta_analysis = jt->second.as<bool>();
      } else if (param_name.compare("lp_analysis")==0) {
      params.lp_analysis = jt->second.as<bool>();
      } else if (param_name.compare("global_order_analysis")==0) {
      params.global_order_analysis = jt->second.as<bool>();
      } else if (param_name.compare("polar_order_analysis")==0) {
      params.polar_order_analysis = jt->second.as<bool>();
      } else if (param_name.compare("polar_order_n_bins")==0) {
      params.polar_order_n_bins = jt->second.as<int>();
      } else if (param_name.compare("polar_order_contact_cutoff")==0) {
      params.polar_order_contact_cutoff = jt->second.as<double>();
      } else if (param_name.compare("polar_order_width")==0) {
      params.polar_order_width = jt->second.as<double>();
      } else if (param_name.compare("msd_analysis")==0) {
      params.msd_analysis = jt->second.as<bool>();
      } else if (param_name.compare("curvature_cluster_analysis")==0) {
      params.curvature_cluster_analysis = jt->second.as<bool>();
      } else if (param_name.compare("cluster_lifetime_min")==0) {
      params.cluster_lifetime_min = jt->second.as<int>();
      } else if (param_name.compare("cluster_by_handedness")==0) {
      params.cluster_by_handedness = jt->second.as<bool>();
      } else if (param_name.compare("spiral_init_flag")==0) {
      params.spiral_init_flag = jt->second.as<bool>();
      } else if (param_name.compare("spiral_analysis")==0) {
      params.spiral_analysis = jt->second.as<bool>();
      } else if (param_name.compare("spiral_number_fail_condition")==0) {
      params.spiral_number_fail_condition = jt->second.as<double>();
      } else if (param_name.compare("orientation_corr_analysis")==0) {
      params.orientation_corr_analysis = jt->second.as<bool>();
      } else if (param_name.compare("orientation_corr_n_steps")==0) {
      params.orientation_corr_n_steps = jt->second.as<int>();
      } else if (param_name.compare("crossing_analysis")==0) {
      params.crossing_analysis = jt->second.as<bool>();
      } else if (param_name.compare("flocking_analysis")==0) {
      params.flocking_analysis = jt->second.as<bool>();
      } else if (param_name.compare("flock_polar_min")==0) {
      params.flock_polar_min = jt->second.as<double>();
      } else if (param_name.compare("flock_contact_min")==0) {
      params.flock_contact_min = jt->second.as<double>();
      } else if (param_name.compare("highlight_flock")==0) {
      params.highlight_flock = jt->second.as<bool>();
      } else if (param_name.compare("flock_color_int")==0) {
      params.flock_color_int = jt->second.as<double>();
      } else if (param_name.compare("flock_color_ext")==0) {
      params.flock_color_ext = jt->second.as<double>();
      } else if (param_name.compare("number_fluctuation_analysis")==0) {
      params.number_fluctuation_analysis = jt->second.as<bool>();
      } else if (param_name.compare("number_fluctuation_boxes")==0) {
      params.number_fluctuation_boxes = jt->second.as<int>();
      } else if (param_name.compare("number_fluctuation_centers")==0) {
      params.number_fluctuation_centers = jt->second.as<int>();
      } else if (param_name.compare("in_out_analysis")==0) {
      params.in_out_analysis = jt->second.as<bool>();
      } else if (param_name.compare("drive_from_bond_center")==0) {
      params.drive_from_bond_center = jt->second.as<bool>();
      } else if (param_name.compare("flagella_flag")==0) {
      params.flagella_flag = jt->second.as<bool>();
      } else if (param_name.compare("flagella_freq")==0) {
      params.flagella_freq = jt->second.as<double>();
      } else if (param_name.compare("flagella_period")==0) {
      params.flagella_period = jt->second.as<double>();
      } else if (param_name.compare("flagella_amplitude")==0) {
      params.flagella_amplitude = jt->second.as<double>();
      } else if (param_name.compare("friction_ratio")==0) {
      params.friction_ratio = jt->second.as<double>();
      } else if (param_name.compare("dynamic_instability_flag")==0) {
      params.dynamic_instability_flag = jt->second.as<bool>();
      } else if (param_name.compare("force_induced_catastrophe_flag")==0) {
      params.force_induced_catastrophe_flag = jt->second.as<bool>();
      } else if (param_name.compare("optical_trap_flag")==0) {
      params.optical_trap_flag = jt->second.as<bool>();
      } else if (param_name.compare("optical_trap_spring")==0) {
      params.optical_trap_spring = jt->second.as<double>();
      } else if (param_name.compare("optical_trap_fixed")==0) {
      params.optical_trap_fixed = jt->second.as<bool>();
      } else if (param_name.compare("cilia_trap_flag")==0) {
      params.cilia_trap_flag = jt->second.as<bool>();
      } else if (param_name.compare("fic_factor")==0) {
      params.fic_factor = jt->second.as<double>();
      } else if (param_name.compare("f_shrink_to_grow")==0) {
      params.f_shrink_to_grow = jt->second.as<double>();
      } else if (param_name.compare("f_shrink_to_pause")==0) {
      params.f_shrink_to_pause = jt->second.as<double>();
      } else if (param_name.compare("f_pause_to_grow")==0) {
      params.f_pause_to_grow = jt->second.as<double>();
      } else if (param_name.compare("f_pause_to_shrink")==0) {
      params.f_pause_to_shrink = jt->second.as<double>();
      } else if (param_name.compare("f_grow_to_pause")==0) {
      params.f_grow_to_pause = jt->second.as<double>();
      } else if (param_name.compare("f_grow_to_shrink")==0) {
      params.f_grow_to_shrink = jt->second.as<double>();
      } else if (param_name.compare("v_poly")==0) {
      params.v_poly = jt->second.as<double>();
      } else if (param_name.compare("v_depoly")==0) {
      params.v_depoly = jt->second.as<double>();
      } else if (param_name.compare("custom_set_tail")==0) {
      params.custom_set_tail = jt->second.as<bool>();
      } else if (param_name.compare("reference_frame_flag")==0) {
      params.reference_frame_flag = jt->second.as<bool>();
      } else {
        Logger::Warning("Unrecognized %s parameter: '%s'", sid.c_str(), param_name.c_str());
      }
    }
    return new filament_parameters(params);
  } else if (sid.compare("br_bead") == 0) {
    br_bead_parameters params;
    parse_species_base_params(params, node);
    for (auto jt = subnode.begin(); jt != subnode.end(); ++jt) {
      std::string param_name = jt->first.as<std::string>();
      if (false) {
      } else if (param_name.compare("name")==0) {
      params.name = jt->second.as<std::string>();
      } else if (param_name.compare("num")==0) {
      params.num = jt->second.as<int>();
      } else if (param_name.compare("diameter")==0) {
      params.diameter = jt->second.as<double>();
      } else if (param_name.compare("length")==0) {
      params.length = jt->second.as<double>();
      } else if (param_name.compare("insertion_type")==0) {
      params.insertion_type = jt->second.as<std::string>();
      } else if (param_name.compare("insert_file")==0) {
      params.insert_file = jt->second.as<std::string>();
      } else if (param_name.compare("overlap")==0) {
      params.overlap = jt->second.as<bool>();
      } else if (param_name.compare("draw_type")==0) {
      params.draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("color")==0) {
      params.color = jt->second.as<double>();
      } else if (param_name.compare("posit_flag")==0) {
      params.posit_flag = jt->second.as<bool>();
      } else if (param_name.compare("spec_flag")==0) {
      params.spec_flag = jt->second.as<bool>();
      } else if (param_name.compare("n_posit")==0) {
      params.n_posit = jt->second.as<int>();
      } else if (param_name.compare("n_spec")==0) {
      params.n_spec = jt->second.as<int>();
      } else if (param_name.compare("driving_factor")==0) {
      params.driving_factor = jt->second.as<double>();
      } else if (param_name.compare("driving_torque")==0) {
      params.driving_torque = jt->second.as<double>();
      } else if (param_name.compare("density")==0) {
      params.density = jt->second.as<double>();
      } else if (param_name.compare("chiral_handedness")==0) {
      params.chiral_handedness = jt->second.as<int>();
      } else if (param_name.compare("rotational_noise")==0) {
      params.rotational_noise = jt->second.as<double>();
      } else if (param_name.compare("translational_noise")==0) {
      params.translational_noise = jt->second.as<double>();
      } else if (param_name.compare("randomize_handedness")==0) {
      params.randomize_handedness = jt->second.as<bool>();
      } else if (param_name.compare("highlight_handedness")==0) {
      params.highlight_handedness = jt->second.as<bool>();
      } else if (param_name.compare("alignment_interaction")==0) {
      params.alignment_interaction = jt->second.as<bool>();
      } else if (param_name.compare("alignment_torque")==0) {
      params.alignment_torque = jt->second.as<double>();
      } else if (param_name.compare("packing_fraction")==0) {
      params.packing_fraction = jt->second.as<double>();
      } else if (param_name.compare("draw_shape")==0) {
      params.draw_shape = jt->second.as<std::string>();
      } else {
        Logger::Warning("Unrecognized %s parameter: '%s'", sid.c_str(), param_name.c_str());
      }
    }
    return new br_bead_parameters(params);
  } else if (sid.compare("spherocylinder") == 0) {
    spherocylinder_parameters params;
    parse_species_base_params(params, node);
    for (auto jt = subnode.begin(); jt != subnode.end(); ++jt) {
      std::string param_name = jt->first.as<std::string>();
      if (false) {
      } else if (param_name.compare("name")==0) {
      params.name = jt->second.as<std::string>();
      } else if (param_name.compare("num")==0) {
      params.num = jt->second.as<int>();
      } else if (param_name.compare("diameter")==0) {
      params.diameter = jt->second.as<double>();
      } else if (param_name.compare("length")==0) {
      params.length = jt->second.as<double>();
      } else if (param_name.compare("insertion_type")==0) {
      params.insertion_type = jt->second.as<std::string>();
      } else if (param_name.compare("insert_file")==0) {
      params.insert_file = jt->second.as<std::string>();
      } else if (param_name.compare("overlap")==0) {
      params.overlap = jt->second.as<bool>();
      } else if (param_name.compare("draw_type")==0) {
      params.draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("color")==0) {
      params.color = jt->second.as<double>();
      } else if (param_name.compare("posit_flag")==0) {
      params.posit_flag = jt->second.as<bool>();
      } else if (param_name.compare("spec_flag")==0) {
      params.spec_flag = jt->second.as<bool>();
      } else if (param_name.compare("n_posit")==0) {
      params.n_posit = jt->second.as<int>();
      } else if (param_name.compare("n_spec")==0) {
      params.n_spec = jt->second.as<int>();
      } else if (param_name.compare("diffusion_analysis")==0) {
      params.diffusion_analysis = jt->second.as<bool>();
      } else if (param_name.compare("n_diffusion_samples")==0) {
      params.n_diffusion_samples = jt->second.as<int>();
      } else if (param_name.compare("midstep")==0) {
      params.midstep = jt->second.as<bool>();
      } else {
        Logger::Warning("Unrecognized %s parameter: '%s'", sid.c_str(), param_name.c_str());
      }
    }
    return new spherocylinder_parameters(params);
  } else if (sid.compare("spindle") == 0) {
    spindle_parameters params;
    parse_species_base_params(params, node);
    for (auto jt = subnode.begin(); jt != subnode.end(); ++jt) {
      std::string param_name = jt->first.as<std::string>();
      if (false) {
      } else if (param_name.compare("name")==0) {
      params.name = jt->second.as<std::string>();
      } else if (param_name.compare("num")==0) {
      params.num = jt->second.as<int>();
      } else if (param_name.compare("diameter")==0) {
      params.diameter = jt->second.as<double>();
      } else if (param_name.compare("length")==0) {
      params.length = jt->second.as<double>();
      } else if (param_name.compare("insertion_type")==0) {
      params.insertion_type = jt->second.as<std::string>();
      } else if (param_name.compare("insert_file")==0) {
      params.insert_file = jt->second.as<std::string>();
      } else if (param_name.compare("overlap")==0) {
      params.overlap = jt->second.as<bool>();
      } else if (param_name.compare("draw_type")==0) {
      params.draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("color")==0) {
      params.color = jt->second.as<double>();
      } else if (param_name.compare("posit_flag")==0) {
      params.posit_flag = jt->second.as<bool>();
      } else if (param_name.compare("spec_flag")==0) {
      params.spec_flag = jt->second.as<bool>();
      } else if (param_name.compare("n_posit")==0) {
      params.n_posit = jt->second.as<int>();
      } else if (param_name.compare("n_spec")==0) {
      params.n_spec = jt->second.as<int>();
      } else if (param_name.compare("n_filaments_bud")==0) {
      params.n_filaments_bud = jt->second.as<int>();
      } else if (param_name.compare("n_filaments_mother")==0) {
      params.n_filaments_mother = jt->second.as<int>();
      } else if (param_name.compare("alignment_potential")==0) {
      params.alignment_potential = jt->second.as<bool>();
      } else if (param_name.compare("k_spring")==0) {
      params.k_spring = jt->second.as<double>();
      } else if (param_name.compare("k_align")==0) {
      params.k_align = jt->second.as<double>();
      } else if (param_name.compare("spring_length")==0) {
      params.spring_length = jt->second.as<double>();
      } else if (param_name.compare("spb_diameter")==0) {
      params.spb_diameter = jt->second.as<double>();
      } else {
        Logger::Warning("Unrecognized %s parameter: '%s'", sid.c_str(), param_name.c_str());
      }
    }
    return new spindle_parameters(params);
  } else if (sid.compare("crosslink") == 0) {
    crosslink_parameters params;
    parse_species_base_params(params, node);
    for (auto jt = subnode.begin(); jt != subnode.end(); ++jt) {
      std::string param_name = jt->first.as<std::string>();
      if (false) {
      } else if (param_name.compare("name")==0) {
      params.name = jt->second.as<std::string>();
      } else if (param_name.compare("num")==0) {
      params.num = jt->second.as<int>();
      } else if (param_name.compare("diameter")==0) {
      params.diameter = jt->second.as<double>();
      } else if (param_name.compare("length")==0) {
      params.length = jt->second.as<double>();
      } else if (param_name.compare("insertion_type")==0) {
      params.insertion_type = jt->second.as<std::string>();
      } else if (param_name.compare("insert_file")==0) {
      params.insert_file = jt->second.as<std::string>();
      } else if (param_name.compare("overlap")==0) {
      params.overlap = jt->second.as<bool>();
      } else if (param_name.compare("draw_type")==0) {
      params.draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("color")==0) {
      params.color = jt->second.as<double>();
      } else if (param_name.compare("posit_flag")==0) {
      params.posit_flag = jt->second.as<bool>();
      } else if (param_name.compare("spec_flag")==0) {
      params.spec_flag = jt->second.as<bool>();
      } else if (param_name.compare("n_posit")==0) {
      params.n_posit = jt->second.as<int>();
      } else if (param_name.compare("n_spec")==0) {
      params.n_spec = jt->second.as<int>();
      } else if (param_name.compare("concentration")==0) {
      params.concentration = jt->second.as<double>();
      } else if (param_name.compare("use_binding_volume")==0) {
      params.use_binding_volume = jt->second.as<bool>();
      } else if (param_name.compare("infinite_reservoir_flag")==0) {
      params.infinite_reservoir_flag = jt->second.as<bool>();
      } else if (param_name.compare("bind_site_density")==0) {
      params.bind_site_density = jt->second.as<double>();
      } else if (param_name.compare("static_flag")==0) {
      params.static_flag = jt->second.as<bool>();
      } else if (param_name.compare("diffusion_s")==0) {
      params.diffusion_s = jt->second.as<double>();
      } else if (param_name.compare("diffusion_d")==0) {
      params.diffusion_d = jt->second.as<double>();
      } else if (param_name.compare("velocity_s")==0) {
      params.velocity_s = jt->second.as<double>();
      } else if (param_name.compare("velocity_d")==0) {
      params.velocity_d = jt->second.as<double>();
      } else if (param_name.compare("k_on_s")==0) {
      params.k_on_s = jt->second.as<double>();
      } else if (param_name.compare("k_off_s")==0) {
      params.k_off_s = jt->second.as<double>();
      } else if (param_name.compare("k_on_d")==0) {
      params.k_on_d = jt->second.as<double>();
      } else if (param_name.compare("k_off_d")==0) {
      params.k_off_d = jt->second.as<double>();
      } else if (param_name.compare("energy_dep_factor")==0) {
      params.energy_dep_factor = jt->second.as<double>();
      } else if (param_name.compare("force_dep_length")==0) {
      params.force_dep_length = jt->second.as<double>();
      } else if (param_name.compare("polar_affinity")==0) {
      params.polar_affinity = jt->second.as<double>();
      } else if (param_name.compare("k_spring")==0) {
      params.k_spring = jt->second.as<double>();
      } else if (param_name.compare("k_spring_compress")==0) {
      params.k_spring_compress = jt->second.as<double>();
      } else if (param_name.compare("f_stall")==0) {
      params.f_stall = jt->second.as<double>();
      } else if (param_name.compare("force_dep_vel_flag")==0) {
      params.force_dep_vel_flag = jt->second.as<bool>();
      } else if (param_name.compare("k_align")==0) {
      params.k_align = jt->second.as<double>();
      } else if (param_name.compare("rest_length")==0) {
      params.rest_length = jt->second.as<double>();
      } else if (param_name.compare("step_direction")==0) {
      params.step_direction = jt->second.as<int>();
      } else if (param_name.compare("tether_draw_type")==0) {
      params.tether_draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("tether_diameter")==0) {
      params.tether_diameter = jt->second.as<double>();
      } else if (param_name.compare("tether_color")==0) {
      params.tether_color = jt->second.as<double>();
      } else if (param_name.compare("end_pausing")==0) {
      params.end_pausing = jt->second.as<bool>();
      } else if (param_name.compare("r_capture")==0) {
      params.r_capture = jt->second.as<double>();
      } else if (param_name.compare("lut_grid_num")==0) {
      params.lut_grid_num = jt->second.as<int>();
      } else {
        Logger::Warning("Unrecognized %s parameter: '%s'", sid.c_str(), param_name.c_str());
      }
    }
    return new crosslink_parameters(params);
  } else if (sid.compare("optical_trap") == 0) {
    optical_trap_parameters params;
    parse_species_base_params(params, node);
    for (auto jt = subnode.begin(); jt != subnode.end(); ++jt) {
      std::string param_name = jt->first.as<std::string>();
      if (false) {
      } else if (param_name.compare("name")==0) {
      params.name = jt->second.as<std::string>();
      } else if (param_name.compare("num")==0) {
      params.num = jt->second.as<int>();
      } else if (param_name.compare("diameter")==0) {
      params.diameter = jt->second.as<double>();
      } else if (param_name.compare("length")==0) {
      params.length = jt->second.as<double>();
      } else if (param_name.compare("insertion_type")==0) {
      params.insertion_type = jt->second.as<std::string>();
      } else if (param_name.compare("insert_file")==0) {
      params.insert_file = jt->second.as<std::string>();
      } else if (param_name.compare("overlap")==0) {
      params.overlap = jt->second.as<bool>();
      } else if (param_name.compare("draw_type")==0) {
      params.draw_type = jt->second.as<std::string>();
      } else if (param_name.compare("color")==0) {
      params.color = jt->second.as<double>();
      } else if (param_name.compare("posit_flag")==0) {
      params.posit_flag = jt->second.as<bool>();
      } else if (param_name.compare("spec_flag")==0) {
      params.spec_flag = jt->second.as<bool>();
      } else if (param_name.compare("n_posit")==0) {
      params.n_posit = jt->second.as<int>();
      } else if (param_name.compare("n_spec")==0) {
      params.n_spec = jt->second.as<int>();
      } else if (param_name.compare("trap_spring")==0) {
      params.trap_spring = jt->second.as<double>();
      } else if (param_name.compare("attach_species")==0) {
      params.attach_species = jt->second.as<std::string>();
      } else if (param_name.compare("trap_diameter")==0) {
      params.trap_diameter = jt->second.as<double>();
      } else if (param_name.compare("bead_diameter")==0) {
      params.bead_diameter = jt->second.as<double>();
      } else if (param_name.compare("trap_color")==0) {
      params.trap_color = jt->second.as<double>();
      } else if (param_name.compare("bead_color")==0) {
      params.bead_color = jt->second.as<double>();
      } else {
        Logger::Warning("Unrecognized %s parameter: '%s'", sid.c_str(), param_name.c_str());
      }
    }
    return new optical_trap_parameters(params);
  } else {
    Logger::Error("Unrecognized SID '%s' in parse_params!", sid.c_str());
  }
  return nullptr;
}

#endif // _CGLASS_PARSE_PARAMS_H_