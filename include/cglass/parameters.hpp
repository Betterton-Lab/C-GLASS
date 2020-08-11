#ifndef _CGLASS_PARAMETERS_H_
#define _CGLASS_PARAMETERS_H_

#include "definitions.hpp"

#include <string>

template <unsigned char S> struct species_parameters {
  std::string name = "species";
  int num = 0;
  double diameter = 1;
  double length = 0;
  std::string insertion_type = "random";
  std::string insert_file = "none";
  bool overlap = false;
  std::string draw_type = "orientation";
  double color = 0;
  bool posit_flag = false;
  bool spec_flag = false;
  int n_posit = 100;
  int n_spec = 100;
  virtual ~species_parameters() {}
};

typedef species_parameters<species_id::none> species_base_parameters;

template <>
struct species_parameters<species_id::rigid_filament>
    : public species_base_parameters {
  double max_length = 500;
  double min_length = 5;
  bool stationary_flag = false;
  bool constrain_motion_flag = false;
  double packing_fraction = -1;
  int n_equil = 0;
};
typedef species_parameters<species_id::rigid_filament> rigid_filament_parameters;

template <>
struct species_parameters<species_id::filament>
    : public species_base_parameters {
  double packing_fraction = -1;
  double persistence_length = 400;
  double perlen_ratio = -1;
  bool polydispersity_flag = false;
  double max_length = 500;
  double min_length = 5;
  double min_bond_length = 1.5;
  double driving_factor = 0;
  int n_equil = 0;
  bool nematic_driving = false;
  double nematic_driving_freq = 0;
  double peclet_number = -1;
  double flexure_number = -1;
  double radius_of_curvature = -1;
  double intrinsic_curvature = 0;
  double intrinsic_curvature_sig = 0;
  bool randomize_intrinsic_curvature_handedness = false;
  double intrinsic_curvature_min = 0;
  bool highlight_handedness = false;
  bool highlight_curvature = false;
  bool draw_center_of_curvature = false;
  bool error_analysis = false;
  bool theta_analysis = false;
  bool lp_analysis = false;
  bool global_order_analysis = false;
  bool polar_order_analysis = false;
  int polar_order_n_bins = 100;
  double polar_order_contact_cutoff = 3;
  double polar_order_width = 10;
  bool msd_analysis = false;
  bool curvature_cluster_analysis = false;
  int cluster_lifetime_min = 100;
  bool cluster_by_handedness = false;
  bool spiral_init_flag = false;
  bool spiral_analysis = false;
  double spiral_number_fail_condition = 0;
  bool orientation_corr_analysis = false;
  int orientation_corr_n_steps = 1000;
  bool crossing_analysis = false;
  bool flocking_analysis = false;
  double flock_polar_min = 0.5;
  double flock_contact_min = 0.5;
  bool highlight_flock = false;
  double flock_color_int = 1.57;
  double flock_color_ext = 4.71;
  bool number_fluctuation_analysis = false;
  int number_fluctuation_boxes = 6;
  int number_fluctuation_centers = 10;
  bool in_out_analysis = false;
  bool drive_from_bond_center = true;
  bool flagella_flag = false;
  double flagella_freq = 1;
  double flagella_period = 2;
  double flagella_amplitude = 1;
  double friction_ratio = 2;
  bool dynamic_instability_flag = false;
  bool force_induced_catastrophe_flag = false;
  bool optical_trap_flag = false;
  double optical_trap_spring = 20;
  bool optical_trap_fixed = false;
  bool cilia_trap_flag = false;
  double fic_factor = 0.828;
  double f_shrink_to_grow = 0.017;
  double f_shrink_to_pause = 0.0;
  double f_pause_to_grow = 0.0;
  double f_pause_to_shrink = 0.0;
  double f_grow_to_pause = 0.0;
  double f_grow_to_shrink = 0.00554;
  double v_poly = 0.44;
  double v_depoly = 0.793;
  bool custom_set_tail = false;
  bool reference_frame_flag = false;
};
typedef species_parameters<species_id::filament> filament_parameters;

template <>
struct species_parameters<species_id::br_bead>
    : public species_base_parameters {
  double driving_factor = 0;
  double driving_torque = 0;
  double density = -1;
  int chiral_handedness = 0;
  double rotational_noise = 1;
  double translational_noise = 1;
  bool randomize_handedness = false;
  bool highlight_handedness = false;
  bool alignment_interaction = false;
  double alignment_torque = 0;
  double packing_fraction = -1;
  std::string draw_shape = "sphere";
};
typedef species_parameters<species_id::br_bead> br_bead_parameters;

template <>
struct species_parameters<species_id::spherocylinder>
    : public species_base_parameters {
  bool diffusion_analysis = false;
  int n_diffusion_samples = 1;
  bool midstep = false;
};
typedef species_parameters<species_id::spherocylinder> spherocylinder_parameters;

template <>
struct species_parameters<species_id::spindle>
    : public species_base_parameters {
  int n_filaments_bud = 0;
  int n_filaments_mother = 0;
  bool alignment_potential = false;
  double k_spring = 1000;
  double k_align = 0;
  double spring_length = 0;
  double spb_diameter = 5;
};
typedef species_parameters<species_id::spindle> spindle_parameters;

template <>
struct species_parameters<species_id::crosslink>
    : public species_base_parameters {
  double concentration = 0;
  bool use_binding_volume = true;
  bool infinite_reservoir_flag = false;
  double bind_site_density = 1;
  bool static_flag = false;
  double diffusion_s = 0;
  double diffusion_d = 0;
  double velocity_s = 0;
  double velocity_d = 0;
  double k_on_s = 10;
  double k_off_s = 2;
  double k_on_d = 10;
  double k_off_d = 2;
  double energy_dep_factor = 0;
  double force_dep_length = 0;
  double polar_affinity = 1;
  double k_spring = 10;
  double k_spring_compress = -1.;
  double f_stall = 100;
  bool force_dep_vel_flag = true;
  double k_align = 0;
  double rest_length = 0;
  int step_direction = 0;
  std::string tether_draw_type = "orientation";
  double tether_diameter = 0.5;
  double tether_color = 3.1416;
  bool end_pausing = false;
  double r_capture = 5;
  int lut_grid_num = 256;
};
typedef species_parameters<species_id::crosslink> crosslink_parameters;

template <>
struct species_parameters<species_id::optical_trap>
    : public species_base_parameters {
  double trap_spring = 1;
  std::string attach_species = "filament";
  double trap_diameter = 4;
  double bead_diameter = 2;
  double trap_color = 3.1416;
  double bead_color = 1.5708;
};
typedef species_parameters<species_id::optical_trap> optical_trap_parameters;

struct system_parameters {
  long seed = 7859459105545;
  int n_runs = 1;
  int n_random = 1;
  std::string run_name = "sc";
  int n_dim = 3;
  int n_periodic = 0;
  int boundary = 0;
  double system_radius = 100;
  int n_steps = 1000000;
  int i_step = 0;
  int prev_step = 0;
  double delta = 0.001;
  bool dynamic_timestep = false;
  double dynamic_timestep_ramp = 0.001;
  bool graph_flag = false;
  int n_graph = 1000;
  double graph_diameter = 0;
  double object_opacity = 1;
  bool invert_background = false;
  bool draw_boundary = true;
  bool load_checkpoint = false;
  std::string checkpoint_run_name = "sc";
  int n_load = 0;
  bool movie_flag = false;
  std::string movie_directory = "frames";
  bool time_analysis = false;
  double bud_height = 680;
  double bud_radius = 300;
  double lj_epsilon = 1;
  double wca_eps = 1;
  double wca_sig = 1;
  double ss_a = 1;
  double ss_rs = 1.5;
  double ss_eps = 1;
  double f_cutoff = 2000;
  bool constant_pressure = false;
  bool constant_volume = false;
  double target_pressure = 0;
  double target_radius = 100;
  int pressure_time = 100;
  double compressibility = 1;
  bool zero_temperature = false;
  bool thermo_flag = false;
  int n_thermo = 1000;
  double insert_radius = -1;
  bool interaction_flag = true;
  bool remove_duplicate_interactions = false;
  bool coarse_grained_mesh_interactions = false;
  int mesh_coarsening = 2;
  int species_insertion_failure_threshold = 10000;
  int species_insertion_reattempt_threshold = 10;
  bool uniform_crystal = false;
  int n_steps_equil = 0;
  int n_steps_target = 100000;
  bool static_particle_number = false;
  bool checkpoint_from_spec = false;
  std::string potential = "wca";
  double soft_potential_mag = 10;
  double soft_potential_mag_target = -1;
  bool like_like_interactions = true;
  bool auto_graph = false;
  bool local_order_analysis = false;
  double local_order_width = 50;
  double local_order_bin_width = 0.5;
  int local_order_n_analysis = 100;
  int density_analysis = 0;
  double density_bin_width = 0.1;
  bool density_com_only = false;
  bool overlap_analysis = false;
  bool highlight_overlaps = false;
  bool reduced = false;
  bool reload_reduce_switch = false;
  bool checkpoint_flag = false;
  int n_checkpoint = 10000;
  bool no_midstep = false;
};

#endif // _CGLASS_PARAMETERS_H_