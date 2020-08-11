
#include "cglass/simulation.hpp"

/* Initialize simulation parameters and run simulation */
void Simulation::Run(YAML::Node sim_params) {
  YAML::Emitter out;
  Logger::Info("Initializing simulation with parameters:\n%s",
               (out << sim_params).c_str());
  // Parse simulation parameters
  parser_.Init(sim_params);
  // Initialize simulation data structures
  InitSimulation();
  // Begin simulation
  RunSimulation();
  // Tear down data structures (e.g. cell lists, etc)
  ClearSimulation();
}

/* Continuously loop over object data structures, solving their equations of
 * motion using symplectic integration schemes for n_steps time steps. */
void Simulation::RunSimulation() {
#ifdef ENABLE_OPENMP
  Logger::Info("Running simulation %s on %d threads", run_name_.c_str(),
               omp_get_max_threads());
#else
  Logger::Info("Running simulation %s", run_name_.c_str());
#endif
  /* Solve equations of motion for all objects */

  /* Note that we are beginning with i_step = 1 here as a convention. This is
   because filaments use a midstep method, so that at the end of step 1 and
   every odd step, they will be in a midstep configuration, whereas at the end
   of every even step, they will be in the fullstep configuration. Since most
   sane step frequency functions (e.g. n_graph, n_thermo, etc) are even numbers,
   we want these to represent filaments in their fullstep configurations. We
   also update other objects' positions on every even step too (e.g. xlinks). */
  double delta_diff = 0;
  bool midstep = true;
  /* With halfstep algorithm, each step only moves the simulation forward
     by 1/2 delta. Defining step_fact prevents unnecessary if-else statements.*/
  double step_fact = params_.no_midstep ? 1. : .5;
  time_ = 0;
  params_.i_step = 0;
  for (i_step_ = 1; params_.i_step < params_.n_steps + 1; ++i_step_) {
    time_ += step_fact * Object::GetDelta();
    params_.i_step = i_step_;
    if (params_.dynamic_timestep) {
      // Calculate nominal timestep
      params_.i_step = (int)round(time_ / (step_fact * params_.delta));
    }
    // Output progress
    PrintComplete();
    /* Zero the force_ array on all objects for bookkeeping */
    ZeroForces();
    /* Calculate forces between objects in system */
    Interact();
    /*Decrease or increase timestep depending on size of forces and torques 
      on objects. Timestep remains the same if equal to initial timestep. */
    if (params_.dynamic_timestep && ix_mgr_.CheckDynamicTimestep()) {
      if (i_step_ <= 1) {
        /* We don't have a previous step to revert back to. */
        Logger::Error(
            "Dynamic timestep triggered on the first simulation step");
      }
      if (midstep && !params_.no_midstep) {
        time_ -= Object::GetDelta();
        i_step_ -= 2;
      } else {
        time_ -= step_fact * Object::GetDelta();
        i_step_ -= 1;
      }
      midstep = true;
      Object::SetDelta(0.5 * Object::GetDelta());
      delta_diff = params_.delta - Object::GetDelta();
      if (Object::GetDelta() < 1e-12) {
        Logger::Warning(
            "Dynamic delta has become ridiculously tiny! (%2.14f) Likely an"
            " interaction error occurred. Triggering an early exit. %d",
            Object::GetDelta(), i_step_);
        return;
      }
      continue;
    } else if (params_.dynamic_timestep && Object::GetDelta() < params_.delta) {
      Object::SetDelta(Object::GetDelta() +
                       params_.dynamic_timestep_ramp * delta_diff);
    }
    /* Generate all output files. The reason that this follows the Interact step
       and precedes the Integrate step is so that we can detect whether the
       timestep used in the Interact step was too large and adjust it
       accordingly before writing the resulting positions from the Integrate
       step. */
    WriteOutputs();
    /* Integrate EOM on all objects to update positions, etc */
    Integrate();
    /* Update system pressure, volume, etc if necessary */
    Statistics();
    /* Draw the particles in graphics window, if using graphics */
    Draw();
    /* Catch soft exceptions and terminate the simulation early */
    if (early_exit) {
      early_exit = false;
      Logger::Info("Early exit triggered. Ending simulation.");
      return;
    }

    midstep = !midstep;
  }
}

/* Print simulation progress as a percentage. If the print_complete parameter
 * evaluates to true, this will print a new line for each progress update,
 * which can be useful for logging. Otherwise, the progress will be refreshed
 * in a single line printed to standard output. */
void Simulation::PrintComplete() {
  long iteration = params_.i_step * 10;
  long steps = params_.n_steps;
  if (iteration % steps == 0) {
    Logger::Info("%d%% complete", 10 * iteration / steps);
    if (params_.dynamic_timestep) {
      Logger::Info("Current simulation time: %2.8f", time_);
      Logger::Info("Current timestep: %2.10f", Object::GetDelta());
    }
  }
  if (i_step_ == log_interval_ + 1) {
    Logger::Info("%d steps completed", log_interval_);
    log_interval_ = (int)floor(2 * log_interval_);
    if (params_.dynamic_timestep) {
      Logger::Info("Current simulation time: %2.8f", time_);
      Logger::Info("Current timestep: %2.10f", Object::GetDelta());
    }
  }
  Logger::Trace("*****Step %d*****", i_step_);
}

/* Update the positions of all objects in the system using numerical
 * integration for one time step defined by the delta parameter. */
void Simulation::Integrate() {
  Logger::Trace("Updating object positions");
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->UpdatePositions();
  }
}

/* Calculate interaction forces between all objects if necessary. */
void Simulation::Interact() {
  Logger::Trace("Calculating object interactions");
  ix_mgr_.Interact();
}

/* Remove forces on all objects. Objects with inertia that use higher order RK
 * schemes with dependence on previous forces need to save their forces in
 * prev_forces_ arrays. */
void Simulation::ZeroForces() {
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->ZeroForces();
  }
}

/* Update system pressure, volume and rescale system size if necessary,
 * handling periodic boundaries in a sane way. */
void Simulation::Statistics() {
  if (params_.i_step % params_.n_thermo == 0 &&
      params_.i_step != params_.prev_step && params_.i_step > 0) {
    Logger::Debug("Calculating system pressure and volume");
    /* Calculate system pressure from stress tensor */
    Logger::Debug("Calculating system thermodynamics");
    ix_mgr_.CalculatePressure();
    if (params_.constant_pressure) {
      space_.ConstantPressure();
    }
    /* Constant volume is used to adiabatically adjust the system size to
     achieve high densities from random insertions. */
    else if (params_.constant_volume) {
      space_.ConstantVolume();
    }
  }
  /* Update object positions, etc and handle periodic BCs */
  if (space_.GetUpdate()) {
    Logger::Debug("Rescaling species scaled positions in dynamic space");
    space_.UpdateSpace();
    ScaleSpeciesPositions();
  }
}

/* Update object positions in periodic space if the system size changed. */
void Simulation::ScaleSpeciesPositions() {
  for (auto spec : species_) {
    spec->ScalePositions();
  }
}

/* Initialize all the data structures in the simulation */
void Simulation::InitSimulation() {
  params_ = parser_.ParseSystemParameters();
  run_name_ = params_.run_name;
  rng_ = new RNG(params_.seed);

#ifdef TRACE
  if (params_.n_steps > 100) {
    Logger::Warning(
        "Simulation run in trace mode with a large number of n_steps (%d).",
        params_.n_steps);
    fprintf(stderr, "Continue anyway? (y/N) ");
    char c;
    if (std::cin.peek() != 'y') {
      Logger::Error("Terminating simulation by user request");
    } else if (!(std::cin >> c)) {
      Logger::Error("Invalid input");
    } else {
      fprintf(stderr, "Resuming simulation\n");
      std::cin.ignore();
    }
  }
#endif
  space_.Init(&params_);
  InitObjects();
  ix_mgr_.Init(&params_, &species_, space_.GetStruct());
  InitSpecies();
  ix_mgr_.InitInteractions();
  InsertSpecies(params_.load_checkpoint, params_.load_checkpoint);
  InitOutputs();
  if (params_.graph_flag) {
    InitGraphics();
  }
  params_.i_step = 0;
  WriteOutputs();
}

/* Initialize static object parameters that are used everywhere */
void Simulation::InitObjects() {
  Object::SetParams(&params_);
  Object::SetNDim(params_.n_dim);
  Object::SetDelta(params_.delta);
  Object::SetSpace(space_.GetStruct());
  SpeciesBase::SetParams(&params_);
  SpeciesBase::SetSpace(space_.GetStruct());
  AnalysisBase::SetParams(&params_);
  AnalysisBase::SetSpace(space_.GetStruct());
}

/* Generate graphics window and draw initial simulation setup */
void Simulation::InitGraphics() {
  GetGraphicsStructure();
  double background_color = (params_.invert_background ? 1 : 0.1);
// If NOGRAPH is defined, skip drawing and grabbing
#ifndef NOGRAPH
  // Initialize graphics structures
  graphics_.Init(&graph_array_, space_.GetStruct(), background_color,
                 params_.draw_boundary, params_.auto_graph,
                 params_.object_opacity);

// This line was interferring with graphics on Windows, and removing it did no
// harm
#ifndef WINGRAPH
  graphics_.DrawLoop();
#endif

#endif
  // Initialize directory for grabbed images
  params_.movie_directory.append("/");
  params_.movie_directory.append(params_.run_name);
#ifndef NOGRAPH
  // Grab first frame
  if (params_.movie_flag) {
    // Record bmp image of frame into movie_directory
    grabber(graphics_.windx_, graphics_.windy_, params_.movie_directory,
            (int)params_.i_step / params_.n_graph);
  }
#endif
}

/* Initialize object types */
void Simulation::InitSpecies() {
  std::vector<sid_label> species_labels = parser_.GetSpeciesLabels();
  species_.reserve(parser_.GetNSpecies());
  SpeciesFactory species_factory;
  for (auto slab = species_labels.begin(); slab != species_labels.end();
       ++slab) {
    const species_id sid = species_id::_from_string(slab->first.c_str());
    if (sid == +species_id::crosslink) {
      ix_mgr_.InitCrosslinkSpecies(*slab, parser_, rng_->GetSeed());
      continue;
    }
    if (sid == +species_id::optical_trap) {
      ix_mgr_.InitOpticalTrapSpecies(*slab, parser_, rng_->GetSeed());
      continue;
    }
    species_.push_back(species_factory.CreateSpecies(sid, rng_->GetSeed()));
    species_.back()->Init(slab->second, parser_);
    if (species_.back()->GetNInsert() > 0) {
#ifdef TRACE
      if (species_.back()->GetNInsert() > 20) {
        Logger::Warning("Simulation run in trace mode with a large number of "
                        "objects in species %s (%d).",
                        sid._to_string(), species_.back()->GetNInsert());
        fprintf(stderr, "Continue anyway? (y/N) ");
        char c;
        if (std::cin.peek() != 'y') {
          Logger::Error("Terminating simulation by user request");
        } else if (!(std::cin >> c)) {
          Logger::Error("Invalid input");
        } else {
          fprintf(stderr, "Resuming simulation\n");
          std::cin.ignore();
        }
      }
#endif
      species_.back()->Reserve();
      double spec_length = species_.back()->GetSpecLength();
      double spec_d = species_.back()->GetSpecDiameter();
      spec_length = (spec_d > spec_length ? spec_d : spec_length);
      CellList::SetMinCellLength(1.5 * spec_length);
    } else {
      species_.back()->CleanUp();
      delete species_.back();
      species_.pop_back();
    }
  }
}

/* Initialize object positions and orientations.*/
void Simulation::InsertSpecies(bool force_overlap, bool processing) {
  Logger::Info("Inserting species");
  for (auto spec = species_.begin(); spec != species_.end(); ++spec) {
    // Check for random insertion
    if (processing ||
        (*spec)->GetInsertionType().find("random") == std::string::npos) {
      /* Insertion is non-random: don't check for overlaps */
      force_overlap = true;
    }
    int num = (*spec)->GetNInsert();
    bool not_done = true;
    int inserted = 0;
    int num_attempts = 0;
    while (num != inserted) {
      inserted = 0;
      int num_failures = 0;
      while (num != inserted) {
        (*spec)->AddMember();
        /* Update the number of particles we need to insert, in case a species
           needs to have a certain packing fraction */
        num = (*spec)->GetNInsert();
        // First check that we are respecting boundary conditions
        std::vector<Object *> last_ixors;
        (*spec)->GetLastInteractors(last_ixors);
        if (params_.boundary != 0 && !processing &&
            ix_mgr_.CheckBoundaryConditions(last_ixors)) {
          (*spec)->PopMember();
          /* We are not counting boundary condition failures in insertion
           failures, since insertion failures are for packing issues */
        }
        // Check if we have an overlap of objects
        else if (!force_overlap && !(*spec)->CanOverlap() && !processing &&
                 ix_mgr_.CheckOverlap(last_ixors)) {
          (*spec)->PopMember();
          num_failures++;
        }
        /* Otherwise update interaction engine to include new interactors and
           update display of percentage of species inserted */
        else {
          inserted++;
          if (!force_overlap && !processing) {
            ix_mgr_.AddInteractors(last_ixors);
          }
        }
        if (num_failures > params_.species_insertion_failure_threshold) {
          Logger::Warning("Too many insertion failures have occurred: managed "
                          "to insert %2.1f%% of objects",
                          100.0 * inserted / (float)num);
          break;
        }
      }
      if (num != inserted && params_.n_dim == 2) {
        // Attempt a lattice-based insertion strategy (only 2d for now)
        Logger::Warning("Attempting lattice-based insertion strategy");
        double pos[3] = {0, 0, 0};
        pos[0] = -params_.system_radius;
        pos[1] = -params_.system_radius;
        double d = 0.5 * (*spec)->GetSpecDiameter();
        double l = 0.25 * (*spec)->GetSpecLength();
        l = (l < d ? d : l);
        int num_x = (int)floor(2 * params_.system_radius / d);
        int num_y = (int)floor(2 * params_.system_radius / l);
        std::vector<std::pair<int, int>> grid_array;
        for (int i = 0; i < num_x; ++i) {
          for (int j = 0; j < num_y; ++j) {
            grid_array.push_back(std::make_pair(i, j));
          }
        }
        int *grid_index = new int[num_x * num_y];
        for (int i = 0; i < num_x * num_y; ++i) {
          grid_index[i] = i;
        }
        rng_->Shuffle<int>(grid_index, num_x * num_y);
        // gsl_ran_shuffle(rng_.r, grid_index, num_x * num_y, sizeof(int));
        for (int i = 0; i < num_x * num_y; ++i) {
          pos[0] = grid_array[grid_index[i]].first * d;
          pos[1] = grid_array[grid_index[i]].second * l;
          (*spec)->AddMember();
          /* Update the number of particles we need to insert, in case a
             species needs to have a certain packing fraction */
          num = (*spec)->GetNInsert();
          (*spec)->SetLastMemberPosition(pos);
          // First check that we are respecting boundary conditions
          std::vector<Object *> last_ixors;
          (*spec)->GetLastInteractors(last_ixors);
          if (params_.boundary != 0 && !processing &&
              ix_mgr_.CheckBoundaryConditions(last_ixors)) {
            (*spec)->PopMember();
            // We are not counting boundary condition failures in insertion
            // failures, since insertion failures are for packing issues
          }
          // Check if we have an overlap of objects
          else if (!force_overlap && !(*spec)->CanOverlap() && !processing &&
                   ix_mgr_.CheckOverlap(last_ixors)) {
            (*spec)->PopMember();
            num_failures++;
          }
          // Otherwise update display of percentage of species inserted
          else {
            inserted++;
            ix_mgr_.AddInteractors(last_ixors);
          }
          if (inserted == num)
            break;
        }
        delete[] grid_index;
      }
      if (num != inserted) {
        Logger::Warning("Species insertion failure threshold of %d reached. "
                        "Reattempting insertion.\n",
                        params_.species_insertion_failure_threshold);
        (*spec)->PopAll();
        ix_mgr_.Reset();
      }
      if (++num_attempts > params_.species_insertion_reattempt_threshold) {
        Logger::Error("Unable to insert species randomly within the reattempt "
                      "threshold of %d.\n",
                      params_.species_insertion_reattempt_threshold);
      }
    }
    if (!processing) {
      (*spec)->InitMembers();
      if ((*spec)->GetInsertionType().find("random") == std::string::npos) {
        (*spec)->ArrangeMembers();
      }
    }
  }
  /* Initialize static crosslink positions */
  ix_mgr_.InsertCrosslinks();
  /* Initialize optical traps attached to species */
  ix_mgr_.InsertOpticalTraps(&species_);
  /* Should do this all the time to force object counting */
  if (params_.load_checkpoint) {
    for (auto spec = species_.begin(); spec != species_.end(); ++spec) {
      (*spec)->LoadFromCheckpoints(run_name_, params_.checkpoint_run_name);
    }
    ix_mgr_.LoadCrosslinksFromCheckpoints(run_name_,
                                          params_.checkpoint_run_name);
    //TODO Load optical traps from checkpoints
  }

  ix_mgr_.ResetCellList(); // Forces rebuild cell list without redundancy
  ix_mgr_.Reset();
  // if (!processing) {
  ix_mgr_.CheckUpdateObjects(); // Forces update as well
  //}
}

/* Tear down data structures, e.g. cell lists, and close graphics window if
 * necessary. */
void Simulation::ClearSimulation() {
  Logger::Debug("Clearing simulation resources");
  output_mgr_.Close();
  ClearSpecies();
  ix_mgr_.Clear();
#ifndef NOGRAPH
  if (params_.graph_flag) {
    graphics_.Clear();
  }
#endif
  delete rng_;
  Logger::Info("Simulation complete");
}

/* Tear down object type data structures */
void Simulation::ClearSpecies() {
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->CleanUp();
    delete (*it);
  }
  species_.clear();
}

/* Update the OpenGL graphics window */
void Simulation::Draw(bool single_frame) {
#ifndef NOGRAPH
  if (params_.graph_flag && params_.i_step % params_.n_graph == 0 &&
      params_.i_step != params_.prev_step) {
    Logger::Trace("Drawing graphable objects");
    /* Get updated object positions and orientations */
    GetGraphicsStructure();
    graphics_.Draw();
    if (params_.movie_flag) {
      /* Record bmp image of frame */
      grabber(graphics_.windx_, graphics_.windy_, params_.movie_directory,
              (single_frame ? 0 : frame_num_++));
    }
  }
#endif
}

/* Loop through objects in simulation to update their positions and
 * orientations in the graphics window */
void Simulation::GetGraphicsStructure() {
  Logger::Trace("Retrieving graphics structures from objects");
  graph_array_.clear();
  for (auto it = species_.begin(); it != species_.end(); ++it) {
    (*it)->Draw(graph_array_);
  }

  /* Visualize interaction forces, crosslinks, etc */
  ix_mgr_.DrawInteractions(graph_array_);
}

/* Initialize output files */
void Simulation::InitOutputs() {
  Logger::Debug("Initializing output files");
  output_mgr_.Init(&params_, &species_, space_.GetStruct());
  // if (!params_.load_checkpoint)
  ix_mgr_.InitOutputs();
  /* If analyzing run time, record cpu time here */
  if (params_.time_analysis) {
    cpu_init_time_ = cpu_time();
  }
}

/* Determine which output files we are reading */
void Simulation::InitInputs(run_options run_opts) {
  output_mgr_.Init(&params_, &species_, space_.GetStruct(), true, &run_opts);
  ix_mgr_.InitOutputs(true, &run_opts);
  /* Initialize object positions from output files if post-processing */
  output_mgr_.ReadInputs();
  ix_mgr_.ReadInputs();
}

/* Write object positions, etc if necessary */
void Simulation::WriteOutputs() {
  output_mgr_.WriteOutputs();
  /* Write interaction information/crosslink positions, etc */
  ix_mgr_.WriteOutputs();
  /* If we are analyzing run time and this is the last step, record final time
   * here. */
  if (params_.time_analysis && params_.i_step == params_.n_steps) {
    double cpu_final_time = cpu_time();
    double cpu_time = cpu_final_time - cpu_init_time_;
    Logger::Info("CPU Time for Initialization: %2.6f", cpu_init_time_);
    Logger::Info("CPU Time: %2.6f", cpu_time);
    Logger::Info("Sim Time: %2.6f", time_);
    Logger::Info("CPU Time/Sim Time: %2.6f\n", cpu_time / time_);
  }
}

/* Run the steps we need for post-processing output files */
void Simulation::ProcessOutputs(YAML::Node sim_params, run_options run_opts) {
  parser_.Init(sim_params);
  InitProcessing(run_opts);
  RunProcessing(run_opts);
  ClearSimulation();
}

// Initialize data structures for post-processing
void Simulation::InitProcessing(run_options run_opts) {
  Logger::Info("Initializing datastructures for post-processing outputs");
  params_ = parser_.ParseSystemParameters();
  run_name_ = params_.run_name;
  rng_ = new RNG(params_.seed);

  // No need to worry about dynamic timestep in post processing
  params_.dynamic_timestep = false;
  // Do not coarse-grain mesh interactions during analysis
  params_.coarse_grained_mesh_interactions = false;
  params_.remove_duplicate_interactions = false;

  /* Ensure that we are not trying to load any checkpoints when processing
     outputs */
  params_.load_checkpoint = 0;

  space_.Init(&params_);
  InitObjects();
  ix_mgr_.Init(&params_, &species_, space_.GetStruct(), true);
  InitSpecies();
  // ix_mgr_.InitInteractions();
  InsertSpecies(true, true);
  InitInputs(run_opts);
  if (run_opts.analysis_flag) {
    for (auto it = species_.begin(); it != species_.end(); ++it) {
      (*it)->InitAnalysis();
      if ((*it)->CheckInteractionAnalysis()) {
        ix_mgr_.SetInteractionAnalysis(true);
      }
    }
  }
  ix_mgr_.InitInteractions();
  ix_mgr_.ResetCellList();
  ix_mgr_.Reset();
  if (run_opts.graphics_flag || run_opts.make_movie) {
    params_.graph_flag = true;
    if (run_opts.use_posits && params_.n_graph < output_mgr_.GetNPosit()) {
      params_.n_graph = output_mgr_.GetNPosit();
    } else if (!run_opts.use_posits &&
               params_.n_graph < output_mgr_.GetNSpec()) {
      params_.n_graph = output_mgr_.GetNSpec();
    }
    if (run_opts.make_movie) {
      params_.movie_flag = true;
    }
    InitGraphics();
  } else {
    params_.graph_flag = false;
  }
}

/* Post-processing on simulation outputs for movie generation, analysis output
 * generation, etc. */
void Simulation::RunProcessing(run_options run_opts) {
  Logger::Info("Processing outputs for %s", run_name_.c_str());
  for (i_step_ = 0; true; ++i_step_) {
    params_.i_step = i_step_;
    time_ = params_.i_step * params_.delta;
    PrintComplete();
    if (early_exit) {
      break;
    }
    if (run_opts.analysis_flag && params_.i_step >= params_.n_steps_equil) {
      bool struct_update = false;
      /* Check if we are running any species analysis to determine whether we
       * run structure analysis */
      for (auto it = species_.begin(); it != species_.end(); ++it) {
        if (((*it)->GetPositFlag() &&
             params_.i_step % (*it)->GetNPosit() == 0) ||
            ((*it)->GetSpecFlag() && params_.i_step % (*it)->GetNSpec() == 0)) {
          struct_update = true;
          break;
        }
      }
      if (struct_update) {
        ix_mgr_.InteractionAnalysis();
        for (auto it = species_.begin(); it != species_.end(); ++it) {
          if (((*it)->GetPositFlag() &&
               params_.i_step % (*it)->GetNPosit() == 0) ||
              ((*it)->GetSpecFlag() &&
               params_.i_step % (*it)->GetNSpec() == 0)) {
            (*it)->RunAnalysis();
          }
        }
      }
    }
    Draw(run_opts.single_frame);
    output_mgr_.ReadInputs();
    ix_mgr_.ReadInputs();
  }
  Draw(run_opts.single_frame);
  early_exit = false;
  Logger::Info("Early exit triggered on step %d. Ending simulation.",
               params_.i_step);
  if (run_opts.analysis_flag && i_step_ > params_.n_steps_equil) {
    for (auto it = species_.begin(); it != species_.end(); ++it) {
      (*it)->FinalizeAnalysis();
    }
  }
}
