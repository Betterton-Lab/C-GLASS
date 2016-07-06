#ifndef _SIMCORE_FORCES_H_
#define _SIMCORE_FORCES_H_

#include "auxiliary.h"
#include "interaction_engine.h"
#include "minimum_distance.h"
#include "potential_manager.h"
#include "particle_tracking.h"
#include "species.h"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

class Forces {
  private:
    int n_dim_,
        n_periodic_,
        max_overlap_,
        draw_flag_,
        nthreads_;
    bool draw_;
    double dr_[3],
           contact1_[3],
           contact2_[3],
           dr_mag_,
           dr_mag2_,
           buffer_mag_,
           buffer_mag2_,
           skin_;
    system_parameters *params_;
    FTYPE force_type_;
    space_struct *space_;
    std::vector<Simple*> simples_; 
    std::vector<SpeciesBase*> *species_;
    PotentialManager potentials_;
    ParticleTracking tracking_;
    InteractionEngine fengine_; //fengine = force engine.  get it?

  public:
    Forces() {}
    ~Forces() {}

    std::vector<graph_struct> draw_array_;
  public:
    void Init(system_parameters *pParams, space_struct *pSpace, std::vector<SpeciesBase*> *pSpecies);
    void InteractMP();
    void DumpAll();
    void InitPotentials();
    void Draw(std::vector<graph_struct*> * graph_array);
};

#endif // _SIMCORE_FORCES_H_
