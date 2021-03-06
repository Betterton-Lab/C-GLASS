#ifndef _CGLASS_R2_POTENTIAL_H_
#define _CGLASS_R2_POTENTIAL_H_

#include "auxiliary.hpp"
#include "interaction.hpp"
#include "potential_base.hpp"

class R2Potential : public PotentialBase {
 protected:
  double eps_, sigma_;

 public:
  R2Potential() {}
  void CalcPotential(Interaction &ix) {
    double rmag = ix.dr_mag2;
    double *dr = ix.dr;
    double ffac = 0;
    double rinv2 = 0;
    rmag = sqrt(rmag);
    double rinv = 1.0 / (rmag);
    rinv2 = rinv * rinv;
    double rinv3 = rinv2 * rinv;
    ffac = -(2.0 * rinv3);
    // Cut off the force at fcut
    if (ABS(ffac) > max_force_) {
      MaxForceViolation();
      ffac = SIGNOF(ffac) * max_force_;
    }
    for (int i = 0; i < n_dim_; ++i) {
      ix.force[i] = ffac * dr[i] / rmag;
    }
    for (int i = 0; i < n_dim_; ++i)
      for (int j = 0; j < n_dim_; ++j)
        ix.stress[n_dim_ * i + j] = -dr[i] * ix.force[j];
    ix.pote = rinv2 - eps_;
  }

  void InitPotentialParams(system_parameters *params) {
    // Initialize potential params
    sigma_ = params->wca_sig;

    // For R2Potential potentials, the rcutoff is
    // restricted to be at 2^(1/6)sigma
    rcut_ = pow(2.0, 1.0 / 6.0) * sigma_;
    rcut2_ = rcut_ * rcut_;
    eps_ = 1.0 / (rcut2_);
  }
};

#endif
