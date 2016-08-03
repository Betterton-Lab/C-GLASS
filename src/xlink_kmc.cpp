// Implementation for simple binding and unbinding

#include <cassert>

#include "xlink_kmc.h"

#include "xlink.h"
#include "xlink_head.h"
#include "xlink_helpers.h"
#include "br_rod.h"

void XlinkKMC::Init(space_struct *pSpace,
                        ParticleTracking *pTracking,
                        SpeciesBase *spec1,
                        SpeciesBase *spec2,
                        int ikmc,
                        YAML::Node &node,
                        long seed) {
  KMCBase::Init(pSpace, pTracking, spec1, spec2, ikmc, node, seed);

  // Grab our specific claims
  // eps effective
  switch (node["kmc"][ikmc]["concentration_0_1"].Type()) {
    case YAML::NodeType::Scalar:
      eps_eff_0_1_[0] = eps_eff_0_1_[1] =
        0.5 * node["kmc"][ikmc]["concentration_0_1"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      eps_eff_0_1_[0] = node["kmc"][ikmc]["concentration_0_1"][0].as<double>();
      eps_eff_0_1_[1] = node["kmc"][ikmc]["concentration_0_1"][1].as<double>();
      break;
  }
  switch (node["kmc"][ikmc]["concentration_1_2"].Type()) {
    case YAML::NodeType::Scalar:
      eps_eff_1_2_[0] = eps_eff_1_2_[1] =
        0.5 * node["kmc"][ikmc]["concentration_1_2"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      eps_eff_1_2_[0] = node["kmc"][ikmc]["concentration_1_2"][0].as<double>();
      eps_eff_1_2_[1] = node["kmc"][ikmc]["concentration_1_2"][1].as<double>();
      break;
  }

  // on rates
  switch (node["kmc"][ikmc]["on_rate_0_1"].Type()) {
    case YAML::NodeType::Scalar:
      on_rate_0_1_[0] = on_rate_0_1_[1] =
        node["kmc"][ikmc]["on_rate_0_1"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      on_rate_0_1_[0] = node["kmc"][ikmc]["on_rate_0_1"][0].as<double>();
      on_rate_0_1_[1] = node["kmc"][ikmc]["on_rate_0_1"][1].as<double>();
      break;
  }
  switch (node["kmc"][ikmc]["on_rate_1_2"].Type()) {
    case YAML::NodeType::Scalar:
      on_rate_1_2_[0] = on_rate_1_2_[1] =
        node["kmc"][ikmc]["on_rate_1_2"].as<double>();
      break;
    case YAML::NodeType::Sequence:
      on_rate_1_2_[0] = node["kmc"][ikmc]["on_rate_1_2"][0].as<double>();
      on_rate_1_2_[1] = node["kmc"][ikmc]["on_rate_1_2"][1].as<double>();
      break;
  }

  alpha_          = node["kmc"][ikmc]["alpha"].as<double>();
  rcutoff_0_1_    = node["kmc"][ikmc]["rcut"].as<double>();
  velocity_       = node["kmc"][ikmc]["velocity"].as<double>();
  barrier_weight_ = node["kmc"][ikmc]["barrier_weight"].as<double>();
  k_stretch_      = node["kmc"][ikmc]["spring_constant"].as<double>();
  r_equil_        = node["kmc"][ikmc]["equilibrium_length"].as<double>();

  CalcCutoff();

  BuildTables();
}

void XlinkKMC::CalcCutoff() {
  BrRodSpecies *prspec = dynamic_cast<BrRodSpecies*>(spec2_);
  // XXX FIXME set back to the correct value
  //max_length_ = 110;
  max_length_ = prspec->GetMaxLength();
  rcutoff_1_2_ = 0.0;
  const double temp = 1.0;
  const double smalleps = 1E-3;
  double eps_eff = eps_eff_1_2_[0] + eps_eff_1_2_[1];
  double rc_0 = sqrt(2.0 / ( (1-barrier_weight_) * k_stretch_) * temp * log(eps_eff * max_length_ / smalleps * sqrt(2.0 * temp / k_stretch_)));
  rcutoff_1_2_ = r_equil_ + rc_0;
}

double XlinkKMC::XKMCErfinv(double x) {
  // See: A handy approximation for the error function and its inverse
  // (Winitzki 2008) (google it).  This isn't a great approximation, but
  // it will do the trick.  It's not programmed for efficiency since it should
  // only be called during program startup once. If you need better prevision
  // or performance, boost apparently has a version (but then we need boot)
  const double a = 0.147;
  double t1 = -2/M_PI/a;
  double t2 = -log(1-x*x)/2;
  double t3 = 2/M_PI/a + log(1-x*x)/2;
  double t4 = -log(1-x*x)/a;

  return sqrt(t1 + t2 + sqrt(t3*t3 + t4));
}

void XlinkKMC::BuildTables() {
  if (r_equil_ != 0.0) {
    std::vector<double> x[2];
    double bin_size = 0.05;
    double alpha = k_stretch_ * (1 - barrier_weight_) / 2;
    double const smalleps = 1E-5;
    double a_cutoff = 1/sqrt(alpha) * XKMCErfinv(1 - 4.0*sqrt(alpha/M_PI)*smalleps) +
      r_equil_;
    double y_cutoff = rcutoff_1_2_;
    
    xlh::xlink_params params;
    params.alpha = alpha;
    params.r0 = r_equil_;

    for (double a = 0.0; a <= a_cutoff; a += bin_size)
      x[0].push_back(a);
    for (double y0 = 0.0; y0 <= y_cutoff; y0 += bin_size)
      x[1].push_back(y0);

    n_exp_lookup_.Init(2, x, &xlh::prob_1_2, &params);
  }
}

void XlinkKMC::Print() {
  printf("Xlink - BR Rod KMC Module\n");
  KMCBase::Print();
  printf("\t {eps_eff 0 -> 1}: [%2.2f, %2.2f]\n", eps_eff_0_1_[0], eps_eff_0_1_[1]);
  printf("\t {eps_eff 1 -> 2}: [%2.2f, %2.2f]\n", eps_eff_1_2_[0], eps_eff_1_2_[1]);
  printf("\t {on_rate 0 -> 1}: [%2.8f, %2.8f]\n", on_rate_0_1_[0], on_rate_0_1_[1]);
  printf("\t {on_rate 1 -> 2}: [%2.8f, %2.8f]\n", on_rate_1_2_[0], on_rate_1_2_[1]);
  printf("\t {barrier_weight: %2.10f}\n", barrier_weight_);
  printf("\t {equilibrium_length: %2.4f}\n", r_equil_);
  printf("\t {k_spring: %2.4f}\n", k_stretch_);
  printf("\t {rcutoff_0_1: %2.8f}\n", rcutoff_0_1_);
  printf("\t {rcutoff_1_2: %2.8f}\n", rcutoff_1_2_);
  printf("\t {alpha: %2.4f}\n", alpha_);
}

void XlinkKMC::PrepKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  // Prepare each composite particle for the upcoming kmc step
  if (!spec1_->IsKMC()) return;
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  double ntot_0_1 = 0.0;
  double ntot_1_2 = 0.0;
  auto xlinks = pxspec->GetXlinks(); 

  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    // Determine what state we're in so that we can do the appropriate thing
    // XXX Do the rest of this
    switch((*xit)->GetBoundState()) {
      case unbound:
        Update_0_1(*xit);
        ntot_0_1 += (*xit)->GetNExp_0_1();
        break;
      case singly:
        Update_1_2(*xit);
        ntot_1_2 += (*xit)->GetNExp_1_2();
        break;
    }
  }

  pxspec->SetNExp_0_1(ntot_0_1);
  pxspec->SetNExp_1_2(ntot_1_2);
}

void XlinkKMC::Update_0_1(Xlink* xit) {
  double nexp_xlink = 0.0;
  auto heads = xit->GetHeads();

  // Final binding affinity is 2x if the eps_eff and on_rate are the same
  for (int i = 0; i < heads->size(); ++i) {
    auto head = &(*heads)[i];
    double nexp = 0.0;
    double binding_affinity = eps_eff_0_1_[i] * on_rate_0_1_[i] * alpha_ * xit->GetDelta();
    auto idx = (*oid_position_map_)[head->GetOID()];
    for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
      nexp += binding_affinity * nldx->kmc_;
    }
    head->SetNExp_0_1(nexp);
    nexp_xlink += nexp;
  }

  xit->SetNExp_0_1(nexp_xlink);
}

void XlinkKMC::Update_1_2(Xlink *xit) {
  auto heads = xit->GetHeads();
  auto head0 = heads->begin();
  auto head1 = heads->begin()+1;
  int free_i, attc_i;
  XlinkHead *attachedhead;
  XlinkHead *freehead;
  if (head0->GetBound()) {
    attc_i = 0;
    free_i = 1;
    attachedhead = &(*head0);
    freehead = &(*head1);
  } else {
    attc_i = 1;
    free_i = 0;
    attachedhead = &(*head1);
    freehead = &(*head0);
  }
  double binding_affinity = eps_eff_1_2_[free_i] * on_rate_1_2_[free_i];
  auto free_idx = (*oid_position_map_)[freehead->GetOID()];
  auto attc_idx = (*oid_position_map_)[attachedhead->GetOID()];
  // Get the attached rod
  auto attach_info = attachedhead->GetAttach();
  auto attach_info_idx = (*oid_position_map_)[attach_info.first];
  auto mrod_attached = (*simples_)[attach_info_idx];
  if (binding_affinity > 0.0) {
    double n_exp = 0.0;

    // We have to look at all of our neighbors withint the mrcut
    for (auto nldx = neighbors_[free_idx].begin(); nldx != neighbors_[free_idx].end(); ++nldx) {
      auto mrod = (*simples_)[nldx->idx_];
      // Check to see if it's really a rod, and if it's the same one we're already attached to
      if (mrod->GetSID() != sid2_) continue;
      if (mrod->GetRID() == mrod_attached->GetRID()) {
        continue;
      }
      // Calculate center to center displacement
      // XXX FIXME CJE possibly don't do this, and store the minimum distance calculation from
      // earlier point point calculation

      // XXX FIXME Polar Affinity
      double polar_affinity = 1.0;
      double r_x[3];
      double s_x[3];
      double r_rod[3];
      double s_rod[3];
      double u_rod[3];
      std::copy(freehead->GetRigidPosition(), freehead->GetRigidPosition()+ndim_, r_x);
      std::copy(freehead->GetRigidScaledPosition(), freehead->GetRigidScaledPosition()+ndim_, s_x);
      std::copy(mrod->GetRigidPosition(), mrod->GetRigidPosition()+ndim_, r_rod);
      std::copy(mrod->GetRigidScaledPosition(), mrod->GetRigidScaledPosition()+ndim_, s_rod);
      std::copy(mrod->GetRigidOrientation(), mrod->GetRigidOrientation()+ndim_, u_rod);
      double l_rod = mrod->GetRigidLength();
      double rcontact[3];
      double dr[3];
      double mu0 = 0.0;
      min_distance_point_carrier_line(ndim_, nperiodic_,
                                      space_->unit_cell, r_x, s_x,
                                      r_rod, s_rod, u_rod, l_rod,
                                      dr, rcontact, &mu0);

      // Now do the integration over the limits on the MT
      // Check the cutoff distance
      double r_min_mag2 = 0.0;
      for (int i = 0; i < ndim_; ++i) {
        r_min_mag2 += SQR(dr[i]);
      }
      if (r_equil_ == 0.0) {
        double kb = k_stretch_ * (1.0 - barrier_weight_);
        double scale_factor = sqrt(0.5 * kb);
        double lim0 = scale_factor * (-mu0 - 0.5*l_rod);
        double term0 = erf(lim0);
        double lim1 = scale_factor * (-mu0 + 0.5*l_rod);
        double term1 = erf(lim1);

        nldx->kmc_ = binding_affinity * sqrt(M_PI_2 / kb) * exp(-0.5*kb*r_min_mag2) * (term1 - term0) * polar_affinity;
        n_exp += nldx->kmc_;
      } else {
        double lim0 = -mu0 - 0.5 * l_rod;
        double lim1 = -mu0 + 0.5 * l_rod;
        double r_min_mag = sqrt(r_min_mag2);
        double x[2] = {fabs(lim0), r_min_mag};
        double term0 = n_exp_lookup_.Lookup(x) * ((lim0 < 0) ? -1.0 : 1.0);
        x[0] = fabs(lim1);
        double term1 = n_exp_lookup_.Lookup(x) * ((lim1 < 0) ? -1.0 : 1.0);
        // OVERRIDE the kmc_ value of this neighbor list
        nldx->kmc_ = binding_affinity * (term1 - term0) * polar_affinity;
        n_exp += nldx->kmc_;
      }
      if (debug_trace)
        printf("[%d] -> neighbor[%d] {kmc: %2.4f}\n", freehead->GetOID(), mrod->GetOID(), nldx->kmc_);
    } // loop over local neighbors of xlink

    freehead->SetNExp_1_2(n_exp);
    xit->SetNExp_1_2(n_exp);
  }
}

void XlinkKMC::StepKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  // Run the bind unbind
  int g[4] = {0, 1, 2, 3};
  for (int i = 0; i < 4; ++i) {
    int j = gsl_rng_uniform_int(rng_.r, 4);
    int swapme = g[i];
    g[i] = g[j];
    g[j] = swapme;
  }

  if (debug_trace)
    printf("XlinkKMC module %d -> %d -> %d -> %d\n", g[0], g[1], g[2], g[3]);

  for (int i = 0; i < 4; ++i) {
    switch (g[i]) {
      case 0:
        KMC_0_1();
        break;
      case 1:
        KMC_1_0();
        break;
      case 2:
        KMC_1_2();
        break;
      case 3:
        KMC_2_1();
        break;
    }
  }
}

void XlinkKMC::KMC_0_1() {
  // Loop over the xlinks to see who binds
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    // Only take free ones
    if ((*xit)->GetBoundState() != unbound) continue;
    auto nexp = (*xit)->GetNExp_0_1(); // Number expected to bind in this timestep, calculated from Update_0_1
    if (nexp <  std::numeric_limits<double>::epsilon() &&
        nexp > -std::numeric_limits<double>::epsilon()) nexp = 0.0;
    // IF we have some probability to fall onto a neighbor, check it
    if (nexp > 0.0) {
      auto mrng = (*xit)->GetRNG();
      double roll = gsl_rng_uniform(mrng->r);
      if (roll < nexp) {
        int head_type = gsl_rng_uniform(mrng->r) < ((eps_eff_0_1_[1])/(eps_eff_0_1_[0]+eps_eff_0_1_[1]));
        auto heads = (*xit)->GetHeads();
        auto head = heads->begin() + head_type;
        double binding_affinity = (eps_eff_0_1_[0] * on_rate_0_1_[0] + eps_eff_0_1_[1] * on_rate_0_1_[1]) *
          alpha_ * head->GetDelta();
        std::ostringstream kmc_event;
        kmc_event << "[" << (*xit)->GetOID() << "] Successful KMC move {0 -> 1}, {nexp: " << nexp;
        kmc_event << "}, {roll: " << roll << "}, {head: " << head_type << "}";
        WriteEvent(kmc_event.str());
        if (debug_trace)
          printf("%s\n", kmc_event.str().c_str());
        double pos = 0.0;
        int idx = (*oid_position_map_)[head->GetOID()];
        // Search through the neighbors of this head to figure out who we want to bind to
        for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
          auto part2 = (*simples_)[nldx->idx_];
          if (part2->GetSID() != sid2_) continue; // Make sure it's what we want to bind to
          pos += binding_affinity * nldx->kmc_;
          if (pos > roll) {
            if (debug_trace)
              printf("[%d] Attaching to [%d]\n", head->GetOID(), part2->GetOID());

            // Here, we do more complicated stuff.  Calculate the coordinate along
            // the rod s.t. the line vector is perpendicular to the separation vec
            // (closest point along carrier line).  In this frame, the position
            // of the xlink should be gaussian distributed
            double r_x[3];
            double s_x[3];
            double r_rod[3];
            double s_rod[3];
            double u_rod[3];
            std::copy(head->GetRigidPosition(), head->GetRigidPosition()+ndim_, r_x);
            std::copy(head->GetRigidScaledPosition(), head->GetRigidScaledPosition()+ndim_, s_x);
            std::copy(part2->GetRigidPosition(), part2->GetRigidPosition()+ndim_, r_rod);
            std::copy(part2->GetRigidScaledPosition(), part2->GetRigidScaledPosition()+ndim_, s_rod);
            std::copy(part2->GetRigidOrientation(), part2->GetRigidOrientation()+ndim_, u_rod);
            double l_rod = part2->GetRigidLength();
            double rcontact[3];
            double dr[3];
            double mu = 0.0;
            min_distance_point_carrier_line(ndim_, nperiodic_,
                                            space_->unit_cell, r_x, s_x,
                                            r_rod, s_rod, u_rod, l_rod,
                                            dr, rcontact, &mu);

            double r_min[3];
            double r_min_mag2 = 0.0;
            for (int i = 0; i < ndim_; ++i) {
              r_min[i] = -mu * u_rod[i] - dr[i];
              r_min_mag2 += SQR(r_min[i]);
            }
            mrcut2_ = rcutoff_0_1_*rcutoff_0_1_;
            double a = sqrt(mrcut2_ - r_min_mag2);
            //double a = sqrt(1.0 - r_min_mag2); //FIXME is this right for 1.0? or mrcut2?
            if (isnan(a))
              a = 0.0;

            double crosspos = 0.0;
            for (int i = 0; i < 100; ++i) {
              double uroll = gsl_rng_uniform(mrng->r);
              crosspos = (uroll - 0.5)*a + mu + 0.5*l_rod;

              if (crosspos >= 0 && crosspos <= l_rod)
                break;
              if (i == 99) {
                crosspos = -mu + 0.5*l_rod;
                if (crosspos < 0)
                  crosspos = 0;
                else if (crosspos > l_rod)
                  crosspos = l_rod;
              }
            }
            // Attach to the OID of the particle, this is done for dynamic instability to work
            head->Attach(part2->GetOID(), crosspos);
            if (debug_trace)
              printf("\t{mu: %2.4f}, {crosspos: %2.4f}\n", mu, crosspos);
            head->SetBound(true);
            (*xit)->CheckBoundState();
            break;
          } // the actual neighbor to fall on
        } // loop over neighbors to figure out which to attach to
      } // successfully bound
    } // found an xlink with some expectation of binding
  } // loop over all xlinks
}

void XlinkKMC::KMC_1_0() {
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  int nbound1[2];
  std::copy(pxspec->GetNBound1(), pxspec->GetNBound1()+2, nbound1);
  double poff_single_a = on_rate_0_1_[0] * alpha_ * pxspec->GetDelta();
  double poff_single_b = on_rate_0_1_[1] * alpha_ * pxspec->GetDelta();

  int noff[2] = {(int)gsl_ran_binomial(rng_.r, poff_single_a, nbound1[0]),
                 (int)gsl_ran_binomial(rng_.r, poff_single_b, nbound1[1])};
  if (debug_trace)
    printf("[Xlink] {poff_single: (%2.8f, %2.8f)}, {noff: (%d, %d)}\n",
        poff_single_a, poff_single_b, noff[0], noff[1]);

  for (int i = 0; i < (noff[0] + noff[1]); ++i) {
    if (debug_trace)
      printf("[KMC_0_1] detaching trial %d/%d\n", i, (noff[0] + noff[1]));
    int head_type = i < noff[1];
    int idxloc = -1;
    int idxoff = gsl_rng_uniform_int(rng_.r, nbound1[head_type]);

    // Find the one to remove
    for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
      if ((*xit)->GetBoundState() != singly) continue;
      auto heads = (*xit)->GetHeads();
      // Check to increment only in the case we have the correctly
      // bound head
      if (!(*heads)[head_type].GetBound()) continue;
      idxloc++;
      if (idxloc == idxoff) {
        auto head0 = heads->begin();
        auto head1 = heads->begin()+1;

        // Figure out which head attached
        // Do some fancy aliasing to make this easier
        XlinkHead *attachedhead;
        XlinkHead *nonattachead;
        if (head_type == 0) {
          attachedhead = &(*head0);
          nonattachead = &(*head1);
        } else {
          attachedhead = &(*head1);
          nonattachead = &(*head0);
        }

        std::ostringstream kmc_event;
        kmc_event << "[x:" << (*xit)->GetOID() << ",head:" << attachedhead->GetOID() << "] Successful KMC move {1 -> 0}, {idxoff=idxloc=";
        kmc_event << idxloc << "}, {head: " << head_type << "}";
        WriteEvent(kmc_event.str());
        if (debug_trace)
          printf("%s\n", kmc_event.str().c_str());
        // Place withint some random distance of the attach point
        double randr[3];
        double mag2 = 0.0;
        mrcut2_ = rcutoff_0_1_*rcutoff_0_1_;
        auto mrng = attachedhead->GetRNG();
        double prevpos[3];
        std::copy(attachedhead->GetRigidPosition(), attachedhead->GetRigidPosition()+ndim_, prevpos);
        do {
          mag2 = 0.0;
          for (int i = 0; i < ndim_; ++i) {
            double mrand = gsl_rng_uniform(mrng->r);
            randr[i] = 2*rcutoff_0_1_*(mrand - 0.5);
            mag2 += SQR(randr[i]);
          }
        } while(mag2 > mrcut2_);
        // Randomly set position based on randr
        for (int i = 0; i < ndim_; ++i) {
          randr[i] = randr[i] + prevpos[i];
        }
        attachedhead->SetPosition(randr);
        attachedhead->SetPrevPosition(prevpos);
        //attachedhead->UpdatePeriodic();
        nonattachead->SetPosition(randr);
        nonattachead->SetPrevPosition(prevpos);
        //nonattachead->UpdatePeriodic();
        attachedhead->SetBound(false);
        (*xit)->CheckBoundState();

        if (debug_trace) {
          auto attachid = attachedhead->GetAttach();
          auto attachidx = (*oid_position_map_)[attachid.first];
          auto idx = (*oid_position_map_)[attachedhead->GetOID()];
          auto part2 = (*simples_)[attachidx];
          printf("[%d] Detached from [%d] (%2.8f, %2.8f) -> (%2.8f, %2.8f)\n",
             attachedhead->GetOID(), part2->GetOID(),
             prevpos[0], prevpos[1],
             attachedhead->GetRigidPosition()[0], attachedhead->GetRigidPosition()[1]);
        }


        attachedhead->Attach(-1, 0.0);
        break;

      } // found it!
    } // find the one to remove

  } // How many to remove
}

void XlinkKMC::KMC_1_2() {
  XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  double nexp_1_2 = pxspec->GetNExp_1_2();
  if (debug_trace)
    printf("[KMC_1_2] ntot: %2.8f\n", nexp_1_2 * pxspec->GetDelta());
  int nattach = gsl_ran_poisson(rng_.r, nexp_1_2 * pxspec->GetDelta());
  for (int itrial = 0; itrial < nattach; ++itrial) {
    if (debug_trace)
      printf("[KMC_1_2] attaching trial %d/%d\n", itrial, nattach);
    double ran_loc = gsl_rng_uniform(rng_.r) * nexp_1_2;
    double loc = 0.0;
    bool foundidx = false;
    for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
      // Only take singly bound
      if ((*xit)->GetBoundState() != singly) continue;
      loc += (*xit)->GetNExp_1_2();
      if (loc > ran_loc) {
        std::ostringstream kmc_event;
        kmc_event << "[" << (*xit)->GetOID() << "] Successful KMC move {1 -> 2}, {nexp_1_2: " << nexp_1_2 << "}, {ran_loc: ";
        kmc_event << ran_loc << "}, {loc: " << loc << "}";
        WriteEvent(kmc_event.str());
        if (debug_trace)
          printf("%s\n", kmc_event.str().c_str());
        // reset loc back to beginning of this xlink
        loc -= (*xit)->GetNExp_1_2();

        // Look at my neighbors and figure out which to fall on
        // Also, get the heads
        auto heads = (*xit)->GetHeads();
        auto head0 = heads->begin();
        auto head1 = heads->begin()+1;

        // Figure out which head attached
        // Do some fancy aliasing to make this easier
        XlinkHead *attachedhead;
        XlinkHead *nonattachead;
        if (head0->GetBound()) {
          attachedhead = &(*head0);
          nonattachead = &(*head1);
        } else {
          attachedhead = &(*head1);
          nonattachead = &(*head0);
        }
        int idx = (*oid_position_map_)[nonattachead->GetOID()];

        // Loop over my neighbors to figure out which to fall on
        for (auto nldx = neighbors_[idx].begin(); nldx != neighbors_[idx].end(); ++nldx) {
          auto mrod = (*simples_)[nldx->idx_];
          // Make sure it's the right species to attach to
          if (mrod->GetSID() != sid2_) continue;

          // Ignore our own rod
          //auto free_idx = (*oid_position_map_)[nonattachead->GetOID()];
          auto free_idx = idx;
          auto attc_idx = (*oid_position_map_)[attachedhead->GetOID()];
          // Get the attached rod
          auto attachinfo = attachedhead->GetAttach();
          auto attachinfoidx = (*oid_position_map_)[attachinfo.first];
          auto mrod_attached = (*simples_)[attachinfoidx];
          if (mrod->GetRID() == mrod_attached->GetRID()) {
            continue;
          }

          loc += nldx->kmc_;
          if (loc > ran_loc) {
            // Found the neighbor to fall onto!!!!
            if (debug_trace) {
              printf("[%d] Attaching to [%d] {loc: %2.4f}\n", nonattachead->GetOID(), mrod->GetOID(), loc);
            }

            // What location (crosspos) do we fall onto?
            double crosspos = 0.0;

            // Here, we do more complicated stuff.  Calculate the coordinate along
            // the rod s.t. the line vector is perpendicular to the separation vec
            // (closest point along carrier line).  In this frame, the position
            // of the xlink should be gaussian distributed
            double r_x[3];
            double s_x[3];
            double r_rod[3];
            double s_rod[3];
            double u_rod[3];
            std::copy(nonattachead->GetRigidPosition(), nonattachead->GetRigidPosition()+ndim_, r_x);
            std::copy(nonattachead->GetRigidScaledPosition(), nonattachead->GetRigidScaledPosition()+ndim_, s_x);
            std::copy(mrod->GetRigidPosition(), mrod->GetRigidPosition()+ndim_, r_rod);
            std::copy(mrod->GetRigidScaledPosition(), mrod->GetRigidScaledPosition()+ndim_, s_rod);
            std::copy(mrod->GetRigidOrientation(), mrod->GetRigidOrientation()+ndim_, u_rod);
            double l_rod = mrod->GetRigidLength();
            double rcontact[3];
            double dr[3];
            double mu = 0.0;
            min_distance_point_carrier_line(ndim_, nperiodic_,
                                            space_->unit_cell, r_x, s_x,
                                            r_rod, s_rod, u_rod, l_rod,
                                            dr, rcontact, &mu);

            // Find it
            auto mrng = (*xit)->GetRNG();
            double y02 = 0.0;
            for (int i = 0; i < ndim_; ++i) {
              y02 += SQR(dr[i]);
            }
            int itrial_loc = 0;
            do {
              itrial_loc++;
              double uroll = gsl_rng_uniform(mrng->r);
              double xvec[2] = {0.0, sqrt(y02)};
              double mpos = ((gsl_rng_uniform(mrng->r) < 0.5) ? -1.0 : 1.0) *
                  n_exp_lookup_.Invert(0, uroll, xvec) + mu + 0.5 * l_rod;
              if (mpos >= 0 && mpos <= l_rod) {
                crosspos = mpos;
                break;
              }
            } while (itrial_loc < 100);

            nonattachead->Attach(mrod->GetOID(), crosspos);
            nonattachead->SetBound(true);
            (*xit)->CheckBoundState();

            foundidx = true;
            break;
          } // got the neighbor to fall onto
        } // check neighbors to see if we need to fall onto this one

      } // found it!

      if (foundidx)
        break;
    } // loop over all xlinks for itrial
  }
}

void XlinkKMC::KMC_2_1() {
  //printf("XlinkKMC::KMC_2_1 begin\n");
  if (barrier_weight_ == 0.0) {
    // All detachments equally probable, do via a poisson distribution
    printf("NOT IMPLEMENTED, EXITING!\n");
    exit(1);
  } else {
    KMC_2_1_ForceDep();
  }
}

void XlinkKMC::KMC_2_1_ForceDep() {
  // Force dependent detachment!
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  int nbound2[2];
  std::copy(pxspec->GetNBound2(), pxspec->GetNBound2()+2, nbound2);
  double poff_base[2] = {
    on_rate_1_2_[0] * pxspec->GetDelta(),
    on_rate_1_2_[1] * pxspec->GetDelta()};
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    if ((*xit)->GetBoundState() != doubly) continue; // only doubly bound
    double kboltzoff = exp(barrier_weight_ * (*xit)->GetInternalU());
    //printf("kboltzoff: %2.8f\n", kboltzoff);
    uint8_t off[2] = {
      gsl_rng_uniform(rng_.r) < poff_base[0] * kboltzoff,
      gsl_rng_uniform(rng_.r) < poff_base[1] * kboltzoff};
    if (off[0] && off[1]) {
      //printf("KMC_1_2 double head removal\n");
      // Easy, set the head to the midpoint of the xlink and fall off both
      //exit(1);
    } else if (off[0] || off[1]) {
      //printf("KMC_1_2 single head removal\n");
      // Randomly insert this head aroudn the detachment area
      //exit(1);
    }
  }
}

void XlinkKMC::UpdateKMC() {
  simples_ = tracking_->GetSimples();
  nsimples_ = tracking_->GetNSimples();
  oid_position_map_ = tracking_->GetOIDPositionMap();
  neighbors_ = tracking_->GetNeighbors();

  nfree_ = 0.0;
  nbound1_[0] = nbound1_[1] = 0.0;
  nbound2_[0] = nbound2_[1] = 0.0;

  // Do a switch on the type that we're examining
  // Loop over the xlinks to see who  does what
  XlinkSpecies* pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  auto xlinks = pxspec->GetXlinks();
  for (auto xit = xlinks->begin(); xit != xlinks->end(); ++xit) {
    switch((*xit)->GetBoundState()) {
      case unbound:
        nfree_++;
        break;
      case singly:
        UpdateStage1(*xit);
        break;
      case doubly:
        UpdateStage2(*xit);
        ApplyStage2Force(*xit);
        break;
    }
  }

  pxspec->SetNFree(nfree_);
  pxspec->SetNBound1(nbound1_[0], nbound1_[1]);
  pxspec->SetNBound2(nbound2_[0], nbound2_[0]);
}

void XlinkKMC::UpdateStage1(Xlink *xit) {

  // Set nexp to zero for all involved
  xit->SetNExp_0_1(0.0);
  auto heads = xit->GetHeads();
  auto head0 = heads->begin();
  auto head1 = heads->begin()+1;
  head0->SetNExp_0_1(0.0);
  head1->SetNExp_0_1(0.0);
  
  // Figure out which head attached
  // Do some fancy aliasing to make this easier
  XlinkHead *attachedhead;
  XlinkHead *nonattachead;
  if (head0->GetBound()) {
    attachedhead = &(*head0);
    nonattachead = &(*head1);
    nbound1_[0]++;
  } else if (head1->GetBound()) {
    attachedhead = &(*head1);
    nonattachead = &(*head0);
    nbound1_[1]++;
  } else {
    printf("Something has gone horribly wrong\n");
    exit(1);
  }
  auto aid = attachedhead->GetAttach().first;
  auto aidx = (*oid_position_map_)[aid];
  auto cross_pos = attachedhead->GetAttach().second; // relative to the -end of the rod!!
  double rx[3];
  std::copy(attachedhead->GetRigidPosition(), attachedhead->GetRigidPosition()+ndim_, rx);
  
  auto part2 = (*simples_)[aidx];
  auto r_rod = part2->GetRigidPosition();
  auto u_rod = part2->GetRigidOrientation();
  auto l_rod = part2->GetRigidLength();

  // If we are moving with some velocity, do that
  // XXX FIXME check for end pausing
  cross_pos += velocity_ * attachedhead->GetDelta();
  if (cross_pos > l_rod) {
    cross_pos = l_rod;
  } else if (cross_pos < 0.0) {
    cross_pos = 0.0;
  }
  attachedhead->Attach(aid, cross_pos);

  double rxnew[3];
  for (int i = 0; i < ndim_; ++i) {
    rxnew[i] = r_rod[i] - 0.5 * u_rod[i] * l_rod + cross_pos * u_rod[i];
  }
  if (debug_trace)
    printf("[%d] attached [%d], (%2.4f, %2.4f) -> setting -> {%2.4f}(%2.4f, %2.4f)\n",
           attachedhead->GetOID(), part2->GetOID(), rx[0], rx[1], cross_pos,
           rxnew[0], rxnew[1]);
  attachedhead->SetPrevPosition(rx);
  attachedhead->SetPosition(rxnew);
  //attachedhead->UpdatePeriodic();
  attachedhead->AddDr();
  nonattachead->SetPrevPosition(rx);
  nonattachead->SetPosition(rxnew);
  //attachedhead->UpdatePeriodic();
  nonattachead->AddDr();
  xit->SetPrevPosition(rx);
  xit->SetPosition(rxnew);
  //xit->UpdatePeriodic();
}

void XlinkKMC::UpdateStage2(Xlink *xit) {
  // Set nexp to zero for all involved
  xit->SetNExp_1_2(0.0);
  double oldxitpos[3];
  std::copy(xit->GetPosition(), xit->GetPosition()+ndim_, oldxitpos);
  double avgpos[3] = {0.0, 0.0, 0.0};
  auto heads = xit->GetHeads();
  int ihead = -1;
  for (auto head = heads->begin(); head != heads->end(); ++head) {
    ihead++;
    nbound2_[ihead]++;
    head->SetNExp_1_2(0.0);

    auto aid = head->GetAttach().first;
    auto aidx = (*oid_position_map_)[aid];
    auto crosspos = head->GetAttach().second;
    auto rx = head->GetRigidPosition();

    auto mrod = (*simples_)[aidx];
    auto rrod = mrod->GetRigidPosition();
    auto urod = mrod->GetRigidOrientation();
    auto lrod = mrod->GetRigidLength();

    // If we are moving at some velocity, do it
    // XXX Check for end pausing
    crosspos += velocity_ * head->GetDelta();
    if (crosspos > lrod) {
      crosspos = lrod;
    } else if (crosspos < 0.0) {
      crosspos = 0.0;
    }
    head->Attach(aid, crosspos);

    double rxnew[3];
    for (int i = 0; i < ndim_; ++i) {
      rxnew[i] = rrod[i] - 0.5 * urod[i] * lrod + crosspos * urod[i];
      avgpos[i] += rxnew[i];
    }
    if (debug_trace) {
      printf("[%d] attached [%d], (%2.4f, %2.4f) -> setting -> {%2.4f}(%2.4f, %2.4f)\n",
          head->GetOID(), mrod->GetOID(), rx[0], rx[1],
          crosspos,
          rxnew[0], rxnew[1]);
    }
    head->SetPrevPosition(rx);
    head->SetPosition(rxnew);
    //head->UpdatePeriodic();
    head->AddDr();
  }

  for (int i = 0; i < ndim_; ++i) {
    avgpos[i] *= 0.5;
  }

  xit->SetPosition(avgpos);
  xit->SetPrevPosition(oldxitpos);
  //xit->UpdatePeriodic();
}

void XlinkKMC::ApplyStage2Force(Xlink *xit) {
  auto heads = xit->GetHeads();
  auto head0 = heads->begin();
  auto head1 = heads->begin()+1;

  auto aid0 = head0->GetAttach().first;
  auto aid1 = head1->GetAttach().first;
  auto aidx0 = (*oid_position_map_)[aid0];
  auto aidx1 = (*oid_position_map_)[aid1];
  auto mrod0 = (*simples_)[aidx0];
  auto mrod1 = (*simples_)[aidx1];
  auto rx0 = head0->GetRigidPosition();
  auto rx1 = head1->GetRigidPosition();
  auto sx0 = head0->GetScaledPosition();
  auto sx1 = head1->GetScaledPosition();

  if (debug_trace)
    printf("[%d:%d] <-> [%d:%d] Applying 2stage force\n", head0->GetOID(), mrod0->GetOID(), head1->GetOID(), mrod1->GetOID());

  double dr[3];
  separation_vector(ndim_, nperiodic_, rx0, sx0, rx1, sx1, space_->unit_cell, dr);
  double rmag2 = 0.0;
  for (int i = 0; i < ndim_; ++i) {
    rmag2 += SQR(dr[i]);
  }
  double rmag = sqrt(rmag2);
  double k = k_stretch_;
  double factor;
  if (r_equil_ == 0.0) {
    factor = k;
  } else {
    factor = k * (1.0 - r_equil_ / sqrt(dot_product(ndim_, dr, dr)));
  }
  double u = 0.0;
  double flink[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < ndim_; ++i) {
    flink[i] = factor * dr[i];
    u += 0.5 * SQR(flink[i]);
  }
  u *= 0.5 / k;
  if (debug_trace)
    printf("\t{uin: %2.4f}\n", u);
  xit->AddPotential(u);
  xit->SetInternalU(u);

  auto rrod0 = mrod0->GetRigidPosition();
  auto rrod1 = mrod1->GetRigidPosition();
  auto crosspos0 = head0->GetAttach().second;
  auto crosspos1 = head1->GetAttach().second;
  auto lrod0 = mrod0->GetRigidLength();
  auto lrod1 = mrod1->GetRigidLength();
  auto urod0 = mrod0->GetRigidOrientation();
  auto urod1 = mrod1->GetRigidOrientation();

  double drmag = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
  /*printf("{rrod0: (%2.4f, %2.4f)}, {rrod1: (%2.4f, %2.4f)}\n", rrod0[0], rrod0[1], rrod1[0], rrod1[1]);
  printf("{urod0: (%2.4f, %2.4f)}, {urod1: (%2.4f, %2.4f)}\n", urod0[0], urod0[1], urod1[0], urod1[1]);
  printf("{rx0: (%2.4f, %2.4f)}, {rx1: (%2.4f, %2.4f)}\n", rx0[0], rx0[1], rx1[0], rx1[1]);
  printf("{dr: (%2.4f, %2.4f)}, {drmag: %2.8f}\n", dr[0], dr[1], drmag);
  printf("{u: %2.4f}, {flink: (%2.4f, %2.4f)}\n", u, flink[0], flink[1]);*/

  // Now, simply add the forces/torques onto the two bonds
  double fbond0[3] = {0.0, 0.0, 0.0};
  double fbond1[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < ndim_; ++i) {
    fbond0[i] += flink[i];
    fbond1[i] -= flink[i];
  }

  double lambda = crosspos0 - 0.5 * lrod0;
  double mu     = crosspos1 - 0.5 * lrod1;

  double rcontact_i[3] = {0.0, 0.0, 0.0};
  double rcontact_j[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < ndim_; ++i) {
    rcontact_i[i] = urod0[i] * lambda;
    rcontact_j[i] = urod1[i] * mu;
  }

  double tau[3];
  double taubond0[3] = {0.0, 0.0, 0.0};
  double taubond1[3] = {0.0, 0.0, 0.0};
  cross_product(rcontact_i, flink, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
    taubond0[i] += tau[i];
  }
  cross_product(rcontact_j, flink, tau, ndim_);
  for (int i = 0; i < 3; ++i) {
    taubond1[i] -= tau[i];
  }

  // Apply the force
  head0->AddForce(fbond0);
  head0->AddTorque(taubond0);
  mrod0->AddForce(fbond0);
  mrod0->AddTorque(taubond0);
  head1->AddForce(fbond1);
  head1->AddTorque(taubond1);
  mrod1->AddForce(fbond1);
  mrod1->AddTorque(taubond1);
}

void XlinkKMC::Dump() {
  // print out the information appropriate to kmc
  if (debug_trace) {
    XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
    printf("XlinkKMC -> dump\n");
    printf("\t{n_exp_0_1: %2.4f, n_exp_1_2: %2.4f}\n", pxspec->GetNExp_0_1(), pxspec->GetNExp_1_2());
    printf("\t{nfree:  %d}\n", pxspec->GetNFree());
    printf("\t{nbound1: %d,%d}\n", pxspec->GetNBound1()[0], pxspec->GetNBound1()[1]);
    printf("\t{nbound2: %d,%d}\n", pxspec->GetNBound2()[0], pxspec->GetNBound2()[1]);
    pxspec->DumpKMC();
  }
}

void XlinkKMC::PrepOutputs() {
  kmc_file_name_ << "sc-kmc-XlinkKMC.log";
  kmc_file.open(kmc_file_name_.str().c_str(), std::ios_base::out);
  kmc_file << "step #ntot #nfree #nbound1[0] #nbound1[1] #nbound2[0] #nbound2[1] nexp01 nexp12\n";
  kmc_file.close();
}

void XlinkKMC::WriteEvent(const std::string &pString) {
  kmc_file.open(kmc_file_name_.str().c_str(), std::ios_base::out | std::ios_base::app);
  kmc_file << pString << std::endl;
  kmc_file.close();
}

void XlinkKMC::WriteOutputs(int istep) {
  XlinkSpecies *pxspec = dynamic_cast<XlinkSpecies*>(spec1_);
  kmc_file.open(kmc_file_name_.str().c_str(), std::ios_base::out | std::ios_base::app);
  kmc_file.precision(16);
  kmc_file.setf(std::ios::fixed, std::ios::floatfield);
  kmc_file << istep << " ";
  kmc_file << pxspec->GetNMembers() << " " << pxspec->GetNFree() << " " <<
    pxspec->GetNBound1()[0] << " " << pxspec->GetNBound1()[1] << " " << 
    pxspec->GetNBound2()[0] << " " << pxspec->GetNBound2()[1] << " ";
  kmc_file << pxspec->GetNExp_0_1() << " "  << pxspec->GetNExp_1_2() << "\n";
  kmc_file.close();
}