#include "cglass/minimum_distance.hpp"

#define SMALL 1.0e-12

int MinimumDistance::n_dim_ = 0;
int MinimumDistance::n_periodic_ = 0;
double *MinimumDistance::unit_cell_ = nullptr;
double MinimumDistance::boundary_cut2_ = 0;
space_struct *MinimumDistance::space_ = nullptr;

void MinimumDistance::Init(space_struct *space, double boundary_cutoff_sq) {
  space_ = space;
  n_dim_ = space_->n_dim;
  n_periodic_ = space_->n_periodic;
  unit_cell_ = space_->unit_cell;
  boundary_cut2_ = boundary_cutoff_sq;
}

/* Find the minimum distance between two particles */
void MinimumDistance::ObjectObject(Interaction &ix) {
  double const *const r1 = ix.obj1->GetInteractorPosition();
  double const *const s1 = ix.obj1->GetInteractorScaledPosition();
  double const *const u1 = ix.obj1->GetInteractorOrientation();
  double const l1 = ix.obj1->GetInteractorLength();
  double const d1 = ix.obj1->GetInteractorDiameter();
  double const *const r2 = ix.obj2->GetInteractorPosition();
  double const *const s2 = ix.obj2->GetInteractorScaledPosition();
  double const *const u2 = ix.obj2->GetInteractorOrientation();
  double const l2 = ix.obj2->GetInteractorLength();
  double const d2 = ix.obj2->GetInteractorDiameter();
  ix.dr_mag2 = 0;
  std::fill(ix.dr, ix.dr + 3, 0.0);
  std::fill(ix.contact1, ix.contact1 + 3, 0.0);
  std::fill(ix.contact2, ix.contact2 + 3, 0.0);
  ix.buffer_mag = 0.5 * (d1 + d2);
  ix.buffer_mag2 = ix.buffer_mag * ix.buffer_mag;
  /* TODO: Right now, we can only find minimum distances between
     point-like particles, and line-like particles. Minimum distance
     algorithms are written for planes, so at some point we can
     incorporate interactions between sides of 3D polygons, etc. */

  if (l1 == 0 && l2 == 0) {
    /* When we have two point-like particles interacting. */
    PointPoint(r1, s1, r2, s2, ix.dr, &ix.dr_mag2, ix.midpoint);
  } else if (l1 == 0 && l2 > 0) {
    /* The case where obj1 is a point-like particle and obj2 is an
       extended, line-like particle */
    SphereSphero(r1, s1, r2, s2, u2, l2, ix.dr, &ix.dr_mag2, ix.contact2);
  } else if (l1 > 0 && l2 == 0) {
    /* Same, but switching the order of obj1 and obj2, so we'll just swap
       the order of obj1 and obj2 in the min distance calculation, then
       reverse the direction of the min distance vector and proceed. */
    SphereSphero(r2, s2, r1, s1, u1, l1, ix.dr, &ix.dr_mag2, ix.contact1);
    for (int i = 0; i < 3; ++i) {
      ix.dr[i] = -ix.dr[i];
    }
  } else if (l1 > 0 && l2 > 0) {
    /* When we have two extended, line-like particles interacting. */
    Sphero(r1, s1, u1, l1, r2, s2, u2, l2, ix.dr, &ix.dr_mag2, ix.contact1,
           ix.contact2);
  }
#ifdef TRACE
  Logger::Trace("Minimum distance between %d and %d is %2.4f",
                ix.obj1->GetOID(), ix.obj2->GetOID(), sqrt(ix.dr_mag2));
#endif
}

/* Returns squared minimum distance (dr_mag2) and minimum distance vector (dr)
 * between two point-like objects centered at r1 and r2 (scaled position of s1
 * and s2 in periodic subspace) */
void MinimumDistance::PointPoint(double const *const r1, double const *const s1,
                                 double const *const r2, double const *const s2,
                                 double *dr, double *dr_mag2,
                                 double *midpoint) {
  // First handle periodic subspace
  double ds[3], mp[3];
  for (int i = 0; i < n_periodic_; ++i) {
    ds[i] = s2[i] - s1[i];
    ds[i] -= NINT(ds[i]);
    mp[i] = s1[i] + 0.5 * ds[i];
    mp[i] -= NINT(mp[i]);
  }
  for (int i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    midpoint[i] = 0.0;
    for (int j = 0; j < n_periodic_; ++j) {
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
      midpoint[i] += unit_cell_[n_dim_ * i + j] * mp[j];
    }
  }
  // Then handle free subspace
  for (int i = n_periodic_; i < n_dim_; ++i) {
    dr[i] = r2[i] - r1[i];
    midpoint[i] = r1[i] + 0.5 * dr[i];
  }
  *dr_mag2 = 0.0;
  for (int i = 0; i < n_dim_; ++i) {
    *dr_mag2 += SQR(dr[i]);
  }
  return;
}

/* Routine to calculate minimum distance between a point and a line of finite
length

output: vector that points from point to line along minimum distance between
        point and line (dr)
        distance from r_line along u_line that indicates point of minimum
        distance (mu) */

void MinimumDistance::PointCarrierLine(double *r_point, double *s_point,
                                       double *r_line, double *s_line,
                                       double *u_line, double length,
                                       double *dr, double *r_contact,
                                       double *mu_ret) {
  int i, j;
  double ds[3], mu;

  /* Compute pair separation vector. */
  /* First handle periodic subspace */
  for (i = 0; i < n_periodic_; ++i) {
    ds[i] = s_line[i] - s_point[i];
    ds[i] -= NINT(ds[i]);
  }
  for (i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    for (j = 0; j < n_periodic_; ++j) {
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
    }
  }
  /* Handle free subspace */
  for (i = n_periodic_; i < n_dim_; ++i) {
    dr[i] = r_line[i] - r_point[i];
  }

  mu = -dot_product(n_dim_, dr, u_line);
  double mu_mag = ABS(mu);
  // Now take into account that the line is finite length
  if (mu_mag > 0.5 * length) {
    mu = SIGNOF(mu) * 0.5 * length;
  }

  for (i = 0; i < n_dim_; ++i) {
    r_contact[i] = mu * u_line[i];
    dr[i] = r_line[i] + r_contact[i] - r_point[i];
  }
  *mu_ret = mu;
}

/* Routine to calculate minimum distance between a point and a line of infinite
length

output: vector that points from point to line along minimum distance between
        point and line (dr)
        distance from r_line along u_line that indicates point of minimum
         distance (mu) */
void MinimumDistance::PointCarrierLineInf(double *r_point, double *s_point,
                                          double *r_line, double *s_line,
                                          double *u_line, double length,
                                          double *dr, double *mu) {
  int i, j;
  double ds[3];

  /* Compute pair separation vector. */
  /* First handle periodic subspace. */
  for (i = 0; i < n_periodic_; ++i) {
    ds[i] = s_line[i] - s_point[i];
    ds[i] -= NINT(ds[i]);
  }
  for (i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    for (j = 0; j < n_periodic_; ++j) {
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
    }
  }
  /* Then handle free subspace. */
  for (i = n_periodic_; i < n_dim_; ++i) {
    dr[i] = r_line[i] - r_point[i];
  }

  *mu = -dot_product(n_dim_, dr, u_line);
}

/* Routine to calculate minimum distance between two spherocylinders, for any
number of spatial dimensions and any type of boundary conditions (free,
periodic, or mixed).

input: real position of first spherocylinder (r_1)
       scaled position of first spherocylinder (s_1)
       director of first spherocylinder (u_1)
       length of first spherocylinder (length_1)
       real position of second spherocylinder (r_2)
       scaled position of second spherocylinder (s_2)
       director of second spherocylinder (u_2)
       length of second spherocylinder (length_2)

output: minimimum separation vector (r_min)
        pointer to squared minimum separation (r_min_mag2)
        vector separating r_1 to point of contact on first sphero (contact1)
        vector separating r_2 to point of contact on second sphero (contact2) */

void MinimumDistance::Sphero(double const *const r_1, double const *const s_1,
                             double const *const u_1, double const length_1,
                             double const *const r_2, double const *const s_2,
                             double const *const u_2, double const length_2,
                             double *r_min, double *r_min_mag2,
                             double *contact_1, double *contact_2) {
  double dr[3], ds[3];

  /* Compute half-length of objects. */
  double half_length_1 = 0.5 * length_1;
  double half_length_2 = 0.5 * length_2;

  /* Compute pair separation vector. */
  /* First handle periodic subspace. */
  for (int i = 0; i < n_periodic_; ++i) {
    ds[i] = s_2[i] - s_1[i];
    ds[i] -= NINT(ds[i]);
  }
  for (int i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    for (int j = 0; j < n_periodic_; ++j) {
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
    }
  }
  /* Then handle free subspace. */
  for (int i = n_periodic_; i < n_dim_; ++i) {
    dr[i] = r_2[i] - r_1[i];
  }

  /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
     First consider two infinitely long lines. */
  double dr_dot_u_1, dr_dot_u_2, u_1_dot_u_2;
  dr_dot_u_1 = dr_dot_u_2 = u_1_dot_u_2 = 0.0;
  for (int i = 0; i < n_dim_; ++i) {
    dr_dot_u_1 += dr[i] * u_1[i];
    dr_dot_u_2 += dr[i] * u_2[i];
    u_1_dot_u_2 += u_1[i] * u_2[i];
  }
  double lambda, mu;
  double denom = 1.0 - SQR(u_1_dot_u_2);
  if (denom < SMALL) {
    lambda = dr_dot_u_1 / 2.0;
    mu = -dr_dot_u_2 / 2.0;
  } else {
    lambda = (dr_dot_u_1 - u_1_dot_u_2 * dr_dot_u_2) / denom;
    mu = (-dr_dot_u_2 + u_1_dot_u_2 * dr_dot_u_1) / denom;
  }
  double lambda_mag = ABS(lambda);
  double mu_mag = ABS(mu);

  /* Now take into account the fact that the two line segments are of finite
   * length. */
  double lambda_a, lambda_b, mu_a, mu_b, r_min_mag2_a, r_min_mag2_b;
  double r_min_a[3], r_min_b[3];
  if (lambda_mag > half_length_1 && mu_mag > half_length_2) {
    /* Calculate first possible case. */
    lambda_a = SIGN(half_length_1, lambda);
    mu_a = -dr_dot_u_2 + lambda_a * u_1_dot_u_2;
    mu_mag = ABS(mu_a);
    if (mu_mag > half_length_2) {
      mu_a = SIGN(half_length_2, mu_a);
    }

    /* Calculate minimum distance between two spherocylinders. */
    r_min_mag2_a = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_min_a[i] = dr[i] - lambda_a * u_1[i] + mu_a * u_2[i];
      r_min_mag2_a += SQR(r_min_a[i]);
    }

    /* Calculate second possible case. */
    mu_b = SIGN(half_length_2, mu);
    lambda_b = dr_dot_u_1 + mu_b * u_1_dot_u_2;
    lambda_mag = ABS(lambda_b);
    if (lambda_mag > half_length_1) {
      lambda_b = SIGN(half_length_1, lambda_b);
    }

    /* Calculate minimum distance between two spherocylinders. */
    r_min_mag2_b = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_min_b[i] = dr[i] - lambda_b * u_1[i] + mu_b * u_2[i];
      r_min_mag2_b += SQR(r_min_b[i]);
    }

    /* Choose the minimum minimum distance. */
    if (r_min_mag2_a < r_min_mag2_b) {
      lambda = lambda_a;
      mu = mu_a;
      *r_min_mag2 = r_min_mag2_a;
      for (int i = 0; i < n_dim_; ++i) {
        r_min[i] = r_min_a[i];
      }
    } else {
      lambda = lambda_b;
      mu = mu_b;
      *r_min_mag2 = r_min_mag2_b;
      for (int i = 0; i < n_dim_; ++i) {
        r_min[i] = r_min_b[i];
      }
    }
  } else if (lambda_mag > half_length_1) {
    /* Adjust lambda and mu. */
    lambda = SIGN(half_length_1, lambda);
    mu = -dr_dot_u_2 + lambda * u_1_dot_u_2;
    mu_mag = ABS(mu);
    if (mu_mag > half_length_2)
      mu = SIGN(half_length_2, mu);

    /* Calculate minimum distance between two spherocylinders. */
    *r_min_mag2 = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
      *r_min_mag2 += SQR(r_min[i]);
    }
  } else if (mu_mag > half_length_2) {
    /* Adjust lambda and mu. */
    mu = SIGN(half_length_2, mu);
    lambda = dr_dot_u_1 + mu * u_1_dot_u_2;
    lambda_mag = ABS(lambda);
    if (lambda_mag > half_length_1)
      lambda = SIGN(half_length_1, lambda);

    /* Calculate minimum distance between two spherocylinders. */
    *r_min_mag2 = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
      *r_min_mag2 += SQR(r_min[i]);
    }
  } else {
    /* Calculate minimum distance between two spherocylinders. */
    *r_min_mag2 = 0.0;
    for (int i = 0; i < n_dim_; ++i) {
      r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
      *r_min_mag2 += SQR(r_min[i]);
    }
  }
  for (int i = 0; i < n_dim_; ++i) {
    contact_1[i] = lambda * u_1[i];
    contact_2[i] = mu * u_2[i];
  }
  return;
}

/* Routine to calculate minimum distance between two spherocylinders and
   center to center separation vector, for any number of
   spatial dimensions and any type of boundary conditions (free, periodic, or
mixed).

input: real position of first spherocylinder (r_1)
       scaled position of first spherocylinder (s_1)
       director of first spherocylinder (u_1)
       length of first spherocylinder (length_1)
       real position of second spherocylinder (r_2)
       scaled position of second spherocylinder (s_2)
       director of second spherocylinder (u_2)
       length of second spherocylinder (length_2)

output: center to center separation vector (dr)
        minimimum separation vector (r_min)
        pointer to squared minimum separation (r_min_mag2)
        pointer to intersection of r_min with axis of first spherocylinder
(lambda) pointer to intersection of r_min with axis of second spherocylinder
(mu). */

void MinimumDistance::SpheroDr(double *r_1, double *s_1, double *u_1,
                               double length_1, double *r_2, double *s_2,
                               double *u_2, double length_2, double *dr,
                               double *r_min, double *r_min_mag2,
                               double *lambda, double *mu) {
  int i, j;
  double half_length_1, half_length_2, dr_dot_u_1, dr_dot_u_2, u_1_dot_u_2,
      denom, lambda_a, lambda_b, mu_a, mu_b, lambda_mag, mu_mag, r_min_mag2_a,
      r_min_mag2_b;
  double ds[3], r_min_a[3], r_min_b[3];

  /* Compute various constants. */
  half_length_1 = 0.5 * length_1;
  half_length_2 = 0.5 * length_2;

  /* Compute pair separation vector. */
  for (i = 0; i < n_periodic_; ++i) { /* First handle periodic subspace. */
    ds[i] = s_2[i] - s_1[i];
    ds[i] -= NINT(ds[i]);
  }
  for (i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    for (int j = 0; j < n_periodic_; ++j)
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
  }
  for (i = n_periodic_; i < n_dim_; ++i) /* Then handle free subspace. */
    dr[i] = r_2[i] - r_1[i];

  /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
     First consider two infinitely long lines. */
  dr_dot_u_1 = dr_dot_u_2 = u_1_dot_u_2 = 0.0;
  for (i = 0; i < n_dim_; ++i) {
    dr_dot_u_1 += dr[i] * u_1[i];
    dr_dot_u_2 += dr[i] * u_2[i];
    u_1_dot_u_2 += u_1[i] * u_2[i];
  }
  denom = 1.0 - SQR(u_1_dot_u_2);
  if (denom < SMALL) {
    *lambda = dr_dot_u_1 / 2.0;
    *mu = -dr_dot_u_2 / 2.0;
  } else {
    *lambda = (dr_dot_u_1 - u_1_dot_u_2 * dr_dot_u_2) / denom;
    *mu = (-dr_dot_u_2 + u_1_dot_u_2 * dr_dot_u_1) / denom;
  }
  lambda_mag = ABS(*lambda);
  mu_mag = ABS(*mu);

  /* Now take into account the fact that the two line segments are of finite
   * length. */
  if (lambda_mag > half_length_1 && mu_mag > half_length_2) {
    /* Calculate first possible case. */
    lambda_a = SIGN(half_length_1, *lambda);
    mu_a = -dr_dot_u_2 + lambda_a * u_1_dot_u_2;
    mu_mag = ABS(mu_a);
    if (mu_mag > half_length_2)
      mu_a = SIGN(half_length_2, mu_a);

    /* Calculate minimum distance between two spherocylinders. */
    r_min_mag2_a = 0.0;
    for (i = 0; i < n_dim_; ++i) {
      r_min_a[i] = dr[i] - lambda_a * u_1[i] + mu_a * u_2[i];
      r_min_mag2_a += SQR(r_min_a[i]);
    }

    /* Calculate second possible case. */
    mu_b = SIGN(half_length_2, *mu);
    lambda_b = dr_dot_u_1 + mu_b * u_1_dot_u_2;
    lambda_mag = ABS(lambda_b);
    if (lambda_mag > half_length_1)
      lambda_b = SIGN(half_length_1, lambda_b);

    /* Calculate minimum distance between two spherocylinders. */
    r_min_mag2_b = 0.0;
    for (i = 0; i < n_dim_; ++i) {
      r_min_b[i] = dr[i] - lambda_b * u_1[i] + mu_b * u_2[i];
      r_min_mag2_b += SQR(r_min_b[i]);
    }

    /* Choose the minimum minimum distance. */
    if (r_min_mag2_a < r_min_mag2_b) {
      *lambda = lambda_a;
      *mu = mu_a;
      *r_min_mag2 = r_min_mag2_a;
      for (i = 0; i < n_dim_; ++i)
        r_min[i] = r_min_a[i];
    } else {
      *lambda = lambda_b;
      *mu = mu_b;
      *r_min_mag2 = r_min_mag2_b;
      for (i = 0; i < n_dim_; ++i)
        r_min[i] = r_min_b[i];
    }
  } else if (lambda_mag > half_length_1) {
    /* Adjust lambda and mu. */
    *lambda = SIGN(half_length_1, *lambda);
    *mu = -dr_dot_u_2 + *lambda * u_1_dot_u_2;
    mu_mag = ABS(*mu);
    if (mu_mag > half_length_2)
      *mu = SIGN(half_length_2, *mu);

    /* Calculate minimum distance between two spherocylinders. */
    *r_min_mag2 = 0.0;
    for (i = 0; i < n_dim_; ++i) {
      r_min[i] = dr[i] - *lambda * u_1[i] + *mu * u_2[i];
      *r_min_mag2 += SQR(r_min[i]);
    }
  } else if (mu_mag > half_length_2) {
    /* Adjust lambda and mu. */
    *mu = SIGN(half_length_2, *mu);
    *lambda = dr_dot_u_1 + *mu * u_1_dot_u_2;
    lambda_mag = ABS(*lambda);
    if (lambda_mag > half_length_1)
      *lambda = SIGN(half_length_1, *lambda);

    /* Calculate minimum distance between two spherocylinders. */
    *r_min_mag2 = 0.0;
    for (i = 0; i < n_dim_; ++i) {
      r_min[i] = dr[i] - *lambda * u_1[i] + *mu * u_2[i];
      *r_min_mag2 += SQR(r_min[i]);
    }
  } else {
    /* Calculate minimum distance between two spherocylinders. */
    *r_min_mag2 = 0.0;
    for (i = 0; i < n_dim_; ++i) {
      r_min[i] = dr[i] - *lambda * u_1[i] + *mu * u_2[i];
      *r_min_mag2 += SQR(r_min[i]);
    }
  }

  return;
}

/* Routine to calculate minimum distance between a sphere and a spherocylinder,
 for any number of spatial dimensions and any type of boundary conditions (free,
 periodic, or mixed).

 input: real position of sphere (r_1)
 scaled position of sphere (s_1)
 real position of spherocylinder (r_2)
 scaled position of spherocylinder (s_2)
 director of spherocylinder (u_2)
 length of spherocylinder (length_2)

 output: minimimum separation vector (r_min)
 pointer to squared minimum separation (r_min_mag2)
 pointer to intersection of r_min with axis of spherocylinder (mu). */

void MinimumDistance::SphereSphero(
    double const *const r_1, double const *const s_1, double const *const r_2,
    double const *const s_2, double const *const u_2, double const length_2,
    double *r_min, double *r_min_mag2, double *contact2) {
  int i, j;
  double half_length_2, dr_dot_u_2, mu_mag, mu;
  double ds[3], dr[3];

  /* Compute various constants. */
  half_length_2 = 0.5 * length_2;

  /* Compute pair separation vector. */
  for (i = 0; i < n_periodic_; ++i) { /* First handle periodic subspace. */
    ds[i] = s_2[i] - s_1[i];
    ds[i] -= NINT(ds[i]);
  }
  for (i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    for (j = 0; j < n_periodic_; ++j)
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
  }
  for (i = n_periodic_; i < n_dim_; ++i) /* Then handle free subspace. */
    dr[i] = r_2[i] - r_1[i];

  /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
     First consider a point and an infinitely long line. */
  dr_dot_u_2 = 0.0;
  for (i = 0; i < n_dim_; ++i)
    dr_dot_u_2 += dr[i] * u_2[i];
  mu = -dr_dot_u_2;
  mu_mag = ABS(mu);

  /* Now take into account the fact that the line segment is of finite length.
   */
  if (mu_mag > half_length_2)
    mu = SIGN(half_length_2, mu);

  /* Calculate minimum distance between sphere and spherocylinder. */
  *r_min_mag2 = 0.0;
  for (i = 0; i < n_dim_; ++i) {
    r_min[i] = dr[i] + mu * u_2[i];
    contact2[i] = mu * u_2[i];
    *r_min_mag2 += SQR(r_min[i]);
  }

  return;
}

void MinimumDistance::SpheroPlane(double *r_mt, double *u_mt, double length,
                                  double *r_plane, double *n_plane,
                                  double *lambda, double *r_min_mag2,
                                  double *r_min) {
  int i;
  double r_1[3], r_2[3], d_1, d_2, d_min;
  double offset = dot_product(3, r_plane, n_plane);
  double costheta = dot_product(3, n_plane, u_mt);

  if (fabs(costheta) > SMALL) {
    for (i = 0; i < 3; ++i)
      r_1[i] = -0.5 * u_mt[i] * length + r_mt[i];
    for (i = 0; i < 3; ++i)
      r_2[i] = 0.5 * u_mt[i] * length + r_mt[i];

    d_1 = dot_product(3, r_1, n_plane) - offset;
    d_2 = dot_product(3, r_2, n_plane) - offset;
    if (fabs(d_1) < fabs(d_2)) {
      *lambda = -0.5 * length;
      d_min = d_1;
    } else {
      *lambda = 0.5 * length;
      d_min = d_2;
    }
  } else {
    *lambda = 0.0;
    d_min = dot_product(3, r_mt, n_plane) - offset;
  }

  for (i = 0; i < 3; ++i)
    r_min[i] = -n_plane[i] * d_min;
  *r_min_mag2 = dot_product(3, r_min, r_min);
}

/* Routine to calculate minimum distance between two liens, for any number of
   spatial dimensions and any type of boundary conditions (free, periodic, or
mixed).

input: real position of first spherocylinder (r_1)
       scaled position of first spherocylinder (s_1)
       director of first spherocylinder (u_1)
       real position of second spherocylinder (r_2)
       scaled position of second spherocylinder (s_2)
       director of second spherocylinder (u_2)

output: minimimum separation vector (r_min)
        pointer to squared minimum separation (r_min_mag2)
        pointer to intersection of r_min with axis of first spherocylinder
(lambda) pointer to intersection of r_min with axis of second spherocylinder
(mu). */

void MinimumDistance::CarrierLines(const double *r_1, const double *s_1,
                                   const double *u_1, const double *r_2,
                                   const double *s_2, const double *u_2,
                                   double *r_min, double *r_min_mag2,
                                   double *lambda, double *mu) {
  int i, j;
  double dr_dot_u_1, dr_dot_u_2, u_1_dot_u_2, denom;
  double ds[3], dr[3];

  /* Compute pair separation vector. */
  for (i = 0; i < n_periodic_; ++i) { /* First handle periodic subspace. */
    ds[i] = s_2[i] - s_1[i];
    ds[i] -= NINT(ds[i]);
  }
  for (i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    for (j = 0; j < n_periodic_; ++j)
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
  }
  for (i = n_periodic_; i < n_dim_; ++i) /* Then handle free subspace. */
    dr[i] = r_2[i] - r_1[i];

  /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
     First consider two infinitely long lines. */
  dr_dot_u_1 = dot_product(n_dim_, dr, u_1);
  dr_dot_u_2 = dot_product(n_dim_, dr, u_2);
  u_1_dot_u_2 = dot_product(n_dim_, u_1, u_2);

  denom = 1.0 - SQR(u_1_dot_u_2);
  if (denom < SMALL) {
    *lambda = dr_dot_u_1 / 2.0;
    *mu = -dr_dot_u_2 / 2.0;
  } else {
    *lambda = (dr_dot_u_1 - u_1_dot_u_2 * dr_dot_u_2) / denom;
    *mu = (-dr_dot_u_2 + u_1_dot_u_2 * dr_dot_u_1) / denom;
  }

  /* Calculate minimum distance between two lines. */
  *r_min_mag2 = 0.0;
  for (i = 0; i < n_dim_; ++i) {
    r_min[i] = dr[i] - *lambda * u_1[i] + *mu * u_2[i];
    *r_min_mag2 += SQR(r_min[i]);
  }

  return;
}

void MinimumDistance::PointSphereBC(double const *const r, double *dr,
                                    double *dr_mag2, double buffer) {
  double r_mag = 0;
  for (int i = 0; i < n_dim_; ++i) {
    r_mag += r[i] * r[i];
  }
  r_mag = sqrt(r_mag);
  for (int i = 0; i < n_dim_; ++i) {
    dr[i] = (space_->radius / r_mag - 1) * r[i];
  }
  *dr_mag2 = 0;
  for (int i = 0; i < n_dim_; ++i) {
    *dr_mag2 += dr[i] * dr[i];
  }
}
void MinimumDistance::SpheroWallBC(double const *const r, double const *const s,
                                   double const *const u, double const length,
                                   double *dr, double *dr_mag2,
                                   double *r_contact, double buffer) {
  /* For a spherocylinder with spherical BCs, the minimum distance will
     always be at one of the endpoints */
  double r_min[3] = {0, 0, 0};
  /* Check which site is closest to the boundary. If the sphero is on the right
     side of the box, check if the orientation is toward or away from the box,
     and vice versa */
  int sign;
  double wall;
  if (s[0] > 0) {
    sign = (u[0] > 0 ? 1 : -1);
    wall = 0.5;
  } else {
    sign = (u[0] < 0 ? 1 : -1);
    wall = -0.5;
  }
  double ds[3] = {0, 0, 0};
  ds[0] = wall - s[0];
  for (int i = 0; i < n_periodic_; ++i) {
    ds[i] -= NINT(ds[i]);
  }
  for (int i = 0; i < n_periodic_; ++i) {
    dr[i] = 0.0;
    for (int j = 0; j < n_periodic_; ++j) {
      dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
    }
  }
  // Then handle free subspace
  if (n_periodic_ == 0) {
    dr[0] = -space_->radius - r[0];
  }
  for (int i = 0; i < n_dim_; ++i) {
    r_contact[i] = sign * 0.5 * length * u[i];
    r_min[i] = r[i] + r_contact[i];
  }
  sign = (s[0] > 0 ? 1 : -1);
  dr[0] -= sign * ABS(r_min[0] - r[0]);
  *dr_mag2 = dr[0] * dr[0];
}

void MinimumDistance::PointWallBC(double const *const r, double const *const s,
                                  double *dr, double *dr_mag2, double buffer) {
  /* Check which site is closest to the boundary. If the sphero is on the right
     side of the box, check if the orientation is toward or away from the box,
     and vice versa */
  if (n_periodic_ > 0) {
    double ds[3] = {0, 0, 0};
    ds[0] = -s[0];
    for (int i = 0; i < n_periodic_; ++i) {
      dr[i] = 0.0;
      for (int j = 0; j < n_periodic_; ++j) {
        dr[i] += unit_cell_[n_dim_ * i + j] * ds[j];
      }
    }
  } else {
    dr[0] = -r[0];
  }
  *dr_mag2 = 0.0;
  for (int i = 0; i < n_dim_; ++i) {
    *dr_mag2 += SQR(dr[i]);
  }
}

void MinimumDistance::SpheroSphereBC(double const *const r,
                                     double const *const u, double const length,
                                     double *dr, double *dr_mag2,
                                     double *r_contact, double buffer) {
  /* For a spherocylinder with spherical BCs, the minimum distance will
     always be at one of the endpoints */
  double r_min[3] = {0, 0, 0};
  /* Check which site is furthest from the origin. This is done by
     looking at the sign of the dot product of the position and
     orientation of the spherocylinder. If it is positive, then it
     is the site in the positive direction of the sphero origin, and
     vice versa */
  int sign = SIGNOF(dot_product(n_dim_, r, u));
  for (int i = 0; i < n_dim_; ++i) {
    r_contact[i] = sign * 0.5 * length * u[i];
    r_min[i] = r[i] + r_contact[i];
  }

  double r_mag = 0;
  for (int i = 0; i < n_dim_; ++i) {
    r_mag += r_min[i] * r_min[i];
  }
  r_mag = sqrt(r_mag);

  for (int i = 0; i < n_dim_; ++i) {
    dr[i] = (space_->radius / r_mag - 1) * r_min[i];
  }
  *dr_mag2 = 0;
  for (int i = 0; i < n_dim_; ++i) {
    *dr_mag2 += dr[i] * dr[i];
  }
}

void MinimumDistance::PointBuddingBC(double const *const r, double *dr,
                                     double *dr_mag2, double buffer) {
  // First see which cell (mother or daughter) we are primarily located in
  bool in_mother = (r[n_dim_ - 1] < space_->bud_neck_height);
  /* There are two regions in which the cusp of the bud neck is always the
     minimum distance region to the object. The first is the cone defined from
     the origin of the mother cell to the bud neck, the second is defined from
     the origin of the daughter cell to the bud neck. */

  /* In 2D, the equation for a (double) cone is: h^2*x^2/r^2 = (y-y0)^2
     In 3D, the equation for a (double) cone is: h^2(x^2+y^2)/r^2 = (z-z0)^2
     where h = bud_neck_height, r = bud_neck_radius, and inside the mother cell
     z0 = 0, inside the daughter cell z0 = bud_height */

  /* First check to see if object has vertical coordinate less than 0, or
     greater than bud_height, in which case, we are definitely not in the cone
     regions */
  bool in_cone_region =
      !(r[n_dim_ - 1] < 0 || r[n_dim_ - 1] > space_->bud_height);

  /* If in_cone_region is false, then we are definitely not in the cone region,
     if in_cone_region is true, we need to check definitively whether we are in
     the cone region. */
  double r_mag = 0.0;
  /* use r_mag to store magnitude of rho^2 for now */
  for (int i = 0; i < n_dim_ - 1; ++i) {
    r_mag += SQR(r[i]);
  }
  double z0 = (in_mother ? 0 : space_->bud_height);
  if (in_cone_region) {
    /* Check to see if distance from vertical axis (rho) has the property
       rho^2 > r^2(z-z0)^2/h^2, if so, we are not in the cone region */
    double cone_rho2 = SQR(space_->bud_neck_radius) * SQR(r[n_dim_ - 1] - z0) /
                       SQR(space_->bud_neck_height);
    in_cone_region = (r_mag < cone_rho2);
  }
  /* Now in_cone_region definitely means what it says it means */
  if (in_cone_region) {
    /* Minimum distance is to the cusp of the bud neck, a little algebra shows
       that in cylindrical coords: dr = { rho ( bud_neck_radius / |rho| - 1) ,
       bud_neck_height - z } */
    double scale_factor = space_->bud_neck_radius / sqrt(r_mag) - 1;
    // Temp test code FIXME
    if (scale_factor < 0)
      Logger::Error("Something went wrong in PointBuddingBC!!\n");
    *dr_mag2 = 0;
    for (int i = 0; i < n_dim_ - 1; ++i) {
      dr[i] = scale_factor * r[i];
      *dr_mag2 += SQR(dr[i]);
    }
    dr[n_dim_ - 1] = space_->bud_neck_height - r[n_dim_ - 1];
    *dr_mag2 += SQR(dr[n_dim_ - 1]);
  }
  /* else we are not in the cone region, we do a typical spherical boundary
     check, with respect to the mother or daughter cell we are inside */
  else {
    r_mag = sqrt(r_mag + SQR(r[n_dim_ - 1] - z0));
    double r_cell = (in_mother ? space_->radius : space_->bud_radius);
    *dr_mag2 = 0;
    for (int i = 0; i < n_dim_ - 1; ++i) {
      dr[i] = (r_cell / r_mag - 1) * r[i];
      *dr_mag2 += SQR(dr[i]);
    }
    dr[n_dim_ - 1] = (r_cell / r_mag - 1) * (r[n_dim_ - 1] - z0);
    *dr_mag2 += SQR(dr[n_dim_ - 1]);
  }
}

void MinimumDistance::SpheroBuddingBC(double const *const r,
                                      double const *const u,
                                      double const length, double *dr,
                                      double *dr_mag2, double *r_contact,
                                      double buffer) {
  /* For spherocylinders, there are two distinct cases we want to consider:
     whether both end sites are above or below the bud neck in the same cell,
     or if they are in different cells. I determine this by determining if
     the z-coordinate of the sites are on the same side of the bud neck */
  double site_z = 0.5 * length * u[n_dim_ - 1];
  double z_offset = space_->bud_neck_height - r[n_dim_ - 1];
  bool same_cell = SIGNOF(z_offset + site_z) == SIGNOF(z_offset - site_z);
  /* If they're in the same cell, determine which site is furthest from the
     origin of their respective cell. This is done by taking a dot product,
     as is done in the sphere boundary, but now the origin can change */
  if (same_cell) {
    bool in_mother = (r[n_dim_ - 1] < space_->bud_neck_height);
    double z0 = (in_mother ? 0 : space_->bud_height);
    double dp = 0.0;
    for (int i = 0; i < n_dim_ - 1; ++i) {
      dp += r[i] * u[i];
    }
    dp += (r[n_dim_ - 1] - z0) * u[n_dim_ - 1];
    int sign = SIGNOF(dp);
    double r_min[3] = {0, 0, 0};
    for (int i = 0; i < n_dim_; ++i) {
      r_contact[i] = sign * 0.5 * length * u[i];
      r_min[i] = r[i] + r_contact[i];
    }
    PointBuddingBC(r_min, dr, dr_mag2, buffer);
  } else {
    // For now just use the centerpoint of the sphero... fix this later FIXME
    for (int i = 0; i < n_dim_; ++i) {
      r_contact[i] = 0.0;
    }
    PointBuddingBC(r, dr, dr_mag2, buffer);
  }
}

bool MinimumDistance::CheckBoundaryInteraction(Interaction &ix) {
  // No interaction with box boundary yet
  if (space_->type == +boundary_type::box ||
      space_->type == +boundary_type::none)
    return false;
  double const *const r1 = ix.obj1->GetInteractorPosition();
  double const *const s1 = ix.obj1->GetInteractorScaledPosition();
  double const *const u1 = ix.obj1->GetInteractorOrientation();
  double const l1 = ix.obj1->GetInteractorLength();
  double const d1 = ix.obj1->GetInteractorDiameter();
  ix.dr_mag2 = 0;
  std::fill(ix.dr, ix.dr + 3, 0.0);
  std::fill(ix.contact1, ix.contact1 + 3, 0.0);
  ix.buffer_mag = 0.5 * d1;
  ix.buffer_mag2 = ix.buffer_mag * ix.buffer_mag;
  if (space_->type == +boundary_type::sphere) {
    if (l1 == 0) {
      PointSphereBC(r1, ix.dr, &(ix.dr_mag2), ix.buffer_mag);
    } else {
      SpheroSphereBC(r1, u1, l1, ix.dr, &(ix.dr_mag2), ix.contact1,
                     ix.buffer_mag);
    }
  } else if (space_->type == +boundary_type::budding) {
    if (l1 == 0) {
      PointBuddingBC(r1, ix.dr, &(ix.dr_mag2), ix.buffer_mag);
    } else {
      SpheroBuddingBC(r1, u1, l1, ix.dr, &(ix.dr_mag2), ix.contact1,
                      ix.buffer_mag);
    }
  } else if (space_->type == +boundary_type::wall) {
    if (l1 == 0) {
      PointWallBC(r1, s1, ix.dr, &(ix.dr_mag2), ix.buffer_mag);
    } else {
      SpheroWallBC(r1, s1, u1, l1, ix.dr, &(ix.dr_mag2), ix.contact1,
                   ix.buffer_mag);
    }
  }
  if (ix.dr_mag2 < boundary_cut2_) {
    return true;
  }
  return false;
}

bool MinimumDistance::CheckOutsideBoundary(Object &obj) {
  if (space_->type == +boundary_type::none)
    return false;
  double const *const r = obj.GetInteractorPosition();
  double const *const u = obj.GetInteractorOrientation();
  double const l = obj.GetInteractorLength();
  double const d = obj.GetInteractorDiameter();
  double r_mag = 0.0;
  double z0 = 0.0;
  double r_boundary = space_->radius;
  int sign = (l > 0 ? SIGNOF(dot_product(n_dim_, r, u)) : 0);
  if (space_->type == +boundary_type::box ||
      space_->type == +boundary_type::wall) {
    for (int j = 0; j < n_dim_; ++j) {
      double r_far = r[j] + sign * 0.5 * l * u[j];
      if (ABS(r_far) > (r_boundary - d))
        return true;
    }
    return false;
  }
  if (space_->type == +boundary_type::budding &&
      r[n_dim_ - 1] > space_->bud_neck_height) {
    z0 = space_->bud_height;
    r_boundary = space_->bud_radius;
  }
  for (int i = 0; i < n_dim_ - 1; ++i) {
    r_mag += SQR(r[i] + sign * 0.5 * l * u[i]);
  }
  r_mag += SQR(r[n_dim_ - 1] + sign * 0.5 * l * u[n_dim_ - 1] - z0);
  return (r_mag > SQR(r_boundary - 0.5 * d));
}

#undef SMALL
