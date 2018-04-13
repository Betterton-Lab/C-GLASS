#include "auxiliary.h"

#define SMALL 1.0e-12

// Returns minimum distance between two points in any given space.
// Along with vector pointing from r1 to r2
// Look around here for additional explanation of variables...
void min_distance_point_point(int n_dim, int n_periodic, double *unit_cell, 
                              double const * const r1, double const * const s1, 
                              double const * const r2, double const * const s2, 
                              double *dr, double *dr_mag2) {
  // First handle periodic subspace
  double ds[3];
  for (int i = 0; i < n_periodic; ++i) {
    ds[i] = s2[i] - s1[i];
    ds[i] -= NINT(ds[i]);
  }
  for (int i = 0; i < n_periodic; ++i) {
    dr[i] = 0.0;
    for (int j = 0; j < n_periodic; ++j)
      dr[i] += unit_cell[n_dim*i+j] * ds[j];
  }
  // Then handle free subspace
  for (int i = n_periodic; i < n_dim; ++i) 
      dr[i] = r2[i] - r1[i];
  *dr_mag2 = 0.0;
  for (int i=0; i<n_dim; ++i) 
      *dr_mag2 += SQR(dr[i]);
  return;
}

/* Routine to calculate minimum distance between a point and a line of finite length

output: vector that points from point to line along minimum distance between
        point and line (dr)
        distance from r_line along u_line that indicates point of minimum
        distance (mu) */

void min_distance_point_carrier_line(int n_dim, int n_periodic, double *h,
                                     double *r_point, double *s_point,
                                     double *r_line, double *s_line, double *u_line,
                                     double length, double *dr, double *r_contact, double *mu_ret) {

    int i, j;
    double ds[3], mu;

    /* Compute pair separation vector. */
    for (i = 0; i < n_periodic; ++i) {  /* First handle periodic subspace. */
        ds[i] = s_line[i] - s_point[i];
        ds[i] -= NINT(ds[i]);
    }
    for (i = 0; i < n_periodic; ++i) {
        dr[i] = 0.0;
        for (j = 0; j < n_periodic; ++j)
            dr[i] += h[n_dim*i+j] * ds[j];
    }
    for (i = n_periodic; i < n_dim; ++i)        /* Then handle free subspace. */
        dr[i] = r_line[i] - r_point[i];

    mu = -dot_product(n_dim, dr, u_line);
    double mu_mag = ABS(mu);
    // Now take into account that the line is finite length
    if (mu_mag > 0.5 * length) {
        mu = SIGNOF(mu) * 0.5 * length;
    }
    
    for (i=0; i<n_dim; ++i) {
        r_contact[i] = mu * u_line[i];
        dr[i] = r_line[i] + r_contact[i] - r_point[i];
    }
    *mu_ret = mu;
}

/* Routine to calculate minimum distance between a point and a line of infinite length

output: vector that points from point to line along minimum distance between
        point and line (dr)
        distance from r_line along u_line that indicates point of minimum
         distance (mu) */
void min_distance_point_carrier_line_inf(int n_dim, int n_periodic, double *h,
                                         double *r_point, double *s_point,
                                         double *r_line, double *s_line, double *u_line,
                                         double length, double *dr, double *mu) {

    int i, j;
    double ds[3];

    /* Compute pair separation vector. */
    for (i = 0; i < n_periodic; ++i) {  /* First handle periodic subspace. */
        ds[i] = s_line[i] - s_point[i];
        ds[i] -= NINT(ds[i]);
    }
    for (i = 0; i < n_periodic; ++i) {
        dr[i] = 0.0;
        for (j = 0; j < n_periodic; ++j)
            dr[i] += h[n_dim*i+j] * ds[j];
    }
    for (i = n_periodic; i < n_dim; ++i)        /* Then handle free subspace. */
        dr[i] = r_line[i] - r_point[i];

    *mu = -dot_product(n_dim, dr, u_line);
}
  
/* Routine to calculate minimum distance between two spherocylinders, for any number of
   spatial dimensions and any type of boundary conditions (free, periodic, or mixed).

input: number of spatial dimensions (n_dim)
       number of periodic dimensions (n_periodic)
       unit cell matrix (h)
       real position of first spherocylinder (r_1)
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

void min_distance_sphero(int n_dim, int n_periodic, double *h,
                         double const * const r_1, 
                         double const * const s_1, 
                         double const * const u_1, 
                         double const length_1,
                         double const * const r_2, 
                         double const * const s_2, 
                         double const * const u_2, 
                         double const length_2,
                         double *r_min, double *r_min_mag2, 
                         double *contact_1, double *contact_2) {
    int i, j;
    double half_length_1, half_length_2,
        dr_dot_u_1, dr_dot_u_2, u_1_dot_u_2, denom, lambda, mu,
        lambda_a, lambda_b, mu_a, mu_b, lambda_mag, mu_mag, r_min_mag2_a, r_min_mag2_b;
    double dr[3], ds[3], r_min_a[3], r_min_b[3];

    /* Compute various constants. */
    half_length_1 = 0.5 * length_1;
    half_length_2 = 0.5 * length_2;
    /* Compute pair separation vector. */
    for (i = 0; i < n_periodic; ++i) {  /* First handle periodic subspace. */
        ds[i] = s_2[i] - s_1[i];
        ds[i] -= NINT(ds[i]);
    }
    for (i = 0; i < n_periodic; ++i) {
        dr[i] = 0.0;
        for (j = 0; j < n_periodic; ++j)
            dr[i] += h[n_dim*i+j] * ds[j];
    }
    for (i = n_periodic; i < n_dim; ++i)        /* Then handle free subspace. */
        dr[i] = r_2[i] - r_1[i];

    /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
       First consider two infinitely long lines. */
    dr_dot_u_1 = dr_dot_u_2 = u_1_dot_u_2 = 0.0;
    for (i = 0; i < n_dim; ++i) {
        dr_dot_u_1 += dr[i] * u_1[i];
        dr_dot_u_2 += dr[i] * u_2[i];
        u_1_dot_u_2 += u_1[i] * u_2[i];
    }
    denom = 1.0 - SQR(u_1_dot_u_2);
    if (denom < SMALL) {
        lambda = dr_dot_u_1 / 2.0;
        mu = -dr_dot_u_2 / 2.0;
    } else {
        lambda = (dr_dot_u_1 - u_1_dot_u_2 * dr_dot_u_2) / denom;
        mu = (-dr_dot_u_2 + u_1_dot_u_2 * dr_dot_u_1) / denom;
    }
    lambda_mag = ABS(lambda);
    mu_mag = ABS(mu);

    /* Now take into account the fact that the two line segments are of finite length. */
    if (lambda_mag > half_length_1 && mu_mag > half_length_2) {

        /* Calculate first possible case. */
        lambda_a = SIGN(half_length_1, lambda);
        mu_a = -dr_dot_u_2 + lambda_a * u_1_dot_u_2;
        mu_mag = ABS(mu_a);
        if (mu_mag > half_length_2)
            mu_a = SIGN(half_length_2, mu_a);

        /* Calculate minimum distance between two spherocylinders. */
        r_min_mag2_a = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min_a[i] = dr[i] - lambda_a * u_1[i] + mu_a * u_2[i];
            r_min_mag2_a += SQR(r_min_a[i]);
        }

        /* Calculate second possible case. */
        mu_b = SIGN(half_length_2, mu);
        lambda_b = dr_dot_u_1 + mu_b * u_1_dot_u_2;
        lambda_mag = ABS(lambda_b);
        if (lambda_mag > half_length_1)
            lambda_b = SIGN(half_length_1, lambda_b);

        /* Calculate minimum distance between two spherocylinders. */
        r_min_mag2_b = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min_b[i] = dr[i] - lambda_b * u_1[i] + mu_b * u_2[i];
            r_min_mag2_b += SQR(r_min_b[i]);
        }

        /* Choose the minimum minimum distance. */
        if (r_min_mag2_a < r_min_mag2_b) {
            lambda = lambda_a;
            mu = mu_a;
            *r_min_mag2 = r_min_mag2_a;
            for (i = 0; i < n_dim; ++i)
                r_min[i] = r_min_a[i];
        } else {
            lambda = lambda_b;
            mu = mu_b;
            *r_min_mag2 = r_min_mag2_b;
            for (i = 0; i < n_dim; ++i)
                r_min[i] = r_min_b[i];
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
        for (i = 0; i < n_dim; ++i) {
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
        for (i = 0; i < n_dim; ++i) {
            r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
            *r_min_mag2 += SQR(r_min[i]);
        }
    } else {
        /* Calculate minimum distance between two spherocylinders. */
        *r_min_mag2 = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min[i] = dr[i] - lambda * u_1[i] + mu * u_2[i];
            *r_min_mag2 += SQR(r_min[i]);
        }
    }
    for (i=0; i<n_dim; ++i) {
      contact_1[i] = lambda * u_1[i];
      contact_2[i] = mu * u_2[i];
    }
    return;
}


/* Routine to calculate minimum distance between two spherocylinders and
   center to center separation vector, for any number of
   spatial dimensions and any type of boundary conditions (free, periodic, or mixed).

input: number of spatial dimensions (n_dim)
       number of periodic dimensions (n_periodic)
       unit cell matrix (h)
       real position of first spherocylinder (r_1)
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
        pointer to intersection of r_min with axis of first spherocylinder (lambda)
        pointer to intersection of r_min with axis of second spherocylinder (mu). */

void min_distance_sphero_dr(int n_dim, int n_periodic, double *h,
                            double *r_1, double *s_1, double *u_1, double length_1,
                            double *r_2, double *s_2, double *u_2, double length_2,
                            double *dr, double *r_min, double *r_min_mag2, 
                            double *lambda, double *mu) {
    int i, j;
    double half_length_1, half_length_2,
        dr_dot_u_1, dr_dot_u_2, u_1_dot_u_2, denom,
        lambda_a, lambda_b, mu_a, mu_b, lambda_mag, mu_mag, r_min_mag2_a, r_min_mag2_b;
    double ds[3], r_min_a[3], r_min_b[3];

    /* Compute various constants. */
    half_length_1 = 0.5 * length_1;
    half_length_2 = 0.5 * length_2;

    /* Compute pair separation vector. */
    for (i = 0; i < n_periodic; ++i) {  /* First handle periodic subspace. */
        ds[i] = s_2[i] - s_1[i];
        ds[i] -= NINT(ds[i]);
    }
    for (i = 0; i < n_periodic; ++i) {
        dr[i] = 0.0;
        for (j = 0; j < n_periodic; ++j)
            dr[i] += h[n_dim*i+j] * ds[j];
    }
    for (i = n_periodic; i < n_dim; ++i)        /* Then handle free subspace. */
        dr[i] = r_2[i] - r_1[i];

    /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
       First consider two infinitely long lines. */
    dr_dot_u_1 = dr_dot_u_2 = u_1_dot_u_2 = 0.0;
    for (i = 0; i < n_dim; ++i) {
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

    /* Now take into account the fact that the two line segments are of finite length. */
    if (lambda_mag > half_length_1 && mu_mag > half_length_2) {

        /* Calculate first possible case. */
        lambda_a = SIGN(half_length_1, *lambda);
        mu_a = -dr_dot_u_2 + lambda_a * u_1_dot_u_2;
        mu_mag = ABS(mu_a);
        if (mu_mag > half_length_2)
            mu_a = SIGN(half_length_2, mu_a);

        /* Calculate minimum distance between two spherocylinders. */
        r_min_mag2_a = 0.0;
        for (i = 0; i < n_dim; ++i) {
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
        for (i = 0; i < n_dim; ++i) {
            r_min_b[i] = dr[i] - lambda_b * u_1[i] + mu_b * u_2[i];
            r_min_mag2_b += SQR(r_min_b[i]);
        }

        /* Choose the minimum minimum distance. */
        if (r_min_mag2_a < r_min_mag2_b) {
            *lambda = lambda_a;
            *mu = mu_a;
            *r_min_mag2 = r_min_mag2_a;
            for (i = 0; i < n_dim; ++i)
                r_min[i] = r_min_a[i];
        } else {
            *lambda = lambda_b;
            *mu = mu_b;
            *r_min_mag2 = r_min_mag2_b;
            for (i = 0; i < n_dim; ++i)
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
        for (i = 0; i < n_dim; ++i) {
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
        for (i = 0; i < n_dim; ++i) {
            r_min[i] = dr[i] - *lambda * u_1[i] + *mu * u_2[i];
            *r_min_mag2 += SQR(r_min[i]);
        }
    } else {

        /* Calculate minimum distance between two spherocylinders. */
        *r_min_mag2 = 0.0;
        for (i = 0; i < n_dim; ++i) {
            r_min[i] = dr[i] - *lambda * u_1[i] + *mu * u_2[i];
            *r_min_mag2 += SQR(r_min[i]);
        }
    }

    return;
}

/* Routine to calculate minimum distance between a sphere and a spherocylinder, for any number of
 spatial dimensions and any type of boundary conditions (free, periodic, or mixed).
 
 input: number of spatial dimensions (n_dim)
 number of periodic dimensions (n_periodic)
 unit cell matrix (h)
 real position of sphere (r_1)
 scaled position of sphere (s_1)
 real position of spherocylinder (r_2)
 scaled position of spherocylinder (s_2)
 director of spherocylinder (u_2)
 length of spherocylinder (length_2)
 
 output: minimimum separation vector (r_min)
 pointer to squared minimum separation (r_min_mag2)
 pointer to intersection of r_min with axis of spherocylinder (mu). */

void min_distance_sphere_sphero(int n_dim, int n_periodic, double *h,
                                double const * const r_1, double const * const s_1,
                                double const * const r_2, double const * const s_2, 
                                double const * const u_2, double const length_2,
                                double *r_min, double *r_min_mag2,
                                double *contact2) {
    int i, j;
    double half_length_2, dr_dot_u_2, mu_mag, mu;
    double ds[3], dr[3];

    /* Compute various constants. */
    half_length_2 = 0.5 * length_2;

    /* Compute pair separation vector. */
    for (i = 0; i < n_periodic; ++i) {  /* First handle periodic subspace. */
        ds[i] = s_2[i] - s_1[i];
        ds[i] -= NINT(ds[i]);
    }
    for (i = 0; i < n_periodic; ++i) {
        dr[i] = 0.0;
        for (j = 0; j < n_periodic; ++j)
            dr[i] += h[n_dim*i+j] * ds[j];
    }
    for (i = n_periodic; i < n_dim; ++i)        /* Then handle free subspace. */
        dr[i] = r_2[i] - r_1[i];

    /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
       First consider a point and an infinitely long line. */
    dr_dot_u_2 = 0.0;
    for (i = 0; i < n_dim; ++i)
        dr_dot_u_2 += dr[i] * u_2[i];
    mu = -dr_dot_u_2;
    mu_mag = ABS(mu);

    /* Now take into account the fact that the line segment is of finite length. */
    if (mu_mag > half_length_2)
        mu = SIGN(half_length_2, mu);

    /* Calculate minimum distance between sphere and spherocylinder. */
    *r_min_mag2 = 0.0;
    for (i = 0; i < n_dim; ++i) {
        r_min[i] = dr[i] + mu * u_2[i];
        contact2[i] = mu * u_2[i];
        *r_min_mag2 += SQR(r_min[i]);
    }

    return;
}

void min_distance_sphero_plane(double *r_mt, double *u_mt, double length,
                               double *r_plane, double *n_plane, double *lambda,
                               double *r_min_mag2, double *r_min) {
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
        }
        else {
            *lambda = 0.5 * length;
            d_min = d_2;
        }
    }
    else {
        *lambda = 0.0;
        d_min = dot_product(3, r_mt, n_plane) - offset;        
    }

    for (i = 0; i < 3; ++i)
        r_min[i] = -n_plane[i] * d_min;
    *r_min_mag2 = dot_product(3, r_min, r_min);
}

/* Routine to calculate minimum distance between two liens, for any number of
   spatial dimensions and any type of boundary conditions (free, periodic, or mixed).

input: number of spatial dimensions (n_dim)
       number of periodic dimensions (n_periodic)
       unit cell matrix (h)
       real position of first spherocylinder (r_1)
       scaled position of first spherocylinder (s_1)
       director of first spherocylinder (u_1)
       real position of second spherocylinder (r_2)
       scaled position of second spherocylinder (s_2)
       director of second spherocylinder (u_2)

output: minimimum separation vector (r_min)
        pointer to squared minimum separation (r_min_mag2)
        pointer to intersection of r_min with axis of first spherocylinder (lambda)
        pointer to intersection of r_min with axis of second spherocylinder (mu). */

void min_distance_carrier_lines(int n_dim, int n_periodic, double *h,
                                double *r_1, double *s_1, double *u_1, 
                                double *r_2, double *s_2, double *u_2, 
                                double *r_min, double *r_min_mag2, 
                                double *lambda, double *mu) {
    int i, j;
    double dr_dot_u_1, dr_dot_u_2, u_1_dot_u_2, denom;
    double ds[3], dr[3];

    /* Compute pair separation vector. */
    for (i = 0; i < n_periodic; ++i) {  /* First handle periodic subspace. */
        ds[i] = s_2[i] - s_1[i];
        ds[i] -= NINT(ds[i]);
    }
    for (i = 0; i < n_periodic; ++i) {
        dr[i] = 0.0;
        for (j = 0; j < n_periodic; ++j)
            dr[i] += h[n_dim*i+j] * ds[j];
    }
    for (i = n_periodic; i < n_dim; ++i)        /* Then handle free subspace. */
        dr[i] = r_2[i] - r_1[i];

    /* Compute minimum distance (see Allen et al., Adv. Chem. Phys. 86, 1 (1993)).
       First consider two infinitely long lines. */
    dr_dot_u_1 = dot_product(n_dim, dr, u_1);
    dr_dot_u_2 = dot_product(n_dim, dr, u_2);
    u_1_dot_u_2 = dot_product(n_dim, u_1, u_2);

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
    for (i = 0; i < n_dim; ++i) {
       r_min[i] = dr[i] - *lambda * u_1[i] + *mu * u_2[i];
       *r_min_mag2 += SQR(r_min[i]);
    }

    return;
}

#undef SMALL