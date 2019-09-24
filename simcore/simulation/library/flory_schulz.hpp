#ifndef _SIMCORE_FLORY_SCHULZ_H_
#define _SIMCORE_FLORY_SCHULZ_H_

#include "auxiliary.hpp"

/* The Flory-Schulz distribution is the distribution used for polydisperse
   polymers. It has the functional form
     w_a(k) = a^2 k (1-a)^(k-1),
   where a is a "polydispersity factor." It has a cumulative distribution
   function
     CDF = 1 - (1-a)^k (1+ak).
   This class is meant to sample randomly from the Flory-Schulz distribution by
   inverting the CDF so that it returns
     F(rand) = {min_k CDF(k) - rand == 0 }
   In order to find the root of the function, the bisection method is used for
   stability.
*/
class FlorySchulz {
private:
  double a_ = 0.03;       // polydispersity parameter, must be between 0 and 1
  double epsilon_ = 1e-3; // level of precision to use in our root-finder
  /* The cumulative distribution function of the Flory-Schulz distribution */
  double CDF(double k) { return 1 - pow(1 - a_, k) * (1 + a_ * k); }
  /* Function we are minimizing in our root-finder */
  double Func(double k, double roll) { return CDF(k) - roll; }
  /* Lower and upper bounds for the Flory-Schulz distribution */
  double upper_bound_ = 1000;
  /* Root finder that effectively finds the inverse of the CDF to map a uniform
     random number between 0 and 1 to the Flory-Schulz distribution */
  double Bisection(double roll) {
    if (roll < 0 || roll > 1) {
      error_exit("Flory-Schulz generator expects a value between 0 and 1");
    }
    double low = 0;
    double high = upper_bound_;
    /* Check that the lower and upper bounds are sane for the given
       polydispersity factor */
    if (Func(low, roll) * Func(high, roll) >= 0) {
      error_exit("Upper and lower bounds are not well-suited for the "
                 "polydispersity factor provided to the Flory-Schulz "
                 "generator");
      return -1;
    }
    double k = low;
    while ((high - low) >= epsilon_) {
      k = (low + high) / 2;
      if (Func(k, roll) == 0.0) {
        break;
      } else if (Func(k, roll) * Func(low, roll) < 0) {
        high = k;
      } else {
        low = k;
      }
    }
    return k;
  }

public:
  /* Initialize the Flory-Schulz distribution generator. The only necessary
     parameter is the polydispersity parameter, which is a sane value in the
     range of 0.001 and 0.5 for the purposes of this simulation, but more
     realistically it should be in the range of 0.01 and 0.1 */
  void Init(double a, double epsilon = 1e-3, double ub = 1000) {
    a_ = a;
    upper_bound_ = ub;
    if (ub < 0) {
      error_exit("The upper bound provided to the Flory-Schulz distribution "
                 "should be a positive number");
    }
    epsilon = epsilon_;
    if (a_ > 0.5 || a < 0.001) {
      error_exit("Flory-Schulz given off-base parameters. If you /really/ want"
                 " these polydispersity parameters, you'll have to remove this "
                 " error\n hint: a reasonable polydispersity parameter is in "
                 "the range of 0.001--0.5");
    }
    if (epsilon_ < 1e-18) {
      error_exit("You have provided the Flory-Schulz generator with an "
                 "unreasonable precision tolerance. Lower your expectations.");
    }
  }
  /* Given a uniform random number in the range of 0 and 1 (roll), return its
     corresponding Flory-Schulz distribution random value */
  double Rand(double roll) { return Bisection(roll); }
};

#endif
