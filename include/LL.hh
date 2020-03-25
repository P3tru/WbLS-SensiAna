//
// Created by zsoldos on 11/26/19.
//

#ifndef _LL_HH_
#define _LL_HH_

// Evaluate Negative log-likelihood
// Add special warning to deal with log(0) issues:
// in that case, NLL = -NPred
double EvalNLL(double Nobs, double Npred);

#endif