#ifndef ACE_H
#define ACE_H

#include <bn-to-cnf/cnf.h>
#include <bnc/timer.h>
#include "evidence.h"

probability_t ace_posterior(const bnmc::Evidence &,Timer &t);

#endif
