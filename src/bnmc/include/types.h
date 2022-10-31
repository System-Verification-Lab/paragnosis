#ifndef BNMC_INCLUDE_TYPES_H_
#define BNMC_INCLUDE_TYPES_H_

#include <bn-to-cnf/cnf.h>
#include "node.h"

namespace bnmc {

enum class ModelType { PWPBDD = 0, WPBDD = 1, MULTIGRAPH = 4, PMULTIGRAPH = 5, TDMULTIGRAPH = 6, PTDMULTIGRAPH = 7, UCLA_ACE = 11 };

typedef std::set<Variable>             VariableSet;

typedef uint8_t TierId;
#define TIER_INIT_VALUE UINT8_MAX

typedef std::vector<TierId>         ConditionTierList;

}
#endif

