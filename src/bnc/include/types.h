#ifndef BNC_TYPES_H
#define BNC_TYPES_H

#include <stdint.h>

typedef double Probability;
typedef uint16_t Variable;
typedef unsigned int Literal;

enum class compilation_t {
    none,
    topdown,
    bottomup,
    topdown_bottomup
};

enum class bdd_t {
    none,
    wpbdd,
    pwpbdd,
    multigraph,
    tdmultigraph
};

#endif
