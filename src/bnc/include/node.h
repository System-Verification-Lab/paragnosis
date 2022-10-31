#ifndef NODE_H
#define NODE_H

#include <bn-to-cnf/cnf.h>
#include <stdint.h>
#include <stdlib.h>
#include <set>
#include <stdio.h>
#include <unordered_set>
#include <queue>
#include <stack>
#include <unordered_map>
#include <numeric>
#include "support.h"

struct weighted_node {
    literal_t l;
    int32_t ref;
    weighted_node *t;
    weighted_node *e;

    #include "nodeweight.h"

    weights *W;

    #include "nodedecl.h"

    typedef weighted_node_map             map;
    typedef weighted_node_set             set;
    typedef weighted_node_table           table;
    typedef weighted_node_stack           stack;
    typedef weighted_node_queue           queue;
    typedef weighted_node_support_map     support_map;
    typedef weighted_node_reference_stack reference_stack;
    typedef weighted_node_computed_table  computed_table;

    #include "nodedef.h"


} __attribute__((aligned(8))) ;


#endif

