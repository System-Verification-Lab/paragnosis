#ifndef NODE_H
#define NODE_H

#include <bn-to-cnf/config.h>
#include <vector>
#include <cstdint>

namespace bnmc {

struct LiteralNode {
    uint32_t l;   // literal
    uint32_t t;   // index of 'then' node
    uint32_t e;   // index of 'else' node
    probability_t w;

    inline bool IsTerminal(){
        return l == 0;
    };
};

typedef uint16_t Variable;
typedef uint16_t VariableValue;
typedef uint32_t NodeIndex;
typedef std::vector< NodeIndex > NodeIdList;

struct VariableNode {
    Variable v;      // variable
    VariableValue i; // i'th value
    NodeIndex t;     // index of 'then' node
    NodeIndex e;     // index of 'else' node
    probability_t w;
    VariableNode& operator=(LiteralNode&);
};

#define TERMINAL_LITERAL_NODE ((LiteralNode) {0,0,0,0})
#define TERMINAL_VARIABLE_NODE ((VariableNode) {0,0,0,0,0})

}

#endif

