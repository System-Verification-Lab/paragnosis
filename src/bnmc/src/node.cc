#include "node.h"
#include "options.h"

namespace bnmc {

VariableNode& VariableNode::operator=(LiteralNode& n){
    if(!n.IsTerminal()){ // if n is not terminal node
        unsigned int variable = manager.mapping.get_literal_to_variable()[n.l];
        unsigned int offset = manager.mapping.get_variable_to_literal()[variable];
        unsigned int value = n.l - offset;

        v = variable;
        i = value;
        t = n.t;
        e = n.e;
        w = n.w;
    } else (*this) = TERMINAL_VARIABLE_NODE;
    return (*this);
}

}
