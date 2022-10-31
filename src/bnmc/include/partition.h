#ifndef BNMC_INCLUDE_PARTITION_H_
#define BNMC_INCLUDE_PARTITION_H_

#include <bnc/ordering.h>
#include <bn-to-cnf/config.h>
#include <type_traits>
#include <string>
#include "node.h"
#include "types.h"
#include "evidence.h"
#include <bnc/dynamicarray.h>
#include "multigraph.h"
namespace bnmc {

template <ModelType> class ModelCounter;

class Partition {
    friend class ModelCounter<ModelType::PWPBDD>;
    friend class ModelCounter<ModelType::WPBDD>;
    friend class ModelCounter<ModelType::MULTIGRAPH>;
    friend class ModelCounter<ModelType::PMULTIGRAPH>;
    friend class ModelCounter<ModelType::TDMULTIGRAPH>;
    friend class ModelCounter<ModelType::PTDMULTIGRAPH>;

    public:

        typedef VariableNode Node;
        typedef DynamicArray<Node> ArithmeticCircuit;
        typedef DynamicArray<LiteralNode> LiteralArithmeticCircuit;

        void Read(const ModelType,unsigned int partition_id = 0);
        const size_t GetCircuitSize() const;
    private:
        void Read(const ModelType, std::string);

        ordering_t ordering_;
        ArithmeticCircuit ac_;
        MultiGraph::Circuit mgc_;
        MultiGraph::Circuit tdmgc_;
};

}

#endif

