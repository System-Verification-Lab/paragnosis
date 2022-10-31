#ifndef BNMC_INCLUDE_DEBUG_H_
#define BNMC_INCLUDE_DEBUG_H_

#include "node.h"
#include "mstack.h"
#include "evidence.h"
#include "partition.h"
#include <vector>
#include "mstack.h"

namespace bnmc {

void WriteDot(std::string filename, const Partition::ArithmeticCircuit&, const EvidenceList &kEvidenceList, std::vector<probability_t> &probabilities, const std::vector<unsigned int>&identity, const unsigned int kTier, const int kCurrentIndex = -1, const int kMaxDepth = -1, mstack::StackTop top = 0, mstack::Stack *stack = NULL);

}


#endif
