#ifndef BNMC_INCLUDE_MODELCOUNTER_H_
#define BNMC_INCLUDE_MODELCOUNTER_H_

#include <vector>
#include <bn-to-cnf/cnf.h>
#include <bnc/partition.h>
#include <bnc/types.h>
#include "architecture.h"
#include "partition.h"
#include "evidence.h"
#include "types.h"
#include "cache.h"
#include <bnc/timer.h>

namespace bnmc {

class Interface;

template <ModelType kModelType>
class ModelCounter {
    friend class Interface;
    public:
        ModelCounter();
        void Read();

        void Init(Evidence&);

        void SetEvidence(const Evidence&);
        void SetEvidence(const Architecture &, const Evidence&);
        void InitCache();
        void PrepareCache();
        void Prepare();

        template <std::size_t N = 1> void SetNodeIds(const Evidence &);
        template <std::size_t N = 1> void SetEvidenceCache(const Evidence &);
        template <std::size_t N = 1> void SetAllocCache();

        template <std::size_t N = 1> std::string Description();
        template <std::size_t N = 1> std::string ParallelDescription();

        template <std::size_t N = 1> probability_t Posterior();
        template <std::size_t N = 1> probability_t Posterior(const Architecture&, Cache&);

        template <std::size_t N = 1> probability_t ParallelPosterior(const unsigned int, Timer *t = NULL);
        template <std::size_t N = 1> probability_t ParallelPosterior(const Architecture&, Cache&);
        template <std::size_t N = 1> probability_t ParallelPosterior(const Architecture&, Cache&,Timer *t);

        template <std::size_t N = 1> void TraverseTiers(const Architecture&, Cache&, const EvidenceList&, const ConditionTierList &);
        template <std::size_t N = 1> void ParallelTraverseTiers(const Architecture&, Cache&, const EvidenceList&, const ConditionTierList &);
        template <std::size_t N = 1> probability_t ParallelGetPartitionResult(const Architecture&, Cache&, const EvidenceList&, const unsigned int kTier);
        template <std::size_t N = 1> void ParallelTraverseArchitecture(const Architecture&, Cache&, EvidenceList&, const ConditionTierList &, const unsigned int, const unsigned int);
        template <std::size_t N = 1> probability_t ParallelTraversePartition(const Architecture&, Cache&, EvidenceList&, const ConditionTierList &, const unsigned int kTier, const unsigned int kNodeId);
        template <std::size_t N = 1> probability_t ParallelTraverse(const Architecture &, const EvidenceList&, const ConditionTierList &, const unsigned int kTier, const unsigned int kNodeId);
        template <std::size_t N = 1> probability_t ParallelTraverse(const Architecture &, Cache&, const EvidenceList&, const ConditionTierList &, const unsigned int kTier, const unsigned int kNodeId);

        template <std::size_t N = 1> probability_t TraverseArchitecture(const Architecture&, Cache&, EvidenceList&, const ConditionTierList &, const unsigned int kTier = 0);
        template <std::size_t N = 1> probability_t TraversePartition(const Architecture&, Cache&, EvidenceList&, const ConditionTierList &, const unsigned int kTier = 0, const unsigned int kNodeId = 0);
        template <std::size_t N = 1> probability_t Traverse(const Architecture&, Cache &, const EvidenceList&, const ConditionTierList &, const unsigned int kTier, const unsigned int kNodeId = 0);
        template <std::size_t N = 1> probability_t Traverse(Cache&,const EvidenceList&, const ConditionTierList &);

        // == composition ==
        template <std::size_t N = 1> probability_t TraverseArchitecture(const Composition::Node * const, Cache&, EvidenceList&, const ConditionTierList &);
        template <std::size_t N = 1> probability_t TraversePartition(const Composition::Node * const, Cache&, EvidenceList&, const ConditionTierList &, const unsigned int kNodeId = 0);
        template <std::size_t N = 1> probability_t Traverse(const Composition::Node * const, Cache &, const EvidenceList&, const ConditionTierList &, const unsigned int kNodeId = 0);
        // =================

        template <std::size_t N = 1> probability_t GetPartitionResult(const Architecture&, Cache&, const EvidenceList&, const unsigned int kTier);

        const bn_partitions_t& GetBnPartitions() const;
        const Architecture &GetArchitecture() const;
    protected:

        Architecture      architecture_;
        ConditionTierList condition_tier_no_evidence_;
        std::vector<bool> is_persistence_;
        bool              has_query_variable_;
        variable_t        query_variable_;

        bn_partitions_t  bn_partitions_;
        std::vector< Partition > partition_;

        Cache             cache_;
        EvidenceList      evidence_list_;
        ConditionTierList condition_tier_;

        Cache             cache2_;
        EvidenceList      evidence_list2_;
        ConditionTierList condition_tier2_;

        const probability_t kNotTraversed;
        const probability_t kTraversed;
        const unsigned int kRootIndex;          // node at index 2 is always the root
        const unsigned int kTrueTerminalIndex;  // node at index 1 is always true terminal
        const unsigned int kFalseTerminalIndex; // node at index 1 is always true terminal
};

}

#endif

