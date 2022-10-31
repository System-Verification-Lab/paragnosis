#ifndef BNMC_INCLUDE_CACHE_H_
#define BNMC_INCLUDE_CACHE_H_

#include "spanning.h"
#include "evidence.h"
#include <vector>
#include <bnc/xary.h>
#include <atomic>
#include <typeinfo>
#include "persistence.h"
#include "architecture.h"
#include "mstack.h"
#include "node.h"
namespace bnmc {

extern Architecture test;
#define AtomicCache
typedef std::atomic<probability_t> AtomicProbability;

#ifdef AtomicCache
typedef AtomicProbability Probability;
#else
typedef probability_t Probability;
#endif

typedef std::vector<Probability>   TierCache;
typedef std::vector<TierCache>     ArchitectureCache;

typedef std::vector< NodeIdList > IdCache;
typedef std::vector< probability_t > ProbabilityList;
typedef std::vector< ProbabilityList > ProbabilityCache;
typedef std::vector< std::vector< EvidenceList > > EvidenceListCache;
typedef std::vector< unsigned int > IdentityList;
typedef std::vector< IdentityList > IdentityCache;

struct TaskGroup {
    TierId tier;
    NodeIdList nodes;
};

typedef std::vector< TaskGroup > TaskCache;
typedef std::atomic<unsigned int> TaskCounter;
typedef std::vector< TaskCounter > TaskDependencyCache;
typedef std::vector< mstack::Stack > StackCache;
typedef std::vector<size_t> SizeList;

class Cache {
    public:
        Cache();

        void InitTierProbabilitiers();
        void Init();

        void Resize(const unsigned int kBufferSize);

        void PrepareEvidenceLists();
        void PrepareTierProbabilities();
        void Prepare();


        void SetAllocSizes();
        void SetStackSizes(const bn_partitions_t&);
        void SetCircuitSizes(const std::vector< Partition >&);
        void SetTierSizes(const Architecture &kArchitecture);
        void SetSizes(const bn_partitions_t&,const std::vector< Partition >&, const Architecture&);

        void SetEvidence(const Architecture&, const Evidence&, bool with_query_variable = true);
        void SetNodeIds(const Architecture&, const Evidence&, bool with_query_variable = true);

        void SetEvidenceList(EvidenceList &dest, const Architecture &kArchitecture, const unsigned int &kTier, const unsigned int &NodeId) const;
        void SetEvidenceList(EvidenceList &dest, const EvidenceList &source, const Architecture &kArchitecture, const unsigned int &kTier, const unsigned int &NodeId) const;

        unsigned int GetNodeId(const Spanning &kInfo, const EvidenceList &kList, const unsigned int kTier) const;
        Probability& GetCacheEntry(const Spanning&, const EvidenceList&, const unsigned int);
        Probability& GetCacheEntry(const unsigned int kNodeId, const unsigned int kTier);

        EvidenceList& GetEvidenceList(const unsigned int kNodeId, const unsigned int kTier);
        const NodeIdList& GetNodeIds(const unsigned int kTier) const;

        IdentityList &GetIdentityList(const unsigned int kId);
        ProbabilityList& GetProbabilityList(const unsigned int kId);
        mstack::Stack& GetStack(const unsigned int kId);

        const unsigned int ObtainIdentityListId();
        const unsigned int ObtainProbabilityListId();
        const unsigned int ObtainStackId();

        void ReleaseIdentityList(const unsigned int);
        void ReleaseProbabilityList(const unsigned int);
        void ReleaseStack(const unsigned int);

        void ShallowCopy(const Cache&);

        const size_t GetMaxAllocSize() const;
        const size_t GetAllocSize(const unsigned int kPartitionId) const;
        const size_t GetStackSize(const unsigned int kPartitionId) const;
        const size_t GetCircuitSize(const unsigned int kPartitionId) const;
        const size_t GetTierSize(const unsigned int kTier) const;
        const size_t GetArchitectureSize() const;

        void SetTasks(const Architecture &kArchitecture, const ConditionTierList &kConditionTierList, const EvidenceList &kEvidenceList);
        size_t GetTaskCount() const;
        const TaskGroup& GetTaskGroup(const unsigned int kTaskId) const;
        TaskCounter& GetTaskCounter(const unsigned int kTaskId);
        const unsigned int GetTaskGroupId(const TierId kTier, const NodeIndex kNodeId) const ;

        static const probability_t kTraversed;
        static const probability_t kNotCached;
        static const unsigned int kMaxBufferSize;
        TaskDependencyCache task_issue_cache_;
        TaskCache task_reverse_cache_;
    private:
        // runtime caches
        StackCache stack_cache_;
        ProbabilityCache probability_cache_;
        IdentityCache identity_cache_;

        TaskCache task_cache_;
        TaskDependencyCache task_dependency_cache_;
        std::vector<unsigned int> task_dependency;
        ArchitectureParentGroups task_id_cache_;

        std::vector<unsigned int> stack_probabilities_;
        std::vector<unsigned int> stack_identity_;
        std::vector<unsigned int> stack_stack_;

        ArchitectureCache cache_;           // store probability of each achitecture node per tier
        EvidenceListCache evidence_cache_;  // store input evidence list for each architecture node per tier
        IdCache id_cache_;                  // store id of each archecture node per tier (given evidence)

        SizeList stack_sizes_;
        SizeList circuit_sizes_;
        SizeList tier_sizes_;
        SizeList alloc_sizes_;
        size_t max_circuit_size_;
        size_t max_stack_size_;
        size_t max_alloc_size_;
        unsigned int buffer_size_;
};

}

#endif

