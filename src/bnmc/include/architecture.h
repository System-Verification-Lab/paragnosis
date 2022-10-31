#ifndef ARCHITECTURE_T
#define ARCHITECTURE_T

#include "evidence.h"
#include "persistence.h"
#include "independence.h"
#include "spanning.h"
#include <bnc/ordering.h>
#include <bnc/partition.h>
#include <bnc/composition.h>

namespace bnmc {


typedef std::vector< std::vector< NodeIdList > > ArchitectureParentLists;
typedef std::vector< std::vector< NodeIndex > > ArchitectureParentGroups;

class Architecture {
    public:

        Architecture();
        ~Architecture();

        void Init(const bn_partitions_t &);
        void Init(const bn_partitions_t &, const Evidence&);

        bool IsPartitioned() const;
        void SetPartitionOrdering(const ordering_t &);
        const ordering_t& GetOrdering() const;
        const ordering_t& GetPartitionOrdering() const;
        const Composition &GetComposition() const;

        size_t Size() const;

        const Spanning& GetSpanning() const;
        const Persistence& GetPersistence() const;
        std::vector<size_t> GetTierWidths() const;

        void CreateGraph();
        const NodeIdList& GetParentList(const TierId, const NodeIndex kGroupId) const;
        const unsigned int GetParentGroup(const TierId, const NodeIndex kNodeId) const;
        const bool HasGraph() const;

        NodeIdList CreateParentList(const TierId, const NodeIndex, const ConditionTierList&, const EvidenceList&) const;
        NodeIdList CreateChildList(const TierId, const NodeIndex, const ConditionTierList&, const EvidenceList&) const;
        NodeIdList CreateSiblingList(const TierId, const NodeIndex, const ConditionTierList&, const EvidenceList&) const;
        NodeIdList CreateNodeList(const TierId, const ConditionTierList&, const EvidenceList&) const;

        ConditionTierList GetConditionTierList() const;
        ConditionTierList GetConditionTierList(const Evidence &, const bool kNoQuery = true) const;
        void ApplyEvidenceToConditionTierList(ConditionTierList&, const Evidence&, const bool kNoQuery = true) const;
        void InitConditionTierList(ConditionTierList &) const;

    private:
        void DeterminePartitionOrdering(const bn_partitions_t&);

        Composition comp_;

        Independence independence_;
        Persistence persistence_;
        Spanning spanning_;
        ArchitectureParentLists parent_lists_;
        ArchitectureParentGroups parent_groups_;
};

}

#endif
