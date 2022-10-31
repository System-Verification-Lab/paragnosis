#ifndef BNMC_INCLUDE_PERSISTENCE_H_
#define BNMC_INCLUDE_PERSISTENCE_H_

#include <map>
#include <vector>
#include <set>
#include <bnc/partition.h>
#include "evidence.h"
#include "types.h"
#include <unordered_map>
#include <stdint.h>

namespace bnmc {

typedef std::vector<unsigned int>   PersistenceVariableCount;
typedef std::vector<Variable>        PersistenceSet;
typedef std::vector<PersistenceSet> PersistenceSets;
typedef std::vector<bool>           PersistenceList;

class Persistence {
    public:
        void Init(const bn_partitions_t&, const ordering_t &);
        void Init(const bn_partitions_t&, const ordering_t &,const Evidence &);

        size_t Size() const;
        const PersistenceSet& GetPersistenceVariables() const;
        const PersistenceSet& GetPersistenceSet(const unsigned int) const;
        const PersistenceSets& GetPersistenceSets() const;
        const PersistenceVariableCount& GetPersistenceVariableCount() const;
        PersistenceList GetPersistenceList() const;
        ConditionTierList GetConditionTierList() const;
        ConditionTierList GetConditionTierList(const Evidence &) const;
        ConditionTierList GetConditionTierListNoQuery(const Evidence &) const;
        void FillConditionTierList(ConditionTierList &) const;
        void ApplyEvidenceToConditionTierList(ConditionTierList&, const Evidence&) const;
        void ApplyEvidenceToConditionTierListNoQuery(ConditionTierList&, const Evidence&) const;

    private:
        void DeterminePersistenceVariables(const bn_partitions_t&, const ordering_t&,const Evidence &kEvidence);

        PersistenceSet           persistence_variables_;
        PersistenceSets          persistence_sets_;
        PersistenceVariableCount persistence_variable_count_;
};

}

#endif

