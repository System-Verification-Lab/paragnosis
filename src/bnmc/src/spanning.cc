#include "spanning.h"
#include "exceptions.h"
#include <algorithm>

namespace bnmc {

void Spanning::Init(const Persistence &kPersistence){
    DetermineSpanningSets(kPersistence);
}

const SpanningSet& Spanning::GetSpanningSet(const unsigned int i) const {
    return spanning_sets_[i];
}

const SpanningSets& Spanning::GetSpanningSets() const {
    return spanning_sets_;
}

void Spanning::DetermineSpanningSets(const Persistence &kPersistence){
    PersistenceVariableCount persistence_variable_count = kPersistence.GetPersistenceVariableCount();
    const unsigned int kTiers = kPersistence.Size();
    spanning_sets_.resize(kTiers);

    Context context;
    for(unsigned int tier = 1; tier < kTiers; tier++)
        DetermineSpanningSet(persistence_variable_count, context, kPersistence, tier);
}

void Spanning::DetermineSpanningSet(PersistenceVariableCount &persistence_variable_count, Context &context, const Persistence &kPersistence, const unsigned int kTier){
    const PersistenceSet  &kPersistenceSet = kPersistence.GetPersistenceSet(kTier-1);
    for(auto it = kPersistenceSet.begin(); it != kPersistenceSet.end(); it++){
        const unsigned int &persistence_variable = *it;
        unsigned int &count = persistence_variable_count[persistence_variable];
        #ifdef DEBUG
        if(count == 0)
            throw SpanningException("[IMPOSSIBLE SITUATION] persistence count got below zero");
        #endif
        count--;
        if(count == 0)
            context.erase(persistence_variable);
        else context.insert(persistence_variable);
    }

    SpanningSet &spanning_set = spanning_sets_[kTier];
    std::copy(context.begin(), context.end(), std::back_inserter(spanning_set));
    std::sort(spanning_set.begin(),spanning_set.end());
    auto end = std::unique(spanning_set.begin(),spanning_set.end());
    spanning_set.resize(end-spanning_set.begin());
}

size_t Spanning::Size() const {
    return spanning_sets_.size();
}

}

