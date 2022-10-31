#ifndef BNMC_INCLUDE_SPANNING_H_
#define BNMC_INCLUDE_SPANNING_H_

#include "types.h"
#include "persistence.h"
#include <bn-to-cnf/cnf.h>
#include <unordered_set>

namespace bnmc {

typedef std::vector<Variable>            SpanningSet;
typedef std::vector<SpanningSet>         SpanningSets;
typedef std::unordered_set<Variable>     Context;


class Spanning {
    public:
        void Init(const Persistence &kPersistence);
        const SpanningSet& GetSpanningSet(const unsigned int) const;
        const SpanningSets& GetSpanningSets() const;
        size_t Size() const;

    private:
        void DetermineSpanningSets(const Persistence&);
        void DetermineSpanningSet(PersistenceVariableCount &, Context &, const Persistence&, const unsigned int);

        SpanningSets        spanning_sets_;
};

}

#endif
