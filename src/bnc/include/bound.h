#ifndef BOUND_H
#define BOUND_H

#include <bn-to-cnf/bayesnet.h>
#include <vector>
#include <gmpxx.h>
#include "ordering.h"
#include "oldbound.h"
#include "types.h"

class Bound {
    public:
        typedef mpz_class BoundSizeBig;
        typedef size_t BoundSize;
        enum Option {
            kNodes     = 0,
            kOperators = 1,
            kWeights   = 2,
            kRatio     = 3,
            kEdges     = 4,
            kScore     = 5
        };

        Bound();
        void Init();
        void Init(const partition_t&);
        void Clear();
        void SetBayesnet(bayesnet*);

        static BoundSize ComputeTree(bayesnet*, const partition_t&, const ordering_t&);
        static BoundSize Compute(bayesnet*, const partition_t&, const ordering_t&);
        static BoundSize Compute(bayesnet*, const ordering_t&);

        template <Bound::Option kOption = Option::kNodes, class T = Bound::BoundSize>
        const T ComputeTree(const ordering_t::value_type* ordering) const;

        template <Bound::Option kOption = Option::kNodes, class T = Bound::BoundSize>
        const T ComputeTree(const ordering_t& ordering) const;

        template <Bound::Option kOption = Option::kNodes,class T = Bound::BoundSize>
        const T Compute(const ordering_t::value_type* ordering) const;

        template <Bound::Option kOption = Option::kNodes,class T = Bound::BoundSize>
        const T Compute(const ordering_t& ordering) const;

        template <Bound::Option kOption = Option::kNodes>
        const double ComputeSafe(const ordering_t&) const ;

        template <Bound::Option kOption = Option::kNodes>
        const double ComputeSafe(const ordering_t::value_type*) const ;

        static inline bool IsOverflow(const BoundSize kSize){ return kSize == std::numeric_limits<BoundSize>::max(); };
    private:

        BoundSize Condition(std::set<Variable> &spanning, const Variable kVariable,unsigned int *, unsigned int *) const;

        std::vector<bool> enabled_;
        std::vector< std::vector<Variable> > variable_to_constraint_;
        std::vector< std::vector<Variable> > constraint_to_variable_;
        bayesnet *bn_;
        size_t size_;
};

#endif

