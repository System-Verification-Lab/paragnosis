#include "bound.h"
#include <cassert>
#include <type_traits>
#include <cmath>
#include <algorithm>
#include <list>
using namespace bnc;

#define MULT_ASSIGN_8(a, b) __builtin_umull_overflow (a, b, &a)
#define ADD_ASSIGN_8(a, b)  __builtin_uaddl_overflow (a, b, &a)

#define MULT_ASSIGN_4(a, b) __builtin_umul_overflow (a, b, &a)
#define ADD_ASSIGN_4(a, b)  __builtin_uadd_overflow (a, b, &a)

template <class T,class U>
inline bool add_assign_overflow(T &a, const U b){
    a += b;

    return false;
}

template <class U>
inline bool add_assign_overflow(size_t &a, const U b){
    return ADD_ASSIGN_8(a, b);
}

template <class U>
inline bool add_assign_overflow(unsigned int &a, const U b){
    return ADD_ASSIGN_4(a, b);
}



template <class T, class U>
inline bool mult_assign_overflow(T &a, const U b){
    a *= b;

    return false;
}

template <class U>
inline bool mult_assign_overflow(size_t &a, const U b){
    return MULT_ASSIGN_8(a, b);
}

template <class U>
inline bool mult_assign_overflow(unsigned int &a, const U b){
    return MULT_ASSIGN_4(a, b);
}

/*  GCC builtin overflow detection

    Built-in Function: bool __builtin_uadd_overflow (unsigned int a, unsigned int b, unsigned int *res)
    Built-in Function: bool __builtin_uaddl_overflow (unsigned long int a, unsigned long int b, unsigned long int *res)
    Built-in Function: bool __builtin_uaddll_overflow (unsigned long long int a, unsigned long long int b, unsigned long long int *res)

    Built-in Function: bool __builtin_umul_overflow (unsigned int a, unsigned int b, unsigned int *res)
    Built-in Function: bool __builtin_umull_overflow (unsigned long int a, unsigned long int b, unsigned long int *res)
    Built-in Function: bool __builtin_umulll_overflow (unsigned long long int a, unsigned long long int b, unsigned long long int *res)
*/


Bound::Bound(){
    bn_ = NULL;
}

void Bound::SetBayesnet(bayesnet *bn){
    bn_ = bn;
    size_ = 0;
}

void Bound::Clear(){
    enabled_.clear();
    variable_to_constraint_.clear();
    constraint_to_variable_.clear();
    size_ = 0;
}

void Bound::Init(){
    assert(bn_ != NULL && "Bayesnet must be set first");
    const unsigned int kNrVariables = bn_->get_nr_variables();
    assert(kNrVariables <= std::numeric_limits<Variable>::max());

    if(size_ == 0)
        size_ = kNrVariables;

    variable_to_constraint_.clear();
    constraint_to_variable_.clear();
    variable_to_constraint_.resize(kNrVariables);
    constraint_to_variable_.resize(kNrVariables);

    if(enabled_.empty()){
        enabled_.resize(kNrVariables);
        std::fill(enabled_.begin(), enabled_.end(), true);
    }

    for(unsigned int cpt = 0; cpt < kNrVariables; cpt++){
        if(enabled_[cpt]){
            variable_to_constraint_[cpt].push_back(cpt);
            constraint_to_variable_[cpt].push_back(cpt);

            const uint32_t* parent = bn_->get_parent(cpt);
            const unsigned int kNrParents = bn_->get_parent_size(cpt);
            for(unsigned int i = 0; i < kNrParents; i++){
                variable_to_constraint_[parent[i]].push_back(cpt);
                constraint_to_variable_[cpt].push_back(parent[i]);
            }
        }
    }
}

void Bound::Init(const partition_t &kPartition){
    assert(bn_ != NULL && "Bayesnet must be set first");
    const unsigned int kNrVariables = bn_->get_nr_variables();

    enabled_.resize(kNrVariables);
    std::fill(enabled_.begin(),enabled_.end(),false);
    for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++)
        enabled_[*it] = true;

    size_ = kPartition.cutset.size() + kPartition.set.size();

    Init();
}

Bound::BoundSize Bound::Condition(std::set<Variable> &spanning, const Variable kVariable, unsigned int *constraint_sat, unsigned int *variable_occurances) const {
    bool insert = false;
    BoundSize satisfied = 0;
    const auto &kConstraints = variable_to_constraint_[kVariable];
    for(auto it = kConstraints.begin(); it != kConstraints.end(); it++){
        const unsigned int &kConstraintId = *it;
        auto &sat = constraint_sat[kConstraintId];;
        #ifdef DEBUG
        assert(sat > 0);
        #endif
        sat--;
        if(sat == 0){ // constraint satisfied
            satisfied++;
            const auto &kConstraint = constraint_to_variable_[kConstraintId];
            for(auto vit = kConstraint.begin(); vit != kConstraint.end(); vit++){
                const auto &kConstraintVariable = *vit;
                auto &occurances = variable_occurances[kConstraintVariable];
                #ifdef BEGIN
                assert(occurances > 0);
                #endif
                occurances--;
                if(occurances == 0)
                    spanning.erase(kConstraintVariable);
            }
        } else insert = true;
    }
    if(insert)
        spanning.insert(kVariable);

    return satisfied;
}

template <>
const double Bound::Compute<Bound::Option::kScore,double>(const ordering_t::value_type* kOrdering) const {
    const unsigned int kNrVariables = bn_->get_nr_variables();
    const uint32_t *kDimension = bn_->get_states();

    assert(constraint_to_variable_.size() == kNrVariables && "Likely not initialized");

    // init
    unsigned int constraint_sat[kNrVariables];
    unsigned int variable_occurances[kNrVariables];
    for(unsigned int i = 0; i < kNrVariables; i++){
        constraint_sat[i] = constraint_to_variable_[i].size();
        variable_occurances[i] = variable_to_constraint_[i].size();
    }

    // compute bound
    double bound = 0;
    std::set<Variable> spanning;
    for(auto it = &(kOrdering[0]); it != &(kOrdering[size_]); it++){
        const auto kVariable = *it;

        // compute level width using spanning set
        double level_bound = 0;
        for(auto it = spanning.begin(); it != spanning.end(); it++){
            const Variable &kSpanningVariable = *it;
            level_bound += log(kDimension[kSpanningVariable]); // analogous to multiplying dimensions
        }

        // condition on kVariable and update spanning set
        Condition(spanning, kVariable,constraint_sat,variable_occurances);

        // compute overall bound
        bound += level_bound;
    }
    #ifdef DEBUG
    assert(spanning.empty());
    #endif
    return bound;
}

template <Bound::Option kOption,class T>
const T Bound::Compute(const ordering_t::value_type* kOrdering) const {
    static_assert(kOption != Bound::Option::kScore,"Use type double as template return type to compute score");
    const unsigned int kNrVariables = bn_->get_nr_variables();
    const uint32_t *kDimension = bn_->get_states();

    assert(constraint_to_variable_.size() == kNrVariables && "Likely not initialized");

    // init
    unsigned int constraint_sat[kNrVariables];
    unsigned int variable_occurances[kNrVariables];
    for(unsigned int i = 0; i < kNrVariables; i++){
        constraint_sat[i] = constraint_to_variable_[i].size();
        variable_occurances[i] = variable_to_constraint_[i].size();
    }

    // compute bound
    T bound = (T) 0;
    std::set<Variable> spanning;
    for(auto it = &(kOrdering[0]); it != &(kOrdering[size_]); it++){
        const auto kVariable = *it;

        // compute level width using spanning set
        T level_bound = (T) 1;
        for(auto it = spanning.begin(); it != spanning.end(); it++){
            const Variable &kSpanningVariable = *it;
            if(mult_assign_overflow(level_bound,kDimension[kSpanningVariable]))
                return std::numeric_limits<BoundSize>::max();
        }

        // condition on kVariable and update spanning set
        BoundSize weights = Condition(spanning, kVariable,constraint_sat,variable_occurances);

        // compute bound per level
        if(kOption != Option::kNodes){
            if(mult_assign_overflow(level_bound, kDimension[kVariable]))
                return std::numeric_limits<BoundSize>::max();

            if(kOption == Option::kWeights || kOption == Option::kOperators) {
                if(mult_assign_overflow(level_bound, (weights+(kOption==Option::kOperators?2:0))))
                    return std::numeric_limits<BoundSize>::max();
            }
        }

        // compute overall bound
        if(add_assign_overflow(bound,level_bound))
            return std::numeric_limits<BoundSize>::max();
    }
    assert(spanning.empty());
    return bound;
}

template <Bound::Option kOption,class T>
const T Bound::Compute(const ordering_t& kOrdering) const {
    assert(kOrdering.size() == size_);
    return Compute<kOption,T>(&(kOrdering[0]));
}

template <Bound::Option kOption>
const double Bound::ComputeSafe(const ordering_t& kOrdering) const {
    return Compute<kOption,BoundSizeBig>(&(kOrdering[0])).get_d();
}

template <>
const double Bound::ComputeSafe<Bound::Option::kRatio>(const ordering_t& kOrdering) const {
    return Compute<Bound::Option::kRatio,BoundSizeBig>(&(kOrdering[0])).get_d() / (double) bn_->get_nr_probabilities();
}

template <>
const double Bound::Compute<Bound::kRatio, double>(const ordering_t::value_type* ordering) const {
    return Bound::Compute<Bound::kWeights, double>(ordering) / (double) bn_->get_nr_probabilities();
};


template <>
const double Bound::Compute<Bound::Option::kRatio,double>(const ordering_t &kOrdering) const {
    return (double) Compute<Option::kRatio,double>(&(kOrdering[0]));
}

Bound::BoundSize Bound::Compute(bayesnet *bn, const partition_t &kPartition, const ordering_t &kOrdering ){
    Bound bound;
    bound.SetBayesnet(bn);
    bound.Init(kPartition);
    return bound.Compute(kOrdering);
}

Bound::BoundSize Bound::Compute(bayesnet *bn, const ordering_t &kOrdering){
    Bound bound;
    bound.SetBayesnet(bn);
    bound.Init();
    return bound.Compute(kOrdering);
}

template <class T>
struct BoundNode {
    Variable variable;
    T bound;
    std::set<Variable> spanning;
};

template <Bound::Option kOption,class T>
const T Bound::ComputeTree(const ordering_t::value_type* kOrdering) const {
    const unsigned int kNrVariables = bn_->get_nr_variables();
    const auto kDimensions = bn_->get_states();

    bool add_message[kNrVariables];
    std::copy(enabled_.begin(), enabled_.end(),&(add_message[0]));
    std::list< BoundNode<T> > roots;


    const auto kEnd = &(kOrdering[-1]);
    auto *it = &(kOrdering[size_-1]);
    while(it != kEnd){
        const Variable kVariable = *it;

        // see to which node to connect
        const auto kRootEnd = roots.end();
        unsigned nr_children = 0;
        auto itr = roots.begin();
        while(itr != kRootEnd){
            if(itr->spanning.find(kVariable) != itr->spanning.end()){
                ++nr_children;
                auto itr2 = itr;
                ++itr2;
                while(itr2 != kRootEnd){
                    if(itr2->spanning.find(kVariable) != itr2->spanning.end()){
                        itr->spanning.insert(itr2->spanning.begin(), itr2->spanning.end());
                        itr->bound += itr2->bound;
                        roots.erase(itr2++);
                        ++nr_children;
                    } else
                        ++itr2;
                }
                break;
            }
            ++itr;
        }

        // create new root if necessary
        if(nr_children == 0){
            roots.emplace_front();
            itr = roots.begin();
            itr->bound = 0;
        }
        itr->variable = kVariable;

        // add new messages
        const auto kConstraintEnd = variable_to_constraint_[kVariable].end();
        auto itc = variable_to_constraint_[kVariable].begin();
        while(itc != kConstraintEnd){
            const Variable kConstraint = *itc;
            if(add_message[kConstraint]){
                itr->spanning.insert(constraint_to_variable_[kConstraint].begin(), constraint_to_variable_[kConstraint].end());
                add_message[kConstraint] = false;
            }
            ++itc;
        }

        // remove current variable
        itr->spanning.erase(kVariable);

        // compute bound
        const auto kSpanningEnd = itr->spanning.end();
        auto its = itr->spanning.begin();
        T level_bound = 1;
        while(its != kSpanningEnd){
            if(kOption == Option::kScore){
                level_bound += log(kDimensions[*its]);
            } else if(kOption == Option::kNodes){
               level_bound *= kDimensions[*its];
                if(!std::is_same<T, double>::value && level_bound >= std::numeric_limits<unsigned int>::max())
                    return std::numeric_limits<unsigned int>::max();
            }
            ++its;
        }

        if(nr_children > 1){
            if(kOption == Option::kScore){
                level_bound += level_bound + log(kDimensions[kVariable]);
            } else if(kOption == Option::kNodes){
                level_bound += level_bound * kDimensions[kVariable];
                if(!std::is_same<T, double>::value && level_bound >= std::numeric_limits<unsigned int>::max())
                    return std::numeric_limits<unsigned int>::max();
            }
        }

        itr->bound += level_bound;
        if(!std::is_same<T, double>::value && itr->bound >= std::numeric_limits<unsigned int>::max())
            return std::numeric_limits<unsigned int>::max()-1;

        --it;
    }

    T bound = 0;
    const auto kRootEnd = roots.end();
    auto itr = roots.begin();
    while(itr != kRootEnd){
        bound += itr->bound;
        if(!std::is_same<T, double>::value && bound >= std::numeric_limits<unsigned int>::max())
            return std::numeric_limits<unsigned int>::max()-2;

        ++itr;
    }
    return bound;
}

template<Bound::Option kOption, class T>
const T Bound::ComputeTree(const ordering_t &kOrdering) const {
    assert(kOrdering.size() == size_);
    return ComputeTree<kOption,T>(&(kOrdering[0]));
}

Bound::BoundSize Bound::ComputeTree(bayesnet *bn, const partition_t &kPartition, const ordering_t &kOrdering ){
    Bound bound;
    bound.SetBayesnet(bn);
    bound.Init(kPartition);
    return bound.ComputeTree<Bound::kNodes, Bound::BoundSize>(kOrdering);
}

template const Bound::BoundSize Bound::ComputeTree<Bound::kNodes,     Bound::BoundSize>(const ordering_t& ordering) const;
template const double Bound::ComputeTree<Bound::kNodes,  double>(const ordering_t& ordering) const;
template const double Bound::ComputeTree<Bound::kScore, double>(const ordering_t& ordering) const;
//template const Bound::BoundSize Bound::ComputeTree<Bound::kOperators, Bound::BoundSize>(const ordering_t& ordering) const;
//template const Bound::BoundSize Bound::ComputeTree<Bound::kWeights,   Bound::BoundSize>(const ordering_t& ordering) const;
//template const Bound::BoundSize Bound::ComputeTree<Bound::kRatio,     Bound::BoundSize>(const ordering_t& ordering) const;
//template const Bound::BoundSize Bound::ComputeTree<Bound::kEdges,     Bound::BoundSize>(const ordering_t& ordering) const;

template const double Bound::ComputeSafe<Bound::kNodes    >(const ordering_t& ordering) const;
template const double Bound::ComputeSafe<Bound::kOperators>(const ordering_t& ordering) const;
template const double Bound::ComputeSafe<Bound::kWeights  >(const ordering_t& ordering) const;
template const double Bound::ComputeSafe<Bound::kRatio    >(const ordering_t& ordering) const;
template const double Bound::ComputeSafe<Bound::kEdges    >(const ordering_t& ordering) const;

template const Bound::BoundSize Bound::Compute<Bound::kNodes,     Bound::BoundSize>(const ordering_t& ordering) const;
template const Bound::BoundSize Bound::Compute<Bound::kOperators, Bound::BoundSize>(const ordering_t& ordering) const;
template const Bound::BoundSize Bound::Compute<Bound::kWeights,   Bound::BoundSize>(const ordering_t& ordering) const;
template const Bound::BoundSize Bound::Compute<Bound::kRatio,     Bound::BoundSize>(const ordering_t& ordering) const;
template const Bound::BoundSize Bound::Compute<Bound::kEdges,     Bound::BoundSize>(const ordering_t& ordering) const;

template const double Bound::Compute<Bound::kNodes,     double>(const ordering_t& ordering) const;
template const double Bound::Compute<Bound::kOperators, double>(const ordering_t& ordering) const;
template const double Bound::Compute<Bound::kWeights,   double>(const ordering_t& ordering) const;
template const double Bound::Compute<Bound::kRatio,     double>(const ordering_t& ordering) const;
template const double Bound::Compute<Bound::kEdges,     double>(const ordering_t& ordering) const;
template const double Bound::Compute<Bound::kScore,     double>(const ordering_t& ordering) const;

template const float Bound::Compute<Bound::kNodes,     float>(const ordering_t& ordering) const;
template const float Bound::Compute<Bound::kOperators, float>(const ordering_t& ordering) const;
template const float Bound::Compute<Bound::kWeights,   float>(const ordering_t& ordering) const;
template const float Bound::Compute<Bound::kRatio,     float>(const ordering_t& ordering) const;
template const float Bound::Compute<Bound::kEdges,     float>(const ordering_t& ordering) const;

