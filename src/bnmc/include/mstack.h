#ifndef MSTACK_H_
#define MSTACK_H_

#include <vector>
#include "types.h"
#include "exceptions.h"
#include "partition.h"
#include "evidence.h"

namespace bnmc {
namespace mstack {

typedef unsigned int StackTop;
typedef std::vector<NodeIndex> Stack;
typedef NodeIndex* CStack;
typedef Partition::ArithmeticCircuit Circuit;
const unsigned int kStackItemSize = 2;

inline void StackTopInit(StackTop &i){
    i = 0;
}

inline bool Empty(const unsigned int &i){
    return i == 0;
}

inline void Increment(StackTop &i){
    i += kStackItemSize;
}

inline void Decrement(StackTop &i){
    i -= kStackItemSize;
}

inline void Pop(unsigned int &i){
    #ifdef DEBUG
    if(!(i >= kStackItemSize))
        throw ModelCounterException("Custom stack is going below 0 elements");
    #endif
    Decrement(i);
}

inline bool Contains(const Stack &s, const unsigned int &from, const unsigned int &to, const unsigned int kStackTop){
    for(StackTop t = kStackItemSize; t <= kStackTop; Increment(t))
       if(s[t-2] == from && s[t-1] == to)
           return true;

    return false;
}

inline void Push(
        Stack &s,
        StackTop &i,
        const NodeIndex kParentIndex,
        const NodeIndex kIndex){

    #ifdef DEBUG
    if(i+kStackItemSize > s.size())
        throw ModelCounterException("Custom stack is going out of bounds");
    #endif

    s[i] = kParentIndex;
    s[i+1] = kIndex;
    Increment(i);
}
inline void Push(
        CStack s,
        StackTop &i,
        const NodeIndex kParentIndex,
        const NodeIndex kIndex){

    s[i] = kParentIndex;
    s[i+1] = kIndex;
    Increment(i);
}
// used by wpbdd traversal
inline void Push(
        Stack &s,
        StackTop &i,
        const NodeIndex &kIndex,
        const Circuit &kCircuit,
        const EvidenceList &kEvidenceList,
        const ConditionTierList &kConditionTier
        ){

    const TierId kTier = 0;
    const Partition::Node &kNode = kCircuit[kIndex];
    const Variable &kVariable = kNode.v;
    const bool kIsConditioned = kConditionTier[kVariable] <= kTier;
    const bool kIsTrue = kEvidenceList[kVariable] == kNode.i;

    if(!kIsConditioned){
        Push(s,i,kIndex,kNode.e);
        Push(s,i,kIndex,kNode.t);
    } else if(kIsTrue) {
        Push(s,i,kIndex,kNode.t);
    } else Push(s,i,kIndex,kNode.e);
}
// used by wpbdd traversal
inline void Push(
        CStack s,
        StackTop &i,
        const NodeIndex &kIndex,
        const Circuit &kCircuit,
        const EvidenceList &kEvidenceList,
        const ConditionTierList &kConditionTier
        ){

    const TierId kTier = 0;
    const Partition::Node &kNode = kCircuit[kIndex];
    const Variable &kVariable = kNode.v;
    const bool kIsConditioned = kConditionTier[kVariable] <= kTier;
    const bool kIsTrue = kEvidenceList[kVariable] == kNode.i;

    if(!kIsConditioned){
        Push(s,i,kIndex,kNode.e);
        Push(s,i,kIndex,kNode.t);
    } else if(kIsTrue) {
        Push(s,i,kIndex,kNode.t);
    } else Push(s,i,kIndex,kNode.e);
}

// used by pwpbdd to traverse last partition
inline void Push(
        Stack &s,
        StackTop &i,
        const NodeIndex &kIndex,
        const Circuit &kCircuit,
        const EvidenceList &kEvidenceList,
        const ConditionTierList &kConditionTier,
        const TierId kTier){

    const Partition::Node &kNode = kCircuit[kIndex];
    const Variable &kVariable = kNode.v;
    const bool kIsConditioned = kConditionTier[kVariable] <= kTier;
    const bool kIsTrue = kEvidenceList[kVariable] == kNode.i;

    if(!kIsConditioned){
        Push(s,i,kIndex,kNode.e);
        Push(s,i,kIndex,kNode.t);
    } else if(kIsTrue) {
        Push(s,i,kIndex,kNode.t);
    } else Push(s,i,kIndex,kNode.e);
}
// used by pwpbdd to traverse last partition
inline void Push(
        CStack s,
        StackTop &i,
        const NodeIndex &kIndex,
        const Circuit &kCircuit,
        const EvidenceList &kEvidenceList,
        const ConditionTierList &kConditionTier,
        const TierId kTier){

    const Partition::Node &kNode = kCircuit[kIndex];
    const Variable &kVariable = kNode.v;
    const bool kIsConditioned = kConditionTier[kVariable] <= kTier;
    const bool kIsTrue = kEvidenceList[kVariable] == kNode.i;

    if(!kIsConditioned){
        Push(s,i,kIndex,kNode.e);
        Push(s,i,kIndex,kNode.t);
    } else if(kIsTrue) {
        Push(s,i,kIndex,kNode.t);
    } else Push(s,i,kIndex,kNode.e);
}

// used by pwpbdd to traverse partition
inline void Push(
        Stack &s,
        StackTop &i,
        const NodeIndex &kIndex,
        const Circuit &kCircuit,
        EvidenceList &evidence_list,
        ConditionTierList &condition_tier,
        const TierId &kTier){

    const Partition::Node &kNode = kCircuit[kIndex];
    const Variable &kVariable = kNode.v;
    TierId &condition_in_tier = condition_tier[kVariable];
    const bool kIsConditioned = condition_in_tier <= kTier;
    VariableValue &value = evidence_list[kVariable];
    const bool kIsTrue = value == kNode.i;

    if(!kIsConditioned){
        Push(s,i,kIndex,kNode.e);
        Push(s,i,kIndex,kNode.t);

        condition_in_tier = kTier+1;
        value = kNode.i;
    } else if(kIsTrue) {
        Push(s,i,kIndex,kNode.t);
    } else Push(s,i,kIndex,kNode.e);
}
// used by pwpbdd to traverse partition
inline void Push(
        CStack s,
        StackTop &i,
        const NodeIndex &kIndex,
        const Circuit &kCircuit,
        EvidenceList &evidence_list,
        ConditionTierList &condition_tier,
        const TierId &kTier){

    const Partition::Node &kNode = kCircuit[kIndex];
    const Variable &kVariable = kNode.v;
    TierId &condition_in_tier = condition_tier[kVariable];
    const bool kIsConditioned = condition_in_tier <= kTier;
    VariableValue &value = evidence_list[kVariable];
    const bool kIsTrue = value == kNode.i;

    if(!kIsConditioned){
        Push(s,i,kIndex,kNode.e);
        Push(s,i,kIndex,kNode.t);

        condition_in_tier = kTier+1;
        value = kNode.i;
    } else if(kIsTrue) {
        Push(s,i,kIndex,kNode.t);
    } else Push(s,i,kIndex,kNode.e);
}
} // namespace mstack

} // namespace bnmc

#endif
