#ifndef BNMC_INCLUDE_EVIDENCE_H_
#define BNMC_INCLUDE_EVIDENCE_H_

#include <string>
#include <bnc/bayesgraph.h>
#include "types.h"
#include <climits>

namespace bnmc {

struct EvidenceVariable {
    EvidenceVariable(Variable variable, VariableValue value) {
        this->variable = variable; this->value = value;
    };
    EvidenceVariable(){};
    Variable variable;
    VariableValue value;
    bool operator<(const EvidenceVariable &other) const;
    bool operator==(const EvidenceVariable &other) const;
};

typedef std::vector<VariableValue>  EvidenceList;
typedef std::vector<literal_t>      EvidenceLiteralList;
typedef std::set<EvidenceVariable>  EvidenceVariableSet;
typedef std::set<Variable>          VariableSet;

class Evidence {
    public:
        Evidence();
        void Clear();

        bool Conflict(const Variable);
        bool Conflict(const EvidenceVariable&);

        int AddPosterior(std::string);
        int Add(const literal_t);
        int Add(const Variable, const unsigned int);
        int Add(const EvidenceVariable&);
        int Add(const Evidence &);

        int Remove(literal_t);
        int Remove(EvidenceVariable);
        int Remove(Variable, unsigned int);
        int Remove(Variable);

        size_t Size() const;

        bool HasEvidence(const Variable&) const;

        void InducePosteriors();
        void ParsePosteriors(std::string);
        void Parse(std::string);
        const bool HaveQueryVariable() const;
        void ClearQueryVariable();
        int SetQueryVariable(const Variable);
        void UnsetQueryVariable();
        int AddQueryLiteral(const literal_t);
        int AddQueryVariable(const Variable, const unsigned int);
        int AddQueryVariable(const EvidenceVariable);
        const EvidenceVariable& GetQueryEvidenceVariable() const;
        const Variable GetQueryVariable() const;
        const VariableValue GetQueryValue() const;

        bool IsQueryVariable(EvidenceVariableSet::iterator) const;
        const literal_t GetQueryLiteral() const ;
        std::string GetQueryString(const bool kRealNames = true) const;
        static std::string GetAssignmentString(const Variable, const unsigned int, const bool kRealNames = true);

        EvidenceLiteralList GetEvidenceLiteralList() const;
        EvidenceList GetEvidenceList() const noexcept;
        EvidenceList GetEvidenceListNoQuery() const;

        const EvidenceVariableSet& GetEvidenceVariableSet() const;
        const VariableSet& GetPosteriorVariableSet() const;
        EvidenceVariable GetEvidenceVariable(literal_t);
    private:
        void InitQueryVariable();
        EvidenceVariableSet::iterator FindEvidenceVariable(const Variable&);
        EvidenceVariableSet::iterator FindEvidenceVariable(const EvidenceVariable&);
        EvidenceVariableSet::iterator query_variable_iterator_;
        EvidenceVariableSet evidence_variable_set_;
        VariableSet posterior_variable_set_;
};

}

#endif
