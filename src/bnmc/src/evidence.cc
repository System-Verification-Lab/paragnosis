#include "options.h"
#include "evidence.h"
#include <bnc/exceptions.h>
#include "exceptions.h"
#include <algorithm>
#include <climits>
#include "io.h"
#include "trim.h"

#define NOT(x) (-1*x)

namespace bnmc {

bool EvidenceVariable::operator<(const EvidenceVariable &other) const {
    return variable < other.variable;
}

bool EvidenceVariable::operator==(const EvidenceVariable &other) const {
    return variable == other.variable;
}

inline const unsigned int GetVariableDimension(const Variable variable){
    return manager.mapping.get_dimension()[variable];
}

inline const Variable GetVariable(const literal_t literal){
    return manager.mapping.get_literal_to_variable()[literal];
}

inline const size_t GetTotalVariables(){
    return manager.mapping.get_nr_variables();
}

inline const size_t GetTotalLiterals(){
    return manager.mapping.get_nr_literals();
}

inline const literal_t GetBaseLiteral(const Variable variable){
    return manager.mapping.get_variable_to_literal()[variable];
}

inline const literal_t GetLiteral(const Variable variable, const unsigned int value){
    const unsigned int offset = GetBaseLiteral(variable);
    return offset + value;
}

inline EvidenceVariable CreateEvidenceVariable(const Variable variable, const unsigned int value){
    return {variable,value};
}

inline const literal_t GetLiteral(const EvidenceVariable &evidence_variable){
    return GetLiteral(evidence_variable.variable, evidence_variable.value);
}

inline void SetEvidenceVariable(EvidenceList& list, const EvidenceVariable &evidence_variable){
    list[evidence_variable.variable] = evidence_variable.value;
}

inline void SetEvidenceVariable(EvidenceLiteralList &literal_list, const EvidenceVariable &evidence_variable){
    const unsigned int kOffset = GetBaseLiteral(evidence_variable.variable);
    const unsigned int kDimension = GetVariableDimension(evidence_variable.variable);
    for(unsigned int value = 0; value < kDimension; value++)
        if(value != evidence_variable.value)
            literal_list[kOffset+value] = 0;
}

inline EvidenceVariable LiteralToEvidenceVariable(const literal_t literal){
    auto variable = GetVariable(literal);
    auto value = literal - GetBaseLiteral(variable);
    return std::move(EvidenceVariable(variable, value));
}

Evidence::Evidence(){
    InitQueryVariable();
}

void Evidence::InitQueryVariable(){
    query_variable_iterator_ = evidence_variable_set_.end();
}

const bool Evidence::HaveQueryVariable() const {
    return query_variable_iterator_ != evidence_variable_set_.end();
}

void Evidence::ClearQueryVariable(){
    if(HaveQueryVariable()){
        evidence_variable_set_.erase(query_variable_iterator_);
        InitQueryVariable();
    }
}

void Evidence::UnsetQueryVariable(){
    InitQueryVariable();
}


EvidenceVariableSet::iterator Evidence::FindEvidenceVariable(const EvidenceVariable &kEvidenceVariable){
    return std::move(evidence_variable_set_.find(kEvidenceVariable));
}

EvidenceVariableSet::iterator Evidence::FindEvidenceVariable(const Variable &kVariable){
    auto lambda = [&kVariable](const EvidenceVariable &kEvidenceVariable) {
        return kEvidenceVariable.variable == kVariable;
    };
    return std::move(std::find_if(evidence_variable_set_.begin(), evidence_variable_set_.end(), lambda));
}

bool Evidence::HasEvidence(const Variable &kVariable) const {
    auto lambda = [&kVariable](const EvidenceVariable &kEvidenceVariable) {
        return kEvidenceVariable.variable == kVariable;
    };
    return std::find_if(evidence_variable_set_.begin(), evidence_variable_set_.end(), lambda) != evidence_variable_set_.end();
}

int Evidence::SetQueryVariable(const Variable kVariable){
    query_variable_iterator_ = FindEvidenceVariable(kVariable);
    return query_variable_iterator_ != evidence_variable_set_.end();
}

int Evidence::AddQueryLiteral(const literal_t literal){
    return AddQueryVariable(LiteralToEvidenceVariable(literal));
}

int Evidence::AddQueryVariable(const Variable variable, const unsigned int value){
    return AddQueryVariable({variable, value});
}

int Evidence::AddQueryVariable(const EvidenceVariable evidence_variable){
    if(Conflict(evidence_variable))
        return 1;
    else {
        auto ret = evidence_variable_set_.insert(evidence_variable);

        if (ret.second==false)
            throw EvidenceException("Could not insert query variable");
        else query_variable_iterator_ = ret.first;
    }
    return 0;

}

const EvidenceVariable& Evidence::GetQueryEvidenceVariable() const {
    return *query_variable_iterator_;
}

const Variable Evidence::GetQueryVariable() const {
    return query_variable_iterator_->variable;
}

const VariableValue Evidence::GetQueryValue() const {
    return query_variable_iterator_->value;
}


bool Evidence::IsQueryVariable(EvidenceVariableSet::iterator it) const {
    return query_variable_iterator_ == it;
}

const literal_t Evidence::GetQueryLiteral() const {
    return GetLiteral(GetQueryEvidenceVariable());
}

int Evidence::Remove(literal_t literal){
    return Remove(LiteralToEvidenceVariable(literal));
}

int Evidence::Remove(Variable variable, unsigned int value){
    return Remove({variable,value});
}

int Evidence::Remove(Variable variable){
    return Remove(variable,0);
}

int Evidence::Remove(const EvidenceVariable kEvidenceVariable){
    evidence_variable_set_.erase(kEvidenceVariable);
    return 0;
}

int Evidence::AddPosterior(std::string variableName){
    InitQueryVariable();

    const unsigned int variableId = manager.bn->get_node_id(variableName);
    if(variableId == -1){
        throw EvidenceException("Unknown variable: %s\n", variableName.c_str());
    }

    if(Conflict(variableId))
        throw EvidenceException("Variable already in evidence: %s\n", variableName.c_str());
    else 
        posterior_variable_set_.insert(variableId);

    return 0;
} 

int Evidence::Add(const Evidence &evidence){
    auto evidenceVariables = evidence.GetEvidenceVariableSet();
    for( auto it = evidenceVariables.begin(); it !=evidenceVariables.end(); it++) {
        int ret = Add(*it);
        if (ret != 0)
            return ret;
    }
    return 0;
}

int Evidence::Add(const literal_t kLiteral){
    return Add(LiteralToEvidenceVariable(kLiteral));
}

int Evidence::Add(const Variable kVariable, const unsigned int kValue){
    return Add({kVariable,kValue});
}

int Evidence::Add(const EvidenceVariable &kEvidenceVariable){
    if(Conflict(kEvidenceVariable))
        return 1;
    else 
      evidence_variable_set_.insert(kEvidenceVariable);

    return 0;
}

EvidenceVariable Evidence::GetEvidenceVariable(literal_t literal){
    return std::move(LiteralToEvidenceVariable(literal));
}

const EvidenceVariableSet& Evidence::GetEvidenceVariableSet() const {
    return evidence_variable_set_;
}

const VariableSet& Evidence::GetPosteriorVariableSet() const {
    return posterior_variable_set_;
}

bool Evidence::Conflict(const EvidenceVariable &kEvidenceVariable){
    return FindEvidenceVariable(kEvidenceVariable) != evidence_variable_set_.end();
}

bool Evidence::Conflict(const Variable kVariable){
    return FindEvidenceVariable(kVariable) != evidence_variable_set_.end();
}

void Evidence::Clear(){
    UnsetQueryVariable();
    evidence_variable_set_.clear();
}

EvidenceLiteralList Evidence::GetEvidenceLiteralList() const {
    EvidenceLiteralList literal_list(GetTotalLiterals(),1);
    for(auto it = evidence_variable_set_.begin(); it != evidence_variable_set_.end(); it++)
        SetEvidenceVariable(literal_list,(*it));

    return std::move(literal_list);
}

EvidenceList Evidence::GetEvidenceListNoQuery() const {
    EvidenceList list(GetTotalVariables());
    for(auto it = evidence_variable_set_.begin(); it != evidence_variable_set_.end(); it++)
        if(!IsQueryVariable(it))
            SetEvidenceVariable(list,(*it));

    return std::move(list);
}

EvidenceList Evidence::GetEvidenceList() const noexcept {
    EvidenceList list(GetTotalVariables());
    for(auto it = evidence_variable_set_.begin(); it != evidence_variable_set_.end(); it++)
        SetEvidenceVariable(list,(*it));

    return std::move(list);
}

size_t Evidence::Size() const {
    return evidence_variable_set_.size();
}

std::string Evidence::GetAssignmentString(const Variable kVariable, const unsigned int  kValue, const bool kRealNames){
    if(kRealNames){
        std::string variable_name = manager.bn->get_node_name(kVariable);
        std::string value_name = manager.bn->get_node_value_name(kVariable,kValue);
        return stringf("%s=%s", variable_name.c_str(), value_name.c_str());
    } else return stringf("%u=%u",kVariable, kValue);
}

std::string Evidence::GetQueryString(const bool kRealNames) const {
    const EvidenceVariableSet &kEvidenceVariables = GetEvidenceVariableSet();

    std::string query;
    if(HaveQueryVariable()){
        const EvidenceVariable &kEvidenceVariable = GetQueryEvidenceVariable();
        query.append(GetAssignmentString(kEvidenceVariable.variable,kEvidenceVariable.value,kRealNames));
        if(kEvidenceVariables.size() > 1)
            query.append(" | ");
    }

    unsigned int nr_evidence = 0;
    for(auto it = kEvidenceVariables.begin(); it != kEvidenceVariables.end(); it++){
        if(!IsQueryVariable(it)){
            const EvidenceVariable &kEvidenceVariable = *it;
            if(nr_evidence > 0)
                query.append(", ");
            query.append(GetAssignmentString(kEvidenceVariable.variable,kEvidenceVariable.value,kRealNames));
            nr_evidence++;
        }
    }

    query.erase(std::remove(query.begin(), query.end(), '\0'), query.end());

    return query;
}

void Evidence::InducePosteriors(){
    InitQueryVariable();

    // fill variable set with all possible variable ids
    posterior_variable_set_.clear();
    for (unsigned int variableId = 0; variableId < manager.mapping.get_nr_variables(); variableId++) 
        posterior_variable_set_.insert(variableId);

    // then remove the variables that are mentioned in the evidence
    auto evidenceSet = GetEvidenceVariableSet();
    for (auto it = evidenceSet.begin(); it != evidenceSet.end(); it++) 
        posterior_variable_set_.erase(it->variable);

    if (posterior_variable_set_.empty())
        throw EvidenceException("There are no posterior variables set. Evidence covers all variables.");
}

void Evidence::ParsePosteriors(std::string query_string){
    posterior_variable_set_.clear();

    // check empty string
    Trim(query_string);
    if(query_string.empty())
        return;

    // parse query
    const unsigned int kNrVariables = std::count(query_string.begin(), query_string.end(), ',');
    while(true){
        size_t cm_pos = query_string.find_first_of(",");
        std::string variable;
        if(cm_pos != std::string::npos)
            variable = query_string.substr(0, query_string.find_first_of(","));
        else variable = query_string;

        std::string variable_name = variable;
        Trim(variable_name);

        if(variable_name.length() == 0)
            throw EvidenceException("Missing a variable name");

        const int a = AddPosterior(variable_name);
        
        if(a != 0)
            throw EvidenceException("Variable already mentioned as evidence: %s\n", variable_name.c_str());

        if(cm_pos != std::string::npos)
            query_string = query_string.substr(cm_pos+1);
        else break;
    }
}

void Evidence::Parse(std::string query_string){
    Clear();

    // check empty string
    Trim(query_string);
    if(query_string.empty())
        return;

    // check string format
    const unsigned int kConditional = std::count(query_string.begin(), query_string.end(), '|');
    if(kConditional > 1)
        throw EvidenceException("Syntax error: only one '|' symbol is allowed");

    // make query_string a comma delimited string
    if(kConditional){
        unsigned int pos = query_string.find_first_of("|");
        if(query_string.find_first_of(",") < pos)
            throw EvidenceException("Syntax error: you can only query one variable at a time");
        query_string.replace(pos,1,",");
    }

    // parse query
    const unsigned int kNrAssignments = std::count(query_string.begin(), query_string.end(), ',') + 1;
    while(true){
        size_t cm_pos = query_string.find_first_of(",");
        std::string assignment;
        if(cm_pos != std::string::npos)
            assignment = query_string.substr(0, query_string.find_first_of(","));
        else assignment = query_string;

        size_t eq_pos = assignment.find_first_of("=");
        if(std::string::npos == eq_pos || std::count(assignment.begin(), assignment.end(), '=') > 1)
            throw EvidenceException("Syntax error: expecting assignment, but got '%s'\n", assignment.c_str());

        std::string variable_name = assignment.substr(0,eq_pos);
        std::string value_name = assignment.substr(eq_pos+1, assignment.length());
        Trim(variable_name);
        Trim(value_name);

        if(variable_name.length() == 0){
            throw EvidenceException("Missing variable at assignment '%s'\n", assignment.c_str());
        }

        if(value_name.length() == 0){
            throw EvidenceException("Missing value at assignment '%s'\n", assignment.c_str());
        }

        const unsigned int l = manager.mapping.get_literal(variable_name,value_name);
        if(l == 0)
            throw EvidenceException("Unknown assignment: %s = %s\n", variable_name.c_str(), value_name.c_str());

        if((kNrAssignments == 1 || kConditional) && !HaveQueryVariable())
            AddQueryLiteral(l);
        else {
            if(Add(l) != 0)
                throw EvidenceException("Conflicting evidence in query at assignment %s = %s, variable %s already has a value\n", variable_name.c_str(), value_name.c_str(), variable_name.c_str());
        }

        if(cm_pos != std::string::npos)
            query_string = query_string.substr(cm_pos+1);
        else break;
    }
}

} // namespace bnmc

