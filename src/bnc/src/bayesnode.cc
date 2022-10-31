#include "bayesnode.h"

using namespace std;

Variable BayesNode::GetVariable() const {
    return variable_;
}

size_t  BayesNode::GetNrVariables() const {
    return variables_.size();
}

const std::vector<Variable>& BayesNode::GetVariables() const {
    return variables_;
}

const std::vector<literal_t>& BayesNode::GetLiterals() const {
    return literals_;
}

unsigned int BayesNode::GetLiteral(unsigned int i) const {
    return literals_[i];
}

unsigned int BayesNode::GetDimension(unsigned int i) const {
    return dimensions_[i];
}

const std::vector<unsigned int>& BayesNode:: GetDimensions() const {
    return dimensions_;
}

const std::vector<unsigned int>& BayesNode::GetWeights() const {
    return weights_;
}

size_t BayesNode::GetNrWeights() const {
    return weights_.size();
}

const std::vector<Probability>& BayesNode::GetProbabilities() const {
    return probabilities_;
}

size_t BayesNode::GetNrProbabilities() const {
    return probabilities_.size();
}

void BayesNode::SetVariable(Variable variable){
   variable_ = variable;
}

void BayesNode::SetVariable(unsigned int i, Variable variable){
    variables_[i] = variable;
}

void BayesNode::SetNrVariables(unsigned int nr_variables) {
    variables_.resize(nr_variables);
    literals_.resize(nr_variables);
    dimensions_.resize(nr_variables);
}

void BayesNode::SetLiteral(unsigned int i, unsigned int literal){
    literals_[i] = literal;
}

void BayesNode::SetDimension(unsigned int i, unsigned int dimension){
    dimensions_[i] = dimension;
}

std::vector<unsigned int>& BayesNode::SetWeights(){
    return weights_;
}


