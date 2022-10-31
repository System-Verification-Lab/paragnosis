#ifndef BAYESNODE_H
#define BAYESNODE_H

#include <bn-to-cnf/bayesnet.h>
#include <bn-to-cnf/cnf.h>
#include <stdarg.h>
#include <vector>
#include "types.h"

class BayesNode {
    friend class bayesgraph;
    public:
        Variable GetVariable() const;
        size_t GetNrVariables() const;
        const std::vector<Variable>& GetVariables() const;

        const std::vector<literal_t>& GetLiterals() const ;
        unsigned int GetLiteral(unsigned int) const;

        const std::vector<unsigned int>& GetDimensions() const;
        unsigned int GetDimension(unsigned int) const;

        const std::vector<unsigned int>& GetWeights() const ;
        const std::vector<Probability>& GetProbabilities() const ;
        size_t GetNrWeights() const;
        size_t GetNrProbabilities() const;
        void SetVariable(Variable);
        void SetVariable(unsigned int, Variable);
        void SetNrVariables(unsigned int);
        void SetLiteral(unsigned int, unsigned int);
        void SetDimension(unsigned int, unsigned int);
        std::vector<unsigned int>& SetWeights();
    private:
        Variable variable_;
        std::vector<Variable> variables_;
        std::vector<literal_t> literals_;
        std::vector<unsigned int> dimensions_;
        std::vector<unsigned int> weights_;
        std::vector<double> probabilities_;
};

#endif
