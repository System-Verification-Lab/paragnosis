#include "exceptions.h"
#include "bayesgraph.h"
#include <map>
#include <vector>
#include <unistd.h>
#include <string.h>
#include <bn-to-cnf/cnf.h>
#include "options.h"
#include <algorithm>
#include "xary.h"
#include <cassert>

using namespace std;
using namespace bnc;

literal_mapping_t::literal_mapping_t(){
    bn = NULL;
    DEFINED = false;
}

probability_t literal_mapping_t::get_probability(unsigned int l) const {
    if(literal_to_variable.size() > l)
        return -1;
    else return weight_to_probability[l - literal_to_variable.size()];
}

bool literal_mapping_t::is_probability(unsigned int l) const {
    return (l >= literal_to_variable.size());
}

void literal_mapping_t::write() const {
    if(DEFINED){
        FILE *file = fopen(files.get_filename_c(MAPPING), "w");
        if(file){
            fprintf(file, "map %u %u %u\n", get_nr_variables(), get_nr_literals(), get_nr_weights());

            // map variables
            for(unsigned int v = 0; v < get_nr_variables(); v++)
                fprintf(file, "%u,%u,%u,%s\n", v, dimension[v], variable_to_literal[v], bn->get_node_name(v).c_str());

            // map literals
            for(unsigned int l = 1; l < literal_to_variable.size(); l++){
                unsigned int variable = literal_to_variable[l];
                unsigned int value = l - variable_to_literal[variable];
                fprintf(file, "%u,%u,%s\n", l, variable, bn->get_node_value_name(variable,value).c_str());
            }

            // map weights
            const unsigned int LITERALS = get_nr_literals();
            for(unsigned int w = 0; w < weight_to_probability.size(); w++){
                fprintf(file, "%u,%lf\n", LITERALS+1+w, weight_to_probability[w]);
            }

            fclose(file);
        } else throw compiler_write_exception("Could not open mapping file '%s'", files.get_filename_c(MAPPING));
    } else throw compiler_mapping_exception("Mapping not defined");
}

void literal_mapping_t::read(const char *filename){
    FILE *file = fopen(filename, "r");
    if(file){
        unsigned int VARIABLES, LITERALS, WEIGHTS;
        if(fscanf(file, "map %u %u %u\n", &VARIABLES, &LITERALS, &WEIGHTS) != 3){
            throw compiler_mapping_exception("Error while reading preemble of '%s'", filename);
        }

        // map variables
        std::vector<std::string> variable_names;
        variable_names.resize(VARIABLES);
        const unsigned int MAX_SIZE = 256;
        char word[MAX_SIZE];
        dimension.resize(VARIABLES);
        variable_to_literal.resize(VARIABLES);
        for(unsigned int i = 0; i < VARIABLES; i++){
            unsigned int v,d,l;
            if(fscanf(file, "%u,%u,%u,%255[^\n]", &v, &d, &l, word) != 4)
                throw compiler_mapping_exception("Error while reading variable %u", i);

            variable_names[v] = std::string(word);
            variable_to_literal[v] = l;
            dimension[v] = d;
        }

        // map literals
        literal_to_variable.resize(LITERALS+1);
        for(unsigned int i = 0; i < LITERALS; i++){
            unsigned int l, variable_id;
            if(fscanf(file, "%u,%u,%255[^\n]", &l, &variable_id, word) != 3)
                throw compiler_mapping_exception("Error while reading nodes in '%s'", filename);

            literal_to_variable[l] = variable_id;
            variable_value_name_to_literal[variable_names[variable_id]][std::string(word)] = l;
        }

        // map weights
        weight_to_probability.resize(WEIGHTS);
        for(unsigned int i = 0; i < WEIGHTS; i++){
            unsigned int l;
            probability_t w;
            if(fscanf(file, "%u,%lf\n", &l, &w) == 2){
                weight_to_probability[l-(LITERALS+1)] = w;
            } else throw compiler_mapping_exception("Error while reading weights in '%s'", filename);
        }
        fclose(file);
    } else throw compiler_mapping_exception("Could not open mapping file '%s'", filename);

}

void literal_mapping_t::read(){
    read(files.get_filename_c(MAPPING));
}

bool literal_mapping_t::defined() const{
    return DEFINED;
}

bayesnet *literal_mapping_t::get_bayesnet() const {
    return bn;
}

const std::vector<probability_t>& literal_mapping_t::get_weight_to_probability() const {
    return weight_to_probability;
}

unsigned int literal_mapping_t::get_nr_literals() const {
    return literal_to_variable.size()-1;
}

unsigned int literal_mapping_t::get_nr_variables() const {
    return variable_to_literal.size();
}

unsigned int literal_mapping_t::get_nr_weights() const {
    return weight_to_probability.size();
}

const std::vector<unsigned int>& literal_mapping_t::get_variable_to_literal() const {
    return variable_to_literal;
}

const variable_literal_map_t& literal_mapping_t::get_variable_value_name_to_literal() const {
    return variable_value_name_to_literal;
}

const std::vector<unsigned int>& literal_mapping_t::get_literal_to_variable() const {
    return literal_to_variable;
}

const std::vector<unsigned int>& literal_mapping_t::get_dimension() const {
    return dimension;
}

unsigned int literal_mapping_t::get_literal(std::string &variable_name, std::string &value_name) const {
    auto &variables = variable_value_name_to_literal;
    auto variable_hit = variables.find(variable_name);
    if(variable_hit == variables.end())
        return 0;

    auto &values = variable_hit->second;
    auto value_hit = values.find(value_name);
    if(value_hit == values.end())
        return 0;

    return value_hit->second;
}

// ========================= bayesgraph ====================================

int bayesgraph::encode(bayesnet *bn){
    this->bn = bn;
    DEFINED = true;

    const unsigned int VARIABLES = bn->get_nr_variables();
    nodes.resize(VARIABLES);
    dimension.resize(VARIABLES);
    memcpy(&dimension.front(), bn->get_states(),sizeof(unsigned int)*VARIABLES);

    unsigned int LITERALS = 0;
    for(unsigned int v = 0; v < VARIABLES; v++)
        LITERALS += dimension[v];

    variable_to_literal.resize(VARIABLES);
    literal_to_variable.resize(LITERALS+1);
    unsigned int l = 1;
    for(unsigned int i = 0; i < VARIABLES; i++){
        variable_to_literal[i] = l;
        for (unsigned int v = 0; v < dimension[i]; v++)
            literal_to_variable[l+v] = i;
        l += dimension[i];
    }

    // encode equal probabilities per CPT
    for(unsigned int v = 0; v < VARIABLES; v++){

        // determine dimensions of variables involved
        vector<unsigned int> dimensions;
        unsigned int PARENTS = bn->get_parent_size(v);
        for(unsigned int i = 0; i < PARENTS; i++)
            dimensions.push_back(dimension[bn->get_parent(v)[i]]);
        dimensions.push_back(dimension[v]);

        // store cpt variables as literals, ordered by variable index
        std::vector<unsigned int> cpt_variables;
        for(unsigned int i = 0; i < PARENTS; i++)
            cpt_variables.push_back(bn->get_parent(v)[i]);
        cpt_variables.push_back(v);
        auto ordered_cpt_variables = cpt_variables;
        std::sort(ordered_cpt_variables.begin(),ordered_cpt_variables.end());

        nodes[v].SetVariable(v);
        nodes[v].SetNrVariables(ordered_cpt_variables.size());
        for(unsigned int i = 0; i < ordered_cpt_variables.size(); i++){
            const auto &kCptVariable = ordered_cpt_variables[i];
            nodes[v].SetVariable(i,kCptVariable);
            nodes[v].SetLiteral(i,variable_to_literal[kCptVariable]);
            nodes[v].SetDimension(i,dimension[kCptVariable]);
        }

        // create HUGIN to natural order mapping
        //
        // hugin node format example:
        // A | B C
        // probabilities: based on [ B C A ], where B is most significant
        // transeform this to
        // [ A B C ] where A is least significant

        rXAry xary_hugin;
        xary_hugin.SetDimension(bn,cpt_variables);

        XAry xary;
        xary.SetDimension(bn,ordered_cpt_variables);
        xary.SetZero();

        XAry::Map map[ordered_cpt_variables.size()];
        XAry::CreateMapUnsorted(cpt_variables,ordered_cpt_variables,map);
        assert(bn->get_cpt_size(v) == xary_hugin.Max());

        // encode probabilities symbolically
        std::map< probability_t, unsigned int> probability_to_weight;
        probability_t *cpt = bn->get_cpt(v);
        auto &weights = nodes[v].weights_;
        auto &probabilities = nodes[v].probabilities_;

        weights.resize(bn->get_cpt_size(v));
        probabilities.resize(bn->get_cpt_size(v));
        for(unsigned int i = 0; i < bn->get_cpt_size(v); i++){
            probability_t p = cpt[i];

            xary_hugin.SetDecimal(i);
            xary.Set(xary_hugin,map);
            const unsigned int ni = xary.GetDecimal();


            // encode equal probabilities as identical symbolic weights
            probabilities[ni] = p;
            if(OPT_DETERMINISM && (p == 1 || p == 0)){
                weights[ni] = (unsigned int) p;
            } else {
                unsigned int symbolic_w;
                if(OPT_ENCODE_STRUCTURE){
                    symbolic_w = weight_to_probability.size();
                    auto hit = probability_to_weight.find(p);
                    if(hit != probability_to_weight.end())
                        symbolic_w = hit->second;
                    else {
                        probability_to_weight[p] = symbolic_w;
                        weight_to_probability.push_back(p);
                    }
                } else {
                    symbolic_w = weight_to_probability.size();
                    weight_to_probability.push_back(p);
                }
                weights[ni] = LITERALS + 1 + symbolic_w;
            }
        }
    }
    return 0;
}

BayesNode* bayesgraph::get_node(unsigned int i){
    return &(nodes[i]);
}

const BayesNode* bayesgraph::get_node(unsigned int i) const {
    return &(nodes[i]);
}


unsigned int bayesgraph::size(){
    return nodes.size();
}

void bayesgraph::WriteUAI(std::string filename) {
    if(filename == "")
        filename = files.get_filename_c(UAI);

    FILE *file = fopen(filename.c_str(), "w");
    if(file){
        fprintf(file, "BAYES\n");

        // print number of variables
        fprintf(file, "%lu\n", get_nr_variables());

        // print dimensions of each variable
        const auto & kDimensions = get_dimension();
        for(auto it = kDimensions.begin(); it != kDimensions.end(); it++)
            fprintf(file, "%lu ", *it);
        fprintf(file, "\n");

        // print nodes and connections
        const auto kNrCpts = get_nr_variables();
        fprintf(file, "%lu\n", kNrCpts);
        for(unsigned int cpt = 0; cpt < kNrCpts; cpt++){
            BayesNode *n = get_node(cpt);

            // print number of connections
            const unsigned int kNrCptVariables = n->GetNrVariables();
            fprintf(file, "%lu", kNrCptVariables);

            const unsigned int kCptVariable = n->GetVariable();
            assert(cpt == kCptVariable);

            // print connections
            const auto kCptVariables = n->GetVariables();
            for(unsigned int i = 0; i < kNrCptVariables; i++)
               fprintf(file, " %lu", kCptVariables[i]);
            fprintf(file, "\n");
        }

        // print cpts
        for(unsigned int cpt = 0; cpt < kNrCpts; cpt++){
            BayesNode *n = get_node(cpt);

            const auto kWeights = n->GetWeights();
            const auto kNrWeights = n->GetNrWeights();
            const auto kNrLiterals = get_nr_literals();
            const auto &kWeightToProbability = get_weight_to_probability();

            // bayesgraph stores probabilities with least significat bit on first position
            // this must be inverted
            const unsigned int kNrCptVariables = n->GetNrVariables();
            const auto kCptVariables = n->GetVariables();

            rXAry xary;
            xary.SetDimension(bn,kCptVariables);

            fprintf(file, "%lu", kNrWeights);
            // for(unsigned int i = 0; i < kNrWeights; i++){
            do {
                const auto kWeightId = xary.GetDecimal();
                assert(kWeightId < kNrWeights);

                const auto kWeight = kWeights[kWeightId];
                if(kWeight == 0)
                    fprintf(file, " 0.0");
                else if(kWeight == 1)
                    fprintf(file, " 1.0");
                else
                    fprintf(file, " %lf", kWeightToProbability[kWeight-(kNrLiterals+1)]);
            } while(xary.Increment());
            fprintf(file, "\n");

        }
        fclose(file);
    } else throw compiler_exception("Could not write UAI to '%s'", filename.c_str());

    filename.append(".map");
    file = fopen(filename.c_str(), "w");
    if(file){
        fprintf(file, "# <id>:<variable name>:<value 1>,...,<value n>\n");
        const auto kNrVariables = bn->get_nr_variables();
        for(unsigned int variable = 0; variable < kNrVariables; variable++){
            std::string variable_name =  bn->get_node_name(variable);
            fprintf(file, "%lu:%s:",variable,variable_name.c_str());

            const auto kNrValues = dimension[variable];
            for(unsigned int value = 0; value < kNrValues; value++){
                std::string value_name = bn->get_node_value_name(variable,value);
                if(value > 0)
                    fprintf(file,",");
                fprintf(file, "%s",value_name.c_str());
            }
            fprintf(file, "\n");
        }
    } else throw compiler_exception("Could not write UAI map to '%s'", filename.c_str());
}

