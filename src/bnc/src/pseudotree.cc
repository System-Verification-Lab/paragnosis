#include "pseudotree.h"
#include "options.h"
#include <vector>
#include <string>
#include <daoopt/Main.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <wordexp.h>
#include <unordered_map>
#include "exceptions.h"
#include <stack>
#include "timer.h"
#include "bound.h"

std::unordered_map<unsigned int,unsigned int> bn_to_daoopt_map;
std::unordered_map<unsigned int,unsigned int> daoopt_to_bn_map;

bool daoopt::Problem::loadBayesnet(::bayesnet *bn, const ::partition_t *kPartition){
    m_name = "problem-name";

    bn_to_daoopt_map.clear();
    daoopt_to_bn_map.clear();
    for(auto it = kPartition->set.begin(); it != kPartition->set.end(); it++){
        unsigned int daoopt_variable = bn_to_daoopt_map.size();
        bn_to_daoopt_map[*it] = daoopt_variable;
        daoopt_to_bn_map[daoopt_variable] = *it;
    }
    for(auto it = kPartition->cutset.begin(); it != kPartition->cutset.end(); it++){
        unsigned int daoopt_variable = bn_to_daoopt_map.size();
        bn_to_daoopt_map[*it] = daoopt_variable;
        daoopt_to_bn_map[daoopt_variable] = *it;
    }
    unsigned int daoopt_variable = bn_to_daoopt_map.size();
    bn_to_daoopt_map[std::numeric_limits<Variable>::max()] = daoopt_variable;
    daoopt_to_bn_map[kPartition->cutset.size() + kPartition->set.size()] = std::numeric_limits<Variable>::max();

    const auto kNrFunctions = kPartition->set.size();
    const auto kNrVariables = kPartition->set.size() + kPartition->cutset.size();

    m_task = TASK_MAX;
    m_prob = PROB_MULT;
    m_n = kNrVariables;   // No. of variables

    const auto *kDomains = bn->get_states();
    for(unsigned int variable = 0; variable < kNrVariables; variable++)
        m_domains.push_back(kDomains[daoopt_to_bn_map[variable]]);

    m_k = *std::max_element(m_domains.begin(), m_domains.end());
    m_c = kNrFunctions;

    // Scope information for functions
    vector< vector<int> > scopes;
    scopes.reserve(m_c);
    m_r = -1;
    for(auto it = kPartition->set.begin(); it != kPartition->set.end(); it++){
        const auto kBnVariable = *it;
        vector<int> scope;
        const auto kNrParents = bn->get_parent_size(kBnVariable);
        const auto *kParents = bn->get_parent(kBnVariable);
        for(unsigned int parent = 0; parent < kNrParents; parent++)
            scope.push_back(bn_to_daoopt_map[kParents[parent]]);
        scope.push_back(bn_to_daoopt_map[kBnVariable]);

        int scope_size = scope.size();
        m_r = std::max(m_r,scope_size);
        scopes.push_back(scope);
    }


    // Read functions
    int i = 0;
    for(auto it = kPartition->set.begin(); it != kPartition->set.end(); it++, i++){
        const auto kBnVariable = *it;

        // create set version of the scope (ordered)
        set<int> scopeSet(scopes[i].begin(), scopes[i].end());

        const auto* kProbabilities = bn->get_cpt(kBnVariable);
        const auto kNrProbabilities = bn->get_cpt_size(kBnVariable);

        double* table = NULL;
        //table = new double[kNrProbabilities];
        //std::copy(kProbabilities, kProbabilities + kNrProbabilities, table);

        Function* f = new FunctionBayes(i,this,scopeSet,table,kNrProbabilities);
        m_functions.push_back(f);
    } // All function tables read


    m_e = 0; // Read evidence? no
    m_m = 0;

    return true;
}


namespace bnc {

class DaooptAPI : private daoopt::Main {
    public:
        bool GenerateTree(bayesnet *bn, const partition_t &);
        bool SetTree(bayesnet *bn, const partition_t &, const ordering_t &, const bool kForceChain = false);
        bool StorePseudoTree(PseudoTreeNode &, ordering_t &, size_t&, size_t&);
        bool StoreMinFillOrdering(ordering_t &);
        bool findOrLoadOrdering(const ordering_t &);
    private:
        // Altered Main functions
        bool loadProblem(bayesnet *bn, const partition_t &);
};

using namespace daoopt;

bool DaooptAPI::loadProblem(bayesnet *bn, const partition_t &kPartition) {
    m_problem.reset(new Problem);

    // load problem file
    m_problem->loadBayesnet(bn,&kPartition);
    //if (!m_problem->parseUAI(m_options->in_problemFile, m_options->in_evidenceFile,
    //                         m_options->in_mmapFile))
    //  return false;
    cout << "Created problem with " << m_problem->getN()
        << " variables and " << m_problem->getC() << " functions." << endl;

    if (!m_options->in_mmapFile.empty()) {
        cout << "Read list of " << m_problem->getM() << " variables for MMAP query." << endl;
    }

    // Remove evidence variables (only works if probabilities are stored)
    // m_problem->removeEvidence();
    //cout << "Removed evidence, now " << m_problem->getN()
    //    << " variables and " << m_problem->getC() << " functions." << endl;

    #if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
    if (!m_options->par_solveLocal) {
        // Re-output problem file for parallel processing
        m_options->out_reducedFile = string("temp_prob.") + m_options->problemName
            + string(".") + m_options->runTag + string(".gz");
        m_problem->writeUAI(m_options->out_reducedFile);
        cout << "Saved problem to file " << m_options->out_reducedFile << endl;
    }
    #else
    // Output reduced network?
    if (!m_options->out_reducedFile.empty()) {
        cout << "Writing reduced network to file " << m_options->out_reducedFile << endl;
        m_problem->writeUAI(m_options->out_reducedFile);
    }
    #endif

    // Some statistics
    cout << "Global constant:\t" << SCALE_LOG(m_problem->globalConstInfo()) << endl;
    cout << "Max. domain size:\t" << (int) m_problem->getK() << endl;
    cout << "Max. function arity:\t" << m_problem->getR() << endl;

    return true;
}

bool DaooptAPI::SetTree(bayesnet *bn, const partition_t &kPartition, const ordering_t &kOrdering, const bool kForceChain){
    // simulate argv argc
    wordexp_t p;
    p.we_offs = 1;
    wordexp("-f foo.uai", &p, WRDE_DOOFFS);
    p.we_wordv[0] = "foo";

    char **argv = p.we_wordv;
    int argc = p.we_wordc+1;

    // disable cout
    #ifndef DEBUG
    std::cout.setstate(std::ios_base::failbit);
    #endif

    // run main
    if (!parseOptions(argc, argv))
        return false;
    if (!loadProblem(bn, kPartition))
        return false;
    if (!runSLS())
        return false;

    Graph g(m_problem->getN());
    const vector<Function*>& fns = m_problem->getFunctions();
    for (vector<Function*>::const_iterator it = fns.begin(); it != fns.end(); ++it) {
        g.addClique((*it)->getScopeVec());
    }
    m_pseudotree.reset(new Pseudotree(m_problem.get(), m_options->subprobOrder));

    vector<int> elim;
    elim.reserve(kOrdering.size());
    for(auto it = kOrdering.begin(); it != kOrdering.end(); it++)
        elim.push_back(bn_to_daoopt_map[*it]);
    std::reverse(elim.begin(), elim.end());

    if(kForceChain)
        m_pseudotree->buildChain(g, elim, m_options->cbound);
    else {
        m_pseudotree->build(g, elim, m_options->cbound);
    }

    wordfree(&p);

    // enable cout
    #ifndef DEBUG
    std::cout.clear();
    #endif

    // check if first variable is indeed dummy variable n+1
    assert(m_pseudotree->getRoot() != NULL);
    assert((size_t)m_pseudotree->getRoot()->getVar() == kPartition.set.size() + kPartition.cutset.size());

    return true;
}


bool DaooptAPI::GenerateTree(bayesnet *bn, const partition_t &kPartition){
    // simulate argv argc
    wordexp_t p;
    p.we_offs = 1;
    wordexp("-f foo.uai --cvo --adaptive --orderIter 500 --orderTime 60", &p, WRDE_DOOFFS);
    p.we_wordv[0] = "foo";

    char **argv = p.we_wordv;
    int argc = p.we_wordc+1;

    // disable cout
    #ifndef DEBUG
    std::cout.setstate(std::ios_base::failbit);
    #endif

    // run main
    if (!parseOptions(argc, argv))
        return false;
    if (!loadProblem(bn, kPartition))
        return false;
    if (!runSLS())
        return false;
    if (!daoopt::Main::findOrLoadOrdering())
        return false;

    wordfree(&p);

    // enable cout
    #ifndef DEBUG
    std::cout.clear();
    #endif

    // check if first variable is indeed dummy variable n+1
    assert(m_pseudotree->getRoot() != NULL);
    assert((size_t)m_pseudotree->getRoot()->getVar() == kPartition.set.size() + kPartition.cutset.size());

    return true;
}


void CopyPseudoTree(PseudoTreeNode &node, PseudotreeNode* daoopt_node, ordering_t &ordering) {
    node.variable = daoopt_to_bn_map[daoopt_node->getVar()];
    if(node.variable < std::numeric_limits<Variable>::max())
        ordering.push_back(node.variable);
    node.children.resize(daoopt_node->getChildren().size());

    auto daoopt_it = daoopt_node->getChildren().begin();
    auto it = node.children.begin();
    while(daoopt_it != daoopt_node->getChildren().end()){
        CopyPseudoTree(*it,*daoopt_it, ordering);
        ++daoopt_it;
        ++it;
    }
}

bool DaooptAPI::StorePseudoTree(PseudoTreeNode &node, ordering_t &ordering, size_t &height, size_t &width){
    height = m_pseudotree->getHeight();
    width = m_pseudotree->getWidth();
    auto *daoopt_node = m_pseudotree->getRoot();
    assert(daoopt_node != NULL);
    assert(daoopt_node->getChildren().size() >= 1);

    CopyPseudoTree(node, daoopt_node, ordering);

    // the first variable is a dummy variable connnecting disconnected components
    node.variable = std::numeric_limits<Variable>::max();


    return true;
}

bool DaooptAPI::StoreMinFillOrdering(ordering_t &ordering){
    const auto &elim_order = m_pseudotree->getElimOrder();

    ordering.clear();
    for(auto it = elim_order.begin(); it != elim_order.end(); it++){
        auto bn_variable = daoopt_to_bn_map[*it];
        if(bn_variable < std::numeric_limits<Variable>::max())
            ordering.push_back(bn_variable);
    }
    std::reverse(ordering.begin(), ordering.end());

    return true;
}


PseudoTree::PseudoTree(){
    Clear();
}

void PseudoTree::SetChain(bayesnet *bn, const partition_t &kPartition, const ordering_t &kOrdering){
    DaooptAPI api;
    const bool kForceChain = true;
    api.SetTree(bn, kPartition, kOrdering, kForceChain);
    api.StorePseudoTree(root_,dfs_ordering_,height_,width_);
    api.StoreMinFillOrdering(minfill_ordering_);
}

void PseudoTree::SetTree(bayesnet *bn, const partition_t &kPartition, const ordering_t &kOrdering){
    DaooptAPI api;
    api.SetTree(bn, kPartition, kOrdering);
    api.StorePseudoTree(root_,dfs_ordering_,height_,width_);
    api.StoreMinFillOrdering(minfill_ordering_);
}

void PseudoTree::GenerateTree(bayesnet *bn, const partition_t &kPartition){
    DaooptAPI api;
    api.GenerateTree(bn, kPartition);
    api.StorePseudoTree(root_,dfs_ordering_,height_,width_);
    api.StoreMinFillOrdering(minfill_ordering_);
}

const ordering_t& PseudoTree::GetPreOrdering() const {
    return dfs_ordering_;
}

const ordering_t& PseudoTree::GetMinFillOrdering() const {
    return minfill_ordering_;
}

const PseudoTreeNode& PseudoTree::GetPseudoTree() const {
    return root_;
}

const size_t PseudoTree::GetSize() const {
    return dfs_ordering_.size();
}

const size_t PseudoTree::GetNrVariables() const {
    return dfs_ordering_.size();
}

bool PseudoTree::Empty() const {
    return root_.children.size() == 0;
}

void PseudoTree::Write() const {
    std::string filename = files.get_filename(PSEUDO_ORDERING);
    Write(filename);
}

void PseudoTree::CreateTreeString(const bnc::PseudoTreeNode* node, std::ostringstream& oss) const {
    oss << "(" << node->variable;
    for (auto it = node->children.begin(); it != node->children.end(); ++it)
        CreateTreeString(&(*it), oss);
    oss << ")";
}

void PseudoTree::Write(std::string filename) const {
    FILE *file = fopen(filename.c_str(), "w");
    if(file){

        std::ostringstream oss;
        CreateTreeString(&root_, oss);  // recursive, call on root node
        fprintf(file, "%s\n", oss.str().c_str());

        fclose(file);
    } else compiler_write_exception("Could not open pseudotree file '%s'", filename.c_str());
}

void PseudoTree::Read(){
    std::string filename = files.get_filename(PSEUDO_ORDERING);
    Read(filename);
}

void PseudoTree::Clear(){
    root_.variable = std::numeric_limits<Variable>::max();
    root_.children.clear();

    dfs_ordering_.clear();
    minfill_ordering_.clear();

    width_ = 0;
    height_ = 0;
}

void PseudoTree::Read(std::string filename){
    ifstream file(filename);
    std::string line;

    std::getline (file,line);

    Clear();
    std::stack<PseudoTreeNode*> s;
    s.push(&root_);

    bool root = true;
    unsigned int state = 1;
    unsigned int i = 0;
    while(i < line.length() && !s.empty()){
        PseudoTreeNode *node = s.top();

        if(line[i] == '(' && state == 1){
            state = 0;
            ++i;
        } else if (line[i] == ')' && state == 1){
            s.pop();
            state = 1;
            ++i;
        } else if(state == 0) {
            if(!std::isdigit(line[i]))
                compiler_read_exception("Read '%c' at pos %u, while expected number", line[i], i);
            unsigned int num = 0;
            while(i < line.length() && std::isdigit(line[i])){
                num = num *10 + (line[i]-'0');
                ++i;
            }

            if(!root){
                dfs_ordering_.push_back(num);
                node->children.resize(node->children.size()+1);
                auto *child = &(node->children.back());
                child->variable = num;
                s.push(child);
            } else {
                node->variable = num;
                root = false;
            }

            state = 1;
        } else {
            compiler_read_exception("Unexpected char '%c' at pos %u", line[i], i);
        }
    }
    if(i < line.length())
        compiler_read_exception("Read pseudo tree but there are %lu characters left in string", line.length() - i);

    if(!s.empty())
        compiler_read_exception("There are %lu missing closing brackets", s.size());

}

size_t PseudoTree::GetWidth() const {
    return width_;
}

size_t PseudoTree::GetHeight() const {
    return height_;
}

ordering_t PseudoTree::MinFill(bayesnet *bn, const partition_t &kPartition){

    DaooptAPI api;
    api.GenerateTree(bn, kPartition);

    ordering_t ordering;
    api.StoreMinFillOrdering(ordering);

    return ordering;
}

}

