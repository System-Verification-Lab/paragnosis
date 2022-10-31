#include "branchnbound.h"
#include <algorithm>
#include <limits>
#include "bnc.h"
#include "options.h"
#include "assert.h"

branchnbound::branchnbound(){
    bn = NULL;
    MAX_UPPER_BOUND = std::numeric_limits<std::size_t>::max();
}

branchnbound::~branchnbound(){
    for(unsigned int i = 0; i < solutions.size(); i++)
        destroy(solutions[i]);
}

void branchnbound::init(bayesnet *bn,const partition_t &kPartition){
    this->bn = bn;
    this->kPartition = &kPartition;
    bound.init(bn);
    const unsigned int VARIABLES = bn->get_nr_variables();
    root = new bbnode;
    root->variable = VARIABLES;
    root->children = 0;
    root->parent = NULL;
    root->upper_bound = 0;
}

void branchnbound::bestguess(){
    ordering_t ordering;
    ordering.generate_variable_ordering(bn,*kPartition);

    // TODO: check for overflow
    MAX_UPPER_BOUND = bnc::get_upper_bound(bn, ordering).get_ui();
}

void branchnbound::print(){
    for(auto it = solutions.begin(); it != solutions.end(); it++){
        bbnode *node = *it;
        printf("%u: ", node->upper_bound);

        std::vector<unsigned int> history;
        while(node->parent != 0){
            history.push_back(node->variable);
            node = node->parent;
        }

        for(auto back_it = history.rbegin(); back_it != history.rend(); back_it++)
            printf("%3u ", *back_it);
        printf("\n");
    }
}

ordering_t branchnbound::get_ordering(){
    ordering_t ordering;

    for(auto it = solutions.begin(); it != solutions.end(); it++){
        bbnode *node = *it;

        while(node->parent != NULL){
            ordering.push_back(node->variable);
            node = node->parent;
        }
        std::reverse(ordering.begin(), ordering.end());

        // TODO select best???
        break;
    }
    return std::move(ordering);
}

void branchnbound::destroy(bbnode *node){
    while(node && node->children == 0){
        bbnode *tmp = node;
        node = node->parent;
        delete tmp;
        counter--;

        if(node)
            node->children--;
    }
}

bbnode* branchnbound::get_ancestor(const unsigned int K, bbnode *node){
    for(unsigned int i = 0; i < K; i++)
        node = node->parent;

    return node;
}

unsigned int branchnbound::get_ancestor_variable(const unsigned int K, bbnode *node){
    return get_ancestor(K, node)->variable;
}

bool branchnbound::is_ancestor(const unsigned int K, const unsigned int variable, bbnode *node){
    for(unsigned int i = 0; i < K; i++)
        node = node->parent;

    return node->variable == variable;
}

size_t branchnbound::get_treewidth(const unsigned int VARIABLES, const unsigned int LEVEL, const unsigned int K){
    size_t width = 1;
    for(unsigned int i = LEVEL; i < LEVEL+(K-1); i++){
        width *= VARIABLES - i;

        if(width < VARIABLES-i)
            width = 0;
    }
    return width;
}

size_t branchnbound::get_treewidth(const unsigned int VARIABLES, const unsigned int K){
    return get_treewidth(VARIABLES, 0, K);
}

void branchnbound::depthfirst(std::queue<bbnode*> &queue, const unsigned int VARIABLES, const unsigned int K){
    std::set<unsigned int> variables;
    for (unsigned int i = 0; i < VARIABLES; i++)
        variables.insert(variables.end(), i);

    MAX_UPPER_BOUND = std::numeric_limits<size_t>::max();
    size_t UPPER_BOUND = std::numeric_limits<size_t>::max();

    bbstate state;
    bound.init(state);
    counter = 1;
    for(unsigned int variable = 0; variable < VARIABLES; variable++){
        root->children++;
        bbnode *node = new bbnode;
        root->child.push_back(node);
        counter++;

        bound.set_value(state, 0);
        bound.condition(state, variable);
        // TODO: check for overflow
        node->upper_bound = bound.get_value(state).get_ui();

        node->children = 0;
        node->variable = variable;
        node->parent = root;

        variables.erase(variable);
        std::queue<bbnode*> local_queue;
        MAX_UPPER_BOUND = std::numeric_limits<size_t>::max();
        depthfirst(local_queue, variables, node, state, K-2);
        if(MAX_UPPER_BOUND <= UPPER_BOUND){
            if(MAX_UPPER_BOUND < UPPER_BOUND){
                while(!queue.empty()){
                    destroy(queue.front());
                    queue.pop();
                }
                queue = std::move(local_queue);
                UPPER_BOUND = MAX_UPPER_BOUND;
            } else {
                while(!local_queue.empty()){
                    queue.push(local_queue.front());
                    local_queue.pop();
                }
            }
        } else {
            while(!local_queue.empty()){
                destroy(local_queue.front());
                local_queue.pop();
            }
        }
        bound.undo(state, variable);
        variables.insert(variable);
    }
}

void branchnbound::depthfirst(std::queue<bbnode*> &queue, std::set<unsigned int> variables, bbnode *root, bbstate &state, const unsigned int K){
    const unsigned int SIZE = variables.size();
    for(unsigned int i = 0; i < SIZE; i++){
        auto it = variables.begin();
        std::advance(it, i);
        unsigned int variable = *it;

        root->children++;
        bbnode *node = new bbnode;
        root->child.push_back(node);
        counter++;

        node->parent = root;
        node->variable = variable;
        node->children = 0;

        bound.set_value(state, root->upper_bound);
        bound.condition(state, variable);
        // TODO: check for overflow
        node->upper_bound = bound.get_value(state).get_ui();
        if(K > 0){
            variables.erase(variable);
            depthfirst(queue, variables, node, state, K - 1);
            variables.insert(variable);
        } else {
            queue.push(node);
            if(node->upper_bound < MAX_UPPER_BOUND)
                MAX_UPPER_BOUND = node->upper_bound;
        }
        bound.undo(state, variable);
    }
}

void branchnbound::lookahead(
        std::queue<bbnode*> &queue,
        const unsigned int VARIABLES,
        const unsigned int K){

    printf("queue size: %u\n", queue.size());
    size_t UPPER_BOUND = std::numeric_limits<size_t>::max();
    std::array< std::set<unsigned int> , 2> ancestors;

    for(unsigned int level = 1; level <= VARIABLES-K; level++){
        const unsigned int CURRENT = level % 2;
        const unsigned int PREVIOUS = (level + 1) % 2;
        const size_t FAMILY_SIZE = get_treewidth(VARIABLES, level, K);
        const size_t WIDTH = get_treewidth(VARIABLES, level+1, K-1);
        const unsigned int ANCESTORS = queue.size() / FAMILY_SIZE;

        assert(queue.size() % FAMILY_SIZE == 0);
        queue.push(NULL);
        // ============================ per level ==============================
        UPPER_BOUND = std::numeric_limits<size_t>::max();
        for(unsigned int a = ANCESTORS; a > 0; a--){

            // skip bad ancestors
            bbnode *check = get_ancestor(K-1, queue.front());
            if(!ancestors[PREVIOUS].empty() && ancestors[PREVIOUS].find(check->variable) == ancestors[PREVIOUS].end()){
                for(unsigned int i = FAMILY_SIZE; i > 0; i--){
                    assert(get_ancestor(K-1, queue.front()) == check);
                    destroy(queue.front());
                    queue.pop();
                }
                continue;
            }

            for(unsigned scnd = 0; scnd < VARIABLES-level; scnd++){

                // determine participating variables by ancestor's predicessors
                std::vector<unsigned int> history;
                std::set<unsigned int> variables;
                for (unsigned int i = 0; i < VARIABLES; i++)
                    variables.insert(variables.end(), i);

                bbnode * const ancestor = get_ancestor(K-2, queue.front());
                bbnode *tmp = ancestor;
                while(tmp->parent != NULL){
                    history.push_back(tmp->variable);
                    variables.erase(tmp->variable);
                    tmp = tmp->parent;
                }

                // restore correct state until 'ancestor'
                bbstate state;
                bound.init(state);

                while(!history.empty()){
                    bound.condition(state, history.back());
                    history.pop_back();
                }

                // ============= per ancestor of Kth generation back ===============
                for(unsigned int i = WIDTH; i > 0; i--){
                    bbnode *descendant = queue.front();
                    queue.pop();

                    // determine participating variables by descendants
                    bbnode *tmp = descendant;
                    while(tmp != ancestor){
                        history.push_back(tmp->variable);
                        variables.erase(tmp->variable);
                        tmp = tmp->parent;
                    }

                    // restore correct state until 'ancestor'
                    for(auto back_it = history.rbegin(); back_it != history.rend(); back_it++)
                        bound.condition(state, *back_it);

                    // =============================================================

                    // expand Kth generation descendants of 'ancestor'
                    for(auto it = variables.begin(); it != variables.end(); it++){
                        unsigned int variable = *it;

                        bound.set_value(state, descendant->upper_bound);
                        bound.condition(state, variable);
                        // TODO: check for overflow
                        size_t upper_bound = bound.get_value(state).get_ui();
                        bound.undo(state, variable);

                        descendant->children++;
                        bbnode *node = new bbnode;
                        descendant->child.push_back(node);
                        counter++;

                        node->upper_bound = upper_bound;
                        node->parent = descendant;
                        node->variable = variable;
                        node->children = 0;
                        queue.push(node);

                        if(upper_bound <= UPPER_BOUND){
                            if(upper_bound < UPPER_BOUND){
                                ancestors[CURRENT].clear();
                                UPPER_BOUND = upper_bound;
                            }
                            ancestors[CURRENT].insert(ancestor->variable);
                        }
                    }


                    // ============================================================

                    // restore correct state until 'ancestor'
                    for(auto it = history.begin(); it != history.end(); it++)
                        bound.undo(state, *it);

                    // restore removed decendants
                    tmp = descendant;
                    while(tmp != ancestor){
                        history.pop_back();
                        variables.insert(tmp->variable);
                        tmp = tmp->parent;
                    }

                }
            }
        }
        assert(queue.front() == NULL);
        queue.pop();

    }

    // filter out worst solutions
    queue.push(NULL);
    bbnode *node;
    while((node = queue.front()) != NULL){
        queue.pop();
        if(node->upper_bound != UPPER_BOUND)
            destroy(node);
        else queue.push(node);
    }
    queue.pop();


        // std::set<unsigned int> ancestors[2];
        // size_t UPPER_BOUND = std::numeric_limits<size_t>::max();

        // const unsigned int i = level%2;
        // queue.push(NULL);
        // bbnode *current_node;
        // while((current_node = queue.front()) != NULL){

        //     unsigned int ancestor = get_ancestor(current_node, K);
        //     for(unsigned int i =
        //         queue.pop();


        //     // determine participating variables by ancestors
        //     std::vector<unsigned int> history;
        //     std::set<unsigned int> variables;
        //     for (unsigned int i = 0; i < VARIABLES; i++)
        //         variables.insert(variables.end(), i);

        //     // perform upper bound computation
        //     if(level == CURRENT_LEVEL+K || variables.size() == 1){
        //     if(get_ancestor(current_node,K)

        //         bbnode *ancestor = current_node;
        //         while(ancestor->parent != NULL){
        //             history.push_back(ancestor->variable);
        //             variables.erase(ancestor->variable);
        //             ancestor = ancestor->parent;
        //         }

        //         // restore correct state
        //         bound.reset();
        //         for(auto back_it = history.rbegin(); back_it != history.rend(); back_it++)
        //             bound.condition(*back_it);

        //         // expand level
        //         for(auto it = variables.begin(); it != variables.end(); it++){
        //             unsigned int variable = *it;

        //             bound.set_value(current_node->upper_bound);
        //             bound.condition(variable);
        //             size_t upper_bound = bound.get_value(state);
        //             bound.undo(variable);

        //             current_node->children++;
        //             bbnode *node = new bbnode;
        //             node->upper_bound = upper_bound;
        //             node->parent = current_node;
        //             node->variable = variable;
        //             node->children = 0;
        //             queue.push(node);

        //             if(upper_bound <= UPPER_BOUND){
        //                 if(upper_bound < UPPER_BOUND){
        //                     ancestors.clear();
        //                     next.clear();
        //                     UPPER_BOUND = upper_bound;
        //                 }

        //                 unsigned int ancestor = get_ancestor(K, node);
        //                 if((ancestors.insert(ancestor)).second)
        //                     next.push_back(ancestor);
        //             }
        //         }
        //     } else {
        //          // expand level
        //         for(auto it = variables.begin(); it != variables.end(); it++){
        //             current_node->children++;
        //             bbnode *node = new bbnode;
        //             node->upper_bound = 0;
        //             node->parent = current_node;
        //             node->variable = *it;
        //             node->children = 0;
        //             queue.push(node);
        //         }
        //     }
        // }
        // queue.pop();
}

// void branchnbound::lookahead(std::queue<bbnode*> &queue, std::vector<unsigned int> &next){
//
//     // determine level
//     bbnode *ancestor = queue.front();
//     unsigned int level = 0;
//     while(ancestor->parent != NULL){
//         ancestor = ancestor->parent;
//         level++;
//     }
//
//     size_t width = get_treewidth(VARIABLES, CURRENT_LEVEL, level);
//     const unsigned int ROOTS = in_queue.size() / width;
//
//     for(; level <= CURRENT_LEVEL+K; level++){
//         width *= VARIABLES - level;
//         for(unsigned int i = width * ROOTS;; i > 0; i++){
//             bbnode *node = in_queue.top();
//             in_queue.pop();
//         }
//     }
//
//
//     queue.push(NULL);
//     bbnode *node;
//     while((node = queue.top()) != NULL){
//         queue.pop();
//
//
//     }
//     queue.pop();
// }

void branchnbound::lookahead(const unsigned int K){
    const unsigned int VARIABLES = bn->get_nr_variables();

    std::queue<bbnode*> queue;
    depthfirst(queue, VARIABLES, K);
    printf("queue size: %u\n", queue.size());

    lookahead(queue, VARIABLES, K);
    printf("queue size: %u\n", queue.size());

    unsigned int SOLUTIONS = 0;
    while(!queue.empty()){
        SOLUTIONS++;
        std::vector<unsigned int> history;
        bbnode *node = queue.front();
        while(node){
            history.push_back(node->variable);
            node = node->parent;
        }

        printf("%3u: ", SOLUTIONS);
        for(auto it = history.rbegin(); it != history.rend(); it++)
            printf("%2u ", *it);
        printf("\n");

        solutions.push_back(queue.front());
        queue.pop();
    }
}

void branchnbound::bestfirst(){
    const unsigned int VARIABLES = bn->get_nr_variables();
    queue.push(root);
    bbstate state;
    bound.init(state);

    while(!queue.empty()){
        bbnode *current_node = queue.top();
        queue.pop();

        // determine participating variables by ancestors
        std::vector<unsigned int> history;
        std::set<unsigned int> variables;
        for (unsigned int i = 0; i < VARIABLES; i++)
            variables.insert(variables.end(), i);

        bbnode *ancestor = current_node;
        while(ancestor->parent != NULL){
            history.push_back(ancestor->variable);
            variables.erase(ancestor->variable);
            ancestor = ancestor->parent;
        }

        // restore correct state
        bound.reset(state);
        for(auto back_it = history.rbegin(); back_it != history.rend(); back_it++)
            bound.condition(state, *back_it);

        for(auto it = variables.begin(); it != variables.end(); it++){
            unsigned int variable = *it;
            bound.set_value(state, current_node->upper_bound);
            bound.condition(state, variable);
            // TODO: check for overflow
            size_t upper_bound = bound.get_value(state).get_ui();
            bound.undo(state, variable);

            if(upper_bound <= MAX_UPPER_BOUND){
                current_node->children++;

                bbnode *node = new bbnode;
                counter++;
                current_node->child.push_back(node);
                node->upper_bound = upper_bound;
                node->parent = current_node;
                node->variable = variable;
                node->children = 0;

                if(variables.size() == 1){
                    if(upper_bound < MAX_UPPER_BOUND)
                        solutions.clear();

                    solutions.push_back(node);
                    MAX_UPPER_BOUND = upper_bound;
                } else
                    queue.push(node);
            }

        }
        destroy(current_node);
    }
}

