#include "bnc.h"
#include "ite.h"
#include "sat.h"
#include "exceptions.h"
#include "options.h"
#include "misc.h"
#include <algorithm>

namespace bnc {

size_t size(node *n){
    size_t size = 0;

    node::queue q;
    node::set done;
    done.insert(NULL);

    node *f_0 = NULL;
    node *f_1 = NULL;

    q.push(n);
    while(!q.empty()){
        node *n = q.front();
        q.pop();
        if(!done.contains(n)){
            done.insert(n);
            if(!node::is_terminal(n)){
                size++;
                q.push(n->e);
                q.push(n->t);
            } else {
                if(node::is_satisfiable(n))
                    f_1 = n;
                else f_0 = n;
            }
        }
    }
    if(f_1)
        size++;
    if(f_0)
        size++;

    return size;
}

size_t operators(node *n){
    size_t OR_OPERATORS = 0;
    size_t AND_OPERATORS = 0;
    size_t EXTRA_AND_OPERATORS = 0;

    node::queue q;
    node::set done;
    done.insert(NULL);

    q.push(n);
    while(!q.empty()){
        node *n = q.front();
        q.pop();
        if(!done.contains(n)){
            done.insert(n);
            if(!node::is_terminal(n)){
                OR_OPERATORS++;
                AND_OPERATORS++;
                if(n->W)
                    EXTRA_AND_OPERATORS += n->W->size();

                q.push(n->e);
                q.push(n->t);
            }
        }
    }
    return OR_OPERATORS + AND_OPERATORS + EXTRA_AND_OPERATORS;
}

node* uncollapse(manager_t *manager, node *n, const unsigned int partition_id){
    // create collapse variable expectation map
    std::map<literal_t,literal_t> expect;
    const support_t support = get_support(manager,n);
    ordering_t ordering = create_ordering(manager, partition_id, support);
    bayesgraph& g = manager->get_bayesgraph();
    const std::vector<unsigned int>& l2v = g.get_literal_to_variable();
    for(unsigned int i = 0; i < ordering.size(); i++){
        literal_t current = ordering[i], next = ordering[i+1];
        if(l2v[current] != l2v[next])
            expect[ordering[i]] = 0;
        else expect[ordering[i]] = ordering[i+1];
    }

    node *root = copy(manager,n);

    if(node::is_terminal(root))
        return root;

    node::stack s;
    node::set done;
    s.push(n);
    while(!s.empty()){
        node *n = s.top();
        if(!done.contains(n)){
            if(!node::is_terminal(n)){
                literal_t l = expect[n->l];
                if(l != 0 && n->e->l != l){
                    // insert deleted node
                    node *tmp = n->e;
                    node *&e = n->e;

                    e = node::reference(create(manager));
                    e->l = l;
                    e->W = bnc::create_weights(manager);
                    *(e->W) = *(n->W);
                    e->t = node::reference(n->t);
                    e->e = tmp;
                }
            }
            s.push(n->e);
            s.push(n->t);
            done.insert(n);
        }
        s.pop();
    }

    return root;
}

inline bool merge(manager_t *manager, node::table &table, node *&n){
    node *canonical = table.find_or_insert(n);
    if(canonical){
        node::dereference(n);
        destroy(manager,n);
        n = node::reference(canonical);
        return true;
    }
    return false;
}

void recursive_merge(manager_t *manager, node::table &table, node *n){
    if(node::is_terminal(n)){
        merge(manager, table, n);
        return;
    }

    node::reference_stack stack;
    node::set done;
    stack.push(&n);
    while(!stack.empty()){
        node *&n = *stack.top();
        if(!done.contains(n->t)){
            done.insert(n->t);
        } else if(!done.contains(n->e)){
            done.insert(n->e);
        } else {
            merge(manager, table, n);
            stack.pop();
        }
    }
}

inline bool collapse_node(manager_t *manager, node *&n){
    #ifdef DEBUG
    if(node::is_terminal(n))
        throw compiler_debug_exception("must not apply collapse on terminal node");
    #endif

    if(node::collapse_candidate(n)){
        node *collapsed = n->e;
        n->e = node::reference(n->e->e);
        node::dereference(collapsed);
        destroy(manager, collapsed);
        return true;
    }
    return false;
}

inline bool collapse(manager_t *manager, node::table &table, node *&n, ite_t *ite = NULL){
    #ifdef DEBUG
    if(node::is_terminal(n))
        throw compiler_debug_exception("must not apply collapse on terminal node");
    #endif

    if(node::collapse_candidate(n)){
        node *collapsed = n->e;
        n->e = node::reference(n->e->e);
        node::dereference(collapsed);
        if(collapsed->ref == 0){
            table.erase(collapsed);
            if(ite)
                ite->remove_cache(collapsed);
        }
        destroy(manager, collapsed);
        return true;
    }
    return false;
}

void collapse(manager_t *manager, node *n){
    if(node::is_terminal(n))
        return;

    node::table table;
    node::reference_stack stack;
    node::set done;
    stack.push(&n);
    while(!stack.empty()){
        node *&n = *stack.top();
        if(!done.contains(n)){
            done.insert(n);
            if(!node::is_terminal(n)){
                stack.push(&(n->e));
                stack.push(&(n->t));
                continue;
            }
        }

        if(!node::is_terminal(n))
            collapse(manager, table, n);

        merge(manager, table, n);
        stack.pop();
    }
}

void write_pdf(std::string aux, bool wait = false){
    std::string filename = files.get_filename_c(DOT,aux);
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-result"
    if(wait)
        system(stringf("dot2pdf %s 2>/dev/null 1>&2",filename.c_str()).c_str());
    else
        system(stringf("dot2pdf %s 2>/dev/null 1>&2 &",filename.c_str()).c_str());
    #pragma GCC diagnostic pop
}

void write_ps(std::string aux){
    std::string filename = files.get_filename_c(DOT,aux);
    exec(stringf("dot -Tps %s -o %s.ps", filename.c_str(), filename.c_str()).c_str());
}

void write_dot_ps(manager_t *manager, std::vector<node*>& n, std::string aux, node *current,bool wait=false){
    if(n.size() == 1)
        write_dot(manager, n[0]);
    else for(unsigned int i = 0; i < n.size(); i++){
        write_dot(manager, n[i], stringf("%s%u",aux.c_str(),i));
        write_pdf(stringf("%s%u",aux.c_str(),i),wait);
    }
}

void write_dot(manager_t *manager, std::vector<node*>& n, node *current){
    if(n.size() == 1)
        write_dot(manager, n[0]);
    else for(unsigned int i = 0; i < n.size(); i++)
        write_dot(manager, n[i], stringf("%u",i), i);
}

void write_dot(manager_t *manager, node *n, std::string aux, unsigned int id, node *current){
    const bayesgraph &g = manager->get_bayesgraph();

    FILE *file = fopen(files.get_filename_c(DOT,aux), "w");
    if(file){
        std::queue<node*> q;
        std::map<node*,bool> done;
        fprintf(file, "digraph %s {\n", "WPBDD");


        fprintf(file, "    center = true;\n");
        fprintf(file, "    edge [dir = none];\n");

        std::map< unsigned int, std::vector<node*> > level;
        if(n){
            q.push(n);
            while(!q.empty()){
                node *n = q.front();
                q.pop();

                if(!done[n]){
                    done[n] = true;
                    if(node::is_terminal(n))
                        level[0].push_back(n);
                    else level[n->l].push_back(n);

                    if(n->t) q.push(n->t);
                    if(n->e) q.push(n->e);
                }
            }
        }

        #ifdef DEBUG
        // ==== build numeric index for each node ====
        // NOTE: index is assigned according to pre-order
        std::unordered_map<node*,unsigned int> index;
        node::stack s;
        s.push(n);
        unsigned int terminals = 0;
        while(!s.empty()){
            node *n = s.top();
            s.pop();

            if(n && index.find(n) == index.end()){
                if(node::is_terminal(n)){
                    index[n] = node::is_satisfiable(n);
                    terminals++;
                } else {
                    index[n] = index.size()+1-terminals;
                    s.push(n->e);
                    s.push(n->t);
                }
            }
        }
        #endif



        #ifdef VERBOSE
        const bool add_literals = true;
        #else
        const bool add_literals = false;
        #endif

        fprintf(file, "\n    // ==================== NODES ====================\n\n");
        for(auto it = level.begin(); it != level.end(); it++){
            if(it->second.size() > 0){
                if(it->first == 0){
                    fprintf(file, "    {   // TERMINALS\n");
                    fprintf(file, "        rank = same;\n");
                    fprintf(file, "        node [shape=box];\n");
                    if(add_literals)
                        fprintf(file, "        \"TERMINALS\";\n");
                    for(auto itn = it->second.begin(); itn != it->second.end(); itn++){
                        node *n = *itn;
                        if(node::is_terminal(n)){
                            if(node::is_satisfiable(n)){
                                #ifdef DEBUG
                                //fprintf(file, "        %lu [label=\"T (%d:%x)\"", n, n->ref, n);
                                fprintf(file, "        %lu [label=\"%u:%p\\nl:%d, ref:%u\\nT\"", (size_t) n, index[n], n, n->l, n->ref);
                                #else
                                fprintf(file, "        %lu [label=\"T\"", (size_t) n);
                                #endif
                                if(current && n == current)
                                    fprintf(file, ",fontcolor=\"red\",color=\"red\"");
                                fprintf(file, "];\n");
                            } else {
                                #ifdef DEBUG
                                //fprintf(file, "        %lu [label=\"F (%d:%x)\"", n, n->ref, n);
                                fprintf(file, "        %lu [label=\"%u:%p\\nl:%d, ref:%u\\nF\"", (size_t) n, index[n], n, n->l, n->ref);
                                #else
                                fprintf(file, "        %lu [label=\"F\"", (size_t) n);
                                #endif
                                if(current && n == current)
                                    fprintf(file, ",fontcolor=\"red\",color=\"red\"");
                                fprintf(file, "];\n");
                            }
                        }
                    }
                    fprintf(file, "    }\n");
                } else {
                    literal_t l = it->first;
                    if(g.defined()){
                        unsigned int variable = g.get_literal_to_variable()[l];
                        unsigned int value = l - g.get_variable_to_literal()[variable];
                        std::string node_name = g.get_bayesnet()->get_node_name(variable);
                        std::string value_name = g.get_bayesnet()->get_node_value_name(variable,value);
                        fprintf(file, "    {   // literal %d (%s=%s)\n", l, node_name.c_str(), value_name.c_str());
                        fprintf(file, "        rank = same;\n");
                        fprintf(file, "        node [ label = \"%s=%s\" ]\n", node_name.c_str(), value_name.c_str());
                        if(add_literals)
                            fprintf(file, "        literal_%d [shape = plaintext];\n", l);
                        for(auto itn = it->second.begin(); itn != it->second.end(); itn++){
                            node* n = *itn;
                            fprintf(file, "        %lu", (size_t) n);
                            #ifdef DEBUG
                            //fprintf(file, " [label=\"%s=%s (%d:%x)\"]", node_name.c_str(), value_name.c_str(), n->ref, n);
                            fprintf(file, " [label=\"%u:%p\\nl:%d, ref:%u\\n%s=%s\"]", index[n], n, n->l, n->ref, node_name.c_str(), value_name.c_str());
                            #endif
                            if(current && n == current)
                                fprintf(file, " [fontcolor=\"red\",color=\"red\"]");
                            fprintf(file, ";\n");
                        }
                        fprintf(file, "    }\n");
                    } else {
                        fprintf(file, "    {   // literal %d\n", l);
                        fprintf(file, "        rank = same;\n");
                        if(add_literals)
                           fprintf(file, "        literal_%d [shape = plaintext];\n", l, l);
                        for(auto itn = it->second.begin(); itn != it->second.end(); itn++){
                            fprintf(file, "        %lu", (size_t) n);
                            #ifdef DEBUG
                            fprintf(file, " [label=\"%u\" (%d:%x)\"]", (size_t) n, l, n->ref, n);
                            #endif
                            if(current && n == current)
                                fprintf(file, "[fontcolor=\"red\",color=\"red\"]");
                            fprintf(file, ";\n");
                        }
                        fprintf(file, "    }\n");
                    }

                }
            }
        }
        // print labels on left side
        if(add_literals){
            fprintf(file, "\n    // ==================== ORDER ====================\n\n");
            fprintf(file, "    { \n");
            fprintf(file, "        node [shape = plaintext];\n");
            fprintf(file, "        edge [style = invis];\n");
            fprintf(file, "        \"TERMINALS\" [style = invis];\n\n");
            const ordering_t &ordering = manager->get_ordering(id);
            for(unsigned int i = 0; i < ordering.size(); i++){
                literal_t l = ordering[i];
                unsigned int variable = g.get_literal_to_variable()[l];
                unsigned int value = l - g.get_variable_to_literal()[variable];
                std::string node_name = g.get_bayesnet()->get_node_name(variable);
                std::string value_name = g.get_bayesnet()->get_node_value_name(variable,value);
                if(level[l].size() > 0)
                    fprintf(file, "        literal_%d [ label = \" Literal %d: %s=%s \\l\"];\n", l, l, node_name.c_str(), value_name.c_str());
            }

            for(unsigned int i = 0; i < ordering.size(); i++){
                literal_t l = ordering[i];
                if(level[l].size() > 0)
                    fprintf(file, "        literal_%d ->\n", l);
                // if(level[l].size() == 0){
                // } else fprintf(file, "        \"(Literal %d)\" ->\n", l);
            }
            fprintf(file, "        \"TERMINALS\";\n");
            fprintf(file, "    }\n");
        }
        fprintf(file, "\n    // ==================== EDGES ====================\n\n");
        done.clear();
        if(n){
            q.push(n);
            while(!q.empty()){
                node *n = q.front();
                q.pop();

                if(!done[n]){
                    done[n] = true;
                    if(n->t){
                        q.push(n->t);
                        fprintf(file, "    %lu -> %lu", (size_t) n, (size_t) n->t);
                        if(n->W){
                            fprintf(file, " [fontsize=8,label=\"");
                            for(auto it = n->W->begin(); it != n->W->end(); it++){
                                weight_t w = *it;

                                if(it != n->W->begin())
                                    fprintf(file, " * ");

                                if(g.defined()){
                                    if(OPT_DETERMINISM && (w == 1 || w == 0))
                                        fprintf(file, "W%u(%.2lf)", w, (double) w);
                                    else
                                        fprintf(file, "W%u(%.2lf)", w, g.get_weight_to_probability()[w-(g.get_nr_literals()+1)]);
                                }
                                else
                                    fprintf(file, "W%u", w);
                            }
                            fprintf(file, "\"]");
                        }
                        fprintf(file, ";\n");
                    }
                    if(n->e){
                        q.push(n->e);
                        fprintf(file, "    %lu -> %lu [style=dotted];\n", (size_t) n, (size_t) n->e);
                    }
                }
            }
        }

        fprintf(file, "}\n");
        fclose(file);
    } else throw compiler_write_exception("Could not open dot file '%s'",files.get_filename_c(DOT,aux));
}

void write_dot_pdf(manager_t *manager, node *n, std::string aux, node *current,bool wait){
    write_dot(manager, n, aux, 0, current);
    write_pdf(aux,wait);
}

void write_state(manager_t *manager, node *root, node *n, node *bdd1, node *bdd2, ite_t &ite,bool wait =false){
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-result"
    system("killall dot 2>/dev/null");
    #pragma GCC diagnostic pop

    write_dot_pdf(manager, root, "bdd", n,wait);
    write_dot_pdf(manager, bdd1, "bdd1", ite.get_current(0),wait);
    write_dot_pdf(manager, bdd2, "bdd2", ite.get_current(1),wait);
}

node* copy(manager_t *manager, node *in){
    if(!in)
        return NULL;

    node::map in_to_out;
    node::stack stack_in;
    node::reference_stack stack_out;

    node *out = NULL;
    stack_in.push(in);
    stack_out.push(&out);
    while(!stack_in.empty()){
        node *n = stack_in.top();
        node *&o = *stack_out.top();
        stack_in.pop(); stack_out.pop();

        node *&hit = in_to_out.find_or_reference(n);
        if(hit){
            o = node::reference(hit);
        } else {
            if(n){
                o = node::reference(create(manager));
                hit = o;

                o->l = n->l;
                if(n->W){
                    o->W = create_weights(manager);
                    *o->W = *n->W;
                }

                stack_in.push(n->e);
                stack_in.push(n->t);
                stack_out.push(&(o->e));
                stack_out.push(&(o->t));
            }
        }
    }
    node::dereference(out);
    return out;
}

inline node* copy(manager_t *manager, node *in, node::table &table){
     if(!in)
        return NULL;

    node::map in_to_out;
    node::stack stack_in;
    node::reference_stack stack_out;

    node *out = NULL;
    stack_in.push(in);
    stack_out.push(&out);
    while(stack_in.size()){
        node *n = stack_in.top();
        node *&o = *stack_out.top();

        node *&hit = in_to_out.find_or_reference(n);
        if(hit){
            o = node::reference(hit);
            stack_in.pop();
            stack_out.pop();
        } else {
            if(!o){
                o = node::reference(create(manager));
                o->l = n->l;
                if(n->W){
                    o->W = create_weights(manager);
                    *o->W = *n->W;
                }

                if(n->t){
                    stack_in.push(n->t);
                    stack_out.push(&(o->t));
                }
            } else if(n->e && !o->e){
                stack_in.push(n->e);
                stack_out.push(&(o->e));
            } else {
                merge(manager, table, o);
                hit = o;
                stack_in.pop();
                stack_out.pop();
            }
        }
    }

    node::dereference(out);
    return out;
}


void print(node *n){
    if(n){
        if(node::is_terminal(n)){
            if(node::is_satisfiable(n))
                printf("%9p: [T] ", n);
            else
                printf("%9p: [F] ", n);
        } else {
            printf("%9p: [%d] ", n, n->l);

            printf("[");
            if(n->t){
                if(node::is_terminal(n->t)){
                    if(node::is_satisfiable(n->t))
                        printf("%9p T", n->t);
                    else
                        printf("%9p F", n->t);
                } else printf("%9p %d", n->t, n->t->l);
            } else printf("NULL");
            printf("] ");

            printf("[");
            if(n->e){
                if(node::is_terminal(n->e)){
                    if(node::is_satisfiable(n->e))
                        printf("%9p T",n->e);
                    else
                        printf("%9p F",n->e);
                } else printf("%9p %d", n->e, n->e->l);
            } else printf("NULL");
            printf("] ");

            printf("[");
            if(n->W){
                for(auto it = n->W->begin(); it != n->W->end(); it++){
                    if(it != n->W->begin())
                        printf(",");
                    printf("W%d", *it);
                }
            } else printf("NULL");
            printf("] ");
        }
        printf("\n");
    }
}

template <bool COLLAPSE>
node* constraints(manager_t *manager, domain_closure_t &closure, support_t &support, ordering_t &ordering){
    return NULL;
}

template <bool COLLAPSE>
void add_clause_context(manager_t *manager, node *&n, domain_closure_t &closure, support_t &support, ordering_t &ordering){
    node *f_0 = create_terminal(manager, false); // unsatisfiable
    node *f_1 = NULL;
    node *root = NULL;

    std::stack<node*> stack;
    int level = 0;
    while(!node::is_terminal(n)){
        literal_t l = ordering[level++];
        if(l == n->l){
            closure.condition(l);
            n->ref = 0;
            stack.push(n);
            n = n->t;
        } else closure.condition(-1*l);
    }
    f_1 = n;
    f_1->ref = 0;
    root = f_1;
    level--;

    while(level >= 0){
        literal_t l = ordering[level--];

        if(!stack.empty() && (n = stack.top()) && n->l == l){
            stack.pop();
            closure.undo(l);
            n->ref = 0;
            n->t = node::reference(root);

            satisfy_t f_l = closure.condition(-1*l);
            if(f_l == unsatisfiable)
                n->e = node::reference(f_0);
            else n->e = node::reference(f_1);
            closure.undo(-1*l);
            #ifdef DEBUG
            if(f_l == redundant)
                throw compiler_debug_exception("impossible domain constraint");
            #endif
            root = n;
        } else {
            if(!closure.is_redundant(l)){
                closure.undo(-1*l);
                n = create(manager);
                n->l = l;
                n->t = node::reference(f_1);
                n->e = node::reference(root);
                root = n;
            }
        }
    }
    destroy(manager, f_0);
    n = node::reference(root);
}

template <bool COLLAPSE>
node* complete(manager_t *manager, node *n, domain_closure_t &closure, support_t &support, ordering_t &ordering){
    return NULL;
}

template <bool COLLAPSE>
void bayesnode_to_wpbdds(manager_t *manager, std::queue<node*> &bdds, BayesNode *n, domain_closure_t &closure, support_t &support, ordering_t &ordering){
    std::map<unsigned int, unsigned int> literal_to_pos;
    for(unsigned int i = 0; i < ordering.size(); i++)
        literal_to_pos[ordering[i]] = i;

    const unsigned int VARIABLES = n->GetNrVariables();

    // get all variable, literals and dimensions involved
    std::vector <literal_t> literals;
    std::vector <unsigned int> dims;
    for(unsigned int i = 0; i < VARIABLES; i++){
        literal_t l = n->GetLiteral(i);
        dims.push_back(n->GetDimension(i));
        literals.push_back(l);
    }

    // determine which literal is present in which clauses
    const auto &kWeights = n->GetWeights();
    const unsigned int M = literals.size() - 1;
    std::vector<unsigned int> ctr(literals.size(),0);
    int clause_nr = 0;
    while(true){
        unsigned int w = kWeights[clause_nr++];
        std::vector<std::pair<unsigned int, literal_t> > ordered_literals;
        for(unsigned int i = 0; i < literals.size(); i++){
            literal_t l = literals[i]+ctr[i];
            ordered_literals.push_back(std::make_pair(literal_to_pos[l],l));
        }
        std::sort(ordered_literals.begin(),ordered_literals.end());

        node *f_1 = create_terminal(manager,true), *clause = NULL, *last = NULL;
        node **parent = &clause;
        for(auto it = ordered_literals.begin(); it != ordered_literals.end(); it++){
            literal_t l = it->second;
            node *n = node::reference(create(manager));
            n->l = l;
            n->e = node::reference(f_1);
            *parent = n;
            parent = &(n->t);
            last = n;
        }
        *parent = node::reference(f_1);
        last->W = create_weights(manager);
        last->W->insert(w);

        // make complete
        add_clause_context<COLLAPSE>(manager, clause, closure, support, ordering);
        bdds.push(clause);

        bool stop = true;
        for(int i = M; i >= 0; i--){
            if(ctr[i] < dims[i]-1){
                stop = false;
                ctr[i]++;
                for(unsigned int j = i+1; j <= M; j++)
                    ctr[j] = 0;
                break;
            }
        }
        if(stop)
            break;
    }
}

template <bool COLLAPSE, bool DETERMINISM>
node* conjoin(manager_t *manager, const unsigned int partition_id, node *bdd1, node *bdd2){
    support_t support = create_support(manager,bdd1, bdd2);
    ordering_t ordering = create_ordering(manager, partition_id, support);
    node *bdd = conjoin<COLLAPSE,DETERMINISM>(manager, partition_id, bdd1, bdd2, ordering);
    if(bdd != NULL)
        add_support(manager,bdd, support);
    return bdd;
}


template <bool COLLAPSE, bool DETERMINISM>
node* conjoin(manager_t *manager, const unsigned int partition_id, node *bdd1, node *bdd2, ordering_t &ordering){
    node::table table;
    bnc::computed_table_clear(manager);
    ite_t ite(manager);

    ite.add_bdd(bdd1);
    ite.add_bdd(bdd2);

    // ASSUMPTION: at least one bdd is non-empty
    node *f_0 = create_terminal(manager, false); // unsatisfiable
    node *f_1 = create_terminal(manager, true);  // satisfiable
    table.insert(f_0);
    table.insert(f_1);

    node::reference_stack stack;
    node *root = NULL;
    stack.push(&root);
    stack.push(&root);
    int level = 0;

    #if defined(INTERMEDIATE_DOT) || defined(DEBUG)
        bool dot = false;
        bool wait = false;
        #define write_intermediate_dot if(dot) write_state(manager, root, n, bdd1, bdd2,ite, wait)
    #else
        #define write_intermediate_dot
    #endif

    try {
        traverse: {
            #ifdef DEBUG
            if((unsigned) level < 0 || (unsigned) level > ordering.size())
                throw compiler_debug_exception("Gone past level threshold. level == %d",level);
            #endif

            satisfy_t f_l;
            node *&n = *stack.top();

            // ==== compute positive cofactor ====
            write_intermediate_dot;
            if(!n){
                positive:
                write_intermediate_dot;
                const literal_t &l = ordering[level];
                f_l = ite.condition<true,COLLAPSE,DETERMINISM>(l);
                write_intermediate_dot;

                switch(f_l){
                    case satisfiable:
                        n = node::reference(create(manager));
                        n->l = l;
                        n->W = ite.get_weights();
                        n->t = node::reference(f_1);
                        ite.undo();
                        goto negative;
                    case unsatisfiable:
                        #ifdef DEBUG
                        if(!DETERMINISM)
                            throw compiler_debug_exception("Negative cofactor at true edge");
                        #endif
                        n = node::reference(create(manager));
                        n->l = l;
                        ite.clear_weight();
                        n->t = node::reference(f_0);
                        ite.undo();
                        goto negative;
                    case trivial:
                        n = node::reference(create(manager));
                        n->l = l;
                        n->W = ite.get_weights();
                        n->t = node::reference(copy(manager, ite.get_trivial(), table));
                        ite.add_cache(n->t);
                        ite.undo();
                        goto negative;
                    case cached:
                        n = node::reference(create(manager));
                        n->l = l;
                        n->t = node::reference(ite.get_cached());
                        n->W = ite.get_weights();
                        ite.undo();
                        goto negative;
                    case unsatisfied:
                        n = node::reference(create(manager));
                        n->l = l;
                        n->W = ite.get_weights();
                        stack.push(&n->t);
                        level++;
                        goto traverse;
                    case redundant:
                        level++;
                        goto positive;
                }
            }

            // === compute negative cofactor ===
            if(!n->e){
                negative:
                write_intermediate_dot;
                const literal_t &l = ordering[level];
                f_l = ite.condition<false,COLLAPSE,DETERMINISM>(l);
                write_intermediate_dot;
                switch(f_l){
                    case satisfiable:
                        n->e = node::reference(f_1);
                        ite.undo();
                        goto canonical;
                    case unsatisfiable:
                        n->e = node::reference(f_0);
                        ite.undo();
                        goto canonical;
                    case trivial:
                        n->e = node::reference(copy(manager, ite.get_trivial(), table));
                        ite.undo();
                        goto canonical;
                    case unsatisfied:
                        stack.push(&n->e);
                        level++;
                        goto traverse;
                    case redundant:
                        level++;
                        goto negative;
                    case cached:
                        throw compiler_debug_exception(
                            "Attempting to use cached entry at negative cofactor");
                        break;
                }
            }

            // ==== make canonical ====
            canonical:
            write_intermediate_dot;

            #ifdef DEBUG
            if(COLLAPSE && !DETERMINISM && node::collapse_candidate(n)){
                // collapse candidates should be impossible, as the collapse should already be
                // applied to the functions representing CPTs. We only encode structure local
                // to CPTs thus conjoining CPT functions cannot cause new collapse candidates to
                // be appear.
                //write_state(manager, root, n, bdd1, bdd2,ite, wait);
                throw compiler_debug_exception("Collapse candidate found");
            }
            #endif

            if(DETERMINISM && COLLAPSE)
                collapse(manager, table, n);

            merge(manager, table, n);
            write_intermediate_dot;

            // ==== pop stack ====
            stack.pop();
            node *&parent = *stack.top();
            while(parent->l != ordering[level]) {level--;}

            // ==== update caching ====
            if(node::is_positive_cofactor(parent,n))
                ite.add_cache(n);
            ite.undo();
            write_intermediate_dot;

            if(stack.size() > 1)
                goto traverse;
        }
    } catch(std::exception &e){
        node* itenode0;
        node* itenode1;
        node* xroot;
        unsigned int context = 4;
        const unsigned int kThreshold = 20;
        if(size(bdd1) > kThreshold)
            itenode0 = ite.get_current(0,context);
        else itenode0 = bdd1;
        if(size(bdd2) > kThreshold)
            itenode1 = ite.get_current(1,context);
        else itenode1 = bdd2;

        if(size(root) > kThreshold){
            for(unsigned int i = 0; i < context && stack.size() > 1; i++)
                stack.pop();
            xroot = *stack.top();
        } else xroot = root;
        write_state(manager, xroot, NULL, itenode0, itenode1, ite, true);
        throw;
    }

    destroy(manager,f_0);
    destroy(manager,f_1);

    #ifdef DEBUG
    if(root == NULL)
        throw compiler_debug_exception("conjoin has NULL as result");
    #endif

    return root;
}

bool equal(node *m, node *n){
    if(m && n){
        node::stack sm;
        node::stack sn;
        node::set done;

        sm.push(m);
        sn.push(n);
        done.insert(NULL);
        while(!sn.empty() && !sm.empty()){
            node *m = sm.top();
            node *n = sn.top();
            sm.pop(); sn.pop();

            bool containsn = done.contains(n);
            bool containsm = done.contains(m);
            if(!containsm && !containsn){
                done.insert(m);
                done.insert(n);

                if(!node::equal(m,n))
                    return false;

                sm.push(m->e);
                sm.push(m->t);
                sn.push(n->e);
                sn.push(n->t);
            } else if(containsm != containsn)
                return false;
        }
        if(sm.empty() != sn.empty())
            return false;
    } else if(m != n)
        return false;

    return true;
}

template <const bool COLLAPSE,const bool DETERMINISM>
node* solve(manager_t *manager, sat_t &sat, const ordering_t &ordering){
    node::table table;
    node *f_0 = create_terminal(manager, false); // contradiction
    node *f_1 = create_terminal(manager, true);  // tautology;
    table.insert(f_0);
    table.insert(f_1);

    #if defined(INTERMEDIATE_DOT) || defined(DEBUG)
        bool dot = false;
        bool wait = false;
        #define write_td_intermediate_dot if(dot) write_dot_pdf(manager, root,"td",n,wait)
    #else
        #define write_td_intermediate_dot
    #endif

    node::reference_stack s;
    node *root = NULL;
    s.push(&root);
    s.push(&root);
    unsigned int level = 0;
    while(s.size() > 1){

        #ifdef DEBUG
        if(level >= ordering.size())
            throw compiler_debug_exception("Gone past level threshold");
        #endif
        literal_t l = ordering[level];

        satisfy_t f_l;
        node *&n = *s.top();
        write_td_intermediate_dot;
        if(!n){
            f_l = sat.condition<DETERMINISM>(l);
            switch(f_l){
                case unsatisfied:
                    n = node::reference(create(manager));
                    n->l = l;
                    n->W = sat.get_weights();
                    s.push(&n->t);
                    level++;
                    break;
                case satisfiable:
                    n = node::reference(create(manager));
                    n->l = l;
                    n->W = sat.get_weights();
                    n->t = node::reference(f_1);
                    sat.undo();
                    break;
                case unsatisfiable:
                    n = node::reference(create(manager));
                    n->l = l;
                    n->W = sat.get_weights();
                    n->t = node::reference(f_0);
                    sat.undo();
                    break;
                case trivial:
                case redundant:
                    level++;
                    break;
                default:
                    break;
            }
            continue;
        } else if(!n->e){
            f_l = sat.condition<DETERMINISM>(-1*l);
            #ifdef DEBUG
            node::weights *w = sat.get_weights();
            if(w != NULL)
                throw compiler_debug_exception("Weights created on negative cofactor");
            #endif
            if(f_l == unsatisfiable){
                n->e = node::reference(f_0);
                sat.undo();
            } else { // unsatisfied
                #ifdef DEBUG
                if(f_l != unsatisfied)
                    throw compiler_debug_exception("Sat concluded impossible negative cofactor\n");
                #endif
                s.push(&n->e);
                level++;
                continue;
            }
        }

        // canonical check
        write_td_intermediate_dot;
        if(COLLAPSE)
            collapse(manager, table, n);
        merge(manager, table, n);
        write_td_intermediate_dot;

        sat.undo();
        s.pop();
        node *m = *s.top();
        while(m && m->l != ordering[level]) level--;
    }

    destroy(manager, f_0); // iff bdd does not have the f_0 node
    destroy(manager, f_1); // iff bdd does not have the f_1 node

    return root;
}

void write_bdd(manager *m, std::vector<node*> &n){
    if(n.size() == 1)
        write_bdd(m, n[0]);
    else for(unsigned int i = 0; i < n.size(); i++)
        write_bdd(m, n[i], stringf("%u",i));
}

void write_bdd(manager *m, node *n, std::string aux){
    FILE *file = fopen(files.get_filename_c(AC,aux), "w");
    if(file){
        fprintf(file, "wpbdd %u\n", size(n));
        if(n){
            // ==== build numeric index for each node ====
            // NOTE: index is assigned according to pre-order
            node *f_0 = NULL;
            node *f_1 = NULL;

            std::unordered_map<node*,unsigned int> index;
            node::stack s;
            s.push(n);
            unsigned int terminals = 0;
            unsigned int non_terminal_index = 2;
            while(!s.empty()){
                node *n = s.top();
                s.pop();

                if(n && index.find(n) == index.end()){
                    if(node::is_terminal(n)){
                        if(node::is_satisfiable(n)){
                            f_1 = n;
                            index[f_1] = 1;
                        } else {
                            f_0 = n;
                            index[f_0] = 0;
                        }

                        terminals++;
                        #ifdef DEBUG
                        if(terminals > 2)
                            throw compiler_debug_exception("indexing error (more than 2 terminals found)");
                        #endif
                    } else {
                        index[n] = non_terminal_index++;
                        s.push(n->e);
                        s.push(n->t);
                    }
                }
            }

            // ==== print WPBDD ====
            fprintf(file, "0 0 0 0\n"); // false
            fprintf(file, "0 0 0 0\n"); // true

            #ifdef DEBUG
            unsigned int check_count = 2;
            #endif
            node::set done;
            done.insert(f_0);
            done.insert(f_1);
            done.insert(NULL);

            s.push(n);
            while(!s.empty()){
                node *n = s.top();
                s.pop();

                if(!done.contains(n)){
                    done.insert(n);
                    s.push(n->e);
                    s.push(n->t);

                    #ifdef DEBUG
                    if(index[n] != check_count++)
                        throw compiler_debug_exception("inconsistency found in node indexing");
                    #endif
                    fprintf(file, "%u %u %u %u", n->l, index[n->t], index[n->e], (n->W?n->W->size():0));
                    if(n->W){
                        #ifdef ENCODE_DETERMINISM
                        if(n->W->size() == 0)
                            fprintf(file, " 0");
                        else
                        #endif
                        for(auto it = n->W->begin(); it != n->W->end(); it++)
                            fprintf(file, " %u", *it);
                    }
                    fprintf(file, "\n");
                }
            }
        }
        fclose(file);
    } else compiler_write_exception("Could not open WPBDD file '%s'", files.get_filename_c(AC,aux));
}

mpz_class get_upper_bound(bound_t<false>& bound, partition_t & partition, ordering_t &o);

mpz_class get_upper_bound(bayesnet *bn, ordering_t &o){
    bound_t<false> bound;
    bound.init(bn);
    return get_upper_bound(bound,o);
}



mpz_class get_upper_bound(bayesnet *bn, partitions_t& partitions, std::vector< ordering_t > &o){
    // we must have an equal number of partitions as orderings
    if(partitions.size() != o.size())
        return 0;

    bound_t<false> bound;
    bound.init(bn);

    // FIXME: simple bound, just sum over all partitions
    mpz_class upperbound = 0;
    for(unsigned int i = 0; i < partitions.size(); i++)
        upperbound += get_upper_bound(bound, partitions[i], o[i]);

    return upperbound;
}

mpz_class get_upper_bound(bayesnet *bn, partition_t & partition, ordering_t &o){
    bound_t<false> bound;
    bound.init(bn);

    return get_upper_bound(bound, partition, o);
}

mpz_class get_upper_bound(bound_t<false>& bound, partition_t & partition, ordering_t &o){
    if(partition.cutset.size() == 0){
        for(auto it = o.begin(); it != o.end(); it++)
            if(partition.set.find(*it) == partition.set.end())
                partition.cutset.insert(*it);
    }

    mpz_class upperbound;

    // set partition
    map_backup_t backup;
    bound.set_partition(partition, backup);

    // compute bound of partition given ordering
    // the following 2 commands can be repeated for other orderings
    bound.reset();
    upperbound = get_upper_bound(bound, o);

    // ordering 2:
    // bound.reset();
    // upperbound = get_upper_bound(bound, o2);

    // ordering 3:
    // bound.reset();
    // upperbound = get_upper_bound(bound, o3);
    // ... etc

    // restore original graph
    bound.unset_partition(partition, backup);

    return upperbound;
}

mpz_class get_upper_bound(bound_t<false>& bound, ordering_t &o){
    for(auto it = o.begin(); it != o.end(); it++){
        const unsigned int variable = *it;
        bound.condition(variable);
    }
    return bound.get_value();
}


template void bayesnode_to_wpbdds<true>(manager_t *manager, std::queue<node*> &bdds, BayesNode *n, domain_closure_t &closure, support_t &support, ordering_t &ordering);
template void bayesnode_to_wpbdds<false>(manager_t *manager, std::queue<node*> &bdds, BayesNode *n, domain_closure_t &closure, support_t &support, ordering_t &ordering);
template node* constraints<true>(manager_t*, domain_closure_t&, support_t&, ordering_t&);
template node* constraints<false>(manager_t*, domain_closure_t&, support_t&, ordering_t&);
template node* complete<true>(manager_t*, node*, domain_closure_t&, support_t&, ordering_t&);
template node* complete<false>(manager_t*, node*, domain_closure_t&, support_t&, ordering_t&);
template void add_clause_context<true>(manager_t*, node*&, domain_closure_t&, support_t&, ordering_t&);
template void add_clause_context<false>(manager_t*, node*&, domain_closure_t&, support_t&, ordering_t&);
template node* conjoin<true,true>(manager_t*, const unsigned int, node*, node*);
template node* conjoin<true,false>(manager_t*, const unsigned int, node*, node*);
template node* conjoin<false,true>(manager_t*, const unsigned int, node*, node*);
template node* conjoin<false,false>(manager_t*, const unsigned int, node*, node*);
template node* conjoin<true,true>(manager_t*, const unsigned int, node*, node*, ordering_t&);
template node* conjoin<true,false>(manager_t*, const unsigned int, node*, node*, ordering_t&);
template node* conjoin<false,true>(manager_t*, const unsigned int, node*, node*, ordering_t&);
template node* conjoin<false,false>(manager_t*, const unsigned int, node*, node*, ordering_t&);
template node* solve<true,true>(manager_t*, sat_t&, const ordering_t&);
template node* solve<true,false>(manager_t*, sat_t&, const ordering_t&);
template node* solve<false,true>(manager_t*, sat_t&, const ordering_t&);
template node* solve<false,false>(manager_t*, sat_t&, const ordering_t&);

} // namespace bnc
