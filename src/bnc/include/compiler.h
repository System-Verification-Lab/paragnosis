#ifndef COMPILER_H
#define COMPILER_H

#include "exceptions.h"
#include <bn-to-cnf/bayesnet.h>
#include "files.h"
#include "timer.h"
#include "types.h"
#include "ordering.h"
#include "sat.h"
#include "ite.h"
#include "queue.h"
#include "bnc.h"
#include "multigraph.h"
#include "multigraphp.h"

struct info_t {
    Timer cpt_timer;
    Timer join_timer;
    Timer total_timer;
    unsigned int CUMULATIVE_CPT_SIZE;
    unsigned int CUMULATIVE_OPT_SIZE;
};

typedef void* (*pthread_function_t)(void*);

class compiler {
    public:
        compiler();
        ~compiler();
        void encode(bayesnet*);
        void load();
        void set_compilation_type(compilation_t);
        void set_bdd_type(bdd_t);

        void print();
        void prepare();
        void compile();

        size_t get_nr_partitions();
        void write(file_t);
        void dot_to_ps();

        void print_result();

        void clear_stats();
        bnc::node* get_wpbdd();
        pthread_function_t get_compile(bool,bool,bool,compilation_t);
    private:
        bdd_t BDD_TYPE;
        compilation_t COMPILATION_TYPE;

        int compile_multigraph();
        std::vector<bnc::MultiGraph> multigraphs;
        int compile_pmultigraph();
        std::vector<bnc::MultiGraphProbability> pmultigraphs;


        int compile_wpbdd();
        std::vector<bnc::node*> wpbdd;

        void info();
        template <bool COLLAPSE, bool DETERMINISM, compilation_t COMPILATION_TYPE> static void* compile_sync(void*);
        template <bool COLLAPSE, bool DETERMINISM, compilation_t COMPILATION_TYPE> static void* compile_async(void*);

        bnc::manager_t manager;
};

#endif

