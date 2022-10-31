#include "bayesgraph.h"
#include "compiler.h"
#include "types.h"
#include "options.h"
#include <unistd.h>
#include <stdio.h>
#include <bn-to-cnf/exceptions.h>
#include <bn-to-cnf/parser.h>
#include <bn-to-cnf/bayesnet.h>
#include <bn-to-cnf/cnf.h>
#include <bn-to-cnf/misc.h>
#include <thread>
#include "memory.h"
#include "threading.h"
#include "stringsep.h"

using namespace bnc;

int isnumber (const char * s){
    if (s == NULL || *s == '\0' || isspace(*s))
        return 0;
    char * p;
    strtod (s, &p);
    return *p == '\0';
}

void help(){
    fprintf(stderr, "Usage:\n   ./bnc [options] <.net file>\n\n");
    fprintf(stderr, "    Options:\n"); // pinvc:o:r:t:w:h
    fprintf(stderr, "        -c <bottomup|topdown|hybrid>      : compilation type (default: hybrid)\n");
    fprintf(stderr, "        -w <filetype>[=<filename>]        : write data to file, of filetype <part|map|dot|elim|var|order|circuit|comp|uai|pseudo|spanning>\n");
    fprintf(stderr, "        -t <wpbdd|mg|tdmg> : decision diagram type (default: tdmg)\n");
    fprintf(stderr, "        -r <part|lit|var|elim>[=<filename>] : read [partition|literal|variable|elimination] ordering\n");
    fprintf(stderr, "        -d <partitioned|monolithic>       : compile as [partitioned|monolithic] representation\n");
    fprintf(stderr, "        -o <option>=<value>               : set option\n");
    fprintf(stderr, "            option:\n");
    fprintf(stderr, "                determinism               (consider determinism, default: %s)\n",(OPT_DETERMINISM?"yes":"no"));
    fprintf(stderr, "                structure                 (encode local structure, default: %s)\n",(OPT_ENCODE_STRUCTURE?"yes":"no"));
    fprintf(stderr, "                use_probability           (use probabilities directly with mg/tdmg, default: %s)\n",(OPT_USE_PROBABILITY?"yes":"no"));
    fprintf(stderr, "                binary                    (write circuit in binary format with mg/tdmg if use_probability is true (default: %s)\n",(OPT_WRITE_BINARY?"yes":"no"));
    fprintf(stderr, "                collapse                  (enable or disable collapse rule, default: %s)\n",(OPT_COLLAPSE?"yes":"no"));
    fprintf(stderr, "                score                     (show score of ordering, default: %s)\n",(OPT_SHOW_SCORE?"yes":"no"));
    fprintf(stderr, "                reserve                   (reserve (pre-allocate) number of nodes before compilation)\n");
    fprintf(stderr, "                time                      (set compilation time limit in seconds)\n");
    fprintf(stderr, "                resources                 (set RAM usage limit by factor, example: 0.8 for 80\% RAM usage)\n");
    fprintf(stderr, "                loadfactor                (set load factor of computed table (hashmap), a value > 0.0 and <= 1.0. Default: %lf)\n",OPT_COMPUTED_TABLE_LOAD_FACTOR);
    fprintf(stderr, "                buckets                   (reserve buckets in computed table (hashmap). Default: %lu)\n", OPT_COMPUTED_TABLE_BUCKETS);
    fprintf(stderr, "                partitions                (number of partitions)\n");
    fprintf(stderr, "                parallel_conjoin          (parallel conjoin, when using parallelism)\n");
    fprintf(stderr, "                parallel_level            (level of parallelism (for multigraphs. Default: %u)\n", OPT_PARALLEL_LEVEL);
    fprintf(stderr, "                order                     (method of determining ordering, default: %u)\n", OPT_ORDER);
    fprintf(stderr, "                    0: induced by topological sort of bn\n");
    fprintf(stderr, "                    1: induced by breath first search of bn\n");
    fprintf(stderr, "                    2: induced by reverse topological sort of bn\n");
    fprintf(stderr, "                    3: branch and bound\n");
    fprintf(stderr, "                    4: brute force all combinations\n");
    fprintf(stderr, "                    5: induced by topological sort of weakly connected components\n");
    fprintf(stderr, "                    6: branch and bound with lookahead\n");
    fprintf(stderr, "                    7: simulated annealing (chain)\n");
    fprintf(stderr, "                    8: induced by all possible topological sorts of bn\n");
    fprintf(stderr, "                    9: simulated annealing (tree)\n");
    fprintf(stderr, "                    10: minhill climbing\n");
    fprintf(stderr, "                simulated annealing options:\n");
    fprintf(stderr, "                    sa_read_elim   : read elim as initialization (default: %s)\n",(OPT_SA_READ_ELIM_ORDERING?"yes":"no") );;
    fprintf(stderr, "                    sa_iterations  : nr of iterations per temperature (default: %ld)\n",OPT_SA_ITERATIONS);
    fprintf(stderr, "                    sa_tries       : nr of tries per step (default: %ld)\n",OPT_SA_TRIES);
    fprintf(stderr, "                    sa_temp_init   : initial temperature (default: %.2lf)\n",OPT_SA_TEMPERATURE_INITIAL);
    fprintf(stderr, "                    sa_score_ratio : score ratio between cutset size and balance (default: %.3lf)\n", OPT_SA_SCORE_RATIO);
    fprintf(stderr, "                    sa_temp_min    : minimal temperature (default: %.3lf)\n",OPT_SA_TEMPERATURE_MIN);
    fprintf(stderr, "                    sa_damp_factor : damping factor (default: %.3lf)\n",OPT_SA_TEMPERATURE_DAMP_FACTOR);
    fprintf(stderr, "                    sa_print       : print current states (default: %s)\n",(OPT_SA_PRINT_ORDERING?"yes":"no") );
    fprintf(stderr, "                    sa_parallel    : run %lu executions in parallel (default: %s)\n", std::thread::hardware_concurrency(), (OPT_SA_PARALLELISM?"yes":"no") );
    fprintf(stderr, "                    sa_runs        : number of sa executions, using previous solution (default: %lu)\n", OPT_SA_RUNS);
    fprintf(stderr, "        -p: parallel compilation\n");
    fprintf(stderr, "        -h: Help\n");
}

void* run(void* data){
    bayesnet *bn = (bayesnet*) data;

    try {
        compiler comp;
        comp.set_compilation_type(OPT_COMPILATION_TYPE);
        comp.set_bdd_type(OPT_BDD_TYPE);
        comp.encode(bn);

        if(OPT_WRITE_UAI){
            auto map = files.get_filename(UAI,"","map");
            auto filename =  files.get_filename(UAI);
            printf("Writting UAI format to %s with mapping %s...\n", filename.c_str(), map.c_str());
            comp.write(UAI);
            printf("Done.\n\n");
            exit(0);
        }

        comp.load();
        if(!OPT_NO_COMPILE)
            comp.compile();

        printf("\n");
        if(OPT_WRITE_DOT){
            if(comp.get_nr_partitions() > 1)
                printf("Writting DOT to %s...\n", files.get_filename_c(DOT,"*"));
            else
                printf("Writting DOT to %s...\n", files.get_filename_c(DOT));
            comp.write(DOT);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_PARTITION){
            printf("Writting partition and ordering to %s...\n", files.get_filename_c(PARTITION));
            comp.write(PARTITION);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_ORDERING){
            printf("Writting ordering to %s...\n", files.get_filename_c(ORDERING));
            comp.write(ORDERING);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_ELIM_ORDERING){
            printf("Writting elim ordering to %s...\n", files.get_filename_c(ELIM_ORDERING, "[0-9]*"));
            comp.write(ELIM_ORDERING);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_VARIABLE_ORDERING){
            printf("Writting var ordering to %s...\n", files.get_filename_c(VARIABLE_ORDERING));
            comp.write(VARIABLE_ORDERING);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_PSEUDO_ORDERING){
            printf("Writting ordering to %s...\n", files.get_filename_c(PSEUDO_ORDERING));
            comp.write(PSEUDO_ORDERING);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_SPANNING){
            printf("Writting spanning tree %s...\n", files.get_filename_c(SPANNING));
            comp.write(SPANNING);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_COMPOSITION_ORDERING){
            printf("Writting composition ordering to %s...\n", files.get_filename_c(COMPOSITION_ORDERING));
            comp.write(COMPOSITION_ORDERING);
            printf("Done.\n\n");
        }
        if(OPT_WRITE_MAPPING){
            printf("Writting mapping to %s...\n", files.get_filename_c(MAPPING));
            comp.write(MAPPING);
            printf("Done.\n\n");
        }

        if(OPT_WRITE_AC){
            std::string filename;
            filename = files.get_filename(AC, "[0-9]*");

            printf("Writting AC(s) to %s...\n", filename.c_str());

            comp.write(AC);
            printf("Done.\n\n");
        }


       // comp.write(COMPOSITION_ORDERING);
        #ifndef DEBUG
        exit(0);
        #endif
    } catch(compiler_exception &e){
        printf("COMPILER ERROR: %s\n", e.what());
    } catch(std::exception &e){
        printf("ERROR: %s\n", e.what());
    }
    return NULL;
}

int main(int argc, char **argv){
    init_options();

    int c;
    while ((c = getopt(argc, argv, "pinvc:o:r:t:w:hd:")) != -1){
        switch (c){
            case 'p': // read ordering file
                OPT_PARALLELISM = true;
                break;
            case 'i':
                OPT_ORDER_POST_WEIGHT = true;
                break;
            case 'd':
                if(strcmp(optarg,"partitioned") == 0){
                    OPT_PARTITION = true;
                } else if (strcmp(optarg,"monolithic") == 0){
                    OPT_PARTITION = false;
                } else {
                    fprintf(stderr, "Unknown option to -d: %s\n", optarg);
                    return 1;
                }
                break;
            case 'n':
                OPT_NODE_LEVEL_COMPILATION = false;
                break;
            case 'v':
                OPT_INFO = true;
                break;
            case 'b':
                OPT_COMPILATION_TYPE = compilation_t::topdown;
                break;
            case'c': // compilation type
                if(strcmp(optarg,"topdown") == 0){
                    OPT_COMPILATION_TYPE = compilation_t::topdown;
                } else if (strcmp(optarg,"bottomup") == 0){
                    OPT_COMPILATION_TYPE = compilation_t::bottomup;
                } else if (strcmp(optarg,"hybrid") == 0){
                    OPT_COMPILATION_TYPE = compilation_t::topdown_bottomup;
                    OPT_TOPDOWN_COMPILATION = false;
                } else {
                    fprintf(stderr, "Unknown option to -c: %s\n", optarg);
                    return 1;
                }
                break;
            case'o':
                {
                    string_sep_t assignment = string_sep(optarg, '=');
                    if(assignment.size() != 2){
                        fprintf(stderr, "Expecting assignment to option -o, but got: %s\n", optarg);
                        return 1;
                    }

                    if(assignment[0] == "reserve"){
                        if(isnumber(assignment[1].c_str())){
                            OPT_RESERVE = std::stoi(assignment[1]);
                            if(OPT_RESERVE < 0){
                                fprintf(stderr, "Must provide positive number to option 'reserve' \n");
                                return 1;
                            }
                        } else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "workers"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_WORKERS = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "structure"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_ENCODE_STRUCTURE = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "determinism"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_DETERMINISM = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    }else if(assignment[0] == "best_composition_ordering"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_BEST_COMPOSITION_ORDERING = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "parallel_level"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_PARALLEL_LEVEL = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "parallel_conjoin"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_PARALLEL_CONJOIN = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "parallel_cube"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_PARALLEL_CUBE = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "parallel_cpt"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_PARALLEL_CPT = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "parallel_partition"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_PARALLEL_PARTITION = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "score"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SHOW_SCORE = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "collapse"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_COLLAPSE = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "binary"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_WRITE_BINARY = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "use_probability"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_USE_PROBABILITY = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "order"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_ORDER = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "lookahead"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_LOOKAHEAD = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_runs"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_RUNS = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_parallel"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_PARALLELISM = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_print"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_PRINT_ORDERING = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "partitions"){
                        if(isnumber(assignment[1].c_str())){
                            OPT_NR_PARTITIONS = std::stoi(assignment[1]);
                            if(OPT_NR_PARTITIONS < 2){
                                fprintf(stderr, "Argument to option 'partitions' must be greater than 1\n");
                                OPT_NR_PARTITIONS = -1;
                                return 1;
                            }
                        } else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_read_elim"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_READ_ELIM_ORDERING = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }

                    } else if(assignment[0] == "sa_iterations"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_ITERATIONS = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_tries"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_TRIES = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_score_ratio"){
                        OPT_SA_SCORE_RATIO = std::stod(assignment[1]);
                        if(!(OPT_SA_SCORE_RATIO >= 0 && OPT_SA_SCORE_RATIO <= 1)){
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number in [0,1]\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_temp_min"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_TEMPERATURE_MIN = atof(assignment[1].c_str());
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_temp_init"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_TEMPERATURE_INITIAL = atof(assignment[1].c_str());
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "sa_damp_factor"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_SA_TEMPERATURE_DAMP_FACTOR = atof(assignment[1].c_str());
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "no_compile"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_NO_COMPILE = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "no_ordering"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_NO_ORDERING = (bool) std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "time"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_TIME_LIMIT = std::stoi(assignment[1]);
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "resources"){
                        double percentage = atof(assignment[1].c_str());
                        if(!(percentage > 0 && percentage < 1)){
                            fprintf(stderr, "Argument to option '%s' (%s) must be in range [0-1]\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                        if(set_memory_limit(get_ram_size(percentage)) != 0)
                            fprintf(stderr, "Unable to set memory limit\n");

                    } else if(assignment[0] == "loadfactor"){
                        float loadfactor = atof(assignment[1].c_str());
                        if(!(loadfactor > 0 && loadfactor <= 1.0)){
                            fprintf(stderr, "Argument to option '%s' (%s) must be in range [0-1]\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else if(assignment[0] == "buckets"){
                        if(isnumber(assignment[1].c_str()))
                            OPT_COMPUTED_TABLE_BUCKETS = std::atoi(assignment[1].c_str());
                        else {
                            fprintf(stderr, "Argument to option '%s' (%s) is not a number\n", assignment[0].c_str(), assignment[1].c_str());
                            return 1;
                        }
                    } else {
                        fprintf(stderr, "Unknown option to -o: %s\n", optarg);
                        return 1;
                    }
                }
                break;
            case 'r':
                {
                    string_sep_t assignment = string_sep(optarg, '=');
                    if(assignment.size() != 2 and assignment.size() != 1){
                        fprintf(stderr, "Invalid argument '%s' passed to option -r\n", optarg);
                        return 1;
                    }

                    if(assignment[0] == "order"){
                        OPT_READ_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::ORDERING,assignment[1]);
                    } else if(assignment[0] == "var"){
                        OPT_READ_VARIABLE_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::VARIABLE_ORDERING,assignment[1]);
                    } else if(assignment[0] == "pseudo"){
                        OPT_READ_PSEUDO_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::PSEUDO_ORDERING,assignment[1]);
                    } else if(assignment[0] == "elim"){
                        OPT_READ_ELIM_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::ELIM_ORDERING,assignment[1]);
                    } else if(assignment[0] == "part"){
                        OPT_READ_PARTITION = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::PARTITION,assignment[1]);
                    } else {
                        fprintf(stderr, "Unknown option to -r: %s\n", optarg);
                        return 1;
                    }
                }
                break;

            case 't':
                if(OPT_BDD_TYPE != bdd_t::none){
                    if(strcmp(optarg,"wpbdd") == 0)
                        OPT_BDD_TYPE = bdd_t::wpbdd;
                    else if(strcmp(optarg,"mg") == 0)
                        OPT_BDD_TYPE = bdd_t::multigraph;
                    else if(strcmp(optarg,"tdmg") == 0)
                        OPT_BDD_TYPE = bdd_t::tdmultigraph;
                    else {
                        fprintf(stderr, "Unknown option to -t: %s\n", optarg);
                        return 1;
                    }
                } else fprintf(stderr, "Warning: ignoring '-t' due to '-m'\n");
                break;
            case 'w':
                {
                    string_sep_t assignment = string_sep(optarg, '=');
                    if(assignment.size() != 2 and assignment.size() != 1){
                        fprintf(stderr, "Invalid argument '%s' passed to option -w\n", optarg);
                        return 1;
                    }

                    if(assignment[0] == "dot"){
                        OPT_WRITE_DOT = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::DOT,assignment[1]);
                    } else if (assignment[0] == "comp"){
                        OPT_WRITE_COMPOSITION_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::COMPOSITION_ORDERING,assignment[1]);
                    } else if (assignment[0] == "spanning"){
                        OPT_WRITE_SPANNING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::SPANNING,assignment[1]);
                    } else if (assignment[0] == "var"){
                        OPT_WRITE_VARIABLE_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::VARIABLE_ORDERING,assignment[1]);
                    } else if (assignment[0] == "pseudo"){
                        OPT_WRITE_PSEUDO_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::PSEUDO_ORDERING,assignment[1]);
                    } else if (assignment[0] == "lit"){
                        OPT_WRITE_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::ORDERING,assignment[1]);
                    } else if (assignment[0] == "elim"){
                        OPT_WRITE_ELIM_ORDERING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::ELIM_ORDERING,assignment[1]);
                    } else if (assignment[0] == "part"){
                        OPT_WRITE_PARTITION = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::PARTITION,assignment[1]);
                    } else if (assignment[0] == "uai"){
                        OPT_WRITE_UAI = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::UAI,assignment[1]);
                    } else if (assignment[0] == "map"){
                        OPT_WRITE_MAPPING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::MAPPING,assignment[1]);
                    } else if (assignment[0] == "circuit"){
                        OPT_WRITE_AC = true;
                        OPT_WRITE_MAPPING = true;
                        if(assignment.size() == 2)
                            files.set_filename(file_t::AC,assignment[1]);
                    } else {
                        fprintf(stderr, "Unknown option to -w: %s\n", optarg);
                        return 1;
                    }
                }
                break;
            case 'h':
                init_options();
                help();
                return 0;
            case '?':
            default:
                help();
                return 1;
        }
    }

    // check if bayesian network is provided
    bool have_bayesian_network = false;
    for (int index = optind; index < argc; index++){
        if(files.get_extension(argv[index]) == "net"){
            if(have_bayesian_network){
                fprintf(stderr, "Only one Hugin input file must be provided.\n");
                return 1;
            }

            have_bayesian_network = true;
            std::string filename = argv[index];
            files.set_filename(file_t::BN,filename);
            if(files.is_set(file_t::AC))
                files.set_basename(files.get_filename(file_t::AC));
            else files.set_basename(filename);
         } else fprintf(stderr, "Warning: option '%s' ignored.\n",argv[index]);
    }
    if(!have_bayesian_network){
        fprintf(stderr, "No Hugin input file provided.\n");
        return 1;
    }

    if(OPT_DETERMINISM && !(
                OPT_BDD_TYPE == bdd_t::wpbdd ||
                OPT_BDD_TYPE == bdd_t::multigraph ||
                OPT_BDD_TYPE == bdd_t::tdmultigraph)){
        fprintf(stdout, "Warning: determinism not supported for output type. Disabling it.\n");
        OPT_DETERMINISM = false;
    }

    if(OPT_PARALLELISM && OPT_ENCODE_STRUCTURE){
        fprintf(stdout, "Warning: local structure not yet supported with parallelism. Disabling it.\n");
        OPT_ENCODE_STRUCTURE = false;
    }

    if((OPT_READ_PARTITION && OPT_WRITE_PARTITION)
        || (OPT_READ_ORDERING && OPT_WRITE_ORDERING)
        || (OPT_READ_ELIM_ORDERING && OPT_WRITE_ELIM_ORDERING)){
        fprintf(stderr, "Option to '-r' also provided to '-w'\n");
        return 1;
    }

    if(OPT_READ_PARTITION
        || OPT_WRITE_PARTITION)
        OPT_PARTITION = true;

    if(OPT_PARTITION){
        if(!(OPT_BDD_TYPE == bdd_t::wpbdd || OPT_BDD_TYPE == bdd_t::multigraph || OPT_BDD_TYPE == bdd_t::tdmultigraph)){
            fprintf(stderr, "Partitioning options are only available for WPBDDs\n");
            return 1;
        }

        if(OPT_WRITE_AC || OPT_WRITE_ORDERING || OPT_WRITE_VARIABLE_ORDERING)
            OPT_WRITE_PARTITION = true;

    }

    bayesnet *bn = NULL;
    try {
        bn = bayesnet::read(files.get_filename_c(BN));
    } catch (bayesnet_exception &e){
        fprintf(stderr, "PARSING ERROR: %s\n", e.what());
        return 1;
    } catch (...){
        fprintf(stderr, "Could not load Bayesian Network\n");
        return 1;
    }

    if(!bn){
        fprintf(stderr, "Could not load Bayesian Network\n");
        return 1;
    }

    if(OPT_TIME_LIMIT > 0){
        timed_thread_t thread;
        thread.run(run,(void*)bn,OPT_TIME_LIMIT);
        if(thread.is_aborted())
            fprintf(stderr, "Time limit exceeded (%d seconds)\n", OPT_TIME_LIMIT);

    } else run((void*)bn);

    delete bn;

    return 0;
}

