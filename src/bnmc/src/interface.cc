#include <limits>
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <bn-to-cnf/misc.h>
#include <bnc/options.h>
#include <string.h>
#include <unistd.h>
#include <unordered_map>
#include <bnc/timer.h>
#include <algorithm>
#include <bn-to-cnf/exceptions.h>
#include <bnc/exceptions.h>
#include <cmath>
#include "exceptions.h"
#include "interface.h"
#include "options.h"
#include "io.h"
#include "trim.h"
#include <thread>
#include "ace.h"
#include <iostream>
#include <sstream>
#include <bits/stdc++.h>

namespace bnmc {

using Fp_info = std::numeric_limits<double>;

inline auto is_ieee754_nan( double const x ) -> bool {
    static constexpr bool   is_claimed_ieee754  = Fp_info::is_iec559;
    static constexpr int    n_bits_per_byte     = CHAR_BIT;
    using Byte = unsigned char;

    static_assert( is_claimed_ieee754, "!" );
    static_assert( n_bits_per_byte == 8, "!" );
    static_assert( sizeof( x ) == sizeof( uint64_t ), "!" );

    #ifdef _MSC_VER
    uint64_t const bits = reinterpret_cast<uint64_t const&>( x );
    #else
    Byte bytes[sizeof(x)];
    memcpy( bytes, &x, sizeof( x ) );
    uint64_t int_value;
    memcpy( &int_value, bytes, sizeof( x ) );
    uint64_t const& bits = int_value;
    #endif

    static constexpr uint64_t   sign_mask       = 0x8000000000000000;
    static constexpr uint64_t   exp_mask        = 0x7FF0000000000000;
    static constexpr uint64_t   mantissa_mask   = 0x000FFFFFFFFFFFFF;

    (void) sign_mask;
    return (bits & exp_mask) == exp_mask and (bits & mantissa_mask) != 0;
}



bool is_nan(probability_t &w){
    return is_ieee754_nan(w) || std::isnan(w) || (w == std::numeric_limits<probability_t>::quiet_NaN()) || (w == std::numeric_limits<probability_t>::infinity()) || (w != w);
}

Interface::Interface(){
    Init();
}

void Interface::Init(){
    try {
        AddCommands();
    } catch (CliException &exception) {
        Print(ERR, "[CLI EXCEPTION] %s\n", exception.what());
    } catch (std::exception &exception) {
        Print(ERR, "[EXCEPTION] %s\n", exception.what());
    }
}

void Interface::Prompt(){
    Cli::Greet();
    while(1) {
        try {
            Cli::Prompt();
            break;
        } catch (CliException &exception) {
            Print(ERR, "[CLI EXCEPTION] %s\n", exception.what());
        } catch (InterfaceException &exception) {
            Print(ERR, "[INTERFACE EXCEPTION] %s\n", exception.what());
        } catch (BnmcException &exception) {
            Print(ERR, "[BNMC EXCEPTION] %s\n", exception.what());
        } catch (std::exception &exception) {
            Print(ERR, "[EXCEPTION] %s\n", exception.what());
        }
    }
}

void Interface::AddCommands(){
    // register functions
    AddCommand("help",          (CommandPtr) &Interface::Help);
    AddCommand("assignments",   (CommandPtr) &Interface::Assignments);
    AddCommand("load",          (CommandPtr) &Interface::Load);
    AddCommand("log",           (CommandPtr) &Interface::Logfile);
    AddCommand("query",         (CommandPtr) &Interface::Query);
    AddCommand("evidence",      (CommandPtr) &Interface::Evid);
    AddCommand("posteriors",    (CommandPtr) &Interface::Posteriors);
    AddCommand("compare",       (CommandPtr) &Interface::Compare);
    #ifdef ACE
    AddCommand("ace",           (CommandPtr) &Interface::InitAce);
    #endif
    AddCommand("verify",        (CommandPtr) &Interface::Verify);

}

void Interface::Read(file_t type){
    switch(type){
        case PWPBDD:
            if(!manager.have_pwpbdd){
                pwpbdd_.Read(); // reads partition file and ACs
                manager.have_pwpbdd = true;
            }
            break;
        case WPBDD:
            if(!manager.have_wpbdd){
                wpbdd_.Read();
                manager.have_wpbdd = true;
            }
            break;
        case MULTIGRAPH:
            if(!manager.have_multigraph){
                multigraph_.Read();
                manager.have_multigraph = true;
            }
            break;
        case TDMULTIGRAPH:
            if(!manager.have_tdmultigraph){
                tdmultigraph_.Read();
                manager.have_tdmultigraph = true;
            }
            break;
        case MAPPING:
            if(!manager.have_mapping){
                manager.Read(MAPPING);
                manager.have_mapping = true;
            }
            break;
        case BN:
            if(!manager.have_bn){
                manager.Read(BN);
                manager.have_bn = true;
            }
            break;
        default:
            throw IoException("Interface: Unknown read option");
    }
}


// ============================= CLI commands ===============================

int Interface::Help(void*){
    if(arguments_.length() > 0){
        Print(ERR, "Unexpected argument: %s\n", arguments_.c_str());
        return 1;
    }

    Print(MSG, "Commands:\n");
    Print(MSG, "    help               : show commands\n");
    Print(MSG, "    quit               : terminate\n");
    Print(MSG, "    load <TYPE> <FILE> : load wpbdd (*.wpbdd), mg (*.mc), (*.tdmc), tdmg (*.tdmg), Bayesian network (*.net (HUGIN), etc\n");
    Print(MSG, "    verify [N]         : exhaustively verify answer of every possible query (of 1 to N instantiations)\n");
    Print(MSG, "    compare [N]        : exhaustively measure time of every possible query (of 1 to N instantiations)\n");
    #ifdef ACE
    Print(MSG, "    ace                : add Ace to comparisons\n");
    #endif

    Print(MSG, "    assignments        : list possible assignments\n");
    Print(MSG, "    query <query type> <query> : compute probability of <query> using <query type>\n");
    Print(MSG, "\n");

    Print(MSG, "Query syntax:\n");
    Print(MSG, "    <query type> : <pwpbdd | wpbdd>\n");
    Print(MSG, "    <query>      : <assignment> [| <assignment> [, <assignment> [, ..]]], or\n");
    Print(MSG, "                 : <assignment> [, <assignment> [, ..]]]\n");
    Print(MSG, "    <assignment> : <variable> = <value>\n");

    return 0;
}

int Interface::Assignments(void*){
    if(!manager.have_mapping){
        Print(MSG, "Assignment list is empty\n");
        return 0;
    }

    if(arguments_.length() > 0){
        Print(ERR, "Unexpected argument: %s\n", arguments_.c_str());
        return 1;
    }

    const variable_literal_map_t &variables = manager.mapping.get_variable_value_name_to_literal();
    unsigned SIZE = 0;
    for(auto it = variables.begin(); it != variables.end(); it++){
        unsigned int length = it->first.length();
        SIZE = (SIZE < length?length:SIZE);
    }

    for(auto it = variables.begin(); it != variables.end(); it++){
        auto &values = it->second;
        std::string value_names;
        value_names += "<";
        bool comma = false;
        for(auto vit = values.begin(); vit != values.end(); vit++){
            if(comma) value_names += " | ";
            value_names += vit->first;
            comma = true;
        }
        value_names += ">";
        Print(MSG, "%.*s = %s\n", SIZE, it->first.c_str(), value_names.c_str());
    }

    return 0;
}

int Interface::Logfile(void*v){
    if(arguments_.empty()){
        Print(ERR, "Provide a filename\n");
        return 1;
    }

    manager.logfile = arguments_;
    Trim(manager.logfile);
    Print(MSG, "Log file set to: %s\n", manager.logfile.c_str());
}

int Interface::Load(void*v){
    if(arguments_.empty()){
        Print(ERR, "Filetype not provided\n");
        return 1;
    }

    // parse first argument
    size_t space_pos = arguments_.find_first_of(' ');
    std::string subcommand = arguments_.substr(0,space_pos);
    if(!(subcommand == "pwpbdd"
            || subcommand == "mg"
            || subcommand == "tdmg"
            || subcommand == "wpbdd"
            || subcommand == "net"
            || subcommand == "map")){
        Print(ERR, "First argument must be <pwpbdd|wpbdd|mg|tdmg|net|map>\n");
        return 1;
    }

    // parse remaining arguments
    std::vector<std::string> arguments;
    while(space_pos != std::string::npos){
        arguments_ = arguments_.substr(space_pos+1);
        space_pos = arguments_.find_first_of(' ');
        if(space_pos == 0)
            continue;

        if(space_pos != std::string::npos)
            arguments.push_back(arguments_.substr(0,space_pos));
        else arguments.push_back(arguments_);
        Trim(arguments.back());
    }

    // check number of provided arguments
    if(subcommand == "pwpbdd"){
        if(arguments.size() != 2 && arguments.size() != 3){
            Print(ERR, "Usage: load pwpbdd <circuit> <partition file> [comp order file].\n");
            return 1;
        }
    } else if(arguments.size() != 1){
        Print(ERR, "Expecting 1 argument (a filename), but got %u arguments.\n",arguments.size()+1);
        return 1;
    }

    // check if files exist
    for(unsigned int i = 0; i < arguments.size(); i++){
        if(!manager.files.exists(arguments[i])){
            Print(ERR,"File '%s' does not exist.\n",arguments[i].c_str());
            return 1;
        }
    }

    // require BN first
    if(!manager.have_bn){
        if(subcommand != "net"){
            Print(ERR, "Must load Bayesian network first\n");
            return 1;
        }
    } else if(!manager.have_mapping){
        if(subcommand != "map"){
            Print(ERR, "Must load mapping first\n");
            return 1;
        }
    }

    // already loaded the input?
    if((subcommand == "pwpbdd" && manager.have_pwpbdd)
            || (subcommand == "wpbdd" && manager.have_wpbdd)
            || (subcommand == "mg" && manager.have_multigraph)
            || (subcommand == "tdmg" && manager.have_tdmultigraph)
            || (subcommand == "net" && manager.have_bn)
            || (subcommand == "map" && manager.have_mapping)){
        Print(ERR, "Already loaded %s\n",subcommand.c_str());
        return 1;
    }

    // load based on first argument, store arguments
    try {
        if(subcommand == "pwpbdd"){
            manager.files.set_filename(file_t::PARTITION,arguments[1]);
            if(arguments.size() == 3)
                manager.files.set_filename(file_t::COMPOSITION_ORDERING,arguments[2]);

            std::string filename = manager.files.remove_extension(
                    manager.files.remove_extension(arguments[0])) +
                    "." + manager.files.get_extension(arguments[0]);

            manager.files.set_filename(file_t::PWPBDD,filename);
            Read(file_t::PWPBDD);
        } else if(subcommand == "wpbdd"){
            manager.files.set_filename(file_t::WPBDD,arguments[0]);
            Read(file_t::WPBDD);
        } else if(subcommand == "mg"){
            manager.files.set_filename(file_t::MULTIGRAPH,arguments[0]);
            Read(file_t::MULTIGRAPH);
        } else if(subcommand == "tdmg"){
            manager.files.set_filename(file_t::TDMULTIGRAPH,arguments[0]);
            Read(file_t::TDMULTIGRAPH);
        } else if(subcommand == "net"){
            manager.files.set_filename(file_t::BN,arguments[0]);
            manager.files.set_basename(arguments[0]);
            manager.have_filename = true;
            Read(file_t::BN);
        } else if(subcommand == "map"){
            manager.files.set_filename(file_t::MAPPING,arguments[0]);
            Read(file_t::MAPPING);
        } else {
            Print(ERR, "Argument '%s' not implemented.\n");
            return 1;
        }
    } catch (IoException &e){
        Print(ERR, "%s\n", e.what());
        return 1;
    } catch (compiler_exception &e){
        Print(ERR, "%s\n", e.what());
        return 1;
    } catch (std::exception &e){
        Print(ERR, "%s\n", e.what());
        return 1;
    } catch (...){
        Print(ERR, "could not load %s\n", arguments_.c_str());
        return 1;
    }
    return 0;
}

int Interface::Evid(void*){
    try {
        manager.evidence.Parse(arguments_);
    } catch (EvidenceException &exception){
        Print(ERR, "%s\n", exception.what());
        return 1;
    }

    if (manager.evidence.HaveQueryVariable() && manager.evidence.GetEvidenceVariableSet().size() != 1){
        Print(ERR, "%s\n", "Evidence is not allowed to be conditional");
        return 1;
    }

    return 0;
}

int Interface::Posteriors(void*){
    size_t argument_length = arguments_.length();
    if(argument_length == 0){
        Print(ERR, "Model not provided\n");
        return 1;
    }

    size_t pos = arguments_.find_first_of(' ');
    std::string subcommand = arguments_.substr(0,pos);
    if(pos != std::string::npos){
        arguments_ = arguments_.substr(pos+1);
        Trim(arguments_);
    } else {
        subcommand = arguments_;
        arguments_ = "";
    }

    if(subcommand == "wpbdd"){
        if(!manager.have_wpbdd){
            Print(ERR, "Must read WPBDD wth 'load' prior to use\n");
            return 1;
        }
        return Posteriors(ModelType::WPBDD);
    } else if(subcommand == "pwpbdd"){
        if(!manager.have_pwpbdd){
            Print(ERR, "Must read PWPBDD with 'load' prior to use\n");
            return 1;
        }
        return Posteriors(ModelType::PWPBDD);
    } else if(subcommand == "mg"){
        if(!manager.have_multigraph){
            Print(ERR, "Must read MULTIGRAPH with 'load' prior to use\n");
            return 1;
        }
        return Posteriors(ModelType::MULTIGRAPH);
    } else if(subcommand == "tdmg"){
        if(!manager.have_tdmultigraph){
            Print(ERR, "Must read Tree-driven WPBDD with 'load' prior to use\n");
            return 1;
        }
        return Posteriors(ModelType::TDMULTIGRAPH);
    } else {
        Print(ERR, "Unknown query type '%s' (supported types: [tdmg|mg|pwpbdd|wpbdd])\n", subcommand.c_str());
        return 1;

    }
    return 1;
}

int Interface::Query(void*){
    size_t argument_length = arguments_.length();
    if(argument_length == 0){
        Print(ERR, "Query not provided\n");
        return 1;
    }

    size_t pos = arguments_.find_first_of(' ');
    std::string subcommand = arguments_.substr(0,pos);
    if(pos != std::string::npos){
        arguments_ = arguments_.substr(pos+1);
        Trim(arguments_);
    } else {
        subcommand = arguments_;
        arguments_ = "";
    }

    if(subcommand == "wpbdd"){
        if(!manager.have_wpbdd){
            Print(ERR, "Must read WPBDD wth 'load' prior to use\n");
            return 1;
        }
        return Query(ModelType::WPBDD);
    } else if(subcommand == "pwpbdd"){
        if(!manager.have_pwpbdd){
            Print(ERR, "Must read PWPBDD with 'load' prior to use\n");
            return 1;
        }
        return Query(ModelType::PWPBDD);
    } else if(subcommand == "mg"){
        if(!manager.have_multigraph){
            Print(ERR, "Must read MULTIGRAPH with 'load' prior to use\n");
            return 1;
        }
        return Query(ModelType::MULTIGRAPH);
    } else if(subcommand == "tdmg"){
        if(!manager.have_tdmultigraph){
            Print(ERR, "Must read Tree-driven WPBDD with 'load' prior to use\n");
            return 1;
        }
        return Query(ModelType::TDMULTIGRAPH);

    } else if(subcommand == "ace"){
        if(!manager.have_ace){
            Print(ERR, "Must run 'ace' first.\n");
            return 1;
        }
        return Query(ModelType::UCLA_ACE);

    } else {
        Print(ERR, "Unknown query type '%s' (supported types: [tdmg|mg|pwpbdd|wpbdd])\n", subcommand.c_str());
        return 1;

    }
    return 1;
}


probability_t Interface::Query(const ModelType kModelType, Evidence &evidence, Timer &t){
    probability_t w;

    switch (kModelType) {
        case ModelType::MULTIGRAPH:
            multigraph_.SetEvidence(evidence);
            t.Start();
            w = multigraph_.Posterior();
            t.Stop();
            t.Add();
            break;

        case ModelType::TDMULTIGRAPH:
            tdmultigraph_.SetEvidence(evidence);
            t.Start();
            w = tdmultigraph_.Posterior();
            t.Stop();
            t.Add();
            break;

        case ModelType::WPBDD:
            wpbdd_.SetEvidence(evidence);
            t.Start();
            w = wpbdd_.Posterior();
            t.Stop();
            t.Add();
            break;

        case ModelType::PWPBDD:
            pwpbdd_.SetEvidence(evidence);
            pwpbdd_.InitCache();
            t.Start();
            w = pwpbdd_.Posterior();
            t.Stop();
            t.Add();
            break;

        default:
            throw InterfaceException("This type of query is not possible with given model type");
            break;
    }

    return w;
}

inline const unsigned int GetVariableDimension(const Variable variable){
    return manager.mapping.get_dimension()[variable];
}


int Interface::Posteriors(const ModelType kModelType){
    try {
        manager.evidence.ParsePosteriors(arguments_);
    } catch (EvidenceException &exception){
        Print(ERR, "%s\n", exception.what());
        return 1;
    }

    auto variables = manager.evidence.GetPosteriorVariableSet();
    if (variables.size() == 0)
        manager.evidence.InducePosteriors();

    // open log file
    std::string filename = manager.logfile;
    FILE *file = !filename.empty() ? fopen(filename.c_str(), "w") : NULL;

    Timer t;
    try {
        Evidence evidence;
        evidence.Add(manager.evidence);

        // compute P(Y)
        probability_t joint = Query(kModelType, evidence, t);

        auto vars = manager.evidence.GetPosteriorVariableSet();
        for(auto posterior_it = vars.begin(); posterior_it != vars.end(); posterior_it++){
            auto variable = *posterior_it;
            auto cardinality = GetVariableDimension(variable);

            // compute P(X | Y) for each instantiation of X
            std::vector<probability_t> P;
            for (auto value = 0; value < cardinality -1; value++) {
                EvidenceVariable v = EvidenceVariable(variable, value);
                evidence.Add(v);

                // compute P(X,Y)
                probability_t p = Query(kModelType, evidence, t);
                evidence.Remove(v);

                // compute P(X | Y) = P(X,Y) / P(Y)
                P.push_back(p/joint);
            }

            // P(X | Y) of the final instantiation of X can be computed without inference, i.e. 1 - Sum( other P(X | Y)'s)
            probability_t sum = 0;
            for (auto p = P.begin(); p != P.end(); p++)
                sum += *p;
            P.push_back(1 - sum);

            // print results
            for (auto value = 0; value < cardinality; value++) {
                std::string assignment = evidence.GetAssignmentString(variable, value);

                Print(MSG, "  %s: %lf\n", assignment.c_str(), P[value]);
                if (file)
                    fprintf(file, "%lf,%s\n", P[value], assignment.c_str());
            }
        }


        Print(MSG,"\n  total time: %.3fms (%.3fus)\n", t.GetTotal<Timer::Milliseconds>(),t.GetTotal<Timer::Microseconds>());


    } catch (ModelCounterException &exception){
        Print(ERR, "%s\n", exception.what());
        return 1;
    }

    // close log file
    if (file) {
        fclose(file);
        Print(MSG,"\n  results are written to: %s\n", filename.c_str());
    }

    return 0;
}

int Interface::Query(const ModelType kModelType){
    try {
        manager.evidence.Parse(arguments_);
    } catch (EvidenceException &exception){
        Print(ERR, "%s\n", exception.what());
        return 1;
    }

    probability_t w;
    Timer t;
    try {
        switch (kModelType) {

            case ModelType::UCLA_ACE:
                #ifdef ACE
                  w = ace_posterior(manager.evidence,t);
                #else
                   throw InterfaceException("ACE was not linked at compilation.");
                #endif
                break;

            case ModelType::MULTIGRAPH:
                multigraph_.SetEvidence(manager.evidence);
                t.Start();
                w = multigraph_.Posterior();
                t.Stop();
                t.Add();
                break;

            case ModelType::TDMULTIGRAPH:
                tdmultigraph_.SetEvidence(manager.evidence);
                t.Start();
                w = tdmultigraph_.Posterior();
                t.Stop();
                t.Add();
                break;

            case ModelType::WPBDD:
                wpbdd_.SetEvidence(manager.evidence);
                t.Start();
                w = wpbdd_.Posterior();
                t.Stop();
                t.Add();
                break;

            case ModelType::PWPBDD:
                pwpbdd_.SetEvidence(manager.evidence);
                pwpbdd_.InitCache();
                t.Start();
                w = pwpbdd_.Posterior();
                t.Stop();
                t.Add();
                break;

            default:
                throw InterfaceException("Undefined query model");
                break;
        }
        Print(MSG,"%lf\n", w);
        Print(MSG,"time: %.3fms (%.3fus)\n", t.GetTotal<Timer::Milliseconds>(),t.GetTotal<Timer::Microseconds>());

        return 0;
    } catch (ModelCounterException &exception){
        Print(ERR, "%s\n", exception.what());
        return 1;
    }

    return 1;
}

void Interface::Exhaustive(const ExhaustiveType kExhaustiveType, const int initial_instantiations, const int kMaxIterations, const int kMaxLocalIterations){

    const unsigned int VARIABLES = manager.bn->get_nr_variables();
    const std::vector<unsigned int> dimension = std::vector<unsigned int>(manager.bn->get_states(), manager.bn->get_states() + VARIABLES);

    Evidence &evidence = manager.evidence;

    std::unordered_map< std::string, probability_t > probabilities;
    std::unordered_map< std::string, Timer > timers;

    std::vector<unsigned int> vars;
    vars.resize(VARIABLES);
    manager.total_query_count = 0;

    if(manager.have_pwpbdd){
        printf("\n");
        printf("Architecture:\n");
        printf("=============\n");
        const Architecture &kArchitecture = pwpbdd_.GetArchitecture();
        const Composition &kComp = kArchitecture.GetComposition();
        kComp.PrintAscii(manager.bn);
    }



    const int next_instance = kMaxIterations >= 0 ? 1 : -1;
    const bool increasing = next_instance >= 0;
    unsigned int instantiations = initial_instantiations >= 0 ? initial_instantiations : (VARIABLES+1) + initial_instantiations;
    int iterations = 0;

    for(; instantiations <= VARIABLES && instantiations > 0 && iterations != kMaxIterations; instantiations += next_instance, iterations += next_instance){
        std::fill(vars.begin(), vars.end(), 0);

        Print(MSG, "\nInference with %u/%u instantiations                                            \n", instantiations, VARIABLES);
        unsigned int count = 0;

        if(!increasing) {
            for(unsigned int v = 0; v < instantiations; v++)
                vars[v] = 1;
        } else {
            for(unsigned int v = instantiations; v-- > 0;)
                vars[v] = 1;
        }

        unsigned int localIterations = 0;
        while(kMaxLocalIterations < 0 || localIterations <= kMaxLocalIterations){

            std::vector<unsigned int> ctr(VARIABLES,0);
            unsigned int M = VARIABLES-1;
            while(M > 0 && !vars[M]) M--;
            while(kMaxLocalIterations < 0 || localIterations <= kMaxLocalIterations){

                // ================== inference ====================
                // add evidence
                evidence.Clear();
                for(unsigned int v = 0; v < VARIABLES; v++)
                    if(vars[v])
                        evidence.Add(v,ctr[v]);

                for(unsigned int ql = 1; ql <= instantiations; ql++){
                    // determine query literal
                    unsigned int counter = 0;
                    for(unsigned int v = 0; v < VARIABLES; v++){
                        if(vars[v]){
                            counter++;
                            if(counter == ql){
                                evidence.SetQueryVariable(v);
                                break;
                            }
                        }
                    }

                    localIterations++;
                    if (kMaxLocalIterations > 0 && localIterations > kMaxLocalIterations)
                        break;

                    manager.total_query_count++;
                    count++;

                    //std::string query = evidence.GetQueryString();
                    //Print(ERR, "Query    : %s\n",query.c_str());

                    // compute probability

                    if(manager.have_tdmultigraph){
                        try {
                            // execute 2x to eliminate cache advantage
                            Timer &time = timers["TDMULTIGRAPH"];
                            probability_t &probability = probabilities["TDMULTIGRAPH"];
                            tdmultigraph_.SetEvidence(evidence);
                            probability = tdmultigraph_.Posterior();
                            time.Start();
                            probability = tdmultigraph_.Posterior();
                            time.Stop();
                            time.Add();
                        } catch (ModelCounterException &exception){
                            Print(ERR, "                                                                      \n");
                            Print(ERR, "TDMULTIGRAPH: %s\n", exception.what());
                            std::string query = evidence.GetQueryString();
                            Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            Print(ERR, "    Query    : %s\n",query.c_str());

                        }
                    }


                    if(manager.have_multigraph){
                        try {
                            // execute 2x to eliminate cache advantage
                            Timer &time = timers["MULTIGRAPH"];
                            probability_t &probability = probabilities["MULTIGRAPH"];
                            multigraph_.SetEvidence(evidence);
                            probability = multigraph_.Posterior();
                            time.Start();
                            probability = multigraph_.Posterior();
                            time.Stop();
                            time.Add();
                        } catch (ModelCounterException &exception){
                            Print(ERR, "                                                                      \n");
                            Print(ERR, "MULTIGRAPH: %s\n", exception.what());
                            std::string query = evidence.GetQueryString();
                            Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            Print(ERR, "    Query    : %s\n",query.c_str());

                        }
                        //try {
                        //    // execute 2x to eliminate cache advantage
                        //    Timer &time = timers["MULTIGRAPH2"];
                        //    probability_t &probability = probabilities["MULTIGRAPH2"];
                        //    multigraph_.SetEvidence(evidence);
                        //    probability = multigraph_.Posterior<2>();
                        //    time.Start();
                        //    probability = multigraph_.Posterior<2>();
                        //    time.Stop();
                        //    time.Add();
                        //} catch (ModelCounterException &exception){
                        //    Print(ERR, "                                                                      \n");
                        //    Print(ERR, "MULTIGRAPH2: %s\n", exception.what());
                        //    std::string query = evidence.GetQueryString();
                        //    Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                        //    Print(ERR, "    Query    : %s\n",query.c_str());
                        //}
                    }

                    if(manager.have_wpbdd){
                        try {
                            // execute 2x to eliminate cache advantage
                            Timer &time = timers["WPBDD"];
                            probability_t &probability = probabilities["WPBDD"];
                            wpbdd_.SetEvidence(evidence);
                            probability = wpbdd_.Posterior();
                            time.Start();
                            probability = wpbdd_.Posterior();
                            time.Stop();
                            time.Add();
                        } catch (ModelCounterException &exception){
                            Print(ERR, "                                                                      \n");
                            Print(ERR, "WPBDD: %s\n", exception.what());
                            std::string query = evidence.GetQueryString();
                            Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            Print(ERR, "    Query    : %s\n",query.c_str());
                        }
                        //try {
                        //    // execute 2x to eliminate cache advantage
                        //    Timer &time = timers["WPBDD2"];
                        //    probability_t &probability = probabilities["WPBDD2"];
                        //    wpbdd_.SetEvidence(evidence);
                        //    probability = wpbdd_.Posterior<2>();
                        //    time.Start();
                        //    probability = wpbdd_.Posterior<2>();
                        //    time.Stop();
                        //    time.Add();
                        //} catch (ModelCounterException &exception){
                        //    Print(ERR, "                                                                      \n");
                        //    Print(ERR, "WPBDD: %s\n", exception.what());
                        //    std::string query = evidence.GetQueryString();
                        //    Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                        //    Print(ERR, "    Query    : %s\n",query.c_str());
                        //}
                        //try {
                        //    // execute 2x to eliminate cache advantage
                        //    Timer &time = timers["WPBDD3"];
                        //    probability_t &probability = probabilities["WPBDD3"];
                        //    wpbdd_.SetEvidence(evidence);
                        //    probability = wpbdd_.Posterior<3>();
                        //    time.Start();
                        //    probability = wpbdd_.Posterior<3>();
                        //    time.Stop();
                        //    time.Add();
                        //} catch (ModelCounterException &exception){
                        //    Print(ERR, "                                                                      \n");
                        //    Print(ERR, "WPBDD: %s\n", exception.what());
                        //    std::string query = evidence.GetQueryString();
                        //    Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                        //    Print(ERR, "    Query    : %s\n",query.c_str());
                        //}
                    }

                    if(manager.have_pwpbdd){
                        try {
                            // execute 2x to eliminate cache advantage
                            Timer &time = timers["PWPBDD"];
                            probability_t &probability = probabilities["PWPBDD"];
                            pwpbdd_.SetEvidence(evidence);
                            pwpbdd_.InitCache();
                            probability = pwpbdd_.Posterior();
                            pwpbdd_.InitCache();
                            time.Start();
                            probability = pwpbdd_.Posterior();
                            time.Stop();
                            time.Add();
                        } catch (ModelCounterException &exception){
                            Print(ERR, "                                                                      \n");
                            Print(ERR, "PWPBDD: %s\n", exception.what());
                            Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                            std::string query = evidence.GetQueryString();
                            Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            Print(ERR, "    Query    : %s\n",query.c_str());
                        }

                        try {
                            // execute 2x to eliminate cache advantage
                            Timer &time = timers["CPWPBDD"];
                            probability_t &probability = probabilities["CPWPBDD"];
                            pwpbdd_.SetEvidence(evidence);
                            pwpbdd_.InitCache();
                            probability = pwpbdd_.Posterior<6>();
                            pwpbdd_.InitCache();
                            time.Start();
                            probability = pwpbdd_.Posterior<6>();
                            time.Stop();
                            time.Add();
                        } catch (ModelCounterException &exception){
                            Print(ERR, "                                                                      \n");
                            Print(ERR, "CPWPBDD: %s\n", exception.what());
                            Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                            std::string query = evidence.GetQueryString();
                            Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            Print(ERR, "    Query    : %s\n",query.c_str());
                        }

                        try {
                            // execute 2x to eliminate cache advantage
                            Timer &time = timers["PWPBDD2"];
                            probability_t &probability = probabilities["PWPBDD2"];
                            pwpbdd_.SetEvidence(evidence);
                            pwpbdd_.InitCache();
                            probability = pwpbdd_.Posterior<2>();
                            pwpbdd_.InitCache();
                            time.Start();
                            probability = pwpbdd_.Posterior<2>();
                            time.Stop();
                            time.Add();
                        } catch (ModelCounterException &exception){
                            Print(ERR, "                                                                      \n");
                            Print(ERR, "PWPBDD: %s\n", exception.what());
                            Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                            std::string query = evidence.GetQueryString();
                            Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            Print(ERR, "    Query    : %s\n",query.c_str());
                        }

                        if(false)
                        for(auto it = manager.workers.begin(); it != manager.workers.end(); it++){
                            const unsigned int &kWorkers = *it;
                            //try {
                            //    // execute 2x to eliminate cache advantage
                            //    const unsigned int imp = 4;
                            //    std::string name = stringf("PPWPBDD%u - %u cores",imp,kWorkers);
                            //    Timer &time = timers[name];
                            //    probability_t &probability = probabilities[name];
                            //    pwpbdd_.SetEvidence(evidence);
                            //    pwpbdd_.SetNodeIds(evidence);
                            //    pwpbdd_.SetEvidenceCache(evidence);
                            //    pwpbdd_.InitCache();
                            //    probability = pwpbdd_.ParallelPosterior<imp>(kWorkers);
                            //    pwpbdd_.InitCache();
                            //    //time.Start();
                            //    probability = pwpbdd_.ParallelPosterior<imp>(kWorkers,&time);
                            //    //time.Stop();
                            //    //time.Add();
                            //} catch (ModelCounterException &exception){
                            //    Print(ERR, "                                                                      \n");
                            //    Print(ERR, "PWPBDD: %s\n", exception.what());
                            //    Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                            //    std::string query = evidence.GetQueryString();
                            //    Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            //    Print(ERR, "    Query    : %s\n",query.c_str());
                            //}
                            try {
                                // execute 2x to eliminate cache advantage
                                const unsigned int imp = 6;
                                std::string name = stringf("PPWPBDD%u - %u cores",imp,kWorkers);
                                Timer &time = timers[name];
                                probability_t &probability = probabilities[name];
                                pwpbdd_.SetEvidence(evidence);
                                pwpbdd_.SetNodeIds(evidence);
                                pwpbdd_.SetEvidenceCache(evidence);
                                pwpbdd_.InitCache();
                                probability = pwpbdd_.ParallelPosterior<imp>(kWorkers);
                                pwpbdd_.InitCache();
                                //time.Start();
                                probability = pwpbdd_.ParallelPosterior<imp>(kWorkers,&time);
                                //time.Stop();
                                //time.Add();
                            } catch (ModelCounterException &exception){
                                Print(ERR, "                                                                      \n");
                                Print(ERR, "PWPBDD: %s\n", exception.what());
                                Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                                std::string query = evidence.GetQueryString();
                                Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                                Print(ERR, "    Query    : %s\n",query.c_str());
                            }
                            //try {
                            //    // execute 2x to eliminate cache advantage
                            //    const unsigned int imp = 8;
                            //    std::string name = stringf("PPWPBDD%u - %u cores",imp,kWorkers);
                            //    Timer &time = timers[name];
                            //    probability_t &probability = probabilities[name];
                            //    pwpbdd_.SetEvidence(evidence);
                            //    pwpbdd_.SetNodeIds(evidence);
                            //    pwpbdd_.SetEvidenceCache(evidence);
                            //    pwpbdd_.InitCache();
                            //    probability = pwpbdd_.ParallelPosterior<imp>(kWorkers);
                            //    pwpbdd_.InitCache();
                            //    //time.Start();
                            //    probability = pwpbdd_.ParallelPosterior<imp>(kWorkers,&time);
                            //    //time.Stop();
                            //    //time.Add();
                            //} catch (ModelCounterException &exception){
                            //    Print(ERR, "                                                                      \n");
                            //    Print(ERR, "PWPBDD: %s\n", exception.what());
                            //    Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                            //    std::string query = evidence.GetQueryString();
                            //    Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            //    Print(ERR, "    Query    : %s\n",query.c_str());
                            //}
                            //try {
                            //    // execute 2x to eliminate cache advantage
                            //    const unsigned int imp = 7;
                            //    std::string name = stringf("PPWPBDD%u - %u cores",imp,kWorkers);
                            //    Timer &time = timers[name];
                            //    probability_t &probability = probabilities[name];
                            //    pwpbdd_.SetEvidence(evidence);
                            //    pwpbdd_.SetNodeIds(evidence);
                            //    pwpbdd_.SetEvidenceCache(evidence);
                            //    pwpbdd_.InitCache();
                            //    probability = pwpbdd_.ParallelPosterior<imp>(kWorkers);
                            //    pwpbdd_.InitCache();
                            //    //time.Start();
                            //    probability = pwpbdd_.ParallelPosterior<imp>(kWorkers,&time);
                            //    //time.Stop();
                            //    //time.Add();
                            //} catch (ModelCounterException &exception){
                            //    Print(ERR, "                                                                      \n");
                            //    Print(ERR, "PWPBDD: %s\n", exception.what());
                            //    Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                            //    std::string query = evidence.GetQueryString();
                            //    Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                            //    Print(ERR, "    Query    : %s\n",query.c_str());
                            //}
                        }
                    }


                    #ifdef ACE
                    if(manager.have_ace){
                        try {
                            Timer &time = timers["ACE"];
                            probability_t &probability = probabilities["ACE"];


                            probability = ace_posterior(evidence,time);
                        } catch (ModelCounterException &exception){
                           Print(ERR, "                                                                      \n");
                           Print(ERR, "ACE: %s\n", exception.what());
                           Print(ERR, "    Query   : %u (%u)\n", count, manager.total_query_count);
                           std::string query = evidence.GetQueryString();
                           Print(ERR, "    Query id : %u (%u)\n", count, manager.total_query_count);
                           Print(ERR, "    Query    : %s\n",query.c_str());
                       }
                    }
                    #endif

                    if(kExhaustiveType == ExhaustiveType::COMPARE){
                        Print(MSG, "  Query %5u (%5u)", count, manager.total_query_count);

                        std::set<std::pair<double,std::string>> sorted_timers;
                        for (auto itr = timers.begin(); itr != timers.end(); ++itr)
                            sorted_timers.insert(std::make_pair(itr->second.GetTotal<Timer::Milliseconds>(),itr->first));

                        //for(auto it = timers.begin(); it != timers.end(); it++)
                        //    Print(MSG, "    %7s (%.4lf)", it->first.c_str(), it->second.GetTotal<Timer::Seconds>());
                        for(auto it = sorted_timers.begin(); it != sorted_timers.end(); it++)
                            Print(MSG, "  %7s (%.4lf)", it->second.c_str(), it->first);
                        Print(MSG, "\r");
                    } else if(kExhaustiveType == ExhaustiveType::VERIFY) {
                        Print(MSG, "  Query %5u (%5u)\r", count, manager.total_query_count);
                    } else {
                        std::string query = evidence.GetQueryString();
                        Print(MSG, "  Query %5u (%5u): %s\n", count, manager.total_query_count,query.c_str());
                    }

                    if(kExhaustiveType == ExhaustiveType::VERIFY){
                        const double kThreshold = 0.01;
                        probability_t probability = -1337;
                        if (probabilities.size() <= 1)
                            throw InterfaceException("Verification requires more that one inference method.");
                        for(auto it = probabilities.begin(); it != probabilities.end(); it++){
                            if(it->first == "pwpbdd+")
                                continue;

                            if(is_nan(it->second) || it->second < -1 || it->second > 1)
                                it->second = -1;

                            bool invalid = false;
                            if(probability == -1337)
                                probability = it->second;
                            else if(std::abs(it->second - probability) > kThreshold)
                                invalid = true;

                            if(invalid){
                                std::string query = evidence.GetQueryString();
                                Print(ERR, "THRESHOLD BREACH\n");
                                Print(ERR, "    On query : %s\n\n",query.c_str());
                                for(auto it = probabilities.begin(); it != probabilities.end(); it++)
                                    Print(ERR, "    %8s : %.4lf  (%.4lf)\n", it->first.c_str(), it->second, std::abs(it->second - probability));
                                printf("\n");

                                FILE *file = fopen("queries.log","a");
                                fprintf(file,"%s\n",query.c_str());
                                fclose(file);
                            }
                        }
                    }
                }

                // ================== inference ====================
                bool stop = true;
                for(int i = M; i >= 0; i--){
                    if(vars[i]){
                        if(ctr[i] < dimension[i]-1){
                            stop = false;
                            ctr[i]++;
                            for(unsigned int j = i+1; j <= M; j++)
                                ctr[j] = 0;
                            break;
                        }
                    }
                }
                if(stop)
                    break;

            }

            if(instantiations < VARIABLES){
                int leftmost = VARIABLES-1;
                while(leftmost >= 0 && vars[leftmost]) leftmost--;
                while(leftmost >= 0 && !vars[leftmost]) leftmost--;
                if(leftmost >= 0){
                    vars[leftmost+1] = 1;
                    vars[leftmost] = 0;
                } else break;
            } else break;
        }
    }

    Print(MSG, "\n\n");
    Print(MSG,     "Total queries   : %lu\n", manager.total_query_count);

    {
        std::set<std::pair<double,std::string>> sorted_timers;
        for (auto itr = timers.begin(); itr != timers.end(); ++itr)
            sorted_timers.insert(std::make_pair(itr->second.GetTotal<Timer::Milliseconds>(),itr->first));

        for(auto it = sorted_timers.begin(); it != sorted_timers.end(); it++)
            Print(MSG, "    %15s time : %10.3lfms\n", it->second.c_str(), it->first);
    }
    Print(MSG, "\n");
}

int Interface::Compare(void*){
    int max_evidence = 9999999;
    int instantiations = 1;
    int maxLocalIterations = -1;

    std::stringstream args(arguments_);
    std::vector<int> ops;
    for (int n; args >> n;)
        ops.push_back(n);

    if (ops.size() >= 1)
        instantiations = ops[0];

    if (ops.size() >= 2)
        max_evidence = ops[1];

    if (ops.size() >= 3)
        maxLocalIterations = ops[2];

    Exhaustive(ExhaustiveType::COMPARE, instantiations, max_evidence, maxLocalIterations);

    return 0;
}

int Interface::Verify(void*){
    int max_evidence = 9999999;
    int instantiations = 1;
    int maxLocalIterations = -1;

    std::stringstream args(arguments_);
    std::vector<int> ops;
    for (int n; args >> n;)
        ops.push_back(n);

    if (ops.size() >= 1)
        instantiations = ops[0];

    if (ops.size() >= 2)
        max_evidence = ops[1];

    if (ops.size() >= 3)
        maxLocalIterations = ops[2];

    Exhaustive(ExhaustiveType::VERIFY, instantiations, max_evidence, maxLocalIterations);

    return 0;
}


#ifdef ACE
int Interface::InitAce(void*){
    manager.have_ace = true;
    return 0;
}
#endif

}
