#include <stdio.h>
#include <signal.h>
#include <setjmp.h>
#include <stddef.h>
#include <readline/readline.h>
#include <readline/history.h>
#include "io.h"
#include "cli.h"
#include "exceptions.h"
#include "trim.h"

#define EXECUTE(function) ((*this).*function)(NULL)

bool sig_quit = false;
sigjmp_buf ctrlc_buf;

void handle_signals(int signo) {
    switch(signo){
        case SIGINT:
            Print(MSG,"\n");
            siglongjmp(ctrlc_buf, 1);
            break;
        case SIGQUIT:
            sig_quit = true;
            break;
    }
}

void* xmalloc (int size){
    void *buf;

    buf = malloc (size);
    if (!buf){
        Print(MSG, "Error: Out of memory. Exiting.'n");
        exit (1);
    }

    return buf;
}

char* dupstr (const char* s) {
    char *r;

    r = (char*) xmalloc ((strlen (s) + 1));
    strcpy (r, s);
    return (r);
}


namespace bnmc {

std::vector< std::string > Cli::command_names_;

Cli::Cli(){
    AddCommand("quit", &Cli::Quit);
}

char* Cli::TabCompletionGenerator(const char* text, int state){
    static unsigned int list_index, len;

    if (!state) {
        list_index = 0;
        len = strlen (text);
    }

    while(list_index < command_names_.size()){
        const char *name = command_names_[list_index++].c_str();
        if (strncmp (name, text, len) == 0)
            return (dupstr(name));
    }

    /* If no names matched, then return NULL. */
    return ((char *)NULL);
}

char** Cli::TabCompletion( const char * text , int start,  int end){
    char **matches;

    matches = (char **)NULL;

    if (start == 0)
        matches = rl_completion_matches ((char*)text, &TabCompletionGenerator);
    else
        rl_bind_key('\t',rl_abort);

    return (matches);
}

void Cli::AddCommand(std::string name, CommandPtr f){
    if(commands_.find(name) == commands_.end()){
        commands_[name] = f;
        command_names_.push_back(name);
    } else throw CliException("Function '%s' is already added", name.c_str());
}

int Cli::RunCommand(std::string &command){
    std::string function;
    arguments_.clear();
    size_t pos = command.find_first_of(" ");
    if(pos != std::string::npos){
        function = command.substr(0,pos);
        arguments_ = command.substr(pos+1);
        Trim(arguments_);
    } else function = command;

    auto hit = commands_.find(function);
    if(hit != commands_.end()){
        CommandPtr command = hit->second;
        return EXECUTE(command);
    } else {
        Print(ERR,"%s: command not found. \n", function.c_str());
        Print(ERR,"Type 'help' to view commands.\n");
    }
    return 1;
}

int Cli::Quit(void*){
    if(arguments_.length() > 0){
        Print(ERR, "Unexpected argument: %s\n", arguments_.c_str());
        return 1;
    }

    sig_quit = true;
    return 0;
}

void Cli::Greet(){
    Print(MSG,"Type 'help' to view commands.\n\n");
    Print(MSG,"Press Ctrl-d to quit...\n\n");
}

void Cli::Prompt(){
    // ctrl-c -> cancel and ctrl-d -> quit
    if ( signal(SIGINT, handle_signals) == SIG_ERR ||
         signal(SIGQUIT, handle_signals) == SIG_ERR) {
        Print(ERR, "failed to register interrupt signals\n");
    }
    while ( sigsetjmp( ctrlc_buf, 1 ) != 0 );

    // register for autocomplete
    // rl_attempted_completion_function = TabCompletion;

    char *line = NULL;
    const char std_prompt[] = "> ";
    const char err_prompt[] = "\033[31;1m>\033[0m ";
    int failed = 0;
    while(!sig_quit) {
        if(failed)
            line = readline(err_prompt);
        else
            line = readline(std_prompt);

        if(!line)
            break;

        //enable auto-complete
        // rl_bind_key('\t',rl_complete);

        std::string command = line;
        Trim(command);
        if(command.length() > 0) {
            add_history(line);
            failed = RunCommand(command);
        }
        free(line);
    }

    Print(MSG,"\n\n\033[32;1mSee you next time!\033[0m\n");
    #ifndef DEBUG
    exit(0);
    #endif
}

}

