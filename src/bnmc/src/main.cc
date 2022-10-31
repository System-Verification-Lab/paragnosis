#include "interface.h"
#include "options.h"
#include <stdio.h>
#include <unistd.h>
#include "memory.h"

using namespace bnmc;

void help(){
    fprintf(stderr, "Usage:\n   ./bnmc [-w <#workers>]\n\n");
}

void limit_memory(){
    double percentage = 0.8;
    const byte_t Gb = 1073741824;
    const byte_t free_bytes = get_free_ram_size();
    const byte_t limit_bytes = ((free_bytes/100)*(percentage*100));
    printf("Free memory : %.2lf Gb\n", (double)free_bytes/(double)Gb);
    printf("Set limit   : %.2lf Gb\n\n", (double)limit_bytes/(double)Gb);
    if(set_memory_limit(get_ram_size(percentage)) != 0)
        fprintf(stderr, "Warning: could not set memory limit.");
}

int main(int argc, char **argv){
    limit_memory();

    int c;
    bool clear_workers = true;
    while ((c = getopt(argc, argv, "w:b:h")) != -1){
        switch (c){
            case'w':
                if(clear_workers){
                    manager.workers.clear();
                    clear_workers = false;
                }
                manager.workers.push_back(std::atoi(optarg));
                break;
            case'b':
                manager.buffer = std::atoi(optarg);
                break;
            case 'h':
                help();
                return 1;
            default:
                help();
                return 1;
        }
    }


    bnmc::Interface cli;
    cli.Prompt();

    return 0;
}

