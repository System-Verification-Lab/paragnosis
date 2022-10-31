#include "memory.h"

#include <unistd.h>
#include <sys/types.h>
#include <sys/resource.h>


namespace bnmc {

#ifdef __APPLE__

#include <sys/sysctl.h>

byte_t get_ram_size(double percentage){
    int mib[4];
    int64_t physical_memory;
    size_t len = sizeof(int64_t);
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    len = sizeof(int64_t);
    sysctl(mib, 2, &physical_memory, &len, NULL, 0);
    return (physical_memory/100)*(percentage*100);
}

byte_t get_free_ram_size(){
    return get_ram_size(50);
}

#else
#include <sys/sysinfo.h>

byte_t get_free_ram_size(){
    struct sysinfo sys_info;

    if(sysinfo(&sys_info) == 0)
        return (sys_info.freeram + sys_info.bufferram) * sys_info.mem_unit;
    else return 0;
}

byte_t get_ram_size(double percentage){

    struct sysinfo sys_info;

    if(sysinfo(&sys_info) == 0)
        return ((sys_info.totalram/100)*(percentage*100)) * sys_info.mem_unit;
    return 0;
}

#endif

int set_memory_limit(byte_t bytes){
    struct rlimit limit;
    if(bytes == 0)
        limit.rlim_cur = get_free_ram_size();
    else limit.rlim_cur = bytes;
    limit.rlim_max = RLIM_INFINITY;

    // Limits total memory size. To limit heap size change RLIMIT_AS to RLIMIT_DATA
    //if (-1 == setrlimit(RLIMIT_AS, &limit))
    //    return 1;
    // else return 0;

    //return prlimit(getpid(),RLIMIT_AS,&limit,NULL);
    return setrlimit(RLIMIT_AS,&limit);
}

}
