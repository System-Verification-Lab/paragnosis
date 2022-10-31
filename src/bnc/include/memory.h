#ifndef MEMORY_H
#define MEMORY_H


typedef unsigned long byte_t;
byte_t get_free_ram_size(); // returns bytes
byte_t get_ram_size(double percentage = 1);      // returns bytes
int set_memory_limit(byte_t);

#endif
