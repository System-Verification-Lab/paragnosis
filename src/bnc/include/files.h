#ifndef FILES_H
#define FILES_H

#include <string>
#include "types.h"
#include <map>

enum file_t { PARTITION, BN, CNF, DOT, OUTPUT, ORDERING, VARIABLE_ORDERING, PSEUDO_ORDERING, ELIM_ORDERING, COMPOSITION_ORDERING, MAPPING, AC, WPBDD, PWPBDD, PMULTIGRAPH, PTDMULTIGRAPH, MULTIGRAPH, TDMULTIGRAPH, PTREE, UAI, SPANNING};

class filename {
    public:
        filename();
        filename& operator=(filename &);

        bool is_set(file_t);
        void clear();
        bool has_basename();
        std::string remove_extension(std::string);
        std::string get_extension(std::string);
        void set_basename(std::string);
        void set_filename(file_t,std::string);
        std::string get_filename(file_t,std::string aux = "",std::string postfix = "");
        const char * get_filename_c(file_t,std::string aux = "",std::string postfix = "");
        void store_filename(file_t,std::string&);
        std::string create_filename(file_t type, bdd_t BDD);
        std::string create_filename(file_t type, std::string aux = "");
        bool exists(std::string);
        bool exists(file_t);
        std::string get_name(file_t, bdd_t);
        std::string get_type(bdd_t);
        std::string get_prefix(std::string aux = "");
    private:
        std::string basename;
        std::string dir;
        std::map<file_t,std::string> filenames;
};

typedef filename filename_t;

#endif
