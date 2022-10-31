#include <stdio.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "exceptions.h"
#include "files.h"

using namespace std;

filename::filename(){
}

void filename::clear(){
    dir.clear();
    basename.clear();
    filenames.clear();
}

filename& filename::operator=(filename &a){
    this->basename = a.basename;
    this->dir = a.dir;
    return *this;
}

std::string filename::get_type(bdd_t BDD){
    if(BDD == bdd_t::wpbdd)
        return std::string("wpbdd");
    else if(BDD == bdd_t::multigraph)
        return std::string("multigraph");
    else if(BDD == bdd_t::tdmultigraph)
        return std::string("tdmultigraph");
    else if(BDD == bdd_t::none)
        return std::string("");
    else
        return std::string("?");
}

bool filename::exists(string filename){
    if( access( filename.c_str(), F_OK ) != -1 )
        return true;
    else return false;
}

bool filename::exists(file_t type){
    return exists(get_filename(type));
}

bool filename::has_basename(){
    if(basename.empty())
        return false;
    else return true;
}

std::string filename::remove_extension(std::string filename){
    size_t found = filename.find_last_of(".");
    if(found != string::npos)
        return filename.substr(0,found);
    else return "";
}

void filename::set_basename(string filename){
    filename = remove_extension(filename);
    dir.clear();
    basename.clear();

    size_t found = filename.find_last_of("/");
    if(found != string::npos){
        dir = filename.substr(0,found);
        filename = filename.substr(found+1);
    }
    basename = filename;
}

std::string filename::get_prefix(string aux){
    string ret;
    if(dir.length() == 0)
        ret = basename;
    else
        ret = dir + "/" + basename;

    if(!aux.empty())
        ret += "." + aux;
    return ret;
}

std::string filename::create_filename(file_t type, bdd_t BDD){
    return create_filename(type, get_type(BDD));
}

std::string filename::get_extension(std::string filename){
    size_t found = filename.find_last_of(".");
    std::string extension = filename.substr(found+1);
    return extension;
}

std::string filename::create_filename(file_t type, std::string aux){
    string ret = get_prefix(aux);
    switch(type){
        case BN:
            ret += ".net";
            break;
        case CNF:
            ret += ".cnf";
            break;
        case SPANNING:
            ret += ".spanning.dot";
            break;
        case DOT:
            ret += ".dot";
            break;
        case OUTPUT:
            ret += ".out";
            break;
        case ORDERING:
            ret += ".lord";
            break;
        case UAI:
            ret += ".uai";
            break;
        case PARTITION:
            ret += ".part";
            break;
        case ELIM_ORDERING:
            ret += ".num";
            break;
        case VARIABLE_ORDERING:
            ret += ".ord";
            break;
        case PSEUDO_ORDERING:
            ret += ".pseudo";
            break;
        case COMPOSITION_ORDERING:
            ret += ".comp";
            break;
        case MULTIGRAPH:
            ret += ".mc.ac";
            break;
        case TDMULTIGRAPH:
            ret += ".tdmc.ac";
            break;
        case PTREE:
            ret += ".ptree";
        case WPBDD:
        case AC:
            ret += ".ac";
            break;
        case MAPPING:
            ret += ".map";
            break;
        //Default:
        //    throw compiler_exception("Non specified file format");
        //    break;
    }
    return ret;
}

void filename::set_filename(file_t filetype, std::string filename){
    filenames[filetype] = filename;
}

const char * filename::get_filename_c(file_t filetype, std::string aux, std::string postfix){
    static std::string filename;
    filename = get_filename(filetype,aux,postfix);
    return filename.c_str();
}

bool filename::is_set(file_t filetype){
    auto hit = filenames.find(filetype);
    if(hit != filenames.end())
        return true;
    else return false;
}

std::string filename::get_filename(file_t filetype,std::string aux, std::string postfix){
    std::string filename;
    auto hit = filenames.find(filetype);
    if(hit != filenames.end())
        filename = hit->second;
    else {
        // create filename when it does not exist
        if(!has_basename())
            basename = "tmp";

        std::string newfilename = create_filename(filetype);
        set_filename(filetype, newfilename);
        filename = filenames[filetype];
    }
    if(aux != ""){
        size_t found = filename.find_last_of(".");
        aux = "." + aux;
        filename.insert(found, aux.c_str());
    }
    if(postfix != ""){
        postfix = "." + postfix;
        filename.append(postfix);
    }
    return filename;
}

