#include "manager.h"
#include "exceptions.h"

namespace bnmc {

Manager::Manager(){
    bn = NULL;

	have_tdmultigraph =
        have_multigraph =
        have_pwpbdd =
        have_partition =
        have_wpbdd =
        have_ace =
        have_ordering =
        have_mapping =
        have_bn =
        have_filename =
        have_evidence = false;

//        have_pwpbdd =
//        have_tdmultigraph =
//        have_multigraph =
//        have_partition =
//        have_wpbdd =
//        have_ace =
//        have_ordering =
//        have_mapping =
//        have_bn =
//        have_filename =
//        have_evidence =
//        false;

    buffer = 0;
    workers.push_back(0); // auto detect
    OPT_PARALLEL_PARENTS = false;
    OPT_PARALLEL_COMPONENTS = true;
    OPT_PARALLEL_JOINTS = false;
}

Manager::~Manager(){
    if(bn)
        delete(bn);

}

void Manager::Read(file_t type){
    switch(type){
        case MAPPING:
            if(!have_mapping){
                mapping.read(files.get_filename_c(file_t::MAPPING));
                have_mapping = true;
            }
            break;
        case BN:
            if(!have_bn){
                try {
                    // store BN
                    bn = bayesnet::read(files.get_filename_c(file_t::BN));
                    have_bn = true;
                } catch(bayesnet_exception &e){
                    throw IoException("%s", e.what());
                }
            }
            break;
        default:
            throw IoException("Unknown read option");
    }
}

}

