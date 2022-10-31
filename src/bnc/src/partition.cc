#include "partition.h"
#include "exceptions.h"
#include <cstring>
#include <set>

void partition_t::determine_cutset(bayesnet *bn){
    // determine cutset
    cutset.clear();
    for(auto it = set.begin(); it != set.end(); it++){
        variable_t variable = *it;
        uint32_t* parents = bn->get_parent(variable);
        for(unsigned int i = 0; i < bn->get_parent_size(variable); i++){
            uint32_t parent = parents[i];
            if(set.find(parent) == set.end())
                cutset.insert(parent);
        }
    }
}

unsigned int variable_name_to_id(bayesnet *bn, char * name){
    int variable_index = bn->get_node_id(name);
    if(variable_index >= 0)
        return variable_index;
    else throw compiler_parse_exception("Error while parsing: variable '%s' unknown", name);
}

void bn_partitions_t::SetOnePartition(bayesnet *bn){
    clear();

    const size_t kNrVariables = bn->get_nr_variables();
    resize(1);

    partition_t &partition = operator[](0).partition;
    for (unsigned int i = 0; i < bn->get_nr_variables(); i++)
        partition.set.insert(partition.set.end(), i);
}

void bn_partitions_t::VerifyOrdering(bayesnet *bn) const{

    for(unsigned int i = 0; i < size(); i++){
        const bn_partition_t &partition = operator[](i);

        // check variable ordering
        if(partition.variable_ordering.size() != 0){
            partition_set_t check, ordering_check;
            check.insert(partition.partition.set.begin(),partition.partition.set.end());
            check.insert(partition.partition.cutset.begin(),partition.partition.cutset.end());
            ordering_check.insert(partition.variable_ordering.begin(),partition.variable_ordering.end());

            if(check != ordering_check){
                std::string full_set = "{";
                for(auto it = partition.partition.set.begin(); it != partition.partition.set.end(); it++){
                    if(it != partition.partition.set.begin())
                        full_set.append(",");
                    variable_t variable = *it;
                    full_set.append(bn->get_node_name(variable));
                }
                for(auto it = partition.partition.cutset.begin(); it != partition.partition.cutset.end(); it++){
                    full_set.append(",");
                    variable_t variable = *it;
                    full_set.append(bn->get_node_name(variable));
                }
                full_set.append("}");

                throw compiler_parse_exception("ordering in partition %u does not consist of partition AND cutset: %s", i+1, full_set.c_str());
            }
        }
    }
}

void bn_partitions_t::verify(bayesnet *bn) const{
    partition_set_t check_cover;
    for (unsigned int i = 0; i < bn->get_nr_variables(); i++)
        check_cover.insert(check_cover.end(), i);

    for(unsigned int i = 0; i < size(); i++){
        const bn_partition_t &partition = operator[](i);
        for(auto it = partition.partition.set.begin(); it != partition.partition.set.end(); it++){
            variable_t variable = *it;
            if(check_cover.find(variable) == check_cover.end())
                throw compiler_parse_exception("variable '%s' occurs in partition %u for the second time", bn->get_node_name(variable).c_str(), i+1);
            check_cover.erase(variable);
        }
    }

    if(check_cover.size() > 0){
        std::string full_set = "{";
        for(auto it = check_cover.begin(); it != check_cover.end(); it++){
            if(it != check_cover.begin())
                full_set.append(",");
            variable_t variable = *it;
            full_set.append(bn->get_node_name(variable));
        }
        full_set.append("}");
        throw compiler_parse_exception("Partitions do not cover all BN variables, missing variable(s): %s", full_set.c_str());
    }
}

void bn_partitions_t::read(std::string filename, bayesnet *bn){
    FILE *file = fopen(filename.c_str(),"r");
    if(file){
        unsigned int nr_partitions;
        if(fscanf(file, "partition %u\n", &nr_partitions) != 1)
            throw compiler_parse_exception("Error while parsing header of partition file");

        if(nr_partitions < 2 || (bn != NULL && nr_partitions > bn->get_nr_variables()))
            throw compiler_parse_exception("Number of partitions must be greater than 1 and smaller than the number of BN variables");

        resize(nr_partitions);
        for(unsigned int i = 0; i < nr_partitions; i++){
            bn_partition_t &partition = operator[](i);
            char delimiter[4];
            if(bn != NULL){
                char word[100] = {0};
                while(true){
                    int ret;
                    ret = fscanf(file,"%99[^,^\n]",word);
                    if(ret == EOF)
                        break;
                    else if(ret == 1)
                        partition.partition.set.insert(variable_name_to_id(bn, word));
                    else throw compiler_parse_exception("Error while parsing: partition variable name expected");

                    ret = fscanf(file,"%1[,\n]", &delimiter);
                    if(ret == EOF)
                        break;
                    else if(strcmp(delimiter,"\n") == 0)
                        break;
                    else if(ret != 1)
                        throw compiler_parse_exception("Error while parsing: delimiter expected");
                }
            } else {
                unsigned int v;
                while(true){
                    int ret;
                    ret = fscanf(file,"%u",&v);
                    if(ret == EOF)
                        break;
                    else if(ret == 1)
                        partition.partition.set.insert(v);
                    else throw compiler_parse_exception("Error while parsing: partition variable number expected");

                    ret = fscanf(file,"%1[,\n]", &delimiter);
                    if(ret == EOF)
                        break;
                    else if(strcmp(delimiter,"\n") == 0)
                        break;
                    else if(ret != 1)
                        throw compiler_parse_exception("Error while parsing: delimiter expected");
                }
            }

            if(partition.partition.set.size() == 0)
                throw compiler_parse_exception("Partition cannot be empty");

            // determine cutset, given partition and bayesian network
            partition.partition.determine_cutset(bn);
        }

        fclose(file);
    } else throw compiler_io_exception("Could not open partition file %s\n", filename.c_str());
}

void bn_partitions_t::write(std::string filename, bayesnet *bn) const {
    FILE *file = fopen(filename.c_str(),"w");
    if(file){
        fprintf(file, "partition %u\n", size());
        for(auto it = begin(); it != end(); it++){
            const partition_t &kPartition = it->partition;
            for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++){
                if(it != kPartition.set.begin())
                    fprintf(file, ",");
                variable_t variable = *it;
                if(bn == NULL)
                    fprintf(file, "%u", variable);
                else
                    fprintf(file, "%s", bn->get_node_name(variable).c_str());
            }
            fprintf(file, "\n");
        }
        fclose(file);
    } else throw compiler_io_exception("Could not open partition file %s\n", filename.c_str());
}

bool bn_partitions_t::monolithic(){
    return size() == 1;
}

unsigned int bn_partitions_t::variable_to_partition_id(variable_t variable){
    for(unsigned int i = 0; i < size(); i++){
        partition_t &partition = operator[](i).partition;
        if(partition.set.find(variable) != partition.set.end())
            return i;
    }
    throw compiler_exception("Could not find partition for variable %u\n", variable);
}

