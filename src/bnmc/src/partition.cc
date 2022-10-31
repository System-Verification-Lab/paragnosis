#include "partition.h"
#include "exceptions.h"
#include "options.h"
#include <climits>
#include <set>
#include <unordered_map>
#include <bnc/multigraphpdef.h>

namespace bnmc {

const size_t Partition::GetCircuitSize() const {
    if(ac_.GetSize())
        return ac_.GetSize();
    else if(mgc_.GetSize())
        return mgc_.GetSize();
    else return tdmgc_.GetSize();
}

void Partition::Read(const ModelType kModelType, std::string filename){
    if(kModelType == ModelType::PWPBDD || kModelType == ModelType::WPBDD){
        FILE *file = fopen(filename.c_str(), "r");
        if(file){
            size_t nr_of_nodes;
            if(fscanf(file, "wpbdd %lu\n", &nr_of_nodes) != 1)
                throw IoException("Error while reading preemble");

            const std::vector<probability_t>& weight_to_probability = manager.mapping.get_weight_to_probability();
            const unsigned int LITERALS = manager.mapping.get_nr_literals();

            if(!ac_.Resize(nr_of_nodes))
                throw IoException("Could not allocate %lu nodes", nr_of_nodes);

            for(unsigned int i = 0; i < nr_of_nodes; i++){
                unsigned int nr_of_weights;
                LiteralNode n;
                if(fscanf(file, "%u %u %u %u", &n.l, &n.t, &n.e, &nr_of_weights) != 4)
                    throw IoException("Error while reading node %u", i);

                probability_t w = 1;
                for(unsigned int j = 0; j < nr_of_weights; j++){
                    unsigned int symbolic_w;
                    if(fscanf(file, " %u", &symbolic_w) != 1)
                        throw IoException("Error while reading weight %u of node %u", j, i);

                    assert((symbolic_w == 0 || symbolic_w == 1) || (symbolic_w >= LITERALS+1));
                    if(symbolic_w == 0 || symbolic_w == 1)
                        w *= (double) symbolic_w;
                    else
                        w *= weight_to_probability[symbolic_w-(LITERALS+1)];
                }
                if(i == 0)
                    w = 0;

                ac_[i] = n;
                ac_[i].w = w;
                if(i == 0 || i == 1)
                    ac_[i].v = manager.mapping.get_nr_variables();
                if(fscanf(file, "\n") == EOF)
                    throw IoException("Error while reading, expecting end-of-line at node %u", i);
            }
            fclose(file);
        } else throw IoException("Could not open bdd file '%s'\n", filename.c_str());
    } else if (kModelType == ModelType::MULTIGRAPH || kModelType == ModelType::PMULTIGRAPH || kModelType == ModelType::TDMULTIGRAPH || kModelType == ModelType::PTDMULTIGRAPH){
        FILE *file = fopen(filename.c_str(), "rb");
        if(file){
            bool is_tree;
            if(fread(&is_tree, sizeof(bool), 1, file) != 1)
                throw IoException("Could not determine WPBDD type");

            decltype(MultiGraph::Size::nodes) nr_of_nodes;
            if(fread(&nr_of_nodes, sizeof(decltype(MultiGraph::Size::nodes)), 1, file) != 1)
                throw IoException("Could not determine nr of nodes");

            decltype(MultiGraph::Size::edges) nr_of_edges;
            if(fread(&nr_of_edges, sizeof(decltype(MultiGraph::Size::edges)), 1, file) != 1)
                throw IoException("Could not determine nr of edges");

            if(kModelType == ModelType::MULTIGRAPH && is_tree)
                throw IoException("Multigraph cannot be tree driven");

            auto &mg = (kModelType == ModelType::TDMULTIGRAPH || kModelType == ModelType::PTDMULTIGRAPH? tdmgc_:mgc_);

            const size_t kNrTerminals = 1;
            if(!mg.Resize(nr_of_nodes+kNrTerminals,nr_of_edges))
                throw IoException("Failed to allocate %lu nodes and %lu edges", nr_of_nodes, nr_of_edges);

            // format
            // weights : #weights weight_1 ... weight_n
            // edge    : index weights
            // edges   : #edges edge_1 ... edge_m
            // type    : < + | * >
            // node    : type variable edges
            // preamble: multigraph <chain|tree> #total_nodes #total_edges
            // line    : <node | preamble>
            // read terminal

            // create terminal(s)
            mg.SetNrTerminals(kNrTerminals);
            mg.SetRootIndex(kNrTerminals);
            for(unsigned int i = 0; i < kNrTerminals; i++){
                auto* terminal = mg.CreateNode(0);
                terminal->variable = 0;
                mg.StoreNode(0);
            }

            // read nodes
            for(size_t i = kNrTerminals; i < nr_of_nodes + kNrTerminals; i++){
                // read variable and type of node
                decltype(MultiGraph::Node::variable) variable;
                if(fread(&variable, sizeof(decltype(MultiGraph::Node::variable)), 1, file) != 1)
                    throw IoException("Could not determine variable of node %lu", i);

                // read nr of edges
                decltype(MultiGraph::Node::size) nr_of_edges;
                if(fread(&nr_of_edges, sizeof(decltype(MultiGraph::Node::size)), 1, file) != 1)
                    throw IoException("Could not determine nr of edges of node %lu", i);

                // create node
                MultiGraph::Node *node = mg.CreateNode(nr_of_edges);
                node->variable = variable;
                mg.StoreNode(nr_of_edges);

                // read edges
                const bool kIsAnd = node->IsAnd();
                const auto kEdgeEnd = node->EdgeEnd();
                auto edge = node->EdgeBegin();
                while(edge != kEdgeEnd){
                    size_t to;
                    if(fread(&to, sizeof(size_t), 1, file) != 1)
                        throw IoException("Could not determine edge %lu of node %lu", node->EdgeBegin()-edge, i);
                    edge->to = mg.GetNode(to);

                    if(kIsAnd)
                        edge->probability = 1;
                    else {
                       if(fread(&(edge->probability), sizeof(decltype(edge->probability)), 1, file) != 1)
                            throw IoException("Could not determine probability of node %lu", i);
                    }

                    ++edge;
                }

            }
            fclose(file);
        }

    } else {
        throw IoException("Read option not supported by model type");
    }
}

void Partition::Read(const ModelType kModelType, const unsigned int kPartitionId){
    std::string filename;
    switch(kModelType){
        case ModelType::WPBDD:
            filename = manager.files.get_filename(file_t::WPBDD);
            Read(kModelType,filename);
            break;
        case ModelType::PWPBDD:
            filename = manager.files.get_filename(file_t::PWPBDD,stringf("%u",kPartitionId));
            Read(kModelType,filename);
            break;
        case ModelType::MULTIGRAPH:
            filename = manager.files.get_filename(file_t::MULTIGRAPH);
            Read(kModelType,filename);
            break;
        case ModelType::TDMULTIGRAPH:
            filename = manager.files.get_filename(file_t::TDMULTIGRAPH);
            Read(kModelType,filename);
            break;
        default:
            throw IoException("Read option not supported for this model type");
            break;
    }
}

}

