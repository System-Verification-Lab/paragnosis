#include "debug.h"
#include "options.h"
#include <queue>
#include <map>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>


namespace bnmc {

void DotToPdf(std::string filename){
    system(stringf("dot -Tpdf %s -o %s.pdf &",filename.c_str(),filename.c_str()).c_str());
}

inline const literal_t GetBaseLiteral(const variable_t variable){
    return manager.mapping.get_variable_to_literal()[variable];
}

inline const literal_t GetLiteral(const variable_t variable, const unsigned int value){
    const unsigned int offset = GetBaseLiteral(variable);
    return offset + value;
}

void WriteDot(std::string filename, const Partition::ArithmeticCircuit &kCircuit, const EvidenceList &kEvidenceList, std::vector<probability_t> &probabilities, const std::vector<unsigned int>&identity, const unsigned int kTier, const int kCurrentIndex, const int kMaxDepth, mstack::StackTop top, mstack::Stack *stack){
    FILE *file = fopen(filename.c_str(), "w");
    if(file){
        const unsigned int kRootIndex = 2;
        const unsigned int kTrueTerminalIndex = 1;
        const unsigned int kFalseTerminalIndex = 0;
        const std::vector<unsigned int>& kVariableToLiteral = manager.mapping.get_variable_to_literal();

        std::queue<unsigned int> q;
        std::map<unsigned int,bool> done;
        fprintf(file, "digraph %s {\n", "WPBDD");
        fprintf(file, "    center = true;\n");
        fprintf(file, "    edge [dir = none];\n");

        // identify which nodes are on which levels
        fprintf(file, "\n    // ==================== EDGES ====================\n\n");
        std::map< unsigned int, std::vector<unsigned int> > level;
        std::vector<int> depth(kCircuit.size(),0);
        depth[kRootIndex] = 1;
        q.push(kRootIndex);
        while(!q.empty()){
            const unsigned int &kIndex = q.front();
            const Partition::Node &kNode = kCircuit[kIndex];
            q.pop();

            if(!done[kIndex]){
                done[kIndex] = true;
                if(kIndex == kTrueTerminalIndex || kIndex == kFalseTerminalIndex)
                    level[0].push_back(kIndex);
                else if(kMaxDepth < 0 || depth[kIndex] <= kMaxDepth){
                    literal_t literal = GetLiteral(kNode.v,kNode.i);
                    level[literal].push_back(kIndex);
                    if(kMaxDepth < 0 || depth[kIndex] < kMaxDepth){
                        if(top > 0 && mstack::Contains(*stack,kIndex,kNode.t,top)){
                            if(kNode.w == 1)
                                fprintf(file, "    %u -> %u [color=green];\n", kIndex, kNode.t);
                            else
                                fprintf(file, "    %u -> %u [color=green,fontsize=8,label=\"%.2f\"];\n", kIndex, kNode.t, kNode.w);
                        } else {
                            if(kNode.w == 1)
                                fprintf(file, "    %u -> %u;\n", kIndex, kNode.t);
                            else
                                fprintf(file, "    %u -> %u [fontsize=8,label=\"%.2f\"];\n", kIndex, kNode.t, kNode.w);
                        }

                        if(top > 0 && mstack::Contains(*stack,kIndex,kNode.e,top))
                            fprintf(file, "    %u -> %u [color=green,style=dotted];\n", kIndex, kNode.e);
                        else
                            fprintf(file, "    %u -> %u [style=dotted];\n", kIndex, kNode.e);

                        q.push(kNode.t);
                        q.push(kNode.e);
                        depth[kNode.t] = depth[kIndex] + 1;
                        depth[kNode.e] = depth[kIndex];
                    }
                }
            }
        }

        fprintf(file, "\n    // ==================== NODES ====================\n\n");
        for(auto it = level.begin(); it != level.end(); it++){
            if(it->second.size() > 0){
                if(it->first == 0){
                    fprintf(file, "    {   // TERMINALS\n");
                    fprintf(file, "        rank = same;\n");
                    for(auto itn = it->second.begin(); itn != it->second.end(); itn++){
                        const unsigned int kIndex = *itn;
                        if(kIndex == kTrueTerminalIndex){
                            fprintf(file, "        %lu [shape=box,label=\"%u\\nT\"", kIndex,kIndex);
                            if(kCurrentIndex >= 0 && kCurrentIndex == kIndex)
                                fprintf(file, ",style=filled,fillcolor=\"grey\"");
                            fprintf(file, "];\n");
                        } else if(kIndex == kFalseTerminalIndex){
                            fprintf(file, "        %lu [shape=box,label=\"%u\\nF\"", kIndex,kIndex);
                            if(kCurrentIndex >= 0 && kCurrentIndex == kIndex)
                                fprintf(file, ",style=filled,fillcolor=\"grey\"");
                            fprintf(file, "];\n");
                        }
                    }
                    fprintf(file, "    }\n");
                } else {
                    const unsigned int &kIndex = it->second[0];

                    if(kMaxDepth > 0 && depth[kIndex] > kMaxDepth)
                        continue;

                    literal_t l = it->first;
                    const Partition::Node &kNode = kCircuit[kIndex];
                    std::string assignment = Evidence::GetAssignmentString(kNode.v,kNode.i);

                    fprintf(file, "    {   // literal %d (%s)\n", l, assignment.c_str());
                    fprintf(file, "        rank = same;\n");
                    for(auto itn = it->second.begin(); itn != it->second.end(); itn++){
                        const unsigned int &kIndex = *itn;
                        fprintf(file, "        %u", kIndex);
                        fprintf(file, " [label=\"%u\\n%u:%.3f\\n%s\\n%u=%u\"", kIndex, identity[kIndex], probabilities[kIndex],assignment.c_str(), kNode.v,kNode.i);


                        if(kCurrentIndex >= 0 && kCurrentIndex == kIndex)
                            fprintf(file, ",fillcolor=\"grey\"");

                        //const unsigned int &kVariable = kNode.v;
                        //const bool kIsTrue = kEvidenceList[kVariable] == kNode.i;

                        //if(kIsPersistent)
                        //    fprintf(file, ",fontcolor=\"blue\"");

                        //if(kIsConditioned){
                        //    if(kIsTrue)
                        //        fprintf(file, ",color=\"green\"");
                        //    else
                        //        fprintf(file, ",color=\"red\"");
                        //}
                        //std::string style;
                        //if(kCurrentIndex >= 0 && kCurrentIndex == kIndex){
                        //    if(kIsConditioned && kVariable.tier > 0)
                        //        style = "style=\"dashed,filled\"";
                        //    else style = "style=\"filled\"";
                        //} else {
                        //    if(kIsConditioned && kVariable.tier > 0)
                        //        style = "style=\"dashed\"";
                        //}
                        //if(!style.empty())
                        //    fprintf(file, ",%s",style.c_str());

                        fprintf(file, "];\n");
                    }
                    fprintf(file, "    }\n");
                }
            }
        }
        fprintf(file, "}\n");
        fclose(file);
        DotToPdf(filename);
    } else fprintf(stderr, "could not open dot file\n");
}

}
