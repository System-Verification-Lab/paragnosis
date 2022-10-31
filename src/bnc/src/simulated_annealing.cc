/* siman/siman_tsp.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Mark Galassi
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

//#include <config.h>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <thread>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>
#include "ordering.h"
#include "bound.h"
#include "options.h"
#include "bnc.h"
#include "partitioner.h"
#include "string.h"
#include <string>
#include "composition.h"
#include "spanning.h"
#include <cstring>
#include <cassert>

using namespace bnc;
/* set up parameters for this simulated annealing run */

#define STEP_SIZE 1.0           /* max step size in random walk */
#define BOLTZMANN 1.0           /* Boltzmann constant */

#ifndef UNUSED
    #define UNUSED __attribute__((unused))
#endif

size_t SA_VARIABLES;
Spanning spanning;

Bound sa_bound;
double sa_ordering_energy(void *xp);
void sa_ordering_step(const gsl_rng * r, void *xp, double step_size);
void sa_ordering_print(void *xp);

template <unsigned int N>
double sa_ordering_energy(void *xp);

template <>
double sa_ordering_energy<2>(void *xp){
    ordering_t::value_type *ordering = (ordering_t::value_type*) xp;
    return sa_bound.ComputeTree<Bound::kNodes, double>(ordering);
}

template <>
double sa_ordering_energy<1>(void *xp){
    ordering_t::value_type *ordering = (ordering_t::value_type*) xp;
    return sa_bound.ComputeTree<Bound::kNodes, Bound::BoundSize>(ordering);
}

template <>
double sa_ordering_energy<0>(void *xp){
    ordering_t::value_type *ordering = (ordering_t::value_type*) xp;
    return sa_bound.Compute<Bound::kRatio,double>(ordering);
}

void sa_ordering_print(void *xp){
    ordering_t::value_type *ordering = (ordering_t::value_type*) xp;

    if(SA_VARIABLES < 10){
        std::string str;
        str = " [";
        for(unsigned int i = 0; i < SA_VARIABLES; i++){
            if(i > 0)
                str += ", ";
            str += std::to_string(ordering[i]);
        }
        str += "]";

        printf("%s", str.c_str());
    } else
        printf(" - ");
}

/* take a step through the TSP space */
void sa_ordering_step(const gsl_rng * r, void *xp, UNUSED double step_size){
    ordering_t::value_type *ordering = (ordering_t::value_type*) xp;

    // /* pick the two variables to swap
    unsigned int x1,x2;
    x1 = gsl_rng_get (r) % SA_VARIABLES;
    do {
        x2 = gsl_rng_get (r) % SA_VARIABLES;
    } while (x2 == x1);

    literal_t l = ordering[x1];
    ordering[x1] = ordering[x2];
    ordering[x2] = l;
}

typedef double (*energy_function_t)(void*);
struct anneal_xp {
    ordering_t ordering;
    unsigned int thread_id;
    energy_function_t energy_f;
};

void* anneal(void *xp){
    anneal_xp *data = (anneal_xp*) xp;
    ordering_t::value_type *basic_ordering = &(data->ordering[0]);
    gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup()) ;
    gsl_ieee_env_setup ();

    gsl_siman_params_t params = {bnc::OPT_SA_TRIES, bnc::OPT_SA_ITERATIONS,
        STEP_SIZE, BOLTZMANN, bnc::OPT_SA_TEMPERATURE_INITIAL,
        bnc::OPT_SA_TEMPERATURE_DAMP_FACTOR, bnc::OPT_SA_TEMPERATURE_MIN};

    double score_original = (data->energy_f)((void*) basic_ordering);
    if(data->thread_id == 0 && bnc::OPT_SA_PRINT_ORDERING){
        gsl_siman_solve(r, basic_ordering, data->energy_f, sa_ordering_step,
            NULL, sa_ordering_print, NULL, NULL, NULL,
            SA_VARIABLES*sizeof(ordering_t::value_type), params);
    } else {
        gsl_siman_solve(r, basic_ordering, data->energy_f, sa_ordering_step,
            NULL, NULL, NULL, NULL, NULL,
            SA_VARIABLES*sizeof(ordering_t::value_type), params);
    }
    double score_new = (data->energy_f)((void*) basic_ordering);

    if(score_original > score_new)
        printf("Anneal score %lf --> %lf (improved)\n", score_original, score_new);
    else
        printf("Anneal score %lf --> %lf\n", score_original, score_new);

    gsl_rng_free(r);

    return 0;
}

int ordering_t::anneal(bayesnet *bn, const partition_t &kPartition, ordering_t &best_ordering, unsigned int score_function){
    // precondition:
    // 0 ordering contains partition and cutset variables
    // 1 partition is defined
    // 2 cutset is defined

    // bad initial ordering ;)
    if(best_ordering.size() == 0){
        for(auto it = kPartition.cutset.begin(); it != kPartition.cutset.end(); it++)
            best_ordering.push_back(*it);
        for(auto it = kPartition.set.begin(); it != kPartition.set.end(); it++)
            best_ordering.push_back(*it);
    }

    if(best_ordering.size() < 2)
        return 0;

    SA_VARIABLES = best_ordering.size();
    sa_bound.SetBayesnet(bn);
    sa_bound.Init(kPartition);

    energy_function_t energy_f;
    switch(score_function){
        case 0: // chain score
            energy_f = ::sa_ordering_energy<0>;
            break;
        case 1: // tree upperbound on nodes
            energy_f = ::sa_ordering_energy<1>;
            break;
        case 2: // tree score
            energy_f = ::sa_ordering_energy<2>;
            break;
        default:
            fprintf(stderr, "Unknown scoring function\n");
            exit(1);
    }

    const unsigned int kNrThreads = std::thread::hardware_concurrency();
    assert(kNrThreads > 0);

    std::string filename = files.get_filename(VARIABLE_ORDERING,"intermediate");

    std::vector<anneal_xp> xp(kNrThreads);

    double best_energy = std::numeric_limits<double>::max();
    unsigned int stable_energy = 0;
    for(unsigned int run = 0; OPT_SA_RUNS == 0 || run < OPT_SA_RUNS; run++){
        if(OPT_SA_PARALLELISM){
            // boot threads
            pthread_t threads[kNrThreads];
            for(unsigned int i = 0; i < kNrThreads; i++){
                xp[i].ordering = best_ordering;
                xp[i].thread_id = i;
                xp[i].energy_f = energy_f;
                if(pthread_create(&(threads[i]), NULL, ::anneal, (void*)&(xp[i]))){
                    fprintf(stderr, "Error creating thread %u of %u\n", i+1, kNrThreads);
                    exit(1);
                }
            }

            // wait until finished
            for(unsigned int i = 0; i < kNrThreads; i++){
                if(pthread_join(threads[i],NULL)){
                    fprintf(stderr, "Error canceling thread %u of %u\n", i+1, kNrThreads);
                    exit(1);
                }
            }

            // select best
            double best_score = std::numeric_limits<double>::max();
            unsigned int best_thread;
            for(unsigned int i = 0; i < kNrThreads; i++){
                auto thread_score = energy_f((void*) &(xp[i].ordering[0]));
                if(thread_score < best_score){
                    best_score = thread_score;
                    best_thread = i;
                }
            }
            best_ordering = xp[best_thread].ordering;

        } else {
            xp[0].thread_id = 0;
            xp[0].energy_f = energy_f;
            xp[0].ordering = best_ordering;
            ::anneal((void*) &(xp[0]));
            best_ordering = xp[0].ordering;
        }

        if(OPT_SA_RUNS != 1){
            best_ordering.write(filename);
            printf("Written intermediate ordering to: %s\n", filename.c_str());
        }

        if(OPT_SA_RUNS == 0){
            double current_energy = energy_f((void*) &(best_ordering[0]));
            if(current_energy < best_energy){
                best_energy = current_energy;
                stable_energy = 1;
            } else if (current_energy == best_energy)
                ++stable_energy;

            printf("Stability: %lu/3\n", stable_energy);
            if(stable_energy >= 3)
                break;
        }
    }

    return 0;
}


/************************************
** partitioning
************************************/
bayesnet *gBn;
bn_partitions_t partitions;
size_t gNrVariables;
size_t gNrPartitions;
size_t gNrRemoveEdges;
const uint32_t *gDimensions;
size_t iterations;
double sa_partition_energy(void *xp);
void sa_partition_step(const gsl_rng * r, void *xp, double step_size);
void sa_partition_print(void *xp);

double sa_max_cutset_score;

/* energy for the travelling salesman problem */
double sa_partition_score_balance(const unsigned int *kPartitionSize){
    unsigned int max_partition_size = 0;
    const unsigned int kNrPartitions = gNrPartitions;
    const unsigned int kNrVariables = gNrVariables;
    for(unsigned int i = 0; i < kNrPartitions; i++)
        if(kPartitionSize[i] > max_partition_size)
            max_partition_size = kPartitionSize[i];

    return max_partition_size / (kNrVariables / (double) kNrPartitions);
}

double sa_partition_score_normalized_stddev_size(const unsigned int *kPartitionSize){
    const unsigned int kNrPartitions = gNrPartitions;
    const unsigned int kNrVariables = gNrVariables;
    double score = 0;
    for(unsigned int i = 0; i < kNrPartitions; i++)
        score += std::pow(kPartitionSize[i] / (kNrVariables/(double)kNrPartitions),2);

    return std::sqrt (score) * (1/(double)kNrVariables);
}


double sa_partition_size_score(const unsigned int *kPartitionSize){
    double score;
    #define SCORE1
    #ifdef SCORE1
    score = sa_partition_score_normalized_stddev_size(kPartitionSize);
    #else
    score = sa_partition_score_balance(kPartitionSize);
    #endif

    return score;

}

double sa_partition_cutset_score(const unsigned int *kPartitionMap, const unsigned int kType = 0){
    double score;
    size_t scores[gNrPartitions];
    std::fill(scores,scores+gNrPartitions,1);
    for(unsigned int variable = 0; variable < gNrVariables; variable++){
        const auto kPartition = kPartitionMap[variable];

        const auto *kParentPtr = gBn->get_parent(variable);
        const auto kEnd = kParentPtr + gBn->get_parent_size(variable);
        while(kParentPtr != kEnd){
            const auto kParentPartition = kPartitionMap[*kParentPtr];
            if(kParentPartition != kPartition){
                if(kType == 0)
                    scores[kParentPartition] *= gDimensions[*kParentPtr];
                else if(kType == 1)
                    scores[kParentPartition] += gDimensions[*kParentPtr];
                else if(kType == 2)
                    scores[kParentPartition] += 1;
            }

            ++kParentPtr;
        }
    }

    score=0;
    for(unsigned int p = 0; p < gNrPartitions; p++)
        if(scores[p] > 1)
            score += scores[p];

    if(score > sa_max_cutset_score)
        sa_max_cutset_score = score;

    return score;
}

inline bool *CastToEdges(void * const kInput, const unsigned int kNrRemoveEdges){
    return (bool*) &(((unsigned int*) kInput)[kNrRemoveEdges]);
}

inline unsigned int *CastToEdgeList(void * const kInput){
    return (unsigned int*) kInput;
}

inline const bool *CastToEdgesR(const void * const kInput, const unsigned int kNrRemoveEdges){
    return (bool*) &(((unsigned int*) kInput)[kNrRemoveEdges]);
}

inline const unsigned int *CastToEdgeListR(const void * const kInput){
    return (unsigned int*) kInput;
}

double sa_partition_energy(void *input){
    unsigned int partition_size[gNrPartitions];
    unsigned int partition_map[gNrVariables];
    const bool *kDeletedEdges = CastToEdges(input, gNrRemoveEdges);
    spanning.ComputePartitionMap(kDeletedEdges, partition_map, partition_size, gNrPartitions);

    double score;
    if(bnc::OPT_SA_SCORE_RATIO == 0){
        //score = sa_partition_cutset_score(partition_map);
        score = sa_partition_cutset_score(partition_map);
    } else if(bnc::OPT_SA_SCORE_RATIO == 1){
        score = sa_partition_size_score(partition_size);
    } else {
        const double kAlphaCutset = bnc::OPT_SA_SCORE_RATIO;
        const double kScoreCutset = sa_partition_cutset_score(partition_map);
        if(kScoreCutset != 0)
            score += kAlphaCutset * (1 / kScoreCutset);

        const double kAlphaBalance = 1 - bnc::OPT_SA_SCORE_RATIO;
        const double kScoreBalance = sa_partition_size_score(partition_size);
        if(kScoreBalance != 0)
            score += kAlphaBalance * (1 / (kScoreBalance*sa_max_cutset_score));

        assert(kAlphaCutset + kAlphaBalance == 1);
        score = 1 / score;
    }

    return score;
}

/* take a step through the TSP space */
void sa_partition_step(const gsl_rng * r, void *input, UNUSED double step_size){
    const unsigned int kNrRemoveEdges = gNrRemoveEdges;
    const unsigned int kNrEdges = spanning.GetSize();
    unsigned int *edge_list = CastToEdgeList(input);
    bool *deleted_edges = CastToEdges(input, kNrRemoveEdges);

    unsigned int edge_to = gsl_rng_get (r) % kNrEdges;
    while(deleted_edges[edge_to])
        edge_to = (edge_to+1) % kNrEdges;

    const unsigned int kEdgeListId = gsl_rng_get (r) % kNrRemoveEdges;
    #ifdef DEBUG
    assert(deleted_edges[edge_list[kEdgeListId]] && !deleted_edges[edge_to]);
    #endif
    deleted_edges[edge_list[kEdgeListId]] = false;
    deleted_edges[edge_to] = true;
    edge_list[kEdgeListId] = edge_to;
}

void partitioner::generate_partitions_sa(bayesnet* bn, bn_partitions_t& result) {
    spanning.SetBayesnet(bn);
    spanning.CreateSpanningTree();
    const unsigned int kNrComponents = spanning.ComputeNrComponents();
    if(kNrComponents > 1)
        fprintf(stderr, "Note: network already consists of %u components\n", kNrComponents);

    // init globals
    const size_t kNrVariables = bn->get_nr_variables();
    gNrVariables = kNrVariables;
    gDimensions = bn->get_states();
    gBn = bn;
    sa_max_cutset_score = 0;
    iterations = 0;
    assert(bnc::OPT_NR_PARTITIONS <= (int)kNrVariables);
    gNrPartitions = 2;
    if(bnc::OPT_NR_PARTITIONS > 1)
        gNrPartitions = bnc::OPT_NR_PARTITIONS;

    // safety check
    if(gNrPartitions < kNrComponents){
        fprintf(stderr, "Partitioning a network into less components than it already consists of is not supported yet (%u partitions, %u components)\n", gNrPartitions, kNrComponents);
        exit(1);
    }

    // init
    const int kNrRemoveEdges = gNrPartitions - kNrComponents;
    const unsigned int kNrEdges = spanning.GetEdges().size();
    const size_t kNrEdgesBytes = kNrEdges*sizeof(bool);
    const size_t kNrRemoveEdgesBytes = kNrRemoveEdges * sizeof(unsigned int);
    const size_t kTotalBytes = kNrRemoveEdgesBytes + kNrEdgesBytes;
    unsigned int input[static_cast<size_t>(std::ceil(kTotalBytes / (double)sizeof(unsigned int)))];
    std::memset(input,0,kTotalBytes);




    if(kNrRemoveEdges > 0){
        gNrRemoveEdges = kNrRemoveEdges;
        unsigned int *edge_list = CastToEdgeList(input);
        bool *deleted_edges = CastToEdges(input, kNrRemoveEdges);
        const unsigned int kComponentSize = kNrEdges/(kNrRemoveEdges+1);
        for(unsigned int i = 0; i < kNrRemoveEdges; i++){
            unsigned int index = (i+1) * kComponentSize;
            edge_list[i] = index;
            deleted_edges[index] = true;
        }

        if(gNrPartitions < gNrVariables){
            // perform simulated annealing
            gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup()) ;
            gsl_ieee_env_setup ();
            gsl_siman_params_t params = {bnc::OPT_SA_TRIES, bnc::OPT_SA_ITERATIONS,
                STEP_SIZE, BOLTZMANN, bnc::OPT_SA_TEMPERATURE_INITIAL,
                bnc::OPT_SA_TEMPERATURE_DAMP_FACTOR, bnc::OPT_SA_TEMPERATURE_MIN};

            gsl_siman_solve(r,
                (void*) &(input[0]),
                sa_partition_energy,
                sa_partition_step,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL,
                kTotalBytes, params);

            gsl_rng_free(r);
        }
    }

    // write result
    unsigned int partition_size[gNrPartitions];
    unsigned int partition_map[kNrVariables];
    const bool *kDeletedEdges = CastToEdges(input, kNrRemoveEdges);
    const unsigned int nr_partitions = spanning.ComputePartitionMap(kDeletedEdges, partition_map, partition_size, gNrPartitions);
    assert(nr_partitions == gNrPartitions);

    result.clear();
    result.resize(gNrPartitions);
    for(unsigned int variable = 0; variable < kNrVariables; variable++){
        const auto kPartitionId = partition_map[variable];
        result[kPartitionId].partition.set.insert(variable);

        const auto *kChildren = bn->get_child(variable);
        const auto kNrChildren = bn->get_child_size(variable);
        for(unsigned int i = 0; i < kNrChildren; i++){
            const unsigned int kChild = kChildren[i];
            const auto kChildPartitionId = partition_map[kChild];
            if(kPartitionId != kChildPartitionId)
                result[kChildPartitionId].partition.cutset.insert(variable);
        }
    }


    const double kScoreCutset = sa_partition_cutset_score(partition_map);
    const double kScoreBalance = sa_partition_size_score(partition_size);
    const double kBestScore = Composition::ComputeBestScore(bn, result);

    printf("%lu ---------------------------\n",++iterations);
    Composition comp;
    comp.SetBayesnet(bn);
    comp.SetPartitions(result);
    ordering_t best_ordering = comp.FindOrdering();
    comp.BuildOrdering(best_ordering);
    printf("Best composition ordering:\n");
    comp.PrintAscii(bn);

    printf("\n");
    printf("sizes             : ");
    for(unsigned int p = 0; p < gNrPartitions; p++){
        double cutset_cardi = 1;
        for(auto it = result[p].partition.cutset.begin(); it != result[p].partition.cutset.end(); it++)
            cutset_cardi *= gDimensions[*it];
        printf("%lu {%lu:%.0lf} ",partition_size[p], result[p].partition.cutset.size(),cutset_cardi);

    }

    printf("\n");
    printf("stddev score      : %lf\n", sa_partition_score_normalized_stddev_size(partition_size));
    printf("balance score     : %lf\n", sa_partition_score_balance(partition_size));
    printf("composition score : %lf\n", comp.ComputeScore(best_ordering));
    printf("cutset score 1    : %lf\n", sa_partition_cutset_score(partition_map,0));
    printf("cutset score 2    : %lf\n", sa_partition_cutset_score(partition_map,1));
    printf("cutset score 3    : %lf\n", sa_partition_cutset_score(partition_map,2));
}

