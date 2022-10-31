#include "partitioner.h"
#include <set>
#include "bnc.h"
#include "options.h"

using namespace std;

// don't try to further split partitions with size bound <= this threshold
mpz_class SCORE_THRESHOLD = 6;
size_t NR_PARTITIONS_THRESHOLD = 9;//900000000;;
bool SPLIT_INTO_CONNECTED_COMPS = false;

void partitioner::generate_partitions(bayesnet* bnet, ordering_t* ordering, bn_partitions_t& partitions) {
    // bnc::OPT_NR_PARTITIONS == -1, iff the option is not set and > than 1 otherwise
    if(bnc::OPT_NR_PARTITIONS > 1){
        SCORE_THRESHOLD = 0;
        NR_PARTITIONS_THRESHOLD = bnc::OPT_NR_PARTITIONS;
        SPLIT_INTO_CONNECTED_COMPS = false;
    }

	nvars = bnet->get_nr_variables();
    bn = bnet;

    partition_t everything;
    for (unsigned int i = 0; i < nvars; i++)
        everything.set.insert(i);
    ordering_t init_ordering = *ordering;
    total_score = bnc::get_upper_bound(bn, everything, init_ordering);

    bn_partition_t bn_part;
    bn_part.partition = everything;
    bn_part.variable_ordering = init_ordering;
    todo.push_back(bn_part);

    while ( (SCORE_THRESHOLD == 0         || total_score > SCORE_THRESHOLD) &&
            (NR_PARTITIONS_THRESHOLD == 0 || todo.size() < NR_PARTITIONS_THRESHOLD) ) {
        printf("current total score: %19s\n", total_score.get_str().c_str());

        mpz_class worst_score = 0;
        list<bn_partition_t>::iterator next_part_iter;

        list<bn_partition_t>::iterator iter = todo.begin();
        for (; iter != todo.end(); ++iter) {
            bn_partition_t& part = *iter;
            const mpz_class score = bnc::get_upper_bound(bn, part.partition, part.variable_ordering);

            if (score > worst_score && part.partition.set.size() > 1) {
                worst_score = score;
                next_part_iter = iter;
            }
        }

        bn_partition_t next_part = *next_part_iter;
        todo.erase(next_part_iter);
        total_score -= worst_score;
	    partition(next_part, todo);
    }

    // add partitions
    printf("partition bounds\n");
    for (bn_partition_t part: todo) {
        partitions.push_back(part);

        for (const variable_t v: part.partition.cutset)
            cutset_union.insert(v);
        mpz_class score = bnc::get_upper_bound(bn, part.partition, part.variable_ordering);
        printf("               %19s\n", score.get_str().c_str());
    }

    printf("nr partitions: %19lu\n", partitions.size());
    printf("#BN vars:      %19lu\n", nvars);
    printf("#cutset union: %19lu\n", cutset_union.size());
    printf("total score:   %19s\n", total_score.get_str().c_str());
    printf("orig network:  %19s\n", bnc::get_upper_bound(bn, init_ordering).get_str().c_str());
}

void partitioner::partition(bn_partition_t& parent, list<bn_partition_t>& todo) {
    const mpz_class score_parent = bnc::get_upper_bound(bn, parent.partition, parent.variable_ordering);

    // determine next split
    two_partitions best_parts;
    mpz_class best_score;
    bool best_found = false;
    for (variable_t start_var = 0; start_var < nvars; start_var++) {
        two_partitions parts;
        parts.right = parent.partition;

        queue<variable_t> queue;
        queue.push(start_var);

        while (!queue.empty()) {
            const variable_t next = queue.front();
            queue.pop();
            if (parent.partition.set.find(next) == parent.partition.set.end()) continue; // not in current parent part
            if (parts.left.set.find(next) != parts.left.set.end())             continue; // already in left part

            parts.left.set.insert(next);
            parts.right.set.erase(next);

            const mpz_class score_left  = compute_score(parent.variable_ordering, parts.left);
            const mpz_class score_right = compute_score(parent.variable_ordering, parts.right);
            const mpz_class score = score_left + score_right;

            if ((!best_found || score < best_score) && parts.right.set.size() > 0) {
                best_found = true;
                best_parts = parts;
                best_score = score;
            }

            for(unsigned int i = 0; i < bn->get_parent_size(next); i++) {
                const variable_t p = bn->get_parent(next)[i];
                queue.push(p);
            }
            for(unsigned int i = 0; i < bn->get_child_size(next); i++) {
                const variable_t c = bn->get_child(next)[i];
                for(unsigned int j = 0; j < bn->get_parent_size(c); j++)
                    queue.push(bn->get_parent(c)[j]);
                queue.push(c);
            }
        }
    }

    // check if improvement was made in this step
    if (score_parent < best_score) {
        fprintf(stderr, "Warning: Could not reduce score below\n");
        // exit(-1);
    }

    // split into connected components in moralised graph
    continue_partition(best_parts.left, parent.variable_ordering, todo);

    if (SPLIT_INTO_CONNECTED_COMPS) {
        partition_set_t& right_part = best_parts.right.set;
        while (!right_part.empty()) {
            partition_t part;
            stack<variable_t> vars_todo;
            vars_todo.push(*right_part.begin());

            while (!vars_todo.empty()) {
                const variable_t v = vars_todo.top();
                vars_todo.pop();
                if (right_part.find(v) == right_part.end()) continue;
                if (part.set.find(v) != part.set.end())     continue;

                part.set.insert(v);
                right_part.erase(v);
                for(unsigned int i = 0; i < bn->get_parent_size(v); i++)
                    vars_todo.push(bn->get_parent(v)[i]);
                for(unsigned int i = 0; i < bn->get_child_size(v); i++) {
                    const variable_t c = bn->get_child(v)[i];
                    for(unsigned int j = 0; j < bn->get_parent_size(c); j++)
                        vars_todo.push(bn->get_parent(c)[j]);
                    vars_todo.push(c);
                }
            }

            continue_partition(part, parent.variable_ordering, todo);
        }
    } else {
        continue_partition(best_parts.right, parent.variable_ordering, todo);
    }
}

void partitioner::continue_partition(partition_t& part, const ordering_t& parent_ordering, list<bn_partition_t>& todo) {
    ordering_t ordering;
    part.determine_cutset(bn);
    partition_ordering(part, parent_ordering, ordering);
    ordering.generate_variable_ordering_simulated_annealing(bn, part);

    const mpz_class score = bnc::get_upper_bound(bn, part, ordering);
    total_score += score;

    bn_partition_t bn_part;
    bn_part.partition = part;
    bn_part.variable_ordering = ordering;
    todo.push_back(bn_part);
}

mpz_class partitioner::compute_score(const ordering_t& parent_ordering, partition_t& part) {
    ordering_t ordering;
    part.determine_cutset(bn);
    partition_ordering(part, parent_ordering, ordering);

    return bnc::get_upper_bound(bn, part, ordering);
}

void partitioner::partition_ordering(partition_t const& partition, const ordering_t& parent_ordering, ordering_t& ordering) {
	const partition_set_t& set = partition.set;
	const partition_set_t& cutset = partition.cutset;

    for (auto const l: parent_ordering)
        if (set.find(l) != set.end() || cutset.find(l) != cutset.end())
            ordering.push_back(l);
}

list<bn_partition_t> partitioner::todo;
mpz_class partitioner::total_score;

size_t partitioner::nvars;
partition_set_t partitioner::cutset_union;
bayesnet* partitioner::bn;

