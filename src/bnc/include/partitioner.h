#ifndef PARTITIONER_H
#define PARTITIONER_H

#include "partition.h"
#include <gmpxx.h>
#include <list>

using namespace std;

struct two_partitions {
	partition_t left;
	partition_t right;
};

class partitioner {
    public:
		static void generate_partitions(bayesnet*, ordering_t*, bn_partitions_t&);
        static void generate_partitions_sa(bayesnet*, bn_partitions_t&);
	private:
        static void partition(bn_partition_t& parent, list<bn_partition_t>& todo);
        static void continue_partition(partition_t& part, const ordering_t& parent_ordering, list<bn_partition_t>& todo);

        static mpz_class compute_score(const ordering_t& parent_ordering, partition_t& part);
        static void partition_ordering(partition_t const&, const ordering_t&, ordering_t&);

        static list<bn_partition_t> todo;
        static mpz_class total_score;

        static size_t nvars;
        static bayesnet* bn;
        static partition_set_t cutset_union;

};

#endif
