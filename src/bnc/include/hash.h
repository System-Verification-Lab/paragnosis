#ifndef BNMC_HASH_H
#define BNMC_HASH_H

#include <cstdint>
#include <functional>


// based on FNV-1a for 64 bits
typedef std::uint64_t Hash;

class Hasher {
    public:
        inline void Seed(Hash seed){
            // or: hash = std::hash<T>()(seed);
            hash = seed;
        }

        template <class T>
        inline void AddHash(T value){
            // or: hash = (hash ^ value) * kPrime;
            hash = (hash ^ std::hash<T>()(value)) * kPrime;
        }

        inline Hash GetHash(){
            // or: return hash ^ (hash>>32);
            return hash;
        }
    private:
        Hash hash;
        const Hash kPrime = 1099511628211;
};

#endif

