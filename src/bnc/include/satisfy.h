#ifndef BNC_SATISFY_H
#define BNC_SATISFY_H

enum satisfy_t {satisfiable = 0, unsatisfied = 1, redundant = 3, trivial = 7, unsatisfiable = 15, cached = 31};

inline satisfy_t operator|(satisfy_t a, satisfy_t b){
    return static_cast<satisfy_t>(static_cast<int>(a) | static_cast<int>(b));
}

inline satisfy_t& operator|=(satisfy_t &a, satisfy_t b){
    return a = a | b;
}

#endif
