#define TRUE (-1)
#define FALSE 0

static inline void init(weighted_node *n){
    *n = (const struct weighted_node){ 0, 0, NULL, NULL, NULL };
}

static inline bool is_terminal(const weighted_node *n){
    return n->l <= FALSE;
}

static inline bool is_satisfiable(const weighted_node *n){
    return n->l == TRUE;
}

static inline bool is_unsatisfiable(const weighted_node *n){
    return n->l == FALSE;
}

static inline bool is_negative_cofactor(const weighted_node *parent,const weighted_node *child){
    return parent->e == child;
}

static inline bool is_positive_cofactor(const weighted_node *parent,const weighted_node *child){
    return parent->t == child;
}

static inline bool is_dead(const weighted_node *n){
    return n->ref == 0;
}

static inline weighted_node* reference(weighted_node *n){
    n->ref++;
    return n;
}

static inline weighted_node* negate(const weighted_node *n){
    return ((weighted_node*)((uintptr_t)(n) ^ (uintptr_t) 01));
}

static inline bool is_negated(const weighted_node *n){
    return ((int) ((uintptr_t) (n) & (uintptr_t) 01));
}

static inline weighted_node* regular(const weighted_node *n){
    return ((weighted_node*)((uintptr_t)(n) & ~(uintptr_t) 01));
}

static inline weighted_node* complement(const weighted_node *n){
    return ((weighted_node*)((uintptr_t)(n) | (uintptr_t) 01));
}

static inline bool collapse_candidate(const weighted_node *n){
    if((n->t == n->e->t)
        && ((n->W == NULL && n->e->W == NULL)
        || (n->W != NULL && n->e->W != NULL && *(n->W) == *(n->e->W))))
        return true;
    else return false;
}

static inline void dereference(weighted_node *n){
    #ifdef debug
    if(n->ref <= 0)
        fprintf(stderr, "dereferencing unreference node (ref count: %d)\n", n->ref);
    #endif
    n->ref--;
}

