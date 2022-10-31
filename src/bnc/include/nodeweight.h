struct weights : public std::set<weight_t> {
    inline unsigned int get_sum() const {
        return std::accumulate(this->begin(),this->end(),0);
    }

    typedef std::stack<weights*> stack;
};

