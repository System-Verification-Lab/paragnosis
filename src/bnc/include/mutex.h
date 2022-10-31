#ifndef MUTEX_H
#define MUTEX_H

#include <pthread.h>

class robust_cond_mutex {
    public:
        robust_cond_mutex();
        ~robust_cond_mutex();

        void reset();
        void lock();
        void unlock();
        void wait();
        void broadcast();
        void signal(); // NOTE a signal() HAS TO be followed by a wait()
        bool condition();
        void condition(bool);

    private:
        pthread_mutex_t mutex;
        pthread_cond_t mutex_cond;
        bool cond;
};


#endif
