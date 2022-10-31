#include "mutex.h"
#include "misc.h"

robust_cond_mutex::robust_cond_mutex(){
    cond = false;
    pthread_mutexattr_t attr;
    pthread_check(pthread_mutexattr_init(&attr));
#ifndef __APPLE__
    pthread_check(pthread_mutexattr_setrobust_np(&attr, PTHREAD_MUTEX_ROBUST_NP));
    pthread_check(pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE_NP));
#endif
    pthread_check(pthread_mutex_init(&mutex, &attr));
    pthread_check(pthread_cond_init(&mutex_cond, NULL));
    pthread_check(pthread_mutexattr_destroy(&attr));
}

robust_cond_mutex::~robust_cond_mutex(){
    lock();
    condition(true);
    broadcast();
    unlock();

    pthread_check(pthread_mutex_destroy(&mutex));
    pthread_check(pthread_cond_destroy(&mutex_cond));
}

void robust_cond_mutex::reset(){
    lock();
    condition(true);
    pthread_check(pthread_cond_broadcast(&mutex_cond));
    unlock();

    lock();
    cond = false;
    unlock();
}

void robust_cond_mutex::lock(){
    int ret;
    if((ret = pthread_mutex_lock(&mutex))){
        if(ret == EOWNERDEAD){
#ifndef __APPLE__
            if((ret = pthread_mutex_consistent_np(&mutex)) != 0) {
                //log(err, "error '%s' could not get lock on robust mutex\n", strerror(ret));
                //log_stacktrace(err);

                //log(dbg, "error '%s' could not get lock on robust mutex\n", strerror(ret));
            }
#endif
        }
    }
}

void robust_cond_mutex::unlock(){
    pthread_check(pthread_mutex_unlock(&mutex));
}

void robust_cond_mutex::broadcast(){
    pthread_check(pthread_cond_broadcast(&mutex_cond));
}

void robust_cond_mutex::signal(){
    pthread_check(pthread_cond_signal(&mutex_cond));
}

void robust_cond_mutex::wait(){
    while (!cond)
        pthread_cond_wait(&mutex_cond, &mutex);
}

bool robust_cond_mutex::condition(){
    return cond;
}

void robust_cond_mutex::condition(bool c){
    cond = c;
}

