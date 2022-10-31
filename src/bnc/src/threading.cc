#include "threading.h"
#include <unistd.h>
#include <stdio.h>
#include "misc.h"

timed_thread_t::timed_thread_t(){
    status = ready;
    delay = 0;
    retval = 0;
    f = NULL;
    args = NULL;
}

bool timed_thread_t::is_aborted(){
    return status == aborted;
}

int timed_thread_t::get_returnvalue(){
    return retval;
}

void* timed_thread_t::timer_thread_function(void *v){
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    timed_thread_t *t = ((timed_thread_t*)v);
    pthread_barrier_wait(&(t->time_barrier));

    sleep(t->delay);
    pthread_cancel(t->timed_thread); // ignore pthread error: "no such process"

    return NULL;
}

void* timed_thread_t::thread_function(void *v){
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    timed_thread_t *t = ((timed_thread_t*)v);
    pthread_barrier_wait(&(t->time_barrier));

    t->retval = (long int) (t->f)(t->args);

    return NULL;
}

void timed_thread_t::run(pthread_function_t f, void *args, unsigned int seconds){
    pthread_barrier_init(&time_barrier,NULL,2);

    this->f = f;
    this->args = args;
    status = ready;
    delay = seconds;
    pthread_create(&timed_thread, NULL, thread_function, (void*) this);
    pthread_create(&timer_thread, NULL, timer_thread_function, (void*) this);

    void *res;
    pthread_join(timed_thread,&res);
    if(res == PTHREAD_CANCELED)
        status = aborted;
    else
        pthread_cancel(timer_thread); // ignore pthread error: "no such process"

    pthread_barrier_destroy(&time_barrier);
}

