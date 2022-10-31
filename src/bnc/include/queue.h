#ifndef QUEUE_H
#define QUEUE_H

#include <queue>
#include <pthread.h>
#include "mutex.h"

template <class T> class synchronous_queue {
    public:
        void push_async(T t){
            q.push(t);
        };

        void push(T t){
            mutex.lock();
            q.push(t);
            mutex.condition(true);
            mutex.signal();
            mutex.unlock();
        }

        void clear_and_free(){
            mutex.lock();
            while (!q.empty()){
                delete q.front();
                q.pop();
            }
            mutex.unlock();
        }

        void clear(){
            mutex.lock();
            while (!q.empty())
                q.pop();
            mutex.unlock();
        }

        unsigned int size_async(){
            return q.size();
        }

        unsigned int size(){
            mutex.lock();
            unsigned int ret = q.size();
            mutex.unlock();

            return ret;
        }

        bool empty(){
            mutex.lock();
            bool ret = q.empty();
            mutex.unlock();

            return ret;
        }

        bool empty_async(){
            return q.empty();
        }

        T pop_ptr(){
            mutex.lock();
            T ret;
            if(!q.empty()){
                ret = q.front();
                q.pop();
            } else ret = NULL;
            mutex.unlock();
            return ret;
        }

        T pop_ptr_async(){
            T ret;
            if(!q.empty()){
                ret = q.front();
                q.pop();
            } else ret = NULL;
            return ret;
        }

        T pop(){
            mutex.lock();
            while (q.empty()){
                mutex.condition(false);
                mutex.wait();
            }
            T ret = q.front();
            q.pop();
            mutex.unlock();
            return ret;
        }

        void clear_and_free_async(){
            while (!q.empty()){
                delete q.front();
                q.pop();
            }
        }

        void clear_async(){
            while (!q.empty())
                q.pop();
        }

        T pop_async(){
            T ret = q.front();
            q.pop();
            return ret;
        }

        void lock(){
            mutex.lock();
        }

        void unlock(){
            mutex.unlock();
        }

        void condition(bool c){
            mutex.condition(c);
        }

        void signal(){
            mutex.signal();
        }

        void wait(){
            mutex.wait();
        }
    private:
        robust_cond_mutex mutex;
        std::queue<T> q;
};


#endif
