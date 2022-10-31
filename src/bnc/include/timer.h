#ifndef TIMER_H
#define TIMER_H

#include <chrono>

class Timer {
    public:
        typedef std::chrono::high_resolution_clock Clock; // alt. std::chrono::system_clock
        typedef std::chrono::time_point<Clock> TimePoint;
        typedef std::chrono::microseconds Microseconds;
        typedef std::chrono::nanoseconds  Nanoseconds;
        typedef std::chrono::milliseconds Milliseconds;
        typedef std::chrono::seconds      Seconds;
        typedef std::chrono::minutes      Minutes;
        typedef std::chrono::hours        Hours;
        typedef Nanoseconds Precision;

    private:
        Precision total_;
        TimePoint  start_, end_;

    public:
        Timer() { Reset(); };
        inline void Reset() { total_ = Precision::zero(); }
        inline void Start() { start_ = Clock::now(); }
        inline void Stop()  { end_   = Clock::now(); }

        inline void Add() {
            total_ += std::chrono::duration_cast<Precision>(end_-start_);
        }

        template <class T>
            inline void Add(Precision::rep t) {
                total_ += std::chrono::duration<Precision::rep,typename T::period>(t);
            }

        template <class T>
            inline double GetDuration() const {
                return std::chrono::duration<double,typename T::period>(end_-start_).count();
            }

        template <class T>
            inline double GetTotal() const {
                return std::chrono::duration<double,typename T::period>(total_).count();
            }
};

#endif
