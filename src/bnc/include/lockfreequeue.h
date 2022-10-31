#ifndef LOCKFREEQUEUE_H
#define LOCKFREEQUEUE_H

#include <array>
#include <atomic>
#include <cassert>
#define SYNC std::memory_order_relaxed

template <typename T, std::uint32_t SIZE>
class LockFreeQueue {
    typedef std::uint32_t Position;

    public:
        //struct Node {
        //    T data;
        //    Position next;
        //} __attribute__ ((aligned (8)));

        LockFreeQueue(){
            std::atomic<Position> test;
            assert(test.is_lock_free());

            read_position_ = 0;
            write_position_ = 0;
            //queue_[SIZE-1].next = 0;
            //for(unsigned int i = 0; i < SIZE -1; i++)
            //    queue_[i].next = i+1;

            atomic_thread_fence(std::memory_order_release);
        }

        virtual ~LockFreeQueue(){

        }

        inline bool empty() const noexcept {
            return read_position_.load(SYNC) == write_position_.load(SYNC);
        }

        inline static constexpr Position Increment(Position pos) noexcept {
            return ((pos + 1 == SIZE) ? 0 : pos + 1);
        }

        inline bool push(const T& element){
            const Position kPos = write_position_.load(SYNC);
            const Position kNewPos = Increment(kPos);
            if (kNewPos != read_position_.load(SYNC)){
                queue_[kPos] = element;
                atomic_thread_fence(std::memory_order_release);
                write_position_.store(kNewPos,SYNC);
                return true;
            }
            return false;
        }

        inline bool pop(T& element){
            const Position kPos = read_position_.load(SYNC);
            if(kPos != write_position_.load(SYNC)){
                std::atomic_thread_fence( std::memory_order_acquire);
                element = queue_[kPos];
                read_position_.store(Increment(kPos),SYNC);
                return true;
            }
            return false;
        }

        inline void clear() noexcept {
            write_position_.store(read_position_.load(SYNC),SYNC);
        }

        inline constexpr Position Max() const noexcept {
            return SIZE;
        }

        inline size_t Size() const noexcept {
            const Position kRead = read_position_.load(SYNC);
            const Position kWrite = write_position_.load(SYNC);
            return (kWrite<kRead?SIZE-kRead+kWrite:kWrite-kRead);
        }

    private:
        std::array<T, SIZE> queue_;
        std::atomic<Position> read_position_;
        std::atomic<Position> write_position_;
};

#endif


