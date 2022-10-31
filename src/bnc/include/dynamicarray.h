#ifndef BNC_DYNAMIC_ARRAY_H
#define BNC_DYNAMIC_ARRAY_H

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */
#include <string.h>
#include <cassert>

template <class T>
class DynamicArray {
    public:
        DynamicArray(){
            Init();
        }

        DynamicArray(size_t size){
            Init();
            bool success = Resize(size);
            assert(success && "Allocation failed");
        }

        ~DynamicArray(){
            Clear();
        }

        template <class U>
        inline T& operator[](U i){
            return data_[i];
        }

        template <class U>
        const T& operator[](U i) const {
            return data_[i];
        }

        inline void Clear() {
            size_ = 0;
            if(data_){
                free(data_);
                data_ = NULL;
            }
        }
        inline size_t GetSize() const {
            return size_;
        }

        inline void SetZero(){
            memset(data_, 0, size_*sizeof(T));
        }

        inline bool Resize(size_t size){
            if(size > 0){
                T *data = (T*) realloc (data_, size * sizeof(T));
                if(data){
                    data_ = data;
                    size_ = size;
                    return true;
                }
                assert(false && "allocation failed");
                return false;
            }
            return true;
        }

        inline T* GetRaw() {
            return &(data_[0]);
        }


        // stl style interface
        inline void clear() { Clear();}
        inline size_t size() const { return GetSize(); }
        inline bool resize(size_t size) { return Resize(size); }
        inline T* begin() { return &(data_[0]);}
        inline const T* begin() const { return &(data_[0]);}
        inline T* end() { return &(data_[size_]);}
        inline const T* end() const { return &(data_[size_]);}
    private:
        inline void Init(){
            data_ = NULL;
            size_ = 0;
        }

        T * data_;
        size_t size_;

};


#endif
