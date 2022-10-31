#ifndef BNMC_INCLUDE_XARY_H_
#define BNMC_INCLUDE_XARY_H_

#include <vector>
#include <set>
#include <exception/exception.h>
#include <bn-to-cnf/bayesnet.h>
#include <bn-to-cnf/cnf.h>
#include <math.h>
#include <algorithm>

create_exception(XAryException);

namespace bnc {


typedef unsigned int XAryBit;
typedef std::vector<XAryBit> XAryBits;

class XAry : public XAryBits {

    // least significal bit placed at index 0
    // most significal bit placed at index n
    public:
        typedef int Map;
        typedef std::vector<Map> Mapping;
        typedef XAryBit Bit;
        typedef unsigned int Dim;
        typedef std::vector<unsigned int> StepsizeList;
        typedef std::vector<unsigned int> Dimension;
        typedef uint16_t                  Variable;
        typedef std::vector<uint16_t>     ValueList;
        typedef std::vector<Variable>     VariableList;
        typedef std::set<Variable>        VariableSet;
        typedef std::vector<char>         RestrictDimension;

    protected:
        Dimension dimension_;

    public:
        inline size_t Size() const {
            return size();
        }

        inline void SetZero(){
            std::fill(begin(),end(), 0);
        }

        inline void Clear(){
            clear();
            dimension_.clear();
        }

        // determine step size per xary bit
        inline void CreateStepsizeList(StepsizeList &stepsize_list){
            stepsize_list.resize(dimension_.size());
            if(stepsize_list.size() > 0){
                stepsize_list[0] = 1;
                for(unsigned int i = 1; i < stepsize_list.size(); i++)
                    stepsize_list[i] = stepsize_list[i-1] * dimension_[i-1];
            }
        }

        // get domain size
        inline size_t Max(const RestrictDimension &kRestrictDimension) const {
            size_t max = 1;
            auto res_it = kRestrictDimension.begin();
            for(auto it = dimension_.begin(); it != dimension_.end(); it++,res_it++)
            if(!*res_it)
            max *= *it;

            return max;
        }

        inline size_t Max() const {
            size_t max = 1;
            for(auto it = dimension_.begin(); it != dimension_.end(); it++)
                max *= *it;

            return max;
        }

        template <class T>
        static size_t Max(bayesnet *bn, const T &kVariableList) {
            const uint32_t *kDimension = bn->get_states();

            size_t max = 1;
            for(auto it = kVariableList.begin(); it != kVariableList.end(); it++)
                max *= kDimension[*it];

            return max;
        }

        // set dimension
        template <class T, class U>
        static inline void SetDimension(T dimension, U variable_begin, const U variable_end, const uint32_t *const kDimensions){
            while(variable_begin != variable_end){
                *dimension = kDimensions[*variable_begin];
                ++dimension;
                ++variable_begin;
            }
        }

        template <class T>
        inline void SetDimension(bayesnet *bn, const T &kVariableList){
            const uint32_t *kDimension = bn->get_states();
            for(auto it = kVariableList.begin(); it != kVariableList.end(); it++)
                dimension_.push_back(kDimension[*it]);

            resize(dimension_.size());
        }

        template <class T>
        inline void SetDimension(bayesnet *bn, const T* kVariableList, const unsigned int kSize){
            const uint32_t *kDimension = bn->get_states();
            const auto kEnd = kVariableList + kSize;
            for(auto it = kVariableList; it != kEnd; it++)
                dimension_.push_back(kDimension[*it]);

            resize(dimension_.size());
        }

        template <class T>
        inline void SetDimension(const T* kDimensionList, const unsigned int kSize){
            const auto kEnd = kDimensionList + kSize;
            for(auto it = kDimensionList; it != kEnd; it++)
                dimension_.push_back(*it);

            resize(kSize);
        }

        inline void SetDimension(const Dimension &kDimension){
            dimension_ = kDimension;
            resize(kDimension.size());
        }

        inline Dimension &GetDimension() {
            return dimension_;
        }

        static inline void Set(Bit *bit, const Bit * const kBitEnd, const Bit * const kOtherBits, const Map* map){
            while(bit != kBitEnd){
                *bit = (*map < 0?0:kOtherBits[*map]);
                ++map;
                ++bit;
            }
        }

        inline void Set(const XAry &kXary, const Map* const &kMap){
            for(unsigned int i = 0; i < Size(); i++)
               (*this)[i] = (kMap[i] < 0?0:kXary[kMap[i]]);
        }

        inline void Set(const XAry &kXary, const Mapping &kMap){
            Set(kXary, &(kMap[0]));
        }

        inline void Set(const XAry &kXary){
            (*this) = kXary;
        }

        // create xary value mapping
        template <class T>
        static inline void CreateMapUnsorted(const T &kFrom, const T &kTo, Map* map){
            unsigned int to_pos = 0;
            for(auto it = kTo.begin(); it != kTo.end(); it++, to_pos++){
                const typename T::value_type &kVariable = *it;
                unsigned int from_pos = std::distance(kFrom.begin(),std::find(kFrom.begin(), kFrom.end(), kVariable));
                map[to_pos] = (from_pos < kFrom.size() ? from_pos:-1);
            }
        }

        template <class T, class U>
        static inline void CreateMap(T first1, const T last1,  U first2, const U last2, Map* result){
            // result must be of size kTo.size()
            unsigned int pos = 0;
            while (first1!=last1 && first2!=last2){
                if (*first1<*first2){
                    ++pos; ++first1;
                } else if (*first2<*first1){
                    *result = -1; ++result; ++first2;
                } else {
                    *result = pos; ++pos; ++result; ++first1; ++first2;
                }
            }
            while(first2!=last2){
                *result = -1; ++result; ++first2;
            }
        }

        template <class T, class U>
        static inline void CreateMap(T first1, const T last1,  U first2, const U last2, Mapping& result){
            result.resize(std::distance(first2, last2));
            CreateMap(first1, last1, first2, last2, &(result[0]));
        }

        template <class T, class U>
        static inline void CreateMap(const T &kFrom, const U &kTo, Map* result){
            CreateMap(kFrom.begin(), kFrom.end(), kTo.begin(), kTo.end(), result);
        }

        template <class T, class U>
        static inline void CreateMap(const T &kFrom, const U &kTo, Mapping &result){
            result.resize(kTo.size());
            CreateMap(kFrom.begin(), kFrom.end(), kTo.begin(), kTo.end(), &(result[0]));
        }

        static inline void CreateRestrictDimension(const Mapping &kMap, RestrictDimension &restrict_dimension){
            restrict_dimension.resize(kMap.size());
            auto kMapEnd = kMap.end();
            auto map_it = kMap.begin();
            auto dim_it = restrict_dimension.begin();
            while(map_it != kMapEnd){
                *dim_it = *map_it >= 0;
                ++dim_it;
                ++map_it;
            }
        }

        // convert xary number to decimal number
        static inline size_t GetDecimal(const Bit * const kBits, const Dim * const kDimensions, const unsigned int kSize){
            if(kSize == 0)
                return 0;

            size_t value = 0;
            auto dimension = &(kDimensions[kSize]);
            auto bit = &(kBits[kSize]);

            do {
                --bit;
                --dimension;

                value = value * (*dimension) + (*bit);
            } while(bit != kBits);

            return value;
        }

        inline size_t GetDecimal() const {
            const size_t MAX = Size()-1;

            size_t value = 0;
            for(int i = MAX; i >= 0; i--)
                value = value * dimension_[i] + (*this)[i];

            return value;
        }

        template <class T>
        static inline size_t GetDecimal(bayesnet *bn, const T &kVariableList, const ValueList &kValueList){
            const uint32_t *kDimension = bn->get_states();

            size_t value = 0;
            for(auto it = kVariableList.rbegin(); it != kVariableList.rend(); it++){
                const variable_t &kVariable = *it;
                value = value * kDimension[kVariable] + kValueList[kVariable];
            }

            return value;
        }

        // convert decimal number to xary number
        static inline void SetDecimal(Bit *const bits, const XAry::Dim* dimension, const unsigned int kSize,  unsigned int value){
            if(kSize == 0)
                return;

            unsigned int steps[kSize];

            auto step = steps;
            *step = 1;
            ++step;
            const auto kStepEnd = &(steps[kSize]);
            while(step != kStepEnd){
                *step = (*(step-1)) * (*dimension);
                ++step;
                ++dimension;
            }

            auto bit = &(bits[kSize]);
            do {
                --bit;
                --step;
                *bit = value / (*step);
                value -= (*bit) * (*step);
            } while(bit != bits);
        }

        inline void SetDecimal(const StepsizeList &kStepsizeList, const unsigned int kValue){
            #ifdef DEBUG
            if(Max() <= kValue)
                XAryException("Decimal exceeds max value for x-ary number");
            #endif

            unsigned int value = kValue;
            for(int i = kStepsizeList.size()-1; i >= 0; i--){
                Bit &bit = (*this)[i];
                bit = value / kStepsizeList[i];
                value -= bit * kStepsizeList[i];
            }
        }

        inline void SetDecimal(const unsigned int kValue){
            StepsizeList stepsize_list;
            CreateStepsizeList(stepsize_list);
            SetDecimal(stepsize_list,kValue);
        }


        // increment xary number by 1
        inline bool Increment(Bit * const kBitBegin, const Bit * const kBitEnd, const Dim *dimension){
            Bit *bit = kBitBegin;
            while(bit != kBitEnd){
                if(*bit < (*dimension)-1){
                    ++(*bit);
                    std::fill(kBitBegin,bit,0);
                    return true;
                }
                ++dimension;
                ++bit;
            }
            return false;
        }

        inline bool Increment(){
            auto dim_itr = dimension_.begin();
            for(auto bit_itr = begin(); bit_itr != end(); bit_itr++,dim_itr++){
                if(*bit_itr < *dim_itr -1){
                    (*bit_itr)++;
                    std::fill(begin(),bit_itr, 0);
                    return true;
                }
            }
            return false;
        }

        static inline bool RestrictedIncrement(Bit* const kBitBegin, const Bit* const kBitEnd, const Dim *dimension, const RestrictDimension::value_type *restricted){
            Bit *bit = kBitBegin;
            while(bit != kBitEnd){
                if(!(*restricted) && (*bit)+1 != (*dimension)){
                    ++(*bit);

                    if(bit != kBitBegin){ // bit > kBitBegin
                        do {
                            --bit;
                            --restricted;
                            if(!(*restricted))
                                *bit = 0;
                        } while(bit != kBitBegin);
                    }
                    return true;
                }
                ++dimension;
                ++bit;
            }
            return false;
        }

        inline bool RestrictedIncrement(const RestrictDimension &kRestrictDimension){
            auto dim_itr = dimension_.begin();
            auto res_itr = kRestrictDimension.begin();
            for(auto bit_itr = begin(); bit_itr != end(); bit_itr++,dim_itr++,res_itr++){
                if(!*res_itr && *bit_itr < *dim_itr -1){
                    (*bit_itr)++;

                    auto res2_itr = kRestrictDimension.begin();
                    for(auto bit2_itr = begin(); bit2_itr != bit_itr; bit2_itr++, res2_itr++)
                        if(!*res2_itr)
                            *bit2_itr = 0;

                    return true;
                }
            }
            return false;
        }
};

class rXAry : public XAry {
    public:

        inline size_t GetDecimal() const {
            const size_t MAX = Size();
            size_t value = 0;

            for(unsigned int i = 0; i < MAX; i++)
                value = value * dimension_[i] + (*this)[i];

            return value;
        }

        inline void CreateStepsizeList(StepsizeList &stepsize_list){
            stepsize_list.resize(dimension_.size());
            if(stepsize_list.size() > 0){
                stepsize_list[stepsize_list.size()-1] = 1;
                for(int i = stepsize_list.size()-2; i >= 0; i--)
                    stepsize_list[i] = stepsize_list[i+1] * dimension_[i+1];
            }
        }

        inline void SetDecimal(const StepsizeList &kStepsizeList, const unsigned int &kValue){
            #ifdef DEBUG
            if(Max() <= kValue)
                XAryException("Decimal exceeds max value for x-ary number");
            #endif

            unsigned int value = kValue;
            for(unsigned int i = 0; i < kStepsizeList.size(); i++){
                Bit &bit = (*this)[i];
                bit = value / kStepsizeList[i];
                value -= bit * kStepsizeList[i];
            }
        }

        inline void SetDecimal(const unsigned int &kValue){
            StepsizeList stepsize_list;
            CreateStepsizeList(stepsize_list);
            SetDecimal(stepsize_list,kValue);
        }

        inline bool Increment(){
            auto dim_itr = dimension_.rbegin();
            for(auto bit_itr = rbegin(); bit_itr != rend(); bit_itr++,dim_itr++){
                if(*bit_itr < *dim_itr -1){
                    (*bit_itr)++;
                    std::fill(rbegin(),bit_itr, 0);
                    return true;
                }
            }
            return false;
        }

        inline bool RestrictedIncrement(const RestrictDimension &kRestrictDimension){
            // not implemented
            return false;
        }
};

class sXAry : public XAry {
    private:
        StepsizeList stepsize_list_;

    public:
        inline void SetDimension(bayesnet *bn, const VariableList &kVariableList){
            XAry::SetDimension(bn, kVariableList);
            CreateStepsizeList(stepsize_list_);
        }

        inline void SetDimension(const Dimension &kDimension){
            XAry::SetDimension(kDimension);
            CreateStepsizeList(stepsize_list_);
        }

        inline void SetDecimal(const unsigned int &kValue){
            XAry::SetDecimal(stepsize_list_,kValue);
        }

        inline void Clear(){
            XAry::Clear();
            stepsize_list_.clear();
        }
};


class srXAry : public rXAry {
    private:
        StepsizeList stepsize_list_;

    public:
        inline void SetDimension(bayesnet *bn, const VariableList &kVariableList){
            rXAry::SetDimension(bn, kVariableList);
            CreateStepsizeList(stepsize_list_);
        }

        inline void SetDimension(const Dimension &kDimension){
            rXAry::SetDimension(kDimension);
            CreateStepsizeList(stepsize_list_);
        }

        template <class T>
        inline void SetDimension(const T* kDimension, const unsigned int kSize){
            rXAry::SetDimension(kDimension,kSize);
            CreateStepsizeList(stepsize_list_);
        }

        inline void SetDecimal(const unsigned int &kValue){
            rXAry::SetDecimal(stepsize_list_,kValue);
        }

        inline void Clear(){
            rXAry::Clear();
            stepsize_list_.clear();
        }
};

}
#endif
