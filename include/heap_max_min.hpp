/***
BSD 2-Clause License

Copyright (c) 2018, Adrián
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/


//
// Created by Adrián on 25/07/2019.
//

#ifndef UTIL_HEAP_MAX_MIN_HPP
#define UTIL_HEAP_MAX_MIN_HPP

#include <vector>
#include <iostream>
#include <algorithm>

namespace util {




    template <class t_object, class t_value>
    class heap_max_min {
    public:
        typedef t_object object_type;
        typedef
        typedef struct {
            object_type object;
            bool deleted;
        } element_type;
        typedef uint64_t ptr_type;
        typedef struct {
            t_value value;
            ptr_type ptr;
        } value_type;
        typedef struct {
            t_value min;
            t_value max;
            ptr_type ptr;
        } value_pair_type;
        typedef struct {
            bool operator ()(const value_type& p1, const value_type& p2) const {
                return p1.value < p2.value;
            }
        } value_compare;
        typedef struct {
            bool operator ()(const value_pair_type& p1, const value_pair_type& p2) const {
                if(p1.max == p2.max) return p1.min < p2.min;
                return p1.max < p2.max;
            }
        } value_pair_compare;
        typedef std::vector<element_type>           container_element_type;
        typedef std::vector<value_type>             container_value_type;
        typedef std::vector<value_pair_type>             container_pair_type;
    private:
        container_pair_type m_max;
        container_value_type m_min;
        container_element_type m_elements;
        value_compare  m_comp;
        value_pair_compare m_pair_comp;

        void copy(const heap_max_min &p){
            m_min = p.m_min;
            m_max = p.m_max;
            m_elements = p.m_elements;
            m_comp = p.m_comp;
            m_pair_comp = p.m_pair_comp;
        }

        bool delete_object(const ptr_type ptr){
            if(m_elements.empty()) return false;
            if(ptr < m_elements.size()-1){
                //std::memmove(&m_elements[ptr], &m_elements[ptr+1], sizeof(element_type) * (m_elements.size()-1));
                m_elements.erase(m_elements.begin() + ptr);
            }else{
                m_elements.pop_back();
            }
            for(uint64_t i = 0; i < m_min.size(); ++i){
                if(m_min[i].ptr>ptr) --m_min[i].ptr;
            }
            for(uint64_t i = 0; i < m_max.size(); ++i){
                if(m_max[i].ptr>ptr) --m_max[i].ptr;
            }
            return true;
        }

    public:


        heap_max_min() = default;

        /*explicit heap(const container_type &values){
            m_values = values;
            std::make_heap(m_values.begin(), m_values.end(), m_comp);
        }*/


        inline void pop(){
           auto ptr =  m_max.front().ptr;
           m_elements[ptr].deleted = true;
           std::pop_heap(m_max.begin(), m_max.end(), m_pair_comp);
           m_max.pop_back();
        }

        inline std::pair<t_value, t_value> top_min_max(){
            auto ptr = m_min.front().ptr;
            while(m_elements[ptr].deleted){
                std::pop_heap(m_min.begin(), m_min.end(), m_comp);
                m_min.pop_back();
                delete_object(ptr);
                ptr = m_min.front().ptr;
            }
            return {m_min.front().value, m_max.front().max};

        }

        inline t_object top(){
            auto ptr = m_max.front().ptr;
            return m_elements[ptr].object;

        }

        inline void push(const t_object &o, const t_value &min, const t_value &max){
            m_elements.push_back(element_type{o, false});
            auto ptr = m_elements.size()-1;
            m_max.push_back(value_pair_type{min, max, ptr});
            m_min.push_back(value_type{min, ptr});
            std::push_heap(m_min.begin(), m_min.end(), m_comp);
            std::push_heap(m_max.begin(), m_max.end(), m_pair_comp);
        }


        inline void clear(){
            m_elements.clear();
            m_max.clear();
            m_min.clear();
        }

        inline bool empty(){
            return m_max.empty();
        }

        inline uint64_t size(){
            return m_max.size();
        }


        //! Copy constructor
        heap_max_min(const heap_max_min& o)
        {
            copy(o);
        }

        //! Move constructor
        heap_max_min(heap_max_min&& o)
        {
            *this = std::move(o);
        }


        heap_max_min &operator=(const heap_max_min &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        heap_max_min &operator=(heap_max_min &&o) {
            if (this != &o) {
                m_elements = std::move(o.m_elements);
                m_min = std::move(o.m_min);
                m_max = std::move(o.m_max);
                m_comp = std::move(o.m_comp);
                m_pair_comp = std::move(o.m_pair_comp);
            }
            return *this;
        }

        void swap(heap_max_min &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            std::swap(m_elements, o.m_elements);
            std::swap(m_min, o.m_min);
            std::swap(m_max, o.m_max);
            std::swap(m_comp, o.m_comp);
            std::swap(m_pair_comp, o.m_pair_comp);
        }

    };

}

#endif
