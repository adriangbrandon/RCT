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
// Created by Adrián on 17/12/2018.
//

#ifndef RCT_RUNS_BITVECTOR_HPP
#define RCT_RUNS_BITVECTOR_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <succ_support_v.hpp>
#include <cstdint>
#include <math_util.hpp>

namespace rct {

    template <uint64_t sample_min = 8, uint64_t sample_max = 4096,
                    class t_sampling = sdsl::bit_vector, class t_values = sdsl::bit_vector,
                    class t_sampling_rank = sdsl::rank_support_v<1>, class t_values_rank = sdsl::rank_support_v<1>,
                    class t_offset = sdsl::int_vector<>, class t_pos = sdsl::int_vector<>>
    class runs_bitvector {

    public:
        typedef t_pos pos_type;
        typedef t_offset offset_type;
        typedef t_sampling samp_bv_type;
        typedef t_sampling_rank rank_samp_bv_type;
        typedef t_values values_type;
        typedef t_values_rank rank_values_type;
        typedef uint64_t size_type;
        typedef uint64_t value_type;
        typedef struct {
            size_type sample_index;
            size_type ones;
            size_type beg;
            size_type pos;
            size_type next_pos;
            size_type length;
        } next_info_type;

    private:

        size_type m_size = 0;
        size_type m_sample;
        samp_bv_type m_sampling;
        rank_samp_bv_type m_rank_sampling;
        offset_type m_offset;
        pos_type m_pos;
        values_type m_values;
        rank_values_type m_rank_values;

        size_type init_tables(const sdsl::bit_vector &c, std::vector<uint8_t> &sampling, std::vector<uint64_t > &start,
                              std::vector<uint64_t > &end, size_type &length){
            const uint64_t* data = c.data();
            uint8_t offset = 0;
            size_type bits = 0, index = 0;
            size_type ones = 0;
            length = 0;
            while(bits < c.size()){
                uint64_t value = sdsl::bits::read_int_and_move(data, offset, sample_min);
                if(value>0){
                    sampling.push_back(1);
                    ++ones;
                    //std::cout << "start: " << sdsl::bits::lo(value) + index * sample_min << std::endl;
                    start.push_back(sdsl::bits::lo(value) + index * sample_min);
                    //std::cout << "end: " << sdsl::bits::hi(value) + index * sample_min << std::endl;
                    end.push_back(sdsl::bits::hi(value) + index * sample_min);
                    length += (end.back() - start.back() + 1);
                }else{
                    sampling.push_back(0);
                    start.push_back(1);
                    end.push_back(0);
                }
                bits += sample_min;
                ++index;
            }
            return 64 + sampling.size() + ones * (sdsl::bits::hi(sample_min)+1) + length + ones * (sdsl::bits::hi(length)+1);
        }

        size_type best_sample_dp_old(const sdsl::bit_vector &c, std::vector<uint8_t> &best_sampling, std::vector<uint64_t> &best_start,
                                  std::vector<uint64_t> &best_end, size_type &length_result){
            std::vector<uint8_t> sampling;
            std::vector<uint64_t > start, end;
            size_type length = 0;
            size_type best_size = init_tables(c, sampling, start, end, length);
            auto sample = sample_min;
            size_type best_sample = sample;
            size_type size = best_size;
            best_sampling = sampling;
            best_start = start;
            best_end = end;
            best_size = size;
            length_result = length;
            while (sample <= sample_max) {
                length = 0;
                size_type ones = 0;
                size_type new_sampling_size = 0;
                for(size_type i = 0; i < sampling.size(); i= i + 2){
                    uint8_t value = sampling[i];
                    if((sampling.size()-i) > 1){
                        value = (uint8_t) (value || sampling[i+1]);
                    }
                    if(value) ++ones;
                    sampling[i/2] = value;
                    ++new_sampling_size;
                }
                sampling.resize(new_sampling_size);

                size_type new_start_end_size = 0;
                for(size_type i = 0; i < start.size(); i = i+ 2){
                    if((start.size()-i) > 1){
                        if(start[i] <= end[i] && start[i+1] <= end[i+1]){
                            start[i/2] = start[i];
                            end[i/2]   = end[i+1];
                        }else if (start[i] > end[i] && start[i+1] <= end[i+1]){
                            start[i/2] = start[i+1];
                            end[i/2]   = end[i+1];
                        }else if (start[i+1] > end[i+1] && start[i] <= end[i]){
                            start[i/2] = start[i];
                            end[i/2]   = end[i];
                        }else{
                            start[i/2] = 1;
                            end[i/2] = 0;
                        }
                    }else{
                        start[i/2] = start[i];
                        end[i/2] = end[i];
                    }
                    if(end[i/2] >= start[i/2]){
                        length += (end[i/2] - start[i/2]+1);
                    }
                    ++new_start_end_size;
                }
                end.resize(new_start_end_size);
                start.resize(new_start_end_size);
                sample *= 2;
                //std::cout << "length: " << length << std::endl;
                size = 64 + sampling.size() + ones * (sdsl::bits::hi(sample)+1) + length + ones * (sdsl::bits::hi(length)+1);
                std::cout << "size: " << size << " with sample: " << sample <<  std::endl;
                if(size < best_size){
                    best_sampling = sampling;
                    best_start = start;
                    best_end = end;
                    best_size = size;
                    best_sample = sample;
                    length_result = length;
                }
                /*for(size_type i = 0; i < start.size(); ++i){
                    std::cout << "[" << start[i] << ", " << end[i] << "], ";
                }
                std::cout << std::endl;*/
            }
            // std::cout << "best_size: " << best_size << " best_sample: " << best_sample << std::endl;
            return best_sample;
        }

        size_type best_sample_dp(const sdsl::bit_vector &c, std::vector<uint8_t> &best_sampling, std::vector<uint64_t> &best_start,
                                 std::vector<uint64_t> &best_end, size_type &length_result){
            std::vector<uint8_t> sampling;
            std::vector<uint64_t > start, end;
            size_type length = 0;
            size_type best_size = init_tables(c, sampling, start, end, length);
            auto sample = sample_min;
            size_type best_sample = sample;
            size_type size = best_size;
            do {
                best_sampling = sampling;
                best_start = start;
                best_end = end;
                best_size = size;
                best_sample = sample;
                length_result = length;
                length = 0;
                size_type ones = 0;
                size_type new_sampling_size = 0;
                for(size_type i = 0; i < sampling.size(); i= i + 2){
                    uint8_t value = sampling[i];
                    if((sampling.size()-i) > 1){
                        value = (uint8_t) (value || sampling[i+1]);
                    }
                    if(value) ++ones;
                    sampling[i/2] = value;
                    ++new_sampling_size;
                }
                sampling.resize(new_sampling_size);

                size_type new_start_end_size = 0;
                for(size_type i = 0; i < start.size(); i = i+ 2){
                    if((start.size()-i) > 1){
                        if(start[i] <= end[i] && start[i+1] <= end[i+1]){
                            start[i/2] = start[i];
                            end[i/2]   = end[i+1];
                        }else if (start[i] > end[i] && start[i+1] <= end[i+1]){
                            start[i/2] = start[i+1];
                            end[i/2]   = end[i+1];
                        }else if (start[i+1] > end[i+1] && start[i] <= end[i]){
                            start[i/2] = start[i];
                            end[i/2]   = end[i];
                        }else{
                            start[i/2] = 1;
                            end[i/2] = 0;
                        }
                    }else{
                        start[i/2] = start[i];
                        end[i/2] = end[i];
                    }
                    if(end[i/2] >= start[i/2]){
                        length += (end[i/2] - start[i/2]+1);
                    }
                    ++new_start_end_size;
                }
                end.resize(new_start_end_size);
                start.resize(new_start_end_size);
                sample *= 2;
                //std::cout << "length: " << length << std::endl;
                size = 64 + sampling.size() + ones * (sdsl::bits::hi(sample)+1) + length + ones * (sdsl::bits::hi(length)+1);
               // std::cout << "size: " << size << " with sample: " << sample <<  std::endl;
            } while(size < best_size);
                /*for(size_type i = 0; i < start.size(); ++i){
                    std::cout << "[" << start[i] << ", " << end[i] << "], ";
                }
                std::cout << std::endl;*/
           // }
            // std::cout << "best_size: " << best_size << " best_sample: " << best_sample << std::endl;
            return best_sample;
        }

        void copy(const runs_bitvector &o){
            m_size = o.m_size;
            m_sample = o.m_sample;
            m_sampling = o.m_sampling;
            m_rank_sampling = o.m_rank_sampling;
            m_rank_sampling.set_vector(&m_sampling);
            m_offset = o.m_offset;
            m_pos = o.m_pos;
            m_values = o.m_values;
            m_rank_values = o.m_rank_values;
            m_rank_values.set_vector(&m_values);
        }

    public:

        const size_type &sample = m_sample;
        const samp_bv_type &sampling = m_sampling;
        const rank_samp_bv_type &rank_sampling = m_rank_sampling;
        const offset_type &offset = m_offset;
        const pos_type &pos = m_pos;
        const values_type &values = m_values;
        const rank_values_type &rank_values = m_rank_values;

        runs_bitvector() = default;

        runs_bitvector(const sdsl::bit_vector &c){
           // std::cout << "size: " << c.size() << std::endl;
            if(c.empty()) return;
            std::vector<uint8_t> best_sampling;
            std::vector<uint64_t > best_start, best_end;
            size_type length_result;
            m_size = c.size();
            m_sample = best_sample_dp(c, best_sampling, best_start, best_end, length_result);

           /* size_type index = 0;
            std::cout << std::endl << "Sampling: " << std::endl;
            for(const auto &v : best_sampling){
                std::cout << "i: " << index << " value: " << (size_type) v << std::endl;
                ++index;
            }

            index =0;
            std::cout << std::endl << "Start: " << std::endl;
            for(const auto &v : best_start){
                std::cout << "i: " << index << " value: " << (size_type) v << std::endl;
                ++index;
            }

            index =0;
            std::cout << std::endl << "End: " << std::endl;
            for(const auto &v : best_end){
                std::cout << "i: " << index << " value: " << (size_type) v << std::endl;
                ++index;
            }*/

            std::vector<size_type> pos, offsets;
            m_sampling.resize(best_sampling.size());
            for(size_type i = 0; i < best_sampling.size(); ++i){
                m_sampling[i] = best_sampling[i];
            }
            size_type j = 0;
            m_values.resize(length_result);
            for(size_type i = 0; i < best_start.size(); ++i){
                if(best_start[i] <= best_end[i]){
                    offsets.push_back(best_start[i] - m_sample*i);
                    pos.push_back(j);
                }
                for(size_type k = best_start[i]; k <= best_end[i]; ++k){
                   // std::cout << "k: " << k << std::endl;
                    m_values[j] = c[k];
                    j++;
                }
            }
            pos.push_back(j);
            sdsl::util::init_support(m_rank_sampling, &m_sampling);
            sdsl::util::init_support(m_rank_values, &m_values);
            m_offset.resize(offsets.size());
            m_pos.resize(pos.size());
            size_type idx = 0;
            for(const auto &v : offsets){
                m_offset[idx++] = v;
            }
            idx = 0;
            for(const auto &v : pos){
                m_pos[idx++] = v;
            }
            sdsl::util::bit_compress(m_offset);
            sdsl::util::bit_compress(m_pos);
        }

        //! Copy constructor
        runs_bitvector(const runs_bitvector& o)
        {
            copy(o);
        }

        //! Move constructor
        runs_bitvector(runs_bitvector&& o)
        {
            *this = std::move(o);
        }


        runs_bitvector &operator=(const runs_bitvector &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }
        runs_bitvector &operator=(runs_bitvector &&o) {
            if (this != &o) {
                m_size = o.m_size;
                m_sample = o.m_sample;
                m_sampling = std::move(o.m_sampling);
                m_rank_sampling = std::move(o.m_rank_sampling);
                m_rank_sampling.set_vector(&m_sampling);
                m_offset = std::move(o.m_offset);
                m_pos = std::move(o.m_pos);
                m_values = std::move(o.m_values);
                m_rank_values = std::move(o.m_rank_values);
                m_rank_values.set_vector(&m_values);
            }
            return *this;
        }

        void swap(runs_bitvector &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            m_size = o.m_size;
            std::swap(m_size, o.m_size);
            std::swap(m_sample, o.m_sample);
            m_sampling.swap(o.m_sampling);
            sdsl::util::swap_support(m_rank_sampling, o.m_rank_sampling, &m_sampling, &o.m_sampling);
            m_offset.swap(o.m_offset);
            m_pos.swap(o.m_pos);
            m_values.swap(o.m_values);
            sdsl::util::swap_support(m_rank_values, o.m_rank_values, &m_values, &o.m_values);
        }

        inline value_type access(const size_type i) const{
            if(!m_sampling[i/m_sample]) return 0;
            auto ones = m_rank_sampling(i/m_sample + 1);
            auto pos = m_pos[ones-1];
            auto beg = i/m_sample* m_sample + m_offset[ones-1];
            auto end = beg + m_pos[ones] - pos - 1;
            return (value_type) (beg <= i && i <= end && m_values[i-beg + pos]);
        }

        inline value_type operator[](const size_type i) const{
            return access(i);
        }

        inline size_type size() const {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_values.serialize(out, child, "values");
            written_bytes += m_rank_values.serialize(out, child, "rank_values");
            written_bytes += m_sampling.serialize(out, child, "sampling");
            written_bytes += m_rank_sampling.serialize(out, child, "rank_sampling");
            written_bytes += m_offset.serialize(out, child, "offset");
            written_bytes += m_pos.serialize(out, child, "m_pos");
            written_bytes += sdsl::serialize(m_size, out, child, "size");
            written_bytes += sdsl::serialize(m_sample, out, child, "sample");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in)
        {
            m_values.load(in);
            m_rank_values.load(in, &m_values);
            m_sampling.load(in);
            m_rank_sampling.load(in, &m_sampling);
            m_offset.load(in);
            m_pos.load(in);
            sdsl::load(m_size, in);
            sdsl::load(m_sample, in);
        }

        void print(){
            std::cout << "sampling " << std::endl;
            for(size_type i = 0; i < m_sampling.size(); ++i){
                std::cout << "i=" << i << " v=" << m_sampling[i] << std::endl;
            }
            std::cout << std::endl;
        }
    };

    template<uint8_t t_b>
    struct rank_support_runs_bitvector_trait {
        typedef runs_bitvector<>::size_type size_type;
        static size_type adjust_rank(size_type r,size_type)
        {
            return r;
        }
    };

    template<>
    struct rank_support_runs_bitvector_trait<0> {
        typedef runs_bitvector<>::size_type size_type;
        static size_type adjust_rank(size_type r, size_type n)
        {
            return n - r;
        }
    };

    template<uint8_t t_b>
    class rank_support_runs_bitvector {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef runs_bitvector<> bit_vector_type;
        typedef typename bit_vector_type::next_info_type next_info_type;
    private:
        const bit_vector_type* m_v;
    public:

        explicit rank_support_runs_bitvector(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }


        //! Returns the position of the i-th occurrence in the bit vector.
        size_type rank(size_type i)const
        {
            assert(i >= 0); assert(i <= m_v->size());
            //auto ones = m_v->rank_sampling(i/m_v->sample + 1);
            size_type r;
            if(i == m_v->size()){
                r = m_v->rank_values(m_v->values.size());
                return rank_support_runs_bitvector_trait<t_b>::adjust_rank(r, i);
            }
            auto ones = m_v->rank_sampling(i/m_v->sample+1);
            //sampling[i] = 0
            if(!m_v->sampling[i/m_v->sample]){
                r = m_v->rank_values(m_v->pos[ones]);
                return rank_support_runs_bitvector_trait<t_b>::adjust_rank(r, i);
            }
            //sampling[i] = 1
            auto pos = m_v->pos[ones-1];
            auto length = m_v->pos[ones] - pos;
            auto beg = i/m_v->sample* m_v->sample + m_v->offset[ones-1];
            if(beg > i){
                r = m_v->rank_values(pos);
            }else if (i >= beg + length){
                r = m_v->rank_values(pos+length);
            }else{
                r = m_v->rank_values(pos + (i-beg));
            }
            return rank_support_runs_bitvector_trait<t_b>::adjust_rank(r, i);
        }


        size_type rank(size_type i, next_info_type &next_info)const
        {
            assert(i >= 0); assert(i <= m_v->size());
            size_type r;
            if(i == m_v->size()){
                r = m_v->rank_values(m_v->values.size());
                return rank_support_runs_bitvector_trait<t_b>::adjust_rank(r, i);
            }
            auto ones = m_v->rank_sampling(i/m_v->sample + 1);
            auto pos = m_v->pos[ones-1];
            auto beg = i/m_v->sample* m_v->sample + m_v->offset[ones-1];
            auto next_pos = m_v->pos[ones];
            auto length = next_pos - pos;
            next_info = {i/m_v->sample, ones, beg, pos, next_pos, length};
            if(!m_v->sampling[i/m_v->sample]){
                r = m_v->rank_values(next_pos);
                return rank_support_runs_bitvector_trait<t_b>::adjust_rank(r, i);
            }
            if(beg > i){
                r = m_v->rank_values(pos);
            }else if (i >= beg + length){
                r = m_v->rank_values(pos+length);
            }else{
                r = m_v->rank_values(pos + (i-beg));
            }
            return rank_support_runs_bitvector_trait<t_b>::adjust_rank(r, i);
        }

        size_type operator()(size_type i)const
        {
            return rank(i);
        }

        size_type operator()(size_type i, size_type &ones, size_type &beg, size_type &length)const
        {
            return rank(i, ones, beg, length);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        rank_support_runs_bitvector& operator=(const rank_support_runs_bitvector& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(rank_support_runs_bitvector&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }


    };


    template <uint8_t t_b>
    class succ_support_runs_bitvector {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef runs_bitvector<> bit_vector_type;
        typedef typename bit_vector_type::next_info_type next_info_type;
    private:
        const bit_vector_type* m_v;
        succ_support_v<t_b> m_succ_sampling;
        succ_support_v<t_b> m_succ_values;
        sdsl::int_vector<> m_delta_next_0;
    public:

        succ_support_runs_bitvector() = default;

        succ_support_runs_bitvector(const bit_vector_type* v);

        size_type succ(size_type i) const;

        size_type succ(size_type i, next_info_type &next_info) const;

        size_type operator()(size_type i)const
        {
            return succ(i);
        }

        size_type operator()(size_type i, size_type &ones, size_type &beg, size_type &length)const
        {
            return succ(i, ones, beg, length);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
            m_succ_sampling.set_vector(&(m_v->sampling));
            m_succ_values.set_vector(&(m_v->values));
        }

        succ_support_runs_bitvector& operator=(const succ_support_runs_bitvector& ss)
        {
            if (this != &ss) {
                m_succ_sampling = ss.m_succ_sampling;
                m_succ_values = ss.m_succ_values;
                m_delta_next_0 = ss.m_delta_next_0;
                set_vector(ss.m_v);
            }
            return *this;
        }

        succ_support_runs_bitvector& operator=(succ_support_runs_bitvector&& ss)
        {
            if (this != &ss) {
                m_succ_sampling = std::move(ss.m_succ_sampling);
                m_succ_values = std::move(ss.m_succ_values);
                m_delta_next_0 = std::move(ss.m_delta_next_0);
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(succ_support_runs_bitvector& ss) {
            if (this != &ss) {
                std::swap(m_v, ss.m_v);
                m_succ_sampling.swap(ss.m_succ_sampling);
                m_succ_values.swap(ss.m_succ_values);
                if(m_v != nullptr){
                    m_succ_sampling.set_vector(&(m_v->sampling));
                    m_succ_values.set_vector(&(m_v->values));
                }
                if(ss.m_v != nullptr){
                    ss.m_succ_sampling.set_vector(&(ss.m_v->sampling));
                    ss.m_succ_values.set_vector(&(ss.m_v->values));
                }
                m_delta_next_0.swap(ss.m_delta_next_0);
            }
        }

        void load(std::istream& in, const bit_vector_type* v=nullptr)
        {
            m_v = v;
            m_succ_sampling.load(in, &(m_v->sampling));
            m_succ_values.load(in, &(m_v->values));
            m_delta_next_0.load(in);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_succ_sampling.serialize(out, child, "succ_sampling");
            written_bytes += m_succ_values.serialize(out, child, "succ_values");
            written_bytes += m_delta_next_0.serialize(out, child, "delta_next");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }


    };

    template<>
    succ_support_runs_bitvector<1>::succ_support_runs_bitvector(const bit_vector_type* v){
        set_vector(v);
        size_type ones = 0;
        std::vector<size_type > deltas;
        for(size_type i = 0; i < m_v->sampling.size(); ++i){
            if(m_v->sampling[i]){
                bool delta = false;
                for(size_type j = (i+1)*m_v->sample; j < m_v->size(); ++j){
                    if(!m_v->access(j)) {
                        auto next = j / m_v->sample;
                        deltas.push_back(next-i);
                        delta = true;
                        break;
                    }
                }
                if(!delta){
                    deltas.push_back(0);
                }
            }
        }
        sdsl::util::init_support(m_succ_sampling, &(m_v->sampling));
        sdsl::util::init_support(m_succ_values, &(m_v->values));
    }

    template<>
    succ_support_runs_bitvector<0>::succ_support_runs_bitvector(const bit_vector_type* v){
        set_vector(v);
        size_type ones = 0;
        std::vector<size_type > deltas;
        for(size_type i = 0; i < m_v->sampling.size(); ++i){
            if(m_v->sampling[i]){
                bool delta = false;
                for(size_type j = (i+1)*m_v->sample; j < m_v->size(); ++j){
                    if(!m_v->access(j)) {
                        auto next = j / m_v->sample;
                        deltas.push_back(next-i);
                        delta = true;
                        break;
                    }
                }
                if(!delta){
                    deltas.push_back(0);
                }
            }
        }

        m_delta_next_0.resize(deltas.size());
        size_type i = 0;
        for(const auto &delta : deltas){
            m_delta_next_0[i] = delta;
            ++i;
        }
        sdsl::util::bit_compress(m_delta_next_0);
        sdsl::util::init_support(m_succ_values, &(m_v->values));
    }

    //! Returns the position of the i-th occurrence in the bit vector.
    template<>
    typename succ_support_runs_bitvector<1>::size_type succ_support_runs_bitvector<1>::succ(size_type i)const
    {
        assert(i >= 0); assert(i < m_v->size());
        auto sample_index = i/m_v->sample;
        auto ones = m_v->rank_sampling(sample_index + 1);
        size_type beg = 0;
        if(m_v->sampling[sample_index]){
            beg = sample_index * m_v->sample + m_v->offset[ones-1];
        }else{
            ++ones;
            ++sample_index;
            beg = m_succ_sampling(sample_index)* m_v->sample + m_v->offset[ones-1];
        }
        auto pos = m_v->pos[ones-1];
        auto length = m_v->pos[ones] - pos;
        size_type r;
        if(beg >= i){
            return beg; //the first value is a one
        }else if (i >= beg + length){ //the first value of the next block is a one
            if(ones == m_v->offset.size()) return m_v->size();
            return m_succ_sampling(sample_index+1) * m_v->sample + m_v->offset[ones];
        }else{
            r = m_succ_values(pos + (i-beg)); //the next one is at the same block
            return beg + (r-pos);
        }
    }



    //! Returns the position of the i-th occurrence in the bit vector.
    template<>
    typename succ_support_runs_bitvector<0>::size_type succ_support_runs_bitvector<0>::succ(size_type i)const
    {
        assert(i >= 0); assert(i < m_v->size());
        auto sample_index = i/m_v->sample;
        if(!m_v->sampling[sample_index]) return i;
        auto ones = m_v->rank_sampling(sample_index + 1);
        auto beg = sample_index * m_v->sample + m_v->offset[ones-1];
        auto pos = m_v->pos[ones-1];
        auto next_pos = m_v->pos[ones];
        auto length = next_pos - pos;
        if(beg > i || beg + length <= i){
            return i;
        }else{
            auto next_0_value = m_succ_values((i-beg) + pos);
            if(next_0_value < next_pos){
                return beg + (next_0_value - pos);
            }
            if((beg + length) % m_v->sample != 0){
                return beg + length;
            }
            if(next_0_value > m_v->values.size()) return m_v->size();
            auto delta = m_delta_next_0[ones-1];
            if(delta == 0) return m_v->size();
            auto new_index = sample_index + delta;
            if (!m_v->sampling[new_index]) return new_index * m_v->sample;
            auto offset = m_v->offset[ones + delta - 1];
            if(offset > 0){
                return new_index * m_v->sample;
            }else{
                pos = m_v->pos[ones + delta - 1];
                next_pos = m_v->pos[ones + delta];
                beg = new_index * m_v->sample + offset;
                if(next_0_value < next_pos){
                    return beg + (next_0_value - pos);
                }else{
                    return beg + (next_pos - pos);
                }
            }
            /*pos = m_v->pos[ones + delta - 1];
            next_pos = m_v->pos[ones + delta];
            beg = new_index * m_v->sample + offset;
            if (pos <= next_0_value && next_0_value < next_pos) {
                return (next_0_value - pos) + beg;
            } else if (next_0_value < pos) {
                return beg;
            } else {
                return beg + (next_pos - pos);
            }*/

        }
    }

    /*template <>
    succ_support_runs_bitvector<1>::size_type succ_support_runs_bitvector<1>::succ(size_type i, next_info_type &next_info)const
    {
        assert(i >= 0); assert(i < m_v->size());
        if(i >= beg + length){
            if(ones == m_v->offset.size()) return m_v->size();
            ++ones;
            beg = m_succ_sampling(beg/m_v->sample+1)* m_v->sample + m_v->offset[ones-1];
            length = m_v->pos[ones] - m_v->pos[ones-1];
            return beg;
        }else if(beg >= i){
            return beg;
        }else{
            auto pos = m_v->pos[ones-1];
            auto r = m_succ_values(pos + (i-beg));
            return beg + (r - pos);
        }
    }*/

    //! Returns the position of the i-th occurrence in the bit vector.
    template<>
    typename succ_support_runs_bitvector<0>::size_type succ_support_runs_bitvector<0>::succ(size_type i, next_info_type &next_info)const
    {
        assert(i >= 0); assert(i < m_v->size());
        auto sample_index = i/m_v->sample;
        if(!m_v->sampling[sample_index]) return i;
        size_type ones, beg, pos, next_pos, length;
        if(sample_index != next_info.sample_index){
            ones = m_v->rank_sampling(sample_index + 1);
            beg = sample_index * m_v->sample + m_v->offset[ones-1];
            pos = m_v->pos[ones-1];
            next_pos = m_v->pos[ones];
            length = next_pos - pos;
            next_info = {i/m_v->sample, ones, beg, pos, next_pos, length};
        }else{
            ones = next_info.ones;
            beg = next_info.beg;
            pos = next_info.pos;
            next_pos = next_info.next_pos;
            length = next_info.length;
        }
        if(beg > i || beg + length <= i){
            return i;
        }else{
            auto next_0_value = m_succ_values((i-beg) + pos);
            if(next_0_value < next_pos){
                return beg + (next_0_value - pos);
            }
            if((beg + length) % m_v->sample != 0){
                return beg + length;
            }
            if(next_0_value > m_v->values.size()) return m_v->size();
            auto delta = m_delta_next_0[ones-1];
            if(delta == 0) return m_v->size();
            auto new_index = sample_index + delta;
            if (!m_v->sampling[new_index]) return new_index * m_v->sample;
            auto offset = m_v->offset[ones + delta - 1];
            if(offset > 0){
                return new_index * m_v->sample;
            }else{
                pos = m_v->pos[ones + delta - 1];
                next_pos = m_v->pos[ones + delta];
                beg = new_index * m_v->sample + offset;
                if(next_0_value < next_pos){
                    return beg + (next_0_value - pos);
                }else{
                    return beg + (next_pos - pos);
                }
            }
        }
    }
}

#endif //RCT_RUNS_BITVECTOR_HPP
