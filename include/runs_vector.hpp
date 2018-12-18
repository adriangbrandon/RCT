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
// Created by Adrián on 16/12/2018.
//

#ifndef RCT_RUNS_VECTOR_HPP
#define RCT_RUNS_VECTOR_HPP

#include <sdsl/vectors.hpp>
#include <sdsl/sd_vector.hpp>
#include <succ_support_v.hpp>
#include <prev_support_v.hpp>

namespace rct {

    template <class t_bitvector = sdsl::sd_vector<>, class t_succ = succ_support_v<>, class t_prev = prev_support_v<>>
    class runs_vector {

    public:

        typedef t_bitvector bitvector_type;
        typedef typename bitvector_type::rank_1_type rank_type;
        typedef typename bitvector_type::value_type value_type;
        typedef typename bitvector_type::size_type size_type;
        typedef t_succ succ_type;
        typedef t_prev prev_type;

    private:

        bitvector_type m_runs;
        rank_type m_rank;
        succ_type m_succ;
        prev_type m_prev;

        void copy(const runs_vector& o){
            m_runs = o.m_runs;
            m_rank = o.m_rank;
            m_rank.set_vector(&m_runs);
            m_succ = o.m_succ;
            m_succ.set_vector(&m_runs.high);
            m_prev = o.m_prev;
        }

    public:

        const bitvector_type& runs = m_runs;

        runs_vector() = default;

        runs_vector(const runs_vector& o)
        {
            copy(o);
        }

        runs_vector(runs_vector&& o)
        {
            *this = std::move(o);
        }

        runs_vector(const sdsl::bit_vector &vb){
            size_type k = 0, j = 1;
            bool ones;
            std::vector<size_type> pos_ones;
            for(const auto &v : vb){
                if(k == 0) {
                    ones = v;
                    pos_ones.push_back(0);
                }else{
                    if((v && ones) || (!v && !ones)){
                        ++j;
                    }else{
                        pos_ones.push_back(j);
                        ++j;
                        ones = !ones;
                    }
                }
                ++k;
            }
            m_runs = bitvector_type(pos_ones.begin(), pos_ones.end());
            sdsl::util::init_support(m_succ, &m_runs.high);
            sdsl::util::init_support(m_prev, &m_runs.high);
            sdsl::util::init_support(m_rank, &m_runs);
        }

        inline value_type operator[](const size_type i)const{
            return (m_rank(i+1) & 0x1);
        }

        inline value_type run(const size_type i, size_type &beg, size_type &end, size_type &r) const{
            r = m_rank(i+1);
            //select_prev
            beg = m_runs.low[r-1] + ((m_prev(i) + 1 - r) << (m_runs.wl));
            //select_next
            end = m_runs.low[r] + ((m_succ(i+1) + 1 - (r+1)) << (m_runs.wl)) - 1;
            return (r & 0x1);
        }

        inline value_type next_run(size_type &beg, size_type &end, size_type &r){
            ++r;
            beg = end +1;
            end = m_runs.low[r] + ((m_succ(beg+1) + 1 - (r+1)) << (m_runs.wl)) - 1;
            return (r & 0x1);
        }

        inline size_type size(){
            return m_runs.size();
        }

        void swap(runs_vector& v)
        {
            if (this != &v) {
                m_runs.swap(v.m_runs);
                sdsl::util::swap_support(m_rank, v.m_rank, &m_runs, &v.m_runs);
            }
        }

        runs_vector& operator=(const runs_vector& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        runs_vector& operator=(runs_vector&& v)
        {
            if (this != &v) {
                m_runs = std::move(v.m_runs);
                m_rank = std::move(v.m_rank);
                m_rank.set_vector(&m_runs);
            }
            return *this;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_runs.serialize(out, child, "runs");
            written_bytes += m_rank.serialize(out, child, "rank");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            m_runs.load(in);
            m_rank.load(in, &m_runs);
        }
    };

}

#endif //RCT_RUNS_VECTOR_HPP
