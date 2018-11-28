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
// Created by Adrián on 27/11/2018.
//

#ifndef RCT_LOG_REFERENCE_HPP
#define RCT_LOG_REFERENCE_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <rmq_succinct_ct.hpp>
#include <sdsl/sd_vector.hpp>

namespace rct {

    template <class t_movements = sdsl::bit_vector>
    class log_reference {

    public:
        typedef uint64_t size_type;
        typedef t_movements move_type;
        typedef sdsl::rank_support_v5<1> rank_move_type;
        typedef typename move_type::select_1_type select_move_type;
        typedef sdsl::sd_vector<> sampling_type;
        typedef typename sampling_type::rank_1_type rank_sampling_type;
        typedef typename sampling_type::select_1_type select_sampling_type;

    private:

        move_type m_x_p;
        rank_move_type m_rank_x_p;
        select_move_type m_select_x_p;

        move_type m_y_p;
        rank_move_type m_rank_y_p;
        select_move_type m_select_y_p;

        move_type m_x_n;
        rank_move_type m_rank_x_n;
        select_move_type m_select_x_n;

        move_type m_y_n;
        rank_move_type m_rank_y_n;
        select_move_type m_select_y_n;

        sampling_type m_rMmq_x_sample;
        sampling_type m_rMmq_y_sample;
        rank_sampling_type m_rank_rMmq_x_sample;
        rank_sampling_type m_rank_rMmq_y_sample;
        select_sampling_type m_select_rMmq_x_sample;
        select_sampling_type m_select_rMmq_y_sample;
        rmq_succinct_ct<true> m_rmq_x;
        rmq_succinct_ct<true> m_rmq_y;
        rmq_succinct_ct<false> m_rMq_x;
        rmq_succinct_ct<false> m_rMq_y;

        void copy(const log_reference &o){
            m_x_p = o.m_x_p;
            m_select_x_p = o.m_select_x_p;
            m_select_x_p.set_vector(&m_x_p);
            m_rank_x_p = o.m_rank_x_p;
            m_rank_x_p.set_vector(&m_x_p);
            m_x_n = o.m_x_n;
            m_select_x_n = o.m_select_x_n;
            m_select_x_n.set_vector(&m_x_n);
            m_rank_x_n = o.m_rank_x_n;
            m_rank_x_n.set_vector(&m_x_n);
            m_y_p = o.m_y_p;
            m_select_y_p = o.m_select_y_p;
            m_select_y_p.set_vector(&m_y_p);
            m_rank_y_p = o.m_rank_y_p;
            m_rank_y_p.set_vector(&m_y_p);
            m_y_n = o.m_y_n;
            m_select_y_n = o.m_select_y_n;
            m_select_y_n.set_vector(&m_y_n);
            m_rank_y_n = o.m_rank_y_n;
            m_rank_y_n.set_vector(&m_y_n);
            m_rMmq_x_sample = o.m_rMmq_x_sample;
            m_rank_rMmq_x_sample = o.m_rank_rMmq_x_sample;
            m_rank_rMmq_x_sample.set_vector(&m_rMmq_x_sample);
            m_select_rMmq_x_sample = o.m_select_rMmq_x_sample;
            m_select_rMmq_x_sample.set_vector(&m_rMmq_x_sample);
            m_rMmq_y_sample = o.m_rMmq_y_sample;
            m_rank_rMmq_y_sample = o.m_rank_rMmq_y_sample;
            m_rank_rMmq_y_sample.set_vector(&m_rMmq_y_sample);
            m_select_rMmq_y_sample = o.m_select_rMmq_y_sample;
            m_select_rMmq_y_sample.set_vector(&m_rMmq_y_sample);
            m_rmq_x = o.m_rmq_x;
            m_rMq_x = o.m_rMq_x;
            m_rmq_y = o.m_rmq_y;
            m_rMq_y = o.m_rMq_y;
        }

        template <class T>
        void compute_movement_repr(const T &movement, std::vector<size_type> &x_p_ones, std::vector<size_type> &x_n_ones,
                            std::vector<size_type> &y_p_ones, std::vector<size_type> &y_n_ones,
                            size_type &x_p_idx, size_type &x_n_idx, size_type &y_p_idx, size_type &y_n_idx){

            if(movement.x < 0){
                x_n_idx -= movement.x;
                x_n_ones.push_back(x_n_idx);
                x_n_idx++;
                x_p_ones.push_back(x_p_idx);
                x_p_idx++;
            }else{
                x_p_idx += movement.x;
                x_p_ones.push_back(x_p_idx);
                x_p_idx++;
                x_n_ones.push_back(x_n_idx);
                x_n_idx++;
            }
            if(movement.y < 0){
                y_n_idx -= movement.y;
                y_n_ones.push_back(y_n_idx);
                y_n_idx++;
                y_p_ones.push_back(y_p_idx);
                y_p_idx++;
            }else{
                y_p_idx += movement.y;
                y_p_ones.push_back(y_p_idx);
                y_p_idx++;
                y_n_ones.push_back(y_n_idx);
                y_n_idx++;
            }
        }

        template <class T>
        void compute_rmMq_repr(const T &diff, const size_type index, std::vector<T> &min, std::vector<T> &max,
                              std::vector<size_type> &rmMq_ones, bool &increase, bool &first, T &value){
            if(index == 0){
                increase = diff >= 0;
            }else{
                if(increase && diff < 0){
                    increase = false;
                    rmMq_ones.push_back(index);
                    max.push_back(value);
                    first = false;
                }else if(!increase && diff > 0){
                    //Flag. Avoiding 2 ranks into a sd_vector
                    if(first) rmMq_ones.push_back(0);
                    increase = true;
                    rmMq_ones.push_back(index);
                    min.push_back(value);
                    first = false;
                }
            }
            value += diff;
        }

    public:

        log_reference() = default;

        template <class ContainerMovements>
        log_reference(const ContainerMovements &movements){

            //1. Bitmap representation of the movements
            {
                std::vector<size_type> x_n_ones, x_p_ones, y_n_ones, y_p_ones;
                size_type x_p_idx = 0, x_n_idx = 0, y_p_idx = 0, y_n_idx = 0;
                for(const auto &movement : movements){
                    compute_movement_repr(movement, x_p_ones, x_n_ones, y_p_ones, y_n_ones, x_p_idx, x_n_idx, y_p_idx, y_n_idx);
                }
                m_x_p = sdsl::bit_vector(x_p_idx, 0);
                m_x_n = sdsl::bit_vector(x_n_idx, 0);
                m_y_p = sdsl::bit_vector(y_p_idx, 0);
                m_y_n = sdsl::bit_vector(y_n_idx, 0);
                for(const auto &one : x_p_ones){
                    m_x_p[one]=1;
                }
                for(const auto &one : x_n_ones){
                    m_x_n[one]=1;
                }
                for(const auto &one : y_p_ones){
                    m_y_p[one]=1;
                }
                for(const auto &one : y_n_ones){
                    m_y_n[one]=1;
                }
                sdsl::util::init_support(m_select_x_p, &m_x_p);
                sdsl::util::init_support(m_rank_x_p, &m_x_p);
                sdsl::util::init_support(m_select_x_n, &m_x_n);
                sdsl::util::init_support(m_rank_x_n, &m_x_n);
                sdsl::util::init_support(m_select_y_p, &m_y_p);
                sdsl::util::init_support(m_rank_y_p, &m_y_p);
                sdsl::util::init_support(m_select_y_n, &m_y_n);
                sdsl::util::init_support(m_rank_y_n, &m_y_n);
            }

            //2. Compute rMmq for x-axis
            {
                std::vector<int32_t> min, max;
                std::vector<size_type> rmMq_ones;
                bool increase = false, first = true;
                int32_t value = 0;
                for(size_type i = 0; i < movements.size(); ++i){
                    compute_rmMq_repr(movements[i].x, i, min, max, rmMq_ones, increase, first, value);
                }
                m_rMmq_x_sample = sampling_type(rmMq_ones.begin(), rmMq_ones.end());
                sdsl::util::init_support(m_rank_rMmq_x_sample, &m_rMmq_x_sample);
                sdsl::util::init_support(m_select_rMmq_x_sample, &m_rMmq_x_sample);
                m_rmq_x = rmq_succinct_ct<true>(min);
                m_rMq_x = rmq_succinct_ct<false>(max);
            }

            //3. Compute rMmq for y-axis
            {
                std::vector<int32_t> min, max;
                std::vector<size_type> rmMq_ones;
                bool increase = false, first = true;
                int32_t value = 0;
                for(size_type i = 0; i < movements.size(); ++i){
                    compute_rmMq_repr(movements[i].y, i, min, max, rmMq_ones, increase, first, value);
                }
                m_rMmq_y_sample = sampling_type(rmMq_ones.begin(), rmMq_ones.end());
                sdsl::util::init_support(m_rank_rMmq_y_sample, &m_rMmq_y_sample);
                sdsl::util::init_support(m_select_rMmq_y_sample, &m_rMmq_y_sample);
                m_rmq_y = rmq_succinct_ct<true>(min);
                m_rMq_y = rmq_succinct_ct<false>(max);
            }
        }

        //! Copy constructor
        log_reference(const log_reference& o)
        {
            copy(o);
        }

        //! Move constructor
        log_reference(log_reference&& o)
        {
            *this = std::move(o);
        }


        log_reference &operator=(const log_reference &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        log_reference &operator=(log_reference &&o) {
            if (this != &o) {
                m_x_p = std::move(o.m_x_p);
                m_select_x_p = std::move(o.m_select_x_p);
                m_select_x_p.set_vector(&m_x_p);
                m_rank_x_p = std::move(o.m_rank_x_p);
                m_rank_x_p.set_vector(&m_x_p);
                m_x_n = std::move(o.m_x_n);
                m_select_x_n = std::move(o.m_select_x_n);
                m_select_x_n.set_vector(&m_x_n);
                m_rank_x_n = std::move(o.m_rank_x_n);
                m_rank_x_n.set_vector(&m_x_n);
                m_y_p = std::move(o.m_y_p);
                m_select_y_p = std::move(o.m_select_y_p);
                m_select_y_p.set_vector(&m_y_p);
                m_rank_y_p = std::move(o.m_rank_y_p);
                m_rank_y_p.set_vector(&m_y_p);
                m_y_n = std::move(o.m_y_n);
                m_select_y_n = std::move(o.m_select_y_n);
                m_select_y_n.set_vector(&m_y_n);
                m_rank_y_n = std::move(o.m_rank_y_n);
                m_rank_y_n.set_vector(&m_y_n);
                m_rMmq_x_sample = std::move(o.m_rMmq_x_sample);
                m_rMmq_y_sample = std::move(o.m_rMmq_y_sample);
                m_rank_rMmq_x_sample = std::move(o.m_rank_rMmq_x_sample);
                m_rank_rMmq_y_sample = std::move(o.m_rank_rMmq_y_sample);
                m_rank_rMmq_x_sample.set_vector(&m_rMmq_x_sample);
                m_rank_rMmq_y_sample.set_vector(&m_rMmq_y_sample);
                m_select_rMmq_x_sample = std::move(o.m_select_rMmq_x_sample);
                m_select_rMmq_y_sample = std::move(o.m_select_rMmq_y_sample);
                m_select_rMmq_x_sample.set_vector(&m_rMmq_x_sample);
                m_select_rMmq_y_sample.set_vector(&m_rMmq_y_sample);
                m_rmq_x = std::move(o.m_rmq_x);
                m_rMq_x = std::move(o.m_rMq_x);
                m_rmq_y = std::move(o.m_rmq_y);
                m_rMq_y = std::move(o.m_rMq_y);
            }
            return *this;
        }

        void swap(log_reference &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            m_x_p.swap(o.m_x_p);
            sdsl::util::swap_support(m_rank_x_p, o.m_rank_x_p, &m_x_p, &o.m_x_p);
            sdsl::util::swap_support(m_select_x_p, o.m_select_x_p, &m_x_p, &o.m_x_p);
            m_x_n.swap(o.m_x_n);
            sdsl::util::swap_support(m_rank_x_n, o.m_rank_x_n, &m_x_n, &o.m_x_n);
            sdsl::util::swap_support(m_select_x_n, o.m_select_x_n, &m_x_n, &o.m_x_n);
            m_y_p.swap(o.m_y_p);
            sdsl::util::swap_support(m_rank_y_p, o.m_rank_y_p, &m_y_p, &o.m_y_p);
            sdsl::util::swap_support(m_select_y_p, o.m_select_y_p, &m_y_p, &o.m_y_p);
            m_y_n.swap(o.m_y_n);
            sdsl::util::swap_support(m_rank_y_n, o.m_rank_y_n, &m_y_n, &o.m_y_n);
            sdsl::util::swap_support(m_select_y_n, o.m_select_y_n, &m_y_n, &o.m_y_n);
            m_rMmq_x_sample.swap(o.m_rMmq_x_sample);
            m_rMmq_y_sample.swap(o.m_rMmq_y_sample);
            sdsl::util::swap_support(m_rank_rMmq_x_sample, o.m_rank_rMmq_x_sample, &m_rMmq_x_sample, &o.m_rMmq_x_sample);
            sdsl::util::swap_support(m_rank_rMmq_y_sample, o.m_rank_rMmq_y_sample, &m_rMmq_y_sample, &o.m_rMmq_y_sample);
            sdsl::util::swap_support(m_select_rMmq_x_sample, o.m_select_rMmq_x_sample, &m_rMmq_x_sample, &o.m_rMmq_x_sample);
            sdsl::util::swap_support(m_select_rMmq_y_sample, o.m_select_rMmq_y_sample, &m_rMmq_y_sample, &o.m_rMmq_y_sample);
            m_rmq_x.swap(o.m_rmq_x);
            m_rmq_y.swap(o.m_rmq_y);
            m_rMq_x.swap(o.m_rMq_x);
            m_rMq_y.swap(o.m_rMq_y);
        }

        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_x_p.serialize(out, child, "x_p");
            written_bytes += m_rank_x_p.serialize(out, child, "rank_x_p");
            written_bytes += m_select_x_p.serialize(out, child, "select_x_p");
            written_bytes += m_x_n.serialize(out, child, "x_n");
            written_bytes += m_rank_x_n.serialize(out, child, "rank_x_n");
            written_bytes += m_select_x_n.serialize(out, child, "select_x_n");
            written_bytes += m_y_p.serialize(out, child, "y_p");
            written_bytes += m_rank_y_p.serialize(out, child, "rank_y_p");
            written_bytes += m_select_y_p.serialize(out, child, "select_y_p");
            written_bytes += m_y_n.serialize(out, child, "y_n");
            written_bytes += m_rank_y_n.serialize(out, child, "rank_y_n");
            written_bytes += m_select_y_n.serialize(out, child, "select_y_n");
            written_bytes += m_rMmq_x_sample.serialize(out, child, "rMmq_x_sample");
            written_bytes += m_rank_rMmq_x_sample.serialize(out, child, "rank_rMmq_x_sample");
            written_bytes += m_select_rMmq_x_sample.serialize(out, child, "select_rMmq_x_sample");
            written_bytes += m_rMmq_y_sample.serialize(out, child, "rMmq_y_sample");
            written_bytes += m_rank_rMmq_y_sample.serialize(out, child, "rank_rMmq_y_sample");
            written_bytes += m_select_rMmq_y_sample.serialize(out, child, "select_rMmq_y_sample");
            written_bytes += m_rmq_x.serialize(out, child, "rmq_x");
            written_bytes += m_rmq_y.serialize(out, child, "rmq_y");
            written_bytes += m_rMq_x.serialize(out, child, "rMq_x");
            written_bytes += m_rMq_y.serialize(out, child, "rMq_y");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

    };
}

#endif //RCT_LOG_REFERENCE_HPP
