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

#ifndef RCT_LOG_REFERENCE_SPARSE_NO_MIN_MAX_HPP
#define RCT_LOG_REFERENCE_SPARSE_NO_MIN_MAX_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <rmq_succinct_ct.hpp>
#include <sdsl/sd_vector.hpp>
#include <succ_support_sd.hpp>
#include <geo_util.hpp>
#include <queue>

namespace rct {

    template <class t_movements = sdsl::sd_vector<>>
    class log_reference_sparse_no_min_max {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type ;
        typedef t_movements move_type;
        typedef sdsl::succ_support_sd<1> succ_move_type;
        typedef typename move_type::select_1_type select_move_type;
        typedef sdsl::sd_vector<> sampling_type;
        typedef typename sampling_type::rank_1_type rank_sampling_type;
        typedef typename sampling_type::select_1_type select_sampling_type;

    private:

        move_type m_x_p;
        succ_move_type m_succ_x_p;
        select_move_type m_select_x_p;

        move_type m_y_p;
        succ_move_type m_succ_y_p;
        select_move_type m_select_y_p;

        move_type m_x_n;
        succ_move_type m_succ_x_n;
        select_move_type m_select_x_n;

        move_type m_y_n;
        succ_move_type m_succ_y_n;
        select_move_type m_select_y_n;


        void copy(const log_reference_sparse_no_min_max &o){
            m_x_p = o.m_x_p;
            m_select_x_p = o.m_select_x_p;
            m_select_x_p.set_vector(&m_x_p);
            m_succ_x_p = o.m_succ_x_p;
            m_succ_x_p.set_vector(&m_x_p);
            m_x_n = o.m_x_n;
            m_select_x_n = o.m_select_x_n;
            m_select_x_n.set_vector(&m_x_n);
            m_succ_x_n = o.m_succ_x_n;
            m_succ_x_n.set_vector(&m_x_n);
            m_y_p = o.m_y_p;
            m_select_y_p = o.m_select_y_p;
            m_select_y_p.set_vector(&m_y_p);
            m_succ_y_p = o.m_succ_y_p;
            m_succ_y_p.set_vector(&m_y_p);
            m_y_n = o.m_y_n;
            m_select_y_n = o.m_select_y_n;
            m_select_y_n.set_vector(&m_y_n);
            m_succ_y_n = o.m_succ_y_n;
            m_succ_y_n.set_vector(&m_y_n);
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


    public:

        log_reference_sparse_no_min_max() = default;

        log_reference_sparse_no_min_max(log_reference<> &ref){
            m_x_p = t_movements(ref.x_p);
            m_x_n = t_movements(ref.x_n);
            m_y_p = t_movements(ref.y_p);
            m_y_n = t_movements(ref.y_n);
            sdsl::util::init_support(m_select_x_p, &m_x_p);
            sdsl::util::init_support(m_succ_x_p, &m_x_p);
            sdsl::util::init_support(m_select_x_n, &m_x_n);
            sdsl::util::init_support(m_succ_x_n, &m_x_n);
            sdsl::util::init_support(m_select_y_p, &m_y_p);
            sdsl::util::init_support(m_succ_y_p, &m_y_p);
            sdsl::util::init_support(m_select_y_n, &m_y_n);
            sdsl::util::init_support(m_succ_y_n, &m_y_n);
        }

        template <class ContainerMovements>
        log_reference_sparse_no_min_max(const ContainerMovements &movements){

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
                sdsl::util::init_support(m_succ_x_p, &m_x_p);
                sdsl::util::init_support(m_select_x_n, &m_x_n);
                sdsl::util::init_support(m_succ_x_n, &m_x_n);
                sdsl::util::init_support(m_select_y_p, &m_y_p);
                sdsl::util::init_support(m_succ_y_p, &m_y_p);
                sdsl::util::init_support(m_select_y_n, &m_y_n);
                sdsl::util::init_support(m_succ_y_n, &m_y_n);
            }

        }

        util::geo::movement compute_movement(const size_type i, const size_type j) const {
            if(i == 0){
                int32_t inc_x = m_select_x_p(j) - m_select_x_n(j);
                int32_t inc_y = m_select_y_p(j) - m_select_y_n(j);
                return util::geo::movement{inc_x, inc_y};
            }else{
                int32_t inc_x = (m_select_x_p(j) - m_select_x_p(i)) - (m_select_x_n(j) - m_select_x_n(i));
                int32_t inc_y = (m_select_y_p(j) - m_select_y_p(i)) - (m_select_y_n(j) - m_select_y_n(i));
                return util::geo::movement{inc_x, inc_y};
            }

        };

        util::geo::movement compute_movement_init(const size_type i, const size_type j, size_type &x_p_prev,
                size_type &x_n_prev, size_type &y_p_prev, size_type &y_n_prev) const {

            x_p_prev = m_select_x_p(j);
            x_n_prev = m_select_x_n(j);
            y_p_prev = m_select_y_p(j);
            y_n_prev = m_select_y_n(j);
            if(i == 0){

                auto inc_x = (int32_t) (x_p_prev - x_n_prev);
                auto inc_y = (int32_t) (y_p_prev - y_n_prev);
                return util::geo::movement{inc_x, inc_y};
            }else{
                int32_t inc_x = (x_p_prev - m_select_x_p(i)) - (x_n_prev - m_select_x_n(i));
                int32_t inc_y = (y_p_prev - m_select_y_p(i)) - (y_n_prev - m_select_y_n(i));
                return util::geo::movement{inc_x, inc_y};
            }
        };

        util::geo::movement compute_movement_init_next(const size_type j, size_type &x_p_prev, size_type &x_n_prev,
                size_type &y_p_prev, size_type &y_n_prev) const {

            x_p_prev = m_select_x_p(j+1);
            x_n_prev = m_select_x_n(j+1);
            y_p_prev = m_select_y_p(j+1);
            y_n_prev = m_select_y_n(j+1);
            if(j == 0){

                auto inc_x = (int32_t) (x_p_prev - x_n_prev);
                auto inc_y = (int32_t) (y_p_prev - y_n_prev);
                return util::geo::movement{inc_x, inc_y};
            }else{
                int32_t inc_x = (x_p_prev - m_select_x_p(j)) - (x_n_prev - m_select_x_n(j));
                int32_t inc_y = (y_p_prev - m_select_y_p(j)) - (y_n_prev - m_select_y_n(j));
                return util::geo::movement{inc_x, inc_y};
            }

        };

        util::geo::movement compute_movement_next(size_type &x_p_prev,
                                                  size_type &x_n_prev, size_type &y_p_prev, size_type &y_n_prev){

            auto x_p = m_succ_x_p(x_p_prev+1);
            auto x_n = m_succ_x_n(x_n_prev+1);
            auto y_p = m_succ_y_p(y_p_prev+1);
            auto y_n = m_succ_y_n(y_n_prev+1);
            auto inc_x = (int32_t) ((x_p - x_p_prev) - (x_n - x_n_prev));
            auto inc_y = (int32_t) ((y_p - y_p_prev) - (y_n - y_n_prev));
            x_p_prev = x_p;
            x_n_prev = x_n;
            y_p_prev = y_p;
            y_n_prev = y_n;
            return util::geo::movement{inc_x, inc_y};
        };

        inline int32_t compute_delta_x(const size_type i) const {
             if(!i) return 0;
             return m_select_x_p(i) - m_select_x_n(i);
        }

        inline int32_t compute_delta_y(const size_type i) const {
            if(!i) return 0;
            return m_select_y_p(i) - m_select_y_n(i);
        }


        util::geo::movement compute_movement_deltas(const size_type i, const int32_t delta_x, const int32_t delta_y) const {

            if(!i) return util::geo::movement{0, 0};
            auto inc_x = (int32_t) (m_select_x_p(i) - m_select_x_n(i) - delta_x);
            auto inc_y = (int32_t) (m_select_y_p(i) - m_select_y_n(i) - delta_y);
            return util::geo::movement{inc_x, inc_y};
        };

        util::geo::movement compute_movement_deltas(const size_type i, const int32_t delta_x, const int32_t delta_y,
                size_type &sel_x_p, size_type &sel_x_n, size_type &sel_y_p, size_type &sel_y_n) const {

            if(!i) return util::geo::movement{0, 0};
            sel_x_p = m_select_x_p(i);
            sel_x_n = m_select_x_n(i);
            sel_y_p = m_select_y_p(i);
            sel_y_n = m_select_y_n(i);
            auto inc_x = (int32_t) (sel_x_p - sel_x_n - delta_x);
            auto inc_y = (int32_t) (sel_y_p - sel_y_n - delta_y);
            return util::geo::movement{inc_x, inc_y};
        };



        //! Copy constructor
        log_reference_sparse_no_min_max(const log_reference_sparse_no_min_max& o)
        {
            copy(o);
        }

        //! Move constructor
        log_reference_sparse_no_min_max(log_reference_sparse_no_min_max&& o)
        {
            *this = std::move(o);
        }


        log_reference_sparse_no_min_max &operator=(const log_reference_sparse_no_min_max &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        log_reference_sparse_no_min_max &operator=(log_reference_sparse_no_min_max &&o) {
            if (this != &o) {
                m_x_p = std::move(o.m_x_p);
                m_select_x_p = std::move(o.m_select_x_p);
                m_select_x_p.set_vector(&m_x_p);
                m_succ_x_p = std::move(o.m_succ_x_p);
                m_succ_x_p.set_vector(&m_x_p);
                m_x_n = std::move(o.m_x_n);
                m_select_x_n = std::move(o.m_select_x_n);
                m_select_x_n.set_vector(&m_x_n);
                m_succ_x_n = std::move(o.m_succ_x_n);
                m_succ_x_n.set_vector(&m_x_n);
                m_y_p = std::move(o.m_y_p);
                m_select_y_p = std::move(o.m_select_y_p);
                m_select_y_p.set_vector(&m_y_p);
                m_succ_y_p = std::move(o.m_succ_y_p);
                m_succ_y_p.set_vector(&m_y_p);
                m_y_n = std::move(o.m_y_n);
                m_select_y_n = std::move(o.m_select_y_n);
                m_select_y_n.set_vector(&m_y_n);
                m_succ_y_n = std::move(o.m_succ_y_n);
                m_succ_y_n.set_vector(&m_y_n);
            }
            return *this;
        }

        void swap(log_reference_sparse_no_min_max &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            m_x_p.swap(o.m_x_p);
            sdsl::util::swap_support(m_succ_x_p, o.m_succ_x_p, &m_x_p, &o.m_x_p);
            sdsl::util::swap_support(m_select_x_p, o.m_select_x_p, &m_x_p, &o.m_x_p);
            m_x_n.swap(o.m_x_n);
            sdsl::util::swap_support(m_succ_x_n, o.m_succ_x_n, &m_x_n, &o.m_x_n);
            sdsl::util::swap_support(m_select_x_n, o.m_select_x_n, &m_x_n, &o.m_x_n);
            m_y_p.swap(o.m_y_p);
            sdsl::util::swap_support(m_succ_y_p, o.m_succ_y_p, &m_y_p, &o.m_y_p);
            sdsl::util::swap_support(m_select_y_p, o.m_select_y_p, &m_y_p, &o.m_y_p);
            m_y_n.swap(o.m_y_n);
            sdsl::util::swap_support(m_succ_y_n, o.m_succ_y_n, &m_y_n, &o.m_y_n);
            sdsl::util::swap_support(m_select_y_n, o.m_select_y_n, &m_y_n, &o.m_y_n);
        }

        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_x_p.serialize(out, child, "x_p");
            written_bytes += m_succ_x_p.serialize(out, child, "succ_x_p");
            written_bytes += m_select_x_p.serialize(out, child, "select_x_p");
            written_bytes += m_x_n.serialize(out, child, "x_n");
            written_bytes += m_succ_x_n.serialize(out, child, "succ_x_n");
            written_bytes += m_select_x_n.serialize(out, child, "select_x_n");
            written_bytes += m_y_p.serialize(out, child, "y_p");
            written_bytes += m_succ_y_p.serialize(out, child, "succ_y_p");
            written_bytes += m_select_y_p.serialize(out, child, "select_y_p");
            written_bytes += m_y_n.serialize(out, child, "y_n");
            written_bytes += m_succ_y_n.serialize(out, child, "succ_y_n");
            written_bytes += m_select_y_n.serialize(out, child, "select_y_n");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        size_type new_space(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_x_p.serialize(out, child, "x_p");
            written_bytes += m_succ_x_p.serialize(out, child, "succ_x_p");
            written_bytes += m_select_x_p.serialize(out, child, "select_x_p");
            written_bytes += m_x_n.serialize(out, child, "x_n");
            written_bytes += m_succ_x_n.serialize(out, child, "succ_x_n");
            written_bytes += m_select_x_n.serialize(out, child, "select_x_n");
            written_bytes += m_y_p.serialize(out, child, "y_p");
            written_bytes += m_succ_y_p.serialize(out, child, "succ_y_p");
            written_bytes += m_select_y_p.serialize(out, child, "select_y_p");
            written_bytes += m_y_n.serialize(out, child, "y_n");
            written_bytes += m_succ_y_n.serialize(out, child, "succ_y_n");
            written_bytes += m_select_y_n.serialize(out, child, "select_y_n");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            m_x_p.load(in);
            m_succ_x_p.load(in, &m_x_p);
            m_select_x_p.load(in, &m_x_p);
            m_x_n.load(in);
            m_succ_x_n.load(in, &m_x_n);
            m_select_x_n.load(in, &m_x_n);
            m_y_p.load(in);
            m_succ_y_p.load(in, &m_y_p);
            m_select_y_p.load(in, &m_y_p);
            m_y_n.load(in);
            m_succ_y_n.load(in, &m_y_n);
            m_select_y_n.load(in, &m_y_n);

        }

    };
}

#endif //RCT_LOG_REFERENCE_HPP
