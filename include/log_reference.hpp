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
#include <succ_support_v.hpp>
#include <geo_util.hpp>
#include <queue>

namespace rct {

    template <class t_movements = sdsl::bit_vector>
    class log_reference {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type ;
        typedef t_movements move_type;
        typedef succ_support_v<1> succ_move_type;
        typedef typename move_type::select_1_type select_move_type;
        typedef sdsl::sd_vector<> sampling_type;
        typedef typename sampling_type::rank_1_type rank_sampling_type;
        typedef typename sampling_type::select_1_type select_sampling_type;

    private:

        typedef struct {
            util::geo::point p_s;
            util::geo::point p_e;
            util::geo::region parent_region;
            size_type beg;
            size_type end;
            size_type sel_p_x;
            size_type sel_n_x;
            size_type sel_p_y;
            size_type sel_n_y;
            size_type min_x_index;
            size_type max_x_index;
            size_type min_y_index;
            size_type max_y_index;
        } mbr_q_type;

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
                sdsl::util::init_support(m_succ_x_p, &m_x_p);
                sdsl::util::init_support(m_select_x_n, &m_x_n);
                sdsl::util::init_support(m_succ_x_n, &m_x_n);
                sdsl::util::init_support(m_select_y_p, &m_y_p);
                sdsl::util::init_support(m_succ_y_p, &m_y_p);
                sdsl::util::init_support(m_select_y_n, &m_y_n);
                sdsl::util::init_support(m_succ_y_n, &m_y_n);
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
                                                  size_type &x_n_prev, size_type &y_p_prev, size_type &y_n_prev) const {

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



        //Time Interval Functions

        inline size_type _rank_min(const size_type rmMq) const{
            return rmMq>>1; //rmMq/2
        };

        inline size_type _rank_max(const size_type rmMq, const sampling_type &rmq_sampling) const{
            return (rmMq>>1) + (rmMq&0x1ULL) - rmq_sampling[0]; //rmMq/2 + rmMq%2 -mM[0]
        };

        inline size_type _select_min(const size_type i, const select_sampling_type &rmq_sampling_select) const{
            return rmq_sampling_select(i<<1) -1;
        }

        inline size_type _select_max(const size_type i, const sampling_type &rmq_sampling,
                                     const select_sampling_type &rmq_sampling_select) const{
            return (rmq_sampling[0]) ? rmq_sampling_select((i<<1)+1)-1 : rmq_sampling_select((i<<1)-1)-1;
            /*if(rmq_sampling[0]){
                return rmq_sampling_select((i<<1)+1)-1;
            }else{
                return rmq_sampling_select((i<<1)-1)-1;
            }*/
        }

        util::geo::region find_MBR(const size_type x, const size_type y, const size_type move_s, const size_type move_e,
                                   const util::geo::point &p_s, const util::geo::point &p_e, const int32_t delta_x,
                                   const int32_t  delta_y, size_type &total_selects) const {

            //1. Init the min and max points
            value_type min_x, max_x, min_y, max_y;
            if(p_s.x < p_e.x){
                min_x = p_s.x;
                max_x = p_e.x;
            }else{
                min_x = p_e.x;
                max_x = p_s.x;
            }
            if(p_s.y < p_e.y){
                min_y = p_s.y;
                max_y = p_e.y;
            }else{
                min_y = p_e.y;
                max_y = p_s.y;
            }

            //Adding extra 1 caused by the flag
            size_type rMmq_x_s = m_rank_rMmq_x_sample(move_s+1);
            size_type rMmq_x_e = m_rank_rMmq_x_sample(move_e+2);
            size_type rMmq_y_s = m_rank_rMmq_y_sample(move_s+1);
            size_type rMmq_y_e = m_rank_rMmq_y_sample(move_e+2);

            //2. Computing the minimum at x-axis
            size_type rmq_x_s = _rank_min(rMmq_x_s)+1;  //rank_x(ms)+1 counts the number of samples before
            size_type rmq_x_e = _rank_min(rMmq_x_e); //rank_x(me+1) counts the number of samples
            if(rmq_x_s > 0 && rmq_x_e > 0 && rmq_x_s <= rmq_x_e){
                /*r = rmq_x(rmq_x_s-1, rmq_x_e-1) returns the sampled position of the minimum
                 *p = select_rmq_x_sample(r+1) returns the position of the movement (we add 1 because it is the (p+1)-th movement)
                 */
                size_type min_x_c = m_rmq_x(rmq_x_s-1, rmq_x_e-1)+1;
                size_type movement_x = _select_min(min_x_c, m_select_rMmq_x_sample)+1;
                auto min_x_aux = (value_type) (compute_delta_x(movement_x) - delta_x + x);
                //total_selects += 2;
                min_x = std::min(min_x, min_x_aux);
            }

            //3. Computing the minimum at y-axis
            size_type rmq_y_s = _rank_min(rMmq_y_s)+1;
            size_type rmq_y_e = _rank_min(rMmq_y_e);
            if(rmq_y_s > 0 && rmq_y_e > 0 && rmq_y_s <= rmq_y_e){
                size_type min_y_c = m_rmq_y(rmq_y_s-1, rmq_y_e-1) + 1;
                size_type movement_y = _select_min(min_y_c, m_select_rMmq_y_sample) + 1;
                auto min_y_aux = (value_type) (compute_delta_y(movement_y) - delta_y + y);
                //total_selects += 2;
                min_y = std::min(min_y, min_y_aux);
            }

            //4. Computing the maximum at x-axis
            size_type rMq_x_s = _rank_max(rMmq_x_s, m_rMmq_x_sample)+1;
            size_type rMq_x_e = _rank_max(rMmq_x_e, m_rMmq_x_sample);
            if(rMq_x_s > 0 && rMq_x_e > 0 && rMq_x_s <= rMq_x_e){
                size_type max_x_c = m_rMq_x(rMq_x_s-1, rMq_x_e-1) + 1;
                size_type movement_x = _select_max(max_x_c, m_rMmq_x_sample, m_select_rMmq_x_sample) + 1;
                auto max_x_aux = (value_type) (compute_delta_x(movement_x) - delta_x + x);
                //total_selects += 2;
                max_x_aux += x;
                max_x = std::max(max_x, max_x_aux);
            }

            //5. Computing the maximum at y-axis
            size_type rMq_y_s = _rank_max(rMmq_y_s, m_rMmq_y_sample)+1;
            size_type rMq_y_e = _rank_max(rMmq_y_e, m_rMmq_y_sample);
            if(rMq_y_s > 0 && rMq_y_e > 0 && rMq_y_s <= rMq_y_e){
                size_type max_y_c = m_rMq_y(rMq_y_s-1, rMq_y_e-1) + 1;
                size_type movement_y = _select_max(max_y_c, m_rMmq_y_sample, m_select_rMmq_y_sample) + 1;
                auto max_y_aux = (value_type) (compute_delta_y(movement_y) - delta_y + y);
                //total_selects += 2;
                max_y = std::max(max_y, max_y_aux);
            }

            return util::geo::region{ util::geo::point{min_x, min_y}, util::geo::point{max_x, max_y}};

        }

        bool contains_region_init(const size_type x, const size_type y, const size_type move_s, const size_type move_e,
                                  const util::geo::point &p_s, const util::geo::point &p_e, const util::geo::region &r_q,
                                  const int32_t delta_x, const int32_t delta_y,
                                  size_type &p_min_x_index, size_type &p_max_x_index,
                                  size_type &p_min_y_index, size_type &p_max_y_index,
                                  util::geo::region &mbr) const {

            //1. Init the min and max points
            value_type min_x, max_x, min_y, max_y;
            if(p_s.x < p_e.x){
                min_x = p_s.x;
                max_x = p_e.x;
            }else{
                min_x = p_e.x;
                max_x = p_s.x;
            }
            if(p_s.y < p_e.y){
                min_y = p_s.y;
                max_y = p_e.y;
            }else{
                min_y = p_e.y;
                max_y = p_s.y;
            }

            p_min_x_index = 0, p_max_x_index = 0, p_min_y_index = 0, p_max_y_index = 0;
            //Adding extra 1 caused by the flag bit
            size_type rMmq_x_s = m_rank_rMmq_x_sample(move_s-1);
            size_type rMmq_x_e = m_rank_rMmq_x_sample(move_e);//move_e + 1 + flag

            //2. Computing the minimum at x-axis
            size_type rmq_x_s = _rank_min(rMmq_x_s)+1;  //rank_x(ms)+1 counts the number of samples before
            size_type rmq_x_e = _rank_min(rMmq_x_e); //rank_x(me+1) counts the number of samples
            if(rmq_x_s > 0 && rmq_x_e > 0 && rmq_x_s <= rmq_x_e){
                /*r = rmq_x(rmq_x_s-1, rmq_x_e-1) returns the sampled position of the minimum
                 *p = select_rmq_x_sample(r+1) returns the position of the movement (we add 1 because it is the (p+1)-th movement)
                 */
                size_type min_x_c = m_rmq_x(rmq_x_s-1, rmq_x_e-1)+1;
                p_min_x_index =  _select_min(min_x_c, m_select_rMmq_x_sample);
                auto min_x_aux = (value_type) (compute_delta_x(p_min_x_index+1) - delta_x + x);
                //total_selects += 2;
                min_x = std::min(min_x, min_x_aux);
            }

            //4. Computing the maximum at x-axis
            size_type rMq_x_s = _rank_max(rMmq_x_s, m_rMmq_x_sample)+1;
            size_type rMq_x_e = _rank_max(rMmq_x_e, m_rMmq_x_sample);
            if(rMq_x_s > 0 && rMq_x_e > 0 && rMq_x_s <= rMq_x_e){
                size_type max_x_c = m_rMq_x(rMq_x_s-1, rMq_x_e-1) + 1;
                p_max_x_index = _select_max(max_x_c, m_rMmq_x_sample, m_select_rMmq_x_sample);
                auto max_x_aux = (value_type) (compute_delta_x(p_max_x_index+1) - delta_x + x);
                //total_selects += 2;
                max_x = std::max(max_x, max_x_aux);
            }

            if(min_x > r_q.max.x || max_x < r_q.min.x) return false;

            size_type rMmq_y_s = m_rank_rMmq_y_sample(move_s-1);
            size_type rMmq_y_e = m_rank_rMmq_y_sample(move_e);

            //3. Computing the minimum at y-axis
            size_type rmq_y_s = _rank_min(rMmq_y_s)+1;
            size_type rmq_y_e = _rank_min(rMmq_y_e);
            if(rmq_y_s > 0 && rmq_y_e > 0 && rmq_y_s <= rmq_y_e){
                size_type min_y_c = m_rmq_y(rmq_y_s-1, rmq_y_e-1) + 1;
                p_min_y_index = _select_min(min_y_c, m_select_rMmq_y_sample);
                auto min_y_aux = (value_type) (compute_delta_y(p_min_y_index+1) - delta_y + y);
                //total_selects += 2;
                min_y = std::min(min_y, min_y_aux);
            }

            //5. Computing the maximum at y-axis
            size_type rMq_y_s = _rank_max(rMmq_y_s, m_rMmq_y_sample)+1;
            size_type rMq_y_e = _rank_max(rMmq_y_e, m_rMmq_y_sample);
            if(rMq_y_s > 0 && rMq_y_e > 0 && rMq_y_s <= rMq_y_e){
                size_type max_y_c = m_rMq_y(rMq_y_s-1, rMq_y_e-1) + 1;
                p_max_y_index = _select_max(max_y_c, m_rMmq_y_sample, m_select_rMmq_y_sample);
                auto max_y_aux = (value_type) (compute_delta_y(p_max_y_index+1) - delta_y + y);
                //total_selects += 2;
                max_y = std::max(max_y, max_y_aux);
            }
            if(min_y > r_q.max.y || max_y < r_q.min.y) return false;

            mbr = {util::geo::point{min_x, min_y}, util::geo::point{max_x, max_y}};
            return true;

        }

        bool contains_region_lazy(const size_type x, const size_type y,
                           const size_type move_s, const size_type move_e,
                           const util::geo::point &p_s, const util::geo::point &p_e, const util::geo::region r_q,
                           const util::geo::region &parent_region, util::geo::region &r,
                           const int32_t delta_x, const int32_t delta_y,
                           size_type &p_min_x_index, size_type &p_max_x_index,
                           size_type &p_min_y_index, size_type &p_max_y_index) const {



            //1. Init the min and max points
            value_type min_x, max_x, min_y, max_y;
            if(p_s.x < p_e.x){
                min_x = p_s.x;
                max_x = p_e.x;
            }else{
                min_x = p_e.x;
                max_x = p_s.x;
            }
            if(p_s.y < p_e.y){
                min_y = p_s.y;
                max_y = p_e.y;
            }else{
                min_y = p_e.y;
                max_y = p_s.y;
            }

            //Check if the local maximum and minimum was computed, during the computation of the parent's MBR
            if(move_s-1 < p_min_x_index && move_e-1 > p_min_x_index && move_s-1 < p_max_x_index && move_e-1 > p_max_x_index){
                r = parent_region;
            }else{
                //Adding extra 1 caused by the flag bit
                size_type rMmq_x_s = m_rank_rMmq_x_sample(move_s-1);
                size_type rMmq_x_e = m_rank_rMmq_x_sample(move_e);//move_e + 1 + flag
                //2. Computing the minimum at x-axis
                if(move_s-1 < p_min_x_index && move_e-1 > p_min_x_index){
                    r.min.x = parent_region.min.x;
                }else{
                    size_type rmq_x_s = _rank_min(rMmq_x_s)+1;  //rank_x(ms)+1 counts the number of samples before
                    size_type rmq_x_e = _rank_min(rMmq_x_e); //rank_x(me+1) counts the number of samples
                    if(rmq_x_s > 0 && rmq_x_e > 0 && rmq_x_s <= rmq_x_e){
                        /*r = rmq_x(rmq_x_s-1, rmq_x_e-1) returns the sampled position of the minimum
                         *p = select_rmq_x_sample(r+1) returns the position of the movement (we add 1 because it is the (p+1)-th movement)
                         */
                        size_type min_x_c = m_rmq_x(rmq_x_s-1, rmq_x_e-1)+1;
                        p_min_x_index = _select_min(min_x_c, m_select_rMmq_x_sample);
                        auto min_x_aux = (value_type) (compute_delta_x(p_min_x_index+1) - delta_x + x);
                        min_x = std::min(min_x, min_x_aux);
                    }
                    r.min.x = min_x;
                }

                if(move_s-1 < p_max_x_index && move_e-1 > p_max_x_index){
                    r.max.x = parent_region.max.x;
                }else{
                    //4. Computing the maximum at x-axis
                    size_type rMq_x_s = _rank_max(rMmq_x_s, m_rMmq_x_sample)+1;
                    size_type rMq_x_e = _rank_max(rMmq_x_e, m_rMmq_x_sample);
                    if(rMq_x_s > 0 && rMq_x_e > 0 && rMq_x_s <= rMq_x_e){
                        size_type max_x_c = m_rMq_x(rMq_x_s-1, rMq_x_e-1) + 1;
                        p_max_x_index = _select_max(max_x_c, m_rMmq_x_sample, m_select_rMmq_x_sample);
                        auto max_x_aux = (value_type) (compute_delta_x(p_max_x_index+1) - delta_x + x);
                        max_x = std::max(max_x, max_x_aux);
                    }
                    r.max.x = max_x;
                }

            }
            if(r.min.x > r_q.max.x || r.max.x < r_q.min.x) return false;

            if(move_s-1 < p_min_y_index && move_e-1 > p_min_y_index && move_s-1 < p_max_y_index && move_e-1 > p_max_y_index){
                r.min.y = parent_region.min.y;
                r.max.y = parent_region.max.y;
            }else{
                size_type rMmq_y_s = m_rank_rMmq_y_sample(move_s-1);
                size_type rMmq_y_e = m_rank_rMmq_y_sample(move_e);//move_e + 1 + flag
                if(move_s-1 < p_min_y_index && move_e-1 > p_min_y_index){
                    r.min.y = parent_region.min.y;
                }else{
                    //3. Computing the minimum at y-axis
                    size_type rmq_y_s = _rank_min(rMmq_y_s)+1;
                    size_type rmq_y_e = _rank_min(rMmq_y_e);
                    if(rmq_y_s > 0 && rmq_y_e > 0 && rmq_y_s <= rmq_y_e){
                        size_type min_y_c = m_rmq_y(rmq_y_s-1, rmq_y_e-1) + 1;
                        p_min_y_index = _select_min(min_y_c, m_select_rMmq_y_sample);
                        auto min_y_aux = (value_type) (compute_delta_y(p_min_y_index+1) - delta_y + y);
                        min_y = std::min(min_y, min_y_aux);
                    }
                    r.min.y = min_y;
                }

                if(move_s-1 < p_max_y_index && move_e-1 > p_max_y_index){
                    r.max.y = parent_region.max.y;
                }else{
                    size_type rMq_y_s = _rank_max(rMmq_y_s, m_rMmq_y_sample)+1;
                    size_type rMq_y_e = _rank_max(rMmq_y_e, m_rMmq_y_sample);
                    if(rMq_y_s > 0 && rMq_y_e > 0 && rMq_y_s <= rMq_y_e){
                        size_type max_y_c = m_rMq_y(rMq_y_s-1, rMq_y_e-1) + 1;
                        p_max_y_index = _select_max(max_y_c, m_rMmq_y_sample, m_select_rMmq_y_sample);
                        auto max_y_aux = (value_type) (compute_delta_y(p_max_y_index+1) - delta_y + y);
                        max_y = std::max(max_y, max_y_aux);
                    }
                    r.max.y = max_y;
                }

            }
            return !(r.min.y > r_q.max.y || r.max.y < r_q.min.y);

        }


        template<size_type leaf_width = 20 >
        bool contains_region_old(const size_type x, const size_type y, const size_type move_s, const size_type move_e,
                             const util::geo::point &p_s, const util::geo::point &p_e, const util::geo::region &r_q,
                             const int32_t delta_x, const int32_t delta_y) const{

            if(util::geo::contains(r_q, p_s) || util::geo::contains(r_q, p_e)) {
                return true;
            }
            util::geo::region mbr;
            std::queue<mbr_q_type> mbrs;
            size_type p_min_x_index = 0, p_max_x_index = 0, p_min_y_index = 0, p_max_y_index = 0;
            size_type sel_p_x, sel_n_x, sel_p_y, sel_n_y;
            mbrs.push(mbr_q_type{p_s, p_e, util::geo::region(), move_s-1, move_e-1, sel_p_x, sel_n_x, sel_p_y, sel_n_y,
                                 p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index});
            size_type i = 0;
            bool intersects;
            while(!mbrs.empty()){
                mbr_q_type q_e = mbrs.front();
                mbrs.pop();
                if(i == 0){
                    intersects = contains_region_init(x, y, move_s, move_e, p_s, p_e, r_q, delta_x, delta_y,
                                                     p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index, mbr);
                }else{
                    p_min_x_index = q_e.min_x_index;
                    p_max_x_index = q_e.max_x_index;
                    p_min_y_index = q_e.min_y_index;
                    p_max_y_index = q_e.max_y_index;
                    intersects = contains_region_lazy(x, y, q_e.beg, q_e.end, q_e.p_s, q_e.p_e, r_q, q_e.parent_region,
                                 mbr, delta_x, delta_y, p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index);

                }

                if(intersects) {

                    if(util::geo::contains(r_q, mbr)){
                        return true;
                    }
                    if((q_e.beg - q_e.end) > leaf_width){
                        size_type move_m_count = (q_e.beg - q_e.end+1) / 2 + q_e.beg + 1; //index -1
                        util::geo::point p_m = compute_movement_init(q_e.beg, q_e.end, sel_p_x, sel_n_x, sel_p_y, sel_n_y);
                        p_m.x += x;
                        p_m.y += y;
#if VERBOSE_AGB
                        std::cout << "m_index = " << move_m_count << std::endl;
                            std::cout << "p_m: (" << p_m.m_x << ", " << p_m.m_y << ")" << std::endl;
#endif
                        if(util::geo::contains(r_q, p_m)) {
#if VERBOSE_AGB
                            std::cout << "the r_q contains p_m" << std::endl;
#endif
                            return true;
                        }
                        mbrs.push(mbr_q_type(q_e.p_s, p_m, q_e.beg, move_m_count-1,
                                             q_e.sel_p_x, q_e.sel_n_x, q_e.sel_p_y, q_e.sel_n_y, mbr,
                                             p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index));
                        mbrs.push(mbr_q_type(p_m, q_e.p_e, move_m_count-1, q_e.end, sel_p_x, sel_n_x, sel_p_y,
                                             sel_n_y,  mbr, p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index));
                    }else{
                        sel_p_x = q_e.sel_p_x;
                        sel_n_x = q_e.sel_n_x;
                        sel_p_y = q_e.sel_p_y;
                        sel_n_y = q_e.sel_n_y;

                        size_type move_index = q_e.beg+1;
                        while(move_index < q_e.end){
#if VERBOSE_AGB
                            std::cout << "(leaf) sel_p_x: " << sel_p_x << " sel_n_x: " << sel_n_x << " sel_p_y: " << sel_p_y << " sel_n_y: " << sel_n_y << std::endl;
                                std::cout << "move_index: " << move_index << std::endl;
#endif
                            util::geo::point p_m = compute_movement_next(sel_p_x, sel_n_x, sel_p_y, sel_n_y);
                            p_m.x += x;
                            p_m.y += y;
#if VERBOSE_AGB
                            std::cout << "p_m leaf: (" << p_m.m_x << ", " << p_m.m_y << ")" << std::endl;
#endif
                            if(util::geo::contains(r_q, p_m)) {
#if VERBOSE_AGB
                                std::cout << "the r_q contains p_m leaf" << std::endl;
#endif
                                return true;
                            }
                            move_index++;
                        }
                    }
                }

                ++i;
            }

            return false;



        }


        template<size_type leaf_width = 20 >
        bool contains_region(const value_type phrase_x, const value_type phrase_y, const size_type phrase_start,
                                const size_type move_s, const size_type move_e, const util::geo::region &r_q) const {

            //1. Compute p_s and p_e
            size_type sel_p_x, sel_n_x, sel_p_y, sel_n_y;
            int32_t delta_x = compute_delta_x(phrase_start);
            int32_t delta_y = compute_delta_y(phrase_start);

            auto m_s = compute_movement_deltas(move_s, delta_x, delta_y, sel_p_x, sel_n_x, sel_p_y, sel_n_y);
            util::geo::point p_s{m_s.x + phrase_x, m_s.y + phrase_y};
            if(util::geo::contains(r_q, p_s)) return true;

            auto m_e = compute_movement_deltas(move_e, delta_x, delta_y);
            util::geo::point p_e{m_e.x + phrase_x, m_e.y + phrase_y};
            if(util::geo::contains(r_q, p_e)) return true;


            util::geo::region mbr;
            std::queue<mbr_q_type> mbrs;
            size_type p_min_x_index = 0, p_max_x_index = 0, p_min_y_index = 0, p_max_y_index = 0;

            //TODO: Check -1 because the first position (-1)
            mbrs.push(mbr_q_type{p_s, p_e, util::geo::region(), move_s, move_e, sel_p_x, sel_n_x, sel_p_y, sel_n_y,
                                 p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index});
            size_type i = 0;
            bool intersects;
            while(!mbrs.empty()){
                mbr_q_type q_e = mbrs.front();
                mbrs.pop();
                if(i == 0){
                    intersects = contains_region_init(phrase_x, phrase_y, move_s, move_e, p_s, p_e, r_q, delta_x, delta_y,
                                                      p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index, mbr);
                }else{
                    p_min_x_index = q_e.min_x_index;
                    p_max_x_index = q_e.max_x_index;
                    p_min_y_index = q_e.min_y_index;
                    p_max_y_index = q_e.max_y_index;
                    intersects = contains_region_lazy(phrase_x, phrase_y, q_e.beg, q_e.end, q_e.p_s, q_e.p_e, r_q, q_e.parent_region,
                                                      mbr, delta_x, delta_y, p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index);

                }

                if(intersects) {

                    if(util::geo::contains(r_q, mbr)){
                        return true;
                    }
                    if((q_e.end - q_e.beg) > leaf_width){
                        size_type move_m_count = (q_e.end - q_e.beg) / 2 + q_e.beg; //count
                        auto movement = compute_movement_deltas(move_m_count, delta_x, delta_y, sel_p_x, sel_n_x, sel_p_y, sel_n_y);
                        util::geo::point p_m {movement.x + phrase_x, movement.y + phrase_y};
#if VERBOSE_AGB
                        std::cout << "m_index = " << move_m_count << std::endl;
                            std::cout << "p_m: (" << p_m.m_x << ", " << p_m.m_y << ")" << std::endl;
#endif
                        if(util::geo::contains(r_q, p_m)) {
#if VERBOSE_AGB
                            std::cout << "the r_q contains p_m" << std::endl;
#endif
                            return true;
                        }
                        mbrs.push(mbr_q_type{q_e.p_s, p_m, mbr, q_e.beg, move_m_count,
                                             q_e.sel_p_x, q_e.sel_n_x, q_e.sel_p_y, q_e.sel_n_y,
                                             p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index});
                        mbrs.push(mbr_q_type{p_m, q_e.p_e, mbr, move_m_count, q_e.end, sel_p_x, sel_n_x, sel_p_y,
                                             sel_n_y, p_min_x_index, p_max_x_index, p_min_y_index, p_max_y_index});
                    }else{
                        sel_p_x = q_e.sel_p_x;
                        sel_n_x = q_e.sel_n_x;
                        sel_p_y = q_e.sel_p_y;
                        sel_n_y = q_e.sel_n_y;

                        size_type move_index = q_e.beg+1;
                        util::geo::point p_i = q_e.p_s;
                        while(move_index <= q_e.end){
#if VERBOSE_AGB
                            std::cout << "(leaf) sel_p_x: " << sel_p_x << " sel_n_x: " << sel_n_x << " sel_p_y: " << sel_p_y << " sel_n_y: " << sel_n_y << std::endl;
                                std::cout << "move_index: " << move_index << std::endl;
#endif
                            util::geo::movement movement = compute_movement_next(sel_p_x, sel_n_x, sel_p_y, sel_n_y);
                            p_i = util::geo::point{movement.x + p_i.x, movement.y + p_i.y};
#if VERBOSE_AGB
                            std::cout << "p_m leaf: (" << p_m.m_x << ", " << p_m.m_y << ")" << std::endl;
#endif
                            if(util::geo::contains(r_q, p_i)) {
#if VERBOSE_AGB
                                std::cout << "the r_q contains p_m leaf" << std::endl;
#endif
                                return true;
                            }
                            move_index++;
                        }
                    }
                }

                ++i;
            }

            return false;



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
            m_rMmq_x_sample.load(in);
            m_rank_rMmq_x_sample.load(in, &m_rMmq_x_sample);
            m_select_rMmq_x_sample.load(in, &m_rMmq_x_sample);
            m_rMmq_y_sample.load(in);
            m_rank_rMmq_y_sample.load(in, &m_rMmq_y_sample);
            m_select_rMmq_y_sample.load(in, &m_rMmq_y_sample);
            m_rmq_x.load(in);
            m_rmq_y.load(in);
            m_rMq_x.load(in);
            m_rMq_y.load(in);

        }

    };
}

#endif //RCT_LOG_REFERENCE_HPP
