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
// Created by Adrián on 26/11/2018.
//

#ifndef RCT_LOG_OBJECT_V1_HPP
#define RCT_LOG_OBJECT_V1_HPP

#include <sdsl/vectors.hpp>
#include <sdsl/sd_vector.hpp>
#include <rmq_succinct_ct.hpp>
#include <rlz_naive.hpp>
#include "alternative_code.hpp"
#include "geo_util.hpp"

namespace rct {

    template <class t_offsets = sdsl::int_vector<>, class t_lengths = sdsl::sd_vector<>,
              class t_values = sdsl::int_vector<>,  class t_disap = sdsl::bit_vector >
    class log_object_v1 {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef t_offsets offsets_type;
        typedef t_lengths lengths_type;
        typedef t_values values_type;
        typedef t_disap disap_type;
        typedef typename t_lengths::rank_1_type rank_1_lengths_type;
        typedef typename t_lengths::select_1_type select_1_lengths_type;
        typedef sdsl::rank_support_v5<1> rank_1_disap_type;
        typedef typename t_disap::select_1_type select_1_disap_type;

    private:
        value_type m_time_start;
        value_type m_x_start;
        value_type m_y_start;

        values_type m_x_values;
        values_type m_y_values;
        rmq_succinct_ct<true> m_rmq_x;
        rmq_succinct_ct<false> m_rMq_x;
        rmq_succinct_ct<true> m_rmq_y;
        rmq_succinct_ct<false> m_rMq_y;

        offsets_type m_offsets;

        lengths_type m_lengths;
        rank_1_lengths_type m_rank_lengths;
        select_1_lengths_type m_select_lengths;

        disap_type m_disap;
        rank_1_disap_type m_rank_disap;
        select_1_disap_type m_select_disap; //TODO: select_next ??? I only need 'select' to compute the next time instant in trajectory

        void copy(const log_object_v1 &o){
            m_time_start = o.m_time_start;
            m_x_start = o.m_x_start;
            m_y_start = o.m_y_start;
            m_x_values = o.m_x_values;
            m_y_values = o.m_y_values;
            m_offsets = o.m_offsets;
            m_lengths = o.m_lengths;
            m_rank_lengths = o.m_rank_lengths;
            m_rank_lengths.set_vector(&m_lengths);
            m_select_lengths = o.m_select_lengths;
            m_select_lengths.set_vector(&m_lengths);
            m_disap = o.m_disap;
            m_rank_disap = o.m_rank_disap;
            m_rank_disap.set_vector(&m_disap);
            m_select_disap = o.m_select_disap;
            m_select_disap.set_vector(&m_disap);
            m_rmq_x = o.m_rmq_x;
            m_rMq_x = o.m_rMq_x;
            m_rmq_y = o.m_rmq_y;
            m_rMq_y = o.m_rMq_y;
        }

    public:

        log_object_v1() = default;


        template<class ContainerTrajectory, class ContainerFactors>
        log_object_v1(const ContainerTrajectory &trajectory, const ContainerFactors &factors) {

            //1. Compute the start values
            m_time_start = trajectory[0].t;
            m_x_start = trajectory[0].x;
            m_y_start = trajectory[0].y;

            //2. Compute the disappearances
            auto last_t = m_time_start;
            size_type factor_i = 0, disap_i = 0;
            m_disap = disap_type(trajectory.back().t - m_time_start + 1, 0);
            for (size_type i = 1; i < trajectory.size(); ++i) {
                const auto &info = trajectory[i];
                for (auto t = last_t + 1; t < info.t; ++t) {
                    m_disap[disap_i++] = 1;
                }
                last_t = info.t;
                disap_i++; //set to zero
            }
            sdsl::util::init_support(m_rank_disap, &m_disap);
            sdsl::util::init_support(m_select_disap, &m_disap);


            //3.1 Prepare the arrays of offsets, x and y
            m_offsets.resize(factors.size());
            m_x_values.resize(factors.size());
            m_y_values.resize(factors.size());

            //3.2 Set offset, length, x and y
            {
                std::vector<size_type> pos_ones_lengths;
                std::vector<value_type> temp_x, temp_y;
                size_type acum_length = 0;
                for (const auto &factor : factors) {
                    acum_length += factor.length;
                    pos_ones_lengths.push_back(acum_length);
                    m_offsets[factor_i] = factor.offset;
                    //Add last value of each phrase
                    temp_x.push_back(trajectory[acum_length].x);
                    temp_y.push_back(trajectory[acum_length].y);
                    ++factor_i;
                }
                sdsl::util::bit_compress(m_offsets);

                //Storing the length in the sd-array
                m_lengths = sdsl::sd_vector<>(pos_ones_lengths.begin(), pos_ones_lengths.end());
                sdsl::util::init_support(m_rank_lengths, &m_lengths);
                sdsl::util::init_support(m_select_lengths, &m_lengths);

                //Construction of structures to solve range minimum queries (rmq) and range maximum queries (rMq)
                m_rmq_x = rmq_succinct_ct<true>(temp_x);
                m_rMq_x = rmq_succinct_ct<false>(temp_x);
                m_rmq_y = rmq_succinct_ct<true>(temp_y);
                m_rMq_y = rmq_succinct_ct<false>(temp_y);

                //The absolute values are stored as differences with respect to the first position
                for (size_type i = 0; i < temp_x.size(); ++i) {
                    m_x_values[i] = alternative_code::encode((int32_t) (temp_x[i] - m_x_start));
                    m_y_values[i] = alternative_code::encode((int32_t) (temp_y[i] - m_y_start));
                }
                sdsl::util::bit_compress(m_x_values);
                sdsl::util::bit_compress(m_y_values);
            }


        }

        //! Copy constructor
        log_object_v1(const log_object_v1& o)
        {
            copy(o);
        }

        //! Move constructor
        log_object_v1(log_object_v1&& o)
        {
            *this = std::move(o);
        }


        log_object_v1 &operator=(const log_object_v1 &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        log_object_v1 &operator=(log_object_v1 &&o) {
            if (this != &o) {
                m_time_start = o.m_time_start;
                m_x_start = o.m_x_start;
                m_y_start = o.m_y_start;
                m_x_values = std::move(o.m_x_values);
                m_y_values = std::move(o.m_y_values);
                m_offsets = std::move(o.m_offsets);
                m_lengths = std::move(o.m_lengths);
                m_rank_lengths = std::move(o.m_rank_lengths);
                m_rank_lengths.set_vector(&m_lengths);
                m_select_lengths = std::move(o.m_select_lengths);
                m_select_lengths.set_vector(&m_lengths);
                m_disap = std::move(o.m_disap);
                m_rank_disap = std::move(o.m_rank_disap);
                m_rank_disap.set_vector(&m_disap);
                m_select_disap = std::move(o.m_select_disap);
                m_select_disap.set_vector(&m_disap);
                m_rmq_x = std::move(o.m_rmq_x);
                m_rMq_x = std::move(o.m_rMq_x);
                m_rmq_y = std::move(o.m_rmq_y);
                m_rMq_y = std::move(o.m_rMq_y);
            }
            return *this;
        }

        void swap(log_object_v1 &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            std::swap(m_time_start, o.m_time_start);
            std::swap(m_x_start, o.m_x_start);
            std::swap(m_y_start, o.m_y_start);
            m_x_values.swap(o.m_x_values);
            m_y_values.swap(o.m_y_values);
            m_offsets.swap(o.m_offsets);
            m_lengths.swap(o.m_lengths);
            sdsl::util::swap_support(m_rank_lengths, o.m_rank_lengths, &m_lengths, &o.m_lengths);
            sdsl::util::swap_support(m_select_lengths, o.m_select_lengths, &m_lengths, &o.m_lengths);
            m_disap.swap(o.m_disap);
            sdsl::util::swap_support(m_rank_disap, o.m_rank_disap, &m_disap, &o.m_disap);
            sdsl::util::swap_support(m_select_disap, o.m_select_disap, &m_disap, &o.m_disap);
            m_rmq_x.swap(o.m_rmq_x);
            m_rmq_y.swap(o.m_rmq_y);
            m_rMq_x.swap(o.m_rMq_x);
            m_rMq_y.swap(o.m_rMq_y);
        }


        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += sdsl::write_member(m_time_start, out, child, "time_start");
            written_bytes += sdsl::write_member(m_x_start, out, child, "x_start");
            written_bytes += sdsl::write_member(m_y_start, out, child, "y_start");
            written_bytes += m_x_values.serialize(out, child, "x_values");
            written_bytes += m_y_values.serialize(out, child, "y_values");
            written_bytes += m_offsets.serialize(out, child, "offsets");
            written_bytes += m_lengths.serialize(out, child, "lengths");
            written_bytes += m_rank_lengths.serialize(out, child, "rank_lengths");
            written_bytes += m_select_lengths.serialize(out, child, "select_lengths");
            written_bytes += m_disap.serialize(out, child, "disap");
            written_bytes += m_rank_disap.serialize(out, child, "rank_disap");
            written_bytes += m_select_disap.serialize(out, child, "select_disap");
            written_bytes += m_rmq_x.serialize(out, child, "rmq_x");
            written_bytes += m_rmq_y.serialize(out, child, "rmq_y");
            written_bytes += m_rMq_x.serialize(out, child, "rMq_x");
            written_bytes += m_rMq_y.serialize(out, child, "rMq_y");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

    };

}

#endif //RCT_LOG_OBJECT_HPP
