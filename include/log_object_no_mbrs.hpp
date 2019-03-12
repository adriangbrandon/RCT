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

#ifndef RCT_LOG_OBJECT_NO_MBRS_HPP
#define RCT_LOG_OBJECT_NO_MBRS_HPP

#include <sdsl/vectors.hpp>
#include <sdsl/sd_vector.hpp>
#include <rmq_succinct_ct.hpp>
#include <rlz_naive.hpp>
#include <succ_support_v.hpp>
#include <queue>
#include "alternative_code.hpp"
#include <oz_vector.hpp>
#include <runs_vector.hpp>
#include "geo_util.hpp"
#include <runs_bitvector.hpp>

#define VERBOSE 0

namespace rct {

    template <class t_offsets = sdsl::dac_vector_dp<sdsl::bit_vector>,
              class t_lengths = sdsl::sd_vector<>,
              class t_values = sdsl::dac_vector_dp<sdsl::bit_vector>,
              class t_values_min_max= sdsl::dac_vector_dp<sdsl::bit_vector>,
              class t_disap = runs_bitvector<>,
              class t_disap_rank_1 = rank_support_runs_bitvector<1>,
              class t_disap_succ_0 = succ_support_runs_bitvector<0>>
    class log_object_no_mbrs {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef t_offsets offsets_type;
        typedef t_lengths lengths_type;
        typedef t_values values_type;
        typedef t_values_min_max values_min_max_type;
        typedef t_disap disap_type;
        typedef typename t_lengths::rank_1_type rank_1_lengths_type;
        typedef typename t_lengths::select_1_type select_1_lengths_type;
        typedef t_disap_rank_1 rank_1_disap_type;
        typedef t_disap_succ_0 succ_0_disap_type;
        typedef typename t_disap::next_info_type next_info_type;

    private:
        value_type m_time_start = 0;
        value_type m_x_start = 0;
        value_type m_y_start = 0;

        values_type m_x_values;
        values_type m_y_values;

        offsets_type m_offsets;

        lengths_type m_lengths;
        rank_1_lengths_type m_rank_lengths;
        select_1_lengths_type m_select_lengths;

        disap_type m_disap;
        rank_1_disap_type m_rank_disap;
        succ_0_disap_type m_succ_0_disap;

        void copy(const log_object_no_mbrs &o){
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
            m_succ_0_disap = o.m_succ_0_disap;
            m_succ_0_disap.set_vector(&m_disap);
        }

        template <class Container, class BitVector>
        void compress_data(sdsl::dac_vector_dp<BitVector> &ref, const Container &data){
            ref = sdsl::dac_vector_dp<BitVector>(data);
        }

        template <class Container>
        void compress_data(sdsl::dac_vector<8> &ref, const Container &data){
            ref = sdsl::dac_vector<8>(data);
        }

        template <class Container>
        void compress_data(sdsl::int_vector<> &ref, const Container &data){
            ref = sdsl::int_vector<>(data.size());
            size_type idx = 0;
            for(const auto &value : data){
                ref[idx++] = value;
            }
            sdsl::util::bit_compress(ref);
        }


    public:

        log_object_no_mbrs() = default;

        const values_type &x_values = m_x_values;
        const values_type &y_values = m_y_values;
        const offsets_type  &offsets = m_offsets;
        const disap_type &disap = m_disap;

        template<class ContainerTrajectory, class ContainerFactors>
        log_object_no_mbrs(const ContainerTrajectory &trajectory, const ContainerFactors &factors) {

            //1. Compute the start values
            m_time_start = trajectory[0].t;
            m_x_start = trajectory[0].x;
            m_y_start = trajectory[0].y;

            //2. Compute the disappearances
            auto last_t = m_time_start;
            size_type disap_i = 1;
            sdsl::bit_vector aux_disap(trajectory.back().t - m_time_start + 1, 0);
            uint64_t d = 0;
            for (size_type i = 1; i < trajectory.size(); ++i) {
                const auto &info = trajectory[i];
                for (auto t = last_t + 1; t < info.t; ++t) {
                    aux_disap[disap_i++] = 1;
                    d++;
                }
                last_t = info.t;
                ++disap_i; //set to zero
            }
            m_disap = disap_type(aux_disap);
            sdsl::util::init_support(m_rank_disap, &m_disap);
            sdsl::util::init_support(m_succ_0_disap, &m_disap);

            //3.2 Set offset, length, x and y
            {
                std::vector<size_type> offset_temp(factors.size());
                std::vector<size_type> x_values_temp(factors.size());
                std::vector<size_type> y_values_temp(factors.size());
                std::vector<size_type> pos_ones_lengths;
                std::vector<value_type> temp_x, temp_y;
                size_type acum_length = 0, start_phrase = 1;
                size_type factor_i = 0;
                for (const auto &factor : factors) {
                    pos_ones_lengths.push_back(acum_length);
                    offset_temp[factor_i] = factor.offset;
                    //Add last value of each phrase (except the last one)
                    temp_x.push_back(trajectory[acum_length].x);
                    temp_y.push_back(trajectory[acum_length].y);
                    start_phrase += factor.length;
                    acum_length += factor.length;
                    ++factor_i;
                }
                if(factor_i > 0) pos_ones_lengths.push_back(acum_length);
                compress_data(m_offsets, offset_temp);

                //Storing the length in the sd-array
                m_lengths = sdsl::sd_vector<>(pos_ones_lengths.begin(), pos_ones_lengths.end());
                sdsl::util::init_support(m_rank_lengths, &m_lengths);
                sdsl::util::init_support(m_select_lengths, &m_lengths);

                //The absolute values are stored as differences with respect to the first position
                for (size_type i = 0; i < temp_x.size(); ++i) {
                    x_values_temp[i] = alternative_code::encode((int32_t) (temp_x[i] - m_x_start));
                    y_values_temp[i] = alternative_code::encode((int32_t) (temp_y[i] - m_y_start));
                }
                compress_data(m_x_values, x_values_temp);
                compress_data(m_y_values, y_values_temp);
            }


        }


        inline util::geo::traj_step start_traj_step() const {
            return util::geo::traj_step{m_time_start, m_x_start, m_y_start};
        };

        inline util::geo::point phrase_point(const size_type phrase) const {
            return util::geo::point{(uint32_t) (alternative_code::decode(x_values[phrase-1]) + m_x_start),
                                    (uint32_t) (alternative_code::decode(y_values[phrase-1]) + m_y_start)};
        };

        inline size_type time_start() const {
            return m_time_start;
        }

        inline size_type time_end() const {
            return m_time_start + m_disap.size()-1;
        }

        inline size_type time_next(const size_type t) const{
            return m_time_start + m_succ_0_disap(t - m_time_start+1);
        }

        inline size_type time_next_fast(const size_type t, next_info_type &next_info) const{
            return m_time_start + m_succ_0_disap.succ(t - m_time_start+1, next_info);
        }


        inline bool time_to_movement(const size_type t_q, size_type &movement_q) const {
            auto idx = t_q - m_time_start;
            if(m_disap[idx]) return false; //Disappeared
            movement_q = idx - m_rank_disap(idx);//index in length
            //idx-1-rank(idx-1+1)+1
            return true;
        }

        inline void time_to_movement(const size_type t_i, const size_type t_j, size_type &movement_i, size_type &movement_j) const{
            auto i = t_i - m_time_start;
            auto j = t_j - m_time_start;
           // movement_i = i - m_rank_disap(i+1); //count
            movement_i = (i-1) - m_rank_disap(i)+1; //count
            movement_j = j - m_rank_disap(j+1); //count
        }

        inline size_type start_movement(const size_type phrase) const {
            return m_select_lengths(phrase)+1;
        }

        inline size_type last_movement(const size_type phrase) const {
            return m_select_lengths(phrase+1);
        }


        inline void interval_phrases(const size_type movement_i, const size_type movement_j,
                                     size_type &c_phrase_i, size_type &c_phrase_j,
                                     size_type &ic_phrase_l, size_type &delta_phrase_l,
                                     size_type &ic_phrase_r, size_type &delta_phrase_r) const {

           // print_length_vector();
            ic_phrase_l = m_rank_lengths(movement_i);
            delta_phrase_l = movement_i - start_movement(ic_phrase_l);
            if(!delta_phrase_l){
                c_phrase_i = ic_phrase_l;
            }else{
                c_phrase_i = ic_phrase_l+1;
            }
            ic_phrase_r = m_rank_lengths(movement_j);
            delta_phrase_r = last_movement(ic_phrase_r) - movement_j;
            if(!delta_phrase_r){
                c_phrase_j = ic_phrase_r;
            }else{
                c_phrase_j = ic_phrase_r-1;
            }
        }


        //Pre: t_q > m_time_start
        inline size_type interval_ref(const size_type movement_q, size_type &idx_beg, size_type &idx_end,
                                      util::geo::movement &r) const {

            auto phrase = m_rank_lengths(movement_q);
            idx_beg = m_offsets[phrase-1];
            idx_end = (movement_q - m_select_lengths(phrase)) + idx_beg;
            r = util::geo::movement{(int32_t) alternative_code::decode(m_x_values[phrase-1]),
                                    (int32_t) alternative_code::decode(m_y_values[phrase-1])};
            return phrase;
        }

        //Pre: t_i > m_time_start && t_j <= m_time_end
        inline size_type interval_ref(const size_type movement_i, size_type &idx_beg, size_type &idx_end,
                size_type &next_phrase_beg, util::geo::movement &r) const {

            auto phrase = m_rank_lengths(movement_i);
            next_phrase_beg = m_select_lengths(phrase+1);
            idx_beg = m_offsets[phrase-1];
            idx_end = (movement_i - m_select_lengths(phrase)) + idx_beg;
            r = util::geo::movement{(int32_t) alternative_code::decode(m_x_values[phrase-1]),
                                    (int32_t) alternative_code::decode(m_y_values[phrase-1])};
            return phrase;
        }

        inline void next_phrase(size_type &phrase, size_type &idx_beg, size_type &next_phrase_beg) const {
            ++phrase;
            next_phrase_beg = m_select_lengths(phrase+1);
            idx_beg = m_offsets[phrase-1];
        }

        //! Copy constructor
        log_object_no_mbrs(const log_object_no_mbrs& o)
        {
            copy(o);
        }

        //! Move constructor
        log_object_no_mbrs(log_object_no_mbrs&& o)
        {
            *this = std::move(o);
        }


        log_object_no_mbrs &operator=(const log_object_no_mbrs &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        log_object_no_mbrs &operator=(log_object_no_mbrs &&o) {
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
                m_succ_0_disap = std::move(o.m_succ_0_disap);
                m_succ_0_disap.set_vector(&m_disap);
            }
            return *this;
        }

        void swap(log_object_no_mbrs &o) {
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
            sdsl::util::swap_support(m_succ_0_disap, o.m_succ_0_disap, &m_disap, &o.m_disap);
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
            written_bytes += m_succ_0_disap.serialize(out, child, "succ_0_disap");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            sdsl::read_member(m_time_start, in);
            sdsl::read_member(m_x_start, in);
            sdsl::read_member(m_y_start, in);
            m_x_values.load(in);
            m_y_values.load(in);
            m_offsets.load(in);
            m_lengths.load(in);
            m_rank_lengths.load(in, &m_lengths);
            m_select_lengths.load(in, &m_lengths);
            m_disap.load(in);
            m_rank_disap.load(in, &m_disap);
            m_succ_0_disap.load(in, &m_disap);

        }

        void print(){
            for(size_type i = 0; i < m_x_values.size(); ++i){
                std::cout << "(x,y)[" <<i << "]=" << alternative_code::decode(m_x_values[i]) << "," << alternative_code::decode(m_y_values[i]) << " ";
            }
            std::cout << std::endl;
            for(size_type i = 0; i < m_lengths.size(); ++i){
                std::cout << "length[" << i <<"]=" << m_lengths[i] << " ";
            }
            std::cout << std::endl;
        }



        template<class Container>
        void stats(size_type &byte1, size_type &byte2, size_type &byte3, size_type &byte4, const Container &cont) const {
            for(const auto &v: cont){
                auto size = sdsl::bits::hi(v);
                if(size < 8){
                    ++byte1;
                }else if (size < 16){
                    ++byte2;
                }else if (size < 24){
                    ++byte3;
                }else{
                    ++byte4;
                }
            }
        }

        template<class Container>
        double diff_dac(const Container &cont) const{
            auto actual_size = sdsl::size_in_mega_bytes(cont);
            sdsl::dac_vector_dp<sdsl::bit_vector> m_dac(cont);
            auto new_size = sdsl::size_in_mega_bytes(m_dac);
            return  new_size - actual_size;
        }

        template<class Container>
        void runs(const Container &cont) const{
            size_type i = 0, j = 1;
            bool ones = false;
            size_type k = 0;
            for(const auto &v : cont){
                if(k == 0) {
                    ones = (v == 1);
                }else{
                    if((v && ones) || (!v && !ones)){
                        ++j;
                    }else{
                        std::cout << "ones: " << ones << " [" << i << ", " << j << "]";
                        i = ++j;
                        ones = !ones;
                    }
                }
                ++k;
            }
            std::cout << std::endl;
        }

        template<class Container>
        double diff_oz_vector(const Container &cont) const{
            auto actual_size = sdsl::size_in_mega_bytes(cont);
            oz_vector<> m_oz(cont);
            auto new_size = sdsl::size_in_mega_bytes(m_oz);
            return  new_size - actual_size;
        }

        template<class Container>
        double diff_runs_vector(const Container &cont) const{
            auto actual_size = sdsl::size_in_mega_bytes(cont);
            runs_vector<> m_oz(cont);
            auto new_size = sdsl::size_in_mega_bytes(m_oz);
            return  new_size - actual_size;
        }

        template<class Container>
        double diff_sd_vector(const Container &cont) const{
            auto actual_size = sdsl::size_in_mega_bytes(cont);
            sdsl::sd_vector<> m_sd(cont);
            auto new_size = sdsl::size_in_mega_bytes(m_sd);
            return  new_size - actual_size;
        }

        template<class Container>
        double diff_runs_bitvector(const Container &cont) const{
            auto actual_size = sdsl::size_in_mega_bytes(cont);
            runs_bitvector<> m_sd(cont);
            auto new_size = sdsl::size_in_mega_bytes(m_sd);
            return  new_size - actual_size;
        }

        template<class Container>
        double size_runs_bitvector(const Container &cont) const{
            auto actual_size = sdsl::size_in_mega_bytes(cont);
            runs_bitvector<> m_sd(cont);
            succ_support_runs_bitvector<0> m_succ;
            sdsl::util::init_support(m_succ, &m_sd);
            auto new_size = sdsl::size_in_mega_bytes(m_sd) + sdsl::size_in_mega_bytes(m_succ);
            return  new_size;
        }

        template<class Container>
        double size_element(const Container &cont) const {
            return sdsl::size_in_mega_bytes(cont);
        }

        template<class Container>
        double diff_offset_vector(const Container &cont) const{
            auto actual_size = sdsl::size_in_mega_bytes(cont);
            sdsl::enc_vector<> enc_v(cont);
            auto new_size = sdsl::size_in_mega_bytes(enc_v);
           // std::cout << "actual_size: " << actual_size << std::endl;
            //std::cout << "new_size: " << new_size << std::endl;
            return  new_size - actual_size;
        }

        void print_offset_vector() const {
            for(const auto &v : m_offsets){
                std::cout << v << ", ";
            }
            std::cout << std::endl;
        }

        void print_length_vector() const {
            size_type  i = 0, one = 1;
            for(const auto &v : m_lengths){
                if(v == 1) {
                    std::cout << "one [" << one << "] at position: " << i << std::endl;
                    ++one;
                }
                ++i;
            }
        }

        void print_lengths() const {
            size_type  i = 0, last_i = 0, one = 1;
            for(const auto &v : m_lengths){
                if(v == 1) {
                    if(i-last_i > 0){
                        std::cerr << i - last_i << std::endl;
                    }
                    last_i = i;
                }
                ++i;
            }
        }


    };

    using log_object_no_mbrs_dac_vector = log_object_no_mbrs< sdsl::dac_vector_dp<sdsl::bit_vector>, sdsl::sd_vector<>,
                                  sdsl::dac_vector_dp<sdsl::bit_vector>, sdsl::dac_vector_dp<sdsl::bit_vector>>;

    using log_object_no_mbrs_int_vector = log_object_no_mbrs< sdsl::int_vector<>, sdsl::sd_vector<>, sdsl::int_vector<>,
                                  sdsl::dac_vector_dp<sdsl::bit_vector> >;

}

#endif //RCT_LOG_OBJECT_NO_MBRS_HPP
