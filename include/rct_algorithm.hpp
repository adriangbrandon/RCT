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
// Created by Adrián on 29/11/2018.
//

#ifndef RCT_ALGORITHM_HPP
#define RCT_ALGORITHM_HPP

#include "geo_util.hpp"
#include <vector>
#include <unordered_map>

namespace rct {


    namespace algorithm {

        /***
         *
         * @param oid
         * @param t_q
         * @param rctIndex
         * @param r
         * @return
         */

        template<class RCTIndex>
        bool search_object(const typename RCTIndex::size_type oid, const typename RCTIndex::value_type t_q,
                           const RCTIndex &rctIndex, util::geo::point &r) {

            auto traj_step = rctIndex.log_objects[oid].start_traj_step();
            if (traj_step.t > t_q) return false;
            if (traj_step.t == t_q) {
                r = util::geo::point{traj_step.x, traj_step.y};
                return true;
            }
            typename RCTIndex::size_type movement_q = 0, idx_beg = 0, idx_end = 0;
            util::geo::movement movement_phrase;
            if (rctIndex.log_objects[oid].time_to_movement(t_q, movement_q)) {
                rctIndex.log_objects[oid].interval_ref(movement_q, idx_beg, idx_end, movement_phrase);
                //std::cout << "movement_phrase: " << movement_phrase.x << ", " << movement_phrase.y << std::endl;
                util::geo::movement movement = rctIndex.log_reference.compute_movement(idx_beg, idx_end);
                //std::cout << "movement: " << movement.x << ", " << movement.y << std::endl;
                r = util::geo::point{traj_step.x + movement_phrase.x + movement.x,
                                     traj_step.y + movement_phrase.y + movement.y};
                return true;
            };
            return false;
        }


        template<class RCTIndex>
        void search_trajectory(const typename RCTIndex::size_type oid, const typename RCTIndex::value_type t_i,
                               const typename RCTIndex::size_type t_j, const RCTIndex &rctIndex,
                               std::vector<util::geo::traj_step> &r) {

            auto t_start = rctIndex.log_objects[oid].time_start();
            auto t_end = rctIndex.log_objects[oid].time_end();


            if (t_end < t_i || t_start > t_j) return;
            auto traj_step = rctIndex.log_objects[oid].start_traj_step();
            typename RCTIndex::size_type t_temp_i = t_i, t_temp_j = t_j;
            if (t_i <= t_start) {
                r.push_back(traj_step);
                t_temp_i = t_start + 1;
            }
            if (t_temp_j > t_end) {
                t_temp_j = t_end;
            }

            typename RCTIndex::size_type movement_i = 0, movement_j = 0, idx_beg = 0, idx_end = 0,
                    phrase_start = 0, phrase_end = 0, x_p_prev, x_n_prev, y_p_prev, y_n_prev, phrase;
            uint32_t t = 0;
            util::geo::movement movement_phrase, movement;
            rctIndex.log_objects[oid].time_to_movement(t_temp_i, t_temp_j, movement_i, movement_j);
            if(movement_i <= movement_j){
                t = rctIndex.log_objects[oid].time_next(t_temp_i-1);
                phrase = rctIndex.log_objects[oid].interval_ref(movement_i, idx_beg, idx_end,
                                                                     phrase_end, movement_phrase);
                movement = rctIndex.log_reference.compute_movement_init(idx_beg, idx_end, x_p_prev, x_n_prev,
                                                                        y_p_prev, y_n_prev);
                r.emplace_back(util::geo::traj_step{t, traj_step.x + movement_phrase.x + movement.x,
                                                    traj_step.y + movement_phrase.y + movement.y});
            }
            auto prev = idx_end-1;
            while (movement_i < movement_j) {
                ++movement_i;
                ++prev;
                t = rctIndex.log_objects[oid].time_next(t);
                if (phrase_end < movement_i) {
                    rctIndex.log_objects[oid].next_phrase(phrase, prev, phrase_end);
                    movement = rctIndex.log_reference.compute_movement_init_next(prev, x_p_prev,
                                                                                 x_n_prev, y_p_prev, y_n_prev);
                } else {
                    movement = rctIndex.log_reference.compute_movement_next(prev, x_p_prev,
                                                                            x_n_prev, y_p_prev, y_n_prev);
                }
                r.emplace_back(util::geo::traj_step{t, r.back().x + movement.x, r.back().y + movement.y});
            }

        }

        template<class RCTIndex>
        void time_interval_reference(const typename RCTIndex::size_type oid,
                                     const util::geo::region& region_q,
                                     const typename RCTIndex::size_type movement_i,
                                     const typename RCTIndex::size_type movement_j,
                                     const std::vector<typename RCTIndex::size_type> &phrases_to_check,
                                     const typename RCTIndex::size_type ic_phrase_l,
                                     const typename RCTIndex::size_type delta_phrase_l, //FIXME: ojo a esto
                                     const typename RCTIndex::size_type ic_phrase_r,
                                     const typename RCTIndex::size_type delta_phrase_r, //FIXME: ojo a esto
                                     const RCTIndex &rctIndex,
                                     std::vector<typename  RCTIndex::value_type> &r){

            auto traj_step = rctIndex.log_objects[oid].start_traj_step();
            typename RCTIndex::size_type phrase_start, phrase_x, phrase_y, last_movement, start_movement;
            if(delta_phrase_l){
                //Check the first incomplete phrase
                //1. Positions at reference
                phrase_start = rctIndex.log_objects[oid].offsets[ic_phrase_l-1];
                phrase_x = rctIndex.log_objects[oid].x_values[ic_phrase_l-1] + traj_step.x;
                phrase_y = rctIndex.log_objects[oid].y_values[ic_phrase_l-1] + traj_step.y;
                last_movement = rctIndex.log_objects[oid].last_movement(ic_phrase_l);
                auto move_ref_start =  phrase_start + delta_phrase_l;
                auto move_ref_end = move_ref_start + (last_movement - movement_i);
                //TODO: revisar
                if(rctIndex.log_reference.contains_region(phrase_x, phrase_y, phrase_start, move_ref_start,
                                                          move_ref_end, region_q)){
                    r.push_back(oid);
                    return;
                };
            }
            for(const auto &phrase : phrases_to_check){
                //Check the rest of phrases
                phrase_start = rctIndex.log_objects[oid].offsets[phrase-1];
                phrase_x = rctIndex.log_objects[oid].x_values[phrase-1] + traj_step.x;
                phrase_y = rctIndex.log_objects[oid].y_values[phrase-1] + traj_step.y;
                start_movement = rctIndex.log_objects[oid].start_movement(phrase);
                last_movement = rctIndex.log_objects[oid].last_movement(phrase);
                auto move_ref_end = phrase_start + (last_movement - start_movement);
                if(rctIndex.log_reference.contains_region(phrase_x, phrase_y, phrase_start, phrase_start,
                                                          move_ref_end, region_q)){
                    r.push_back(oid);
                    return;
                };
            }
            if(delta_phrase_r && ic_phrase_l < ic_phrase_r){
                //Check the last incomplete phrase
                phrase_start = rctIndex.log_objects[oid].offsets[ic_phrase_r-1];
                phrase_x = rctIndex.log_objects[oid].x_values[ic_phrase_r-1] + traj_step.x;
                phrase_y = rctIndex.log_objects[oid].y_values[ic_phrase_r-1] + traj_step.y;
                start_movement = rctIndex.log_objects[oid].start_movement(ic_phrase_r);
                auto move_ref_end = phrase_start + (movement_j - start_movement);
                if(rctIndex.log_reference.contains_region(phrase_x, phrase_y, phrase_start, phrase_start,
                                                          move_ref_end, region_q)){
                    r.push_back(oid);
                    return;
                };
            }

        }

        template<class RCTIndex>
        void time_interval(const util::geo::region& region_q, const typename RCTIndex::value_type t_i,
                           const typename RCTIndex::size_type t_j, const RCTIndex &rctIndex,
                           std::vector<typename  RCTIndex::value_type> &r) {


            std::vector<typename RCTIndex::value_type> ids;
            ids.push_back(0);
            std::unordered_map<typename RCTIndex::value_type, char> processed_ids;
            typename RCTIndex::size_type movement_i = 0, movement_j = 0, c_phrase_i = 0, c_phrase_j = 0,
            ic_phrase_l = 0, ic_phrase_r = 0, delta_phrase_l = 0, delta_phrase_r = 0;
            for(const auto &oid: ids){
                if(processed_ids.count(oid) == 0){
                    rctIndex.log_objects[oid].time_to_movement(t_i, t_j, movement_i, movement_j);
                    if(movement_i <= movement_j){
                        rctIndex.log_objects[oid].interval_phrases(movement_i, movement_j, c_phrase_i, c_phrase_j,
                                ic_phrase_l, delta_phrase_l, ic_phrase_r, delta_phrase_r);
                        std::vector<typename RCTIndex::size_type> phrases_to_check;
                        if(!rctIndex.log_objects[oid].contains_region(c_phrase_i, c_phrase_j, region_q, phrases_to_check)){
                            time_interval_reference(oid, region_q, movement_i, movement_j, phrases_to_check, ic_phrase_l,
                                                    delta_phrase_l, ic_phrase_r, delta_phrase_r, rctIndex, r);
                        }else{
                            r.push_back(oid);
                        }
                    }
                    processed_ids[oid]=1;
                }
            }
        }
    };
}

#endif //RCT_SEARCH_OBJECT_SUPPORT_HPP
