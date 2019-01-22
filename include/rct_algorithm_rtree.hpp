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
#include <math_util.hpp>
#include <cmath>
#include <rct_index.hpp>

namespace rct {


    class algorithm {

    private:

        template<class RCTIndex>
        static void time_interval_reference(const typename RCTIndex::size_type oid,
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

            typename RCTIndex::size_type phrase_start, last_movement, start_movement;


            if(delta_phrase_l || ic_phrase_l == ic_phrase_r){
                //Check the first incomplete phrase
                //1. Positions at reference
                phrase_start = rctIndex.log_objects[oid].offsets[ic_phrase_l-1];
                auto point_phrase = rctIndex.log_objects[oid].phrase_point(ic_phrase_l);
                last_movement = std::min(rctIndex.log_objects[oid].last_movement(ic_phrase_l), movement_j);
                auto move_ref_start =  phrase_start + delta_phrase_l+1; //count
                auto move_ref_end = move_ref_start + (last_movement - movement_i); //count
                //TODO: revisar
                if(rctIndex.log_reference.contains_region(point_phrase.x, point_phrase.y, phrase_start, move_ref_start,
                                                          move_ref_end, region_q)){
                    r.push_back(oid);
                    return;
                };
            }
            for(const auto &phrase : phrases_to_check){
                //Check the rest of phrases
                phrase_start = rctIndex.log_objects[oid].offsets[phrase-1];
                auto point_phrase = rctIndex.log_objects[oid].phrase_point(phrase);
                start_movement = rctIndex.log_objects[oid].start_movement(phrase);
                last_movement = rctIndex.log_objects[oid].last_movement(phrase);
                auto move_ref_end = phrase_start+1 + (last_movement - start_movement); //count
                if(rctIndex.log_reference.contains_region(point_phrase.x, point_phrase.y, phrase_start, phrase_start+1,
                                                          move_ref_end, region_q)){
                    r.push_back(oid);
                    return;
                };
            }
            if(delta_phrase_r && ic_phrase_l < ic_phrase_r){
                //Check the last incomplete phrase
                phrase_start = rctIndex.log_objects[oid].offsets[ic_phrase_r-1];
                auto point_phrase = rctIndex.log_objects[oid].phrase_point(ic_phrase_r);
                start_movement = rctIndex.log_objects[oid].start_movement(ic_phrase_r);
                auto move_ref_end = phrase_start+1 + (movement_j - start_movement); //count
                if(rctIndex.log_reference.contains_region(point_phrase.x, point_phrase.y, phrase_start, phrase_start+1,
                                                          move_ref_end, region_q)){
                    r.push_back(oid);
                    return;
                };
            }

        }


    public:

        /***
         *
         * @param oid
         * @param t_q
         * @param rctIndex
         * @param r
         * @return
         */

        template<class RCTIndex>
        static bool search_object(const typename RCTIndex::size_type oid, const typename RCTIndex::value_type t_q,
                           const RCTIndex &rctIndex, util::geo::point &r) {

            auto traj_step = rctIndex.log_objects[oid].start_traj_step();
            if (traj_step.t > t_q || t_q > rctIndex.log_objects[oid].time_end()) return false;
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
        static void search_trajectory(const typename RCTIndex::size_type oid, const typename RCTIndex::value_type t_i,
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
                    phrase_start = 0, next_phrase_beg = 0, x_p_prev, x_n_prev, y_p_prev, y_n_prev, phrase;
            uint32_t t = 0;
            util::geo::movement movement_phrase, movement;
            if(t_temp_i <= t_temp_j){
                rctIndex.log_objects[oid].time_to_movement(t_temp_i, t_temp_j, movement_i, movement_j);
                if(movement_i <= movement_j){
                    t = rctIndex.log_objects[oid].time_next(t_temp_i-1);
                    phrase = rctIndex.log_objects[oid].interval_ref(movement_i, idx_beg, idx_end,
                                                                    next_phrase_beg, movement_phrase);
                    movement = rctIndex.log_reference.compute_movement_init(idx_beg, idx_end, x_p_prev, x_n_prev,
                                                                            y_p_prev, y_n_prev);
                    r.emplace_back(util::geo::traj_step{t, traj_step.x + movement_phrase.x + movement.x,
                                                        traj_step.y + movement_phrase.y + movement.y});
                    ++movement_i; //next movement
                }
                auto prev = idx_end;
                while (movement_i <= movement_j) {
                    t = rctIndex.log_objects[oid].time_next(t);
                    if (next_phrase_beg < movement_i) { //next_phrase_beg is index, and movement_i count
                        rctIndex.log_objects[oid].next_phrase(phrase, prev, next_phrase_beg);
                        movement = rctIndex.log_reference.compute_movement_init_next(prev, x_p_prev,
                                                                                     x_n_prev, y_p_prev, y_n_prev);
                    } else {
                        movement = rctIndex.log_reference.compute_movement_next(x_p_prev,
                                                                                x_n_prev, y_p_prev, y_n_prev);
                    }
                    r.emplace_back(util::geo::traj_step{t, r.back().x + movement.x, r.back().y + movement.y});
                    ++movement_i;
                    ++prev;
                }
            }

        }


        template<class RCTIndex>
        static void search_trajectory_fast(const typename RCTIndex::size_type oid, const typename RCTIndex::value_type t_i,
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
                    phrase_start = 0, next_phrase_beg = 0, x_p_prev, x_n_prev, y_p_prev, y_n_prev, phrase;
            uint32_t t = 0;
            util::geo::movement movement_phrase, movement;
            typename RCTIndex::next_info_type next_info = {false};
            if(t_temp_i <= t_temp_j){
                rctIndex.log_objects[oid].time_to_movement(t_temp_i, t_temp_j, movement_i, movement_j);
                if(movement_i <= movement_j){
                    t = rctIndex.log_objects[oid].time_next_fast(t_temp_i-1, next_info);
                    phrase = rctIndex.log_objects[oid].interval_ref(movement_i, idx_beg, idx_end,
                                                                    next_phrase_beg, movement_phrase);
                    movement = rctIndex.log_reference.compute_movement_init(idx_beg, idx_end, x_p_prev, x_n_prev,
                                                                            y_p_prev, y_n_prev);
                    r.emplace_back(util::geo::traj_step{t, traj_step.x + movement_phrase.x + movement.x,
                                                        traj_step.y + movement_phrase.y + movement.y});
                    ++movement_i; //next movement
                }
                auto prev = idx_end;
                while (movement_i <= movement_j) {
                    t = rctIndex.log_objects[oid].time_next_fast(t, next_info);
                    if (next_phrase_beg < movement_i) { //next_phrase_beg is index, and movement_i count
                        rctIndex.log_objects[oid].next_phrase(phrase, prev, next_phrase_beg);
                        movement = rctIndex.log_reference.compute_movement_init_next(prev, x_p_prev,
                                                                                     x_n_prev, y_p_prev, y_n_prev);
                    } else {
                        movement = rctIndex.log_reference.compute_movement_next(x_p_prev,
                                                                                x_n_prev, y_p_prev, y_n_prev);
                    }
                    r.emplace_back(util::geo::traj_step{t, r.back().x + movement.x, r.back().y + movement.y});
                    ++movement_i;
                    ++prev;
                }
            }

        }

        template<class RCTIndex>
        static void time_slice(const util::geo::region& region_q, const typename RCTIndex::value_type t_q,
                                  const RCTIndex &rctIndex,
                                  std::vector<util::geo::id_point> &r) {

            auto snap_q = t_q / rctIndex.period_snapshot;
            std::vector<uint32_t > data;
            util::geo::point p;
            rctIndex.snapshots[snap_q].intersection(region_q, data);
            for(const auto &id : data){
                if(search_object(id, t_q, rctIndex, p) && util::geo::contains(region_q, p)){
                    r.emplace_back(util::geo::id_point{id, p.x, p.y});
                }
            }
        }


        template<class RCTIndex>
        static void time_interval(const util::geo::region& region_q, const typename RCTIndex::size_type t_i,
                                  const typename RCTIndex::size_type t_j, const RCTIndex &rctIndex,
                                  std::vector<typename  RCTIndex::value_type> &r) {

            auto snap_start = t_i / rctIndex.period_snapshot;
            auto snap_end = t_j /  rctIndex.period_snapshot+1;

            //std::unordered_map<typename RCTIndex::value_type, char> processed_ids;
            sdsl::bit_vector processed_ids(rctIndex.total_objects, 0);
            typename RCTIndex::size_type movement_i = 0, movement_j = 0, c_phrase_i = 0, c_phrase_j = 0,
                    ic_phrase_l = 0, ic_phrase_r = 0, delta_phrase_l = 0, delta_phrase_r = 0;
            typename RCTIndex::size_type t_beg = t_i, t_end = t_j;

            auto time_interval_object = [&] (typename RCTIndex::value_type oid){
                auto beg = std::max(t_beg, rctIndex.log_objects[oid].time_start());
                auto end = std::min(t_end, rctIndex.log_objects[oid].time_end());
                auto traj_step = rctIndex.log_objects[oid].start_traj_step();
                if(end >= traj_step.t) {
                    if(beg <= traj_step.t){
                        if(util::geo::contains(region_q, util::geo::point{traj_step.x, traj_step.y})){
                            r.push_back(oid);
                            processed_ids[oid]=1;
                            return;
                        };
                        beg = traj_step.t+1;
                    }
                    if(beg <= end){
                        rctIndex.log_objects[oid].time_to_movement(beg, end, movement_i, movement_j);
                        if (movement_i > 0 && movement_i <= movement_j) {
                            rctIndex.log_objects[oid].interval_phrases(movement_i, movement_j, c_phrase_i, c_phrase_j,
                                                                       ic_phrase_l, delta_phrase_l,
                                                                       ic_phrase_r, delta_phrase_r);
                            std::vector<typename RCTIndex::size_type> phrases_to_check;
                            if (!rctIndex.log_objects[oid].contains_region(c_phrase_i, c_phrase_j, region_q,
                                                                           phrases_to_check)) {
                                time_interval_reference(oid, region_q, movement_i, movement_j, phrases_to_check,
                                                        ic_phrase_l, delta_phrase_l,
                                                        ic_phrase_r, delta_phrase_r, rctIndex, r);
                            } else {
                                r.push_back(oid);
                            }
                        }
                    }

                }
                processed_ids[oid]=1;
            };

            auto time_interval_left = [&] (typename RCTIndex::value_type snap_id) {
                std::vector<uint32_t > data;
                rctIndex.snapshots[snap_id].intersection(region_q, data);
                //std::cerr << "snapshot: " << data.size() << std::endl;
                for(const auto &id: data){
                    if(!processed_ids[id]){
                        // std::cout << "succ: " << id << std::endl;
                        time_interval_object(id);
                    }
                }
            };

            int64_t snap_id = snap_start;
            while(snap_id < snap_end){
                time_interval_left(snap_id);
                t_beg = ++snap_id * rctIndex.period_snapshot;
            }

        }
    };
}

#endif //RCT_SEARCH_OBJECT_SUPPORT_HPP