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
#include <algorithm>
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
                                            const typename RCTIndex::size_type delta_phrase_l,
                                            const typename RCTIndex::size_type ic_phrase_r,
                                            const typename RCTIndex::size_type delta_phrase_r,
                                            const RCTIndex &rctIndex,
                                            std::vector<typename  RCTIndex::value_type> &r){


            //std::ofstream log("log_ti_l.log", std::ios_base::app | std::ios_base::out);
            //log << "Phrases size: " << phrases_to_check.size() << std::endl;
            typename RCTIndex::size_type phrase_start, last_movement, start_movement;
            if(delta_phrase_l || ic_phrase_l == ic_phrase_r){
                //log << "Incomplete left phrase" << std::endl;
                //Check the first incomplete phrase
                const auto &mbr = rctIndex.log_objects[oid].MBR(ic_phrase_l);
                if(util::geo::contains(region_q, mbr)){
                    r.push_back(oid);
                    return;
                }

                //std::cout << "Oid: " << oid << " touches left: " << util::geo::touches(region_q, mbr) << std::endl;
                //std::cout << "MBR: " << mbr << std::endl;
                if(util::geo::touches(region_q, mbr)){
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
                //log << "Incomplete right phrase" << std::endl;
                const auto &mbr = rctIndex.log_objects[oid].MBR(ic_phrase_r);
                if(util::geo::contains(region_q, mbr)){
                    r.push_back(oid);
                    return;
                }

                if(util::geo::touches(region_q, mbr)) {
                    phrase_start = rctIndex.log_objects[oid].offsets[ic_phrase_r - 1];
                    auto point_phrase = rctIndex.log_objects[oid].phrase_point(ic_phrase_r);
                    start_movement = rctIndex.log_objects[oid].start_movement(ic_phrase_r);
                    auto move_ref_end = phrase_start + 1 + (movement_j - start_movement); //count
                    if (rctIndex.log_reference.contains_region(point_phrase.x, point_phrase.y, phrase_start,
                                                               phrase_start + 1,
                                                               move_ref_end, region_q)) {
                        r.push_back(oid);
                        return;
                    };
                }
            }

        }

        template<class RCTIndex>
        static void time_interval_reference(const typename RCTIndex::size_type oid,
                                            const util::geo::region& region_q,
                                            const typename RCTIndex::size_type movement_i,
                                            const typename RCTIndex::size_type movement_j,
                                            const typename RCTIndex::size_type phrase_i,
                                            const typename RCTIndex::size_type phrase_j,
                                            const typename RCTIndex::size_type ic_phrase_l,
                                            const typename RCTIndex::size_type delta_phrase_l,
                                            const typename RCTIndex::size_type ic_phrase_r,
                                            const typename RCTIndex::size_type delta_phrase_r,
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
            for(auto phrase = phrase_i; phrase <= phrase_j; ++phrase){
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

        template<class RCTIndex>
        static void time_interval_object (typename RCTIndex::value_type oid,
                                          const util::geo::region& region_q,
                                          const typename RCTIndex::size_type t_beg,
                                          const typename RCTIndex::size_type t_end,
                                          sdsl::bit_vector &processed_ids,
                                          const RCTIndex &rctIndex,
                                          std::vector<typename  RCTIndex::value_type> &r){

            typename RCTIndex::size_type movement_i = 0, movement_j = 0, c_phrase_i = 0, c_phrase_j = 0,
                    ic_phrase_l = 0, ic_phrase_r = 0, delta_phrase_l = 0, delta_phrase_r = 0;
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

        template<class RCTIndex>
        static void time_interval_object_brute_force (typename RCTIndex::value_type oid,
                                          const util::geo::region& region_q,
                                          const typename RCTIndex::size_type t_beg,
                                          const typename RCTIndex::size_type t_end,
                                          sdsl::bit_vector &processed_ids,
                                          const RCTIndex &rctIndex,
                                          std::vector<typename  RCTIndex::value_type> &r){

            typename RCTIndex::size_type movement_i = 0, movement_j = 0, c_phrase_i = 0, c_phrase_j = 0,
                    ic_phrase_l = 0, ic_phrase_r = 0, delta_phrase_l = 0, delta_phrase_r = 0;
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
                        time_interval_reference(oid, region_q, movement_i, movement_j, c_phrase_i, c_phrase_j,
                                                ic_phrase_l, delta_phrase_l,
                                                ic_phrase_r, delta_phrase_r, rctIndex, r);
                    }
                }
            }
            processed_ids[oid]=1;
        };

        template<class RCTIndex>
        static void time_interval_left(typename RCTIndex::value_type snap_id, const util::geo::region& region_q,
                                       const typename RCTIndex::size_type t_i,
                                       const typename RCTIndex::size_type t_j,
                                       sdsl::bit_vector &processed_ids,
                                       const RCTIndex &rctIndex,
                                       std::vector<typename  RCTIndex::value_type> &r) {

            auto t_left = std::min(t_j, (snap_id+1)*rctIndex.period_snapshot);
            auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max, t_left - snap_id * rctIndex.period_snapshot,
                                                     rctIndex.x_max, rctIndex.y_max);
            // std::cerr << "Expanded region: " << region_expanded << std::endl;
            auto data = rctIndex.snapshots[snap_id].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                                                                           region_expanded.min.y, region_expanded.max.y);
            //std::cerr << "snapshot: " << data.size() << std::endl;
            for(const auto &e: data){
                //std::cout << "normal: " << e.id << std::endl;
                if(!processed_ids[e.id]){
                    if(t_i <= snap_id*rctIndex.period_snapshot && snap_id*rctIndex.period_snapshot <= t_j
                       && util::geo::contains(region_q, util::geo::point{e.x, e.y})){
                        r.push_back(e.id);
                        processed_ids[e.id]=1;
                        continue;
                    }
                    time_interval_object(e.id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
            //std::cerr << "reap: " << rctIndex.reap[snap_id].size() << std::endl;
            for(const uint32_t &id : rctIndex.reap[snap_id]){
                if(!processed_ids[id]){
                    // std::cout << "succ: " << id << std::endl;
                    time_interval_object(id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
        };

        template<class RCTIndex>
        static void time_interval_left_brute_force(typename RCTIndex::value_type snap_id, const util::geo::region& region_q,
                                       const typename RCTIndex::size_type t_i,
                                       const typename RCTIndex::size_type t_j,
                                       sdsl::bit_vector &processed_ids,
                                       const RCTIndex &rctIndex,
                                       std::vector<typename  RCTIndex::value_type> &r) {

            auto t_left = std::min(t_j, (snap_id+1)*rctIndex.period_snapshot);
            auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max, t_left - snap_id * rctIndex.period_snapshot,
                                                     rctIndex.x_max, rctIndex.y_max);
            // std::cerr << "Expanded region: " << region_expanded << std::endl;
            auto data = rctIndex.snapshots[snap_id].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                                                                           region_expanded.min.y, region_expanded.max.y);
            //std::cerr << "snapshot: " << data.size() << std::endl;
            for(const auto &e: data){
                //std::cout << "normal: " << e.id << std::endl;
                if(!processed_ids[e.id]){
                    if(t_i <= snap_id*rctIndex.period_snapshot && snap_id*rctIndex.period_snapshot <= t_j
                       && util::geo::contains(region_q, util::geo::point{e.x, e.y})){
                        r.push_back(e.id);
                        processed_ids[e.id]=1;
                        continue;
                    }
                    time_interval_object_brute_force(e.id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
            //std::cerr << "reap: " << rctIndex.reap[snap_id].size() << std::endl;
            for(const uint32_t &id : rctIndex.reap[snap_id]){
                if(!processed_ids[id]){
                    // std::cout << "succ: " << id << std::endl;
                    time_interval_object_brute_force(id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
        };

        template<class RCTIndex>
        static void time_interval_right(typename RCTIndex::value_type snap_id, const util::geo::region& region_q,
                                        const typename RCTIndex::size_type t_i,
                                        const typename RCTIndex::size_type t_j,
                                        sdsl::bit_vector &processed_ids,
                                        const RCTIndex &rctIndex,
                                        std::vector<typename  RCTIndex::value_type> &r) {
            auto t_right = std::max(t_i, (snap_id-1)*rctIndex.period_snapshot);
            auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max, snap_id * rctIndex.period_snapshot - t_right,
                                                     rctIndex.x_max, rctIndex.y_max);
            //std::cerr << "Expanded region: " << region_expanded << std::endl;
            auto data = rctIndex.snapshots[snap_id].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                                                                           region_expanded.min.y, region_expanded.max.y);
            //std::cerr << "snapshot: " << data.size() << std::endl;
            for(const auto &e: data){
                //std::cout << "normal: " << e.id << std::endl;
                if(!processed_ids[e.id]){
                    if(t_i <= snap_id*rctIndex.period_snapshot && snap_id*rctIndex.period_snapshot <= t_j
                       && util::geo::contains(region_q, util::geo::point{e.x, e.y})){
                        r.push_back(e.id);
                        processed_ids[e.id]=1;
                        continue;
                    }
                    time_interval_object(e.id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
            //std::cerr << "disap: " << rctIndex.disap[snap_id-1].size() << std::endl;
            for(const uint32_t &id : rctIndex.disap[snap_id-1]){
                if(!processed_ids[id]){
                    // std::cout << "succ: " << id << std::endl;
                    time_interval_object(id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
        };


        template<class RCTIndex>
        static void time_interval_right_brute_force(typename RCTIndex::value_type snap_id, const util::geo::region& region_q,
                                        const typename RCTIndex::size_type t_i,
                                        const typename RCTIndex::size_type t_j,
                                        sdsl::bit_vector &processed_ids,
                                        const RCTIndex &rctIndex,
                                        std::vector<typename  RCTIndex::value_type> &r) {
            auto t_right = std::max(t_i, (snap_id-1)*rctIndex.period_snapshot);
            auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max, snap_id * rctIndex.period_snapshot - t_right,
                                                     rctIndex.x_max, rctIndex.y_max);
            //std::cerr << "Expanded region: " << region_expanded << std::endl;
            auto data = rctIndex.snapshots[snap_id].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                                                                           region_expanded.min.y, region_expanded.max.y);
            //std::cerr << "snapshot: " << data.size() << std::endl;
            for(const auto &e: data){
                //std::cout << "normal: " << e.id << std::endl;
                if(!processed_ids[e.id]){
                    if(t_i <= snap_id*rctIndex.period_snapshot && snap_id*rctIndex.period_snapshot <= t_j
                       && util::geo::contains(region_q, util::geo::point{e.x, e.y})){
                        r.push_back(e.id);
                        processed_ids[e.id]=1;
                        continue;
                    }
                    time_interval_object_brute_force(e.id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
            //std::cerr << "disap: " << rctIndex.disap[snap_id-1].size() << std::endl;
            for(const uint32_t &id : rctIndex.disap[snap_id-1]){
                if(!processed_ids[id]){
                    // std::cout << "succ: " << id << std::endl;
                    time_interval_object_brute_force(id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
        };

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
            util::geo::movement info_phrase;
            if (rctIndex.log_objects[oid].time_to_movement(t_q, movement_q)) {
                rctIndex.log_objects[oid].interval_ref(movement_q, idx_beg, idx_end, info_phrase);
                //std::cout << "info_phrase: " << info_phrase.x << ", " << info_phrase.y << std::endl;
                util::geo::movement movement = rctIndex.log_reference.compute_deltas(idx_end);
                //std::cout << "movement: " << movement.x << ", " << movement.y << std::endl;
                //std::cout << "traj_step.x: " << traj_step.x << " info_phrase.x: " << info_phrase.x << " movement.x: " << movement.x << std::endl;
                r = util::geo::point{traj_step.x + info_phrase.x + movement.x,
                                     traj_step.y + info_phrase.y + movement.y};
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
            util::geo::movement info_phrase, movement;
            if(t_temp_i <= t_temp_j){
                rctIndex.log_objects[oid].time_to_movement(t_temp_i, t_temp_j, movement_i, movement_j);
                if(movement_i <= movement_j){
                    t = rctIndex.log_objects[oid].time_next(t_temp_i-1);
                    phrase = rctIndex.log_objects[oid].interval_ref(movement_i, idx_beg, idx_end,
                                                                    next_phrase_beg, info_phrase);
                    movement = rctIndex.log_reference.compute_deltas_init(idx_beg, idx_end, x_p_prev, x_n_prev,
                                                                            y_p_prev, y_n_prev);
                    r.emplace_back(util::geo::traj_step{t, traj_step.x + info_phrase.x + movement.x,
                                                        traj_step.y + info_phrase.y + movement.y});
                    ++movement_i; //next movement
                }
                auto prev = idx_end;
                while (movement_i <= movement_j) {
                    t = rctIndex.log_objects[oid].time_next(t);
                    if (next_phrase_beg < movement_i) { //next_phrase_beg is index, and movement_i count
                        rctIndex.log_objects[oid].next_phrase(phrase, prev, next_phrase_beg);
                        movement = rctIndex.log_reference.compute_deltas_init_next(prev, x_p_prev,
                                                                                     x_n_prev, y_p_prev, y_n_prev);
                    } else {
                        movement = rctIndex.log_reference.compute_deltas_next(x_p_prev,
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
                    movement = rctIndex.log_reference.compute_deltas_init(idx_end, x_p_prev, x_n_prev,
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
                       // std::cout << phrase << std::endl;
                        auto point = rctIndex.log_objects[oid].phrase_point(phrase);
                        movement = rctIndex.log_reference.compute_deltas_init(prev+1, x_p_prev, x_n_prev, y_p_prev, y_n_prev);
                        //std::cout << "Adding t=" << t << " x-> " << point.x << " " << movement.x << std::endl;
                        r.emplace_back(util::geo::traj_step{t, (uint32_t) (point.x + movement.x),
                                                            (uint32_t) (point.y + movement.y)});
                    } else {
                        movement = rctIndex.log_reference.compute_deltas_next(x_p_prev,
                                                                                x_n_prev, y_p_prev, y_n_prev);
                        r.emplace_back(util::geo::traj_step{t, r.back().x + movement.x, r.back().y + movement.y});
                    }
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
            if(t_q == snap_q * rctIndex.period_snapshot){
                r = rctIndex.snapshots[snap_q].find_objects_in_region(region_q.min.x, region_q.max.x,
                                                                      region_q.min.y, region_q.max.y);
                return;
            }

            std::vector<util::geo::id_point> data;
            uint64_t total = 0;
            if(t_q - snap_q * rctIndex.period_snapshot > rctIndex.period_snapshot / 2
               && snap_q < rctIndex.last_snapshot() ){
                //std::cout << "right" << std::endl;
                auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max,
                        (snap_q+1) * rctIndex.period_snapshot - t_q,
                        rctIndex.x_max, rctIndex.y_max);
                //std::cout << "Expanded region: " << region_expanded << std::endl;
                data = rctIndex.snapshots[snap_q+1].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                        region_expanded.min.y, region_expanded.max.y);
                util::geo::point p;
                //std::cerr << " [candidates_snap: " << data.size() << " ] " << std::endl;
                //sdsl::bit_vector processed_ids(rctIndex.total_objects, 0);;
                for(const auto &d : data){
                    if(search_object(d.id, t_q, rctIndex, p) && util::geo::contains(region_q, p)){
                        r.emplace_back(util::geo::id_point{d.id, p.x, p.y});
                    }
                    //processed_ids[d.id]=1;
                }

                //std::cerr << "snap done" << std::endl;
                /*typename RCTIndex::value_type id = rctIndex.succs_disap[snap_q](0);
                uint64_t succs = 0;
                while(id < rctIndex.total_objects){
                    if(search_object(id, t_q, rctIndex, p) && util::geo::contains(region_q, p)){
                        r.emplace_back(util::geo::id_point{id, p.x, p.y});
                    }
                    ++total;
                    id = rctIndex.succs_disap[snap_q](++id);
                    ++succs;
                }*/

                for(const uint32_t &id : rctIndex.disap[snap_q]){
                    if(search_object(id, t_q, rctIndex, p) && util::geo::contains(region_q, p)){
                        r.emplace_back(util::geo::id_point{id, p.x, p.y});
                    }
                }
                //std::cerr << "[candidates_disap: " << succs << "]" << std::endl;
                //std::cerr << "time_slice_right: " << total << std::endl;
            }else{
                //std::cout << "left" << std::endl;
                auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max, t_q - snap_q * rctIndex.period_snapshot,
                                                         rctIndex.x_max, rctIndex.y_max);
                //std::cout << "Expanded region: " << region_expanded << std::endl;
                data = rctIndex.snapshots[snap_q].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                                                                              region_expanded.min.y, region_expanded.max.y);

                util::geo::point p;
                //std::cerr << " [candidates_snap: " << data.size() << " ] " << std::endl;
                //std::unordered_map<typename RCTIndex::value_type, char> processed_ids;
                for(const auto &d : data){
                    ++total;
                    if(search_object(d.id, t_q, rctIndex, p) && util::geo::contains(region_q, p)){
                        r.emplace_back(util::geo::id_point{d.id, p.x, p.y});
                    }
                    //processed_ids[d.id]=1;
                }

                //std::cerr << "snap done" << std::endl;
               /* typename RCTIndex::value_type id = rctIndex.succs_reap[snap_q](0);
                uint64_t succs = 0;
                while(id < rctIndex.total_objects){
                    ++total;
                    if(search_object(id, t_q, rctIndex, p) && util::geo::contains(region_q, p)){
                        r.emplace_back(util::geo::id_point{id, p.x, p.y});
                    }
                    id = rctIndex.succs_reap[snap_q](++id);
                    ++succs;
                }*/
                for(const uint32_t &id : rctIndex.reap[snap_q]){
                    if(search_object(id, t_q, rctIndex, p) && util::geo::contains(region_q, p)){
                        r.emplace_back(util::geo::id_point{id, p.x, p.y});
                    }
                }
                //std::cerr << "[candidates_reap: " << succs << "]" << std::endl;
                //std::cerr << "time_slice_right: " << total << std::ednl;
            }

            //std::cerr << " [candidates_reap: " << succs << " ] " << std::endl;
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

            /*auto time_interval_object = [&] (typename RCTIndex::value_type oid){
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
                auto t_left = std::min(t_j, (snap_id+1)*rctIndex.period_snapshot);
                auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max, t_left - snap_id * rctIndex.period_snapshot,
                                                    rctIndex.x_max, rctIndex.y_max);
               // std::cerr << "Expanded region: " << region_expanded << std::endl;
                auto data = rctIndex.snapshots[snap_id].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                                                                               region_expanded.min.y, region_expanded.max.y);
                //std::cerr << "snapshot: " << data.size() << std::endl;
                for(const auto &e: data){
                    //std::cout << "normal: " << e.id << std::endl;
                    if(!processed_ids[e.id]){
                        if(t_i <= snap_id*rctIndex.period_snapshot && snap_id*rctIndex.period_snapshot <= t_j
                           && util::geo::contains(region_q, util::geo::point{e.x, e.y})){
                            r.push_back(e.id);
                            processed_ids[e.id]=1;
                            continue;
                        }
                        time_interval_object(e.id);
                    }
                }
                //std::cerr << "reap: " << rctIndex.reap[snap_id].size() << std::endl;
                for(const uint32_t &id : rctIndex.reap[snap_id]){
                    if(!processed_ids[id]){
                        // std::cout << "succ: " << id << std::endl;
                        time_interval_object(id);
                    }
                }
            };

            auto time_interval_right = [&] (typename RCTIndex::value_type snap_id) {
                auto t_right = std::max(t_i, (snap_id-1)*rctIndex.period_snapshot);
                auto region_expanded = util::geo::expand(region_q, rctIndex.speed_max, snap_id * rctIndex.period_snapshot - t_right,
                                                    rctIndex.x_max, rctIndex.y_max);
                //std::cerr << "Expanded region: " << region_expanded << std::endl;
                auto data = rctIndex.snapshots[snap_id].find_objects_in_region(region_expanded.min.x, region_expanded.max.x,
                                                                               region_expanded.min.y, region_expanded.max.y);
                //std::cerr << "snapshot: " << data.size() << std::endl;
                for(const auto &e: data){
                    //std::cout << "normal: " << e.id << std::endl;
                    if(!processed_ids[e.id]){
                        if(t_i <= snap_id*rctIndex.period_snapshot && snap_id*rctIndex.period_snapshot <= t_j
                           && util::geo::contains(region_q, util::geo::point{e.x, e.y})){
                            r.push_back(e.id);
                            processed_ids[e.id]=1;
                            continue;
                        }
                        time_interval_object(e.id);
                    }
                }
                //std::cerr << "disap: " << rctIndex.disap[snap_id-1].size() << std::endl;
                for(const uint32_t &id : rctIndex.disap[snap_id-1]){
                    if(!processed_ids[id]){
                        // std::cout << "succ: " << id << std::endl;
                        time_interval_object(id);
                    }
                }
            };*/

            if(t_i - snap_start * rctIndex.period_snapshot > snap_end * rctIndex.period_snapshot - t_j
               && snap_end <= rctIndex.last_snapshot()){
                int64_t snap_id = snap_end;
                while(snap_id > snap_start){
                    time_interval_right(snap_id, region_q, t_beg, t_end, processed_ids, rctIndex, r);
                    t_end = --snap_id * rctIndex.period_snapshot;
                }
            }else{
                int64_t snap_id = snap_start;
                while(snap_id < snap_end){
                    time_interval_left(snap_id, region_q, t_beg, t_end, processed_ids, rctIndex, r);
                    t_beg = ++snap_id * rctIndex.period_snapshot;
                }
            }
            //std::cerr << std::endl;


        }

        template<class RCTIndex>
        static void time_interval_brute_force(const util::geo::region& region_q, const typename RCTIndex::size_type t_i,
                                  const typename RCTIndex::size_type t_j, const RCTIndex &rctIndex,
                                  std::vector<typename  RCTIndex::value_type> &r) {

            auto snap_start = t_i / rctIndex.period_snapshot;
            auto snap_end = t_j /  rctIndex.period_snapshot+1;

            //std::unordered_map<typename RCTIndex::value_type, char> processed_ids;
            sdsl::bit_vector processed_ids(rctIndex.total_objects, 0);
            typename RCTIndex::size_type movement_i = 0, movement_j = 0, c_phrase_i = 0, c_phrase_j = 0,
                    ic_phrase_l = 0, ic_phrase_r = 0, delta_phrase_l = 0, delta_phrase_r = 0;
            typename RCTIndex::size_type t_beg = t_i, t_end = t_j;


            if(t_i - snap_start * rctIndex.period_snapshot > snap_end * rctIndex.period_snapshot - t_j
               && snap_end <= rctIndex.last_snapshot()){
                int64_t snap_id = snap_end;
                while(snap_id > snap_start){
                    time_interval_right_brute_force(snap_id, region_q, t_beg, t_end, processed_ids, rctIndex, r);
                    t_end = --snap_id * rctIndex.period_snapshot;
                }
            }else{
                int64_t snap_id = snap_start;
                while(snap_id < snap_end){
                    time_interval_left_brute_force(snap_id, region_q, t_beg, t_end, processed_ids, rctIndex, r);
                    t_beg = ++snap_id * rctIndex.period_snapshot;
                }
            }
            //std::cerr << std::endl;


        }
    };

}

#endif //RCT_SEARCH_OBJECT_SUPPORT_HPP
