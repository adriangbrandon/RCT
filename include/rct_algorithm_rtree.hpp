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
#include <knn_support_helper_rtree.hpp>
#include <heap_max_min.hpp>
#include <rmMq.hpp>

using namespace rct::knn_support_helper_rtree;

namespace rct {


    class algorithm {

    private:

        static void insert_result(knn_element &candidate, pq_knn_result_type &pq_results, pq_knn_movement_type &pq_distances, uint64_t k){
            /*if(pq_distances.empty()){
                std::cout << "Insert: " << candidate.id << " weight=" << candidate.weight << " distance=" << candidate.distance << " top= <none>"  << std::endl;
            }else{
                std::cout << "Insert: " << candidate.id << " weight=" << candidate.weight << " distance=" << candidate.distance << " top=" <<  pq_distances.top() << std::endl;
            }*/
            pq_results.push(candidate);
            if (pq_distances.size() < k) {
                pq_distances.push(candidate.distance);
            } else if (candidate.distance < pq_distances.top()) {
                pq_distances.pop();
                pq_distances.push(candidate.distance);
            }
        }

        static void insert_result2(knn_element &candidate, std::vector<knn_element> &results, pq_knn_movement_type &pq_distances, uint64_t k){
            /*if(pq_distances.empty()){
                std::cout << "Insert: " << candidate.id << " weight=" << candidate.weight << " distance=" << candidate.distance << " top= <none>"  << std::endl;
            }else{
                std::cout << "Insert: " << candidate.id << " weight=" << candidate.weight << " distance=" << candidate.distance << " top=" <<  pq_distances.top() << std::endl;
            }*/
            results.push_back(candidate);
            if (pq_distances.size() < k) {
                pq_distances.push(candidate.distance);
            } else if (candidate.distance < pq_distances.top()) {
                pq_distances.pop();
                pq_distances.push(candidate.distance);
            }
        }

        template<class RCTIndex, class t_array>
        static void init_knn(const typename RCTIndex::size_type t_q,
                             const util::geo::point &p_q,
                             const typename RCTIndex::size_type k,
                             pq_knn_result_type &knn_pq,
                             pq_knn_result_type &pq_results,
                             pq_knn_movement_type &pq_distances,
                             const t_array &disap,
                             const uint64_t snap_id, const RCTIndex &rctIndex){
            std::vector<knn_element> elements;
            for(const auto &oid : disap){
                util::geo::point p;
                if(search_object(oid, t_q, rctIndex, p)){
                    //double_t distance = util::geo::distance(p, p_q);
                    knn_element knn_e(oid, p_q, p, t_q);
                    insert_result2(knn_e, elements, pq_distances, k);
                };
            }
            pq_results = pq_knn_result_type(elements.begin(), elements.end());

            knn_element knn_e(0, util::geo::point{0,0}, util::geo::point{0,0}, 0, snap_id * rctIndex.period_snapshot, 0);
            knn_pq.push(knn_e);
        };

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

        template<class RCTIndex>
        static void time_interval_reference(const typename RCTIndex::size_type oid,
                                            const util::geo::region& region_q,
                                            const typename RCTIndex::size_type movement_i,
                                            const typename RCTIndex::size_type movement_j,
                                            const typename RCTIndex::size_type phrase_i,
                                            const typename RCTIndex::size_type phrase_j,
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

            std::vector<uint32_t > data;
            rctIndex.snapshots[snap_id].intersection(region_q, data);
            //std::cerr << "snapshot: " << data.size() << std::endl;
            for(const auto &id: data){
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

            std::vector<uint32_t > data;
            rctIndex.snapshots[snap_id].intersection(region_q, data);
            //std::cerr << "snapshot: " << data.size() << std::endl;
            for(const auto &id: data){
                if(!processed_ids[id]){
                    // std::cout << "succ: " << id << std::endl;
                    time_interval_object_brute_force(id, region_q, t_i, t_j, processed_ids, rctIndex, r);
                }
            }
        };

        static void update_MBR(util::geo::region &mbr_result, util::geo::region &mbr, bool &init_mbr){
            if(init_mbr){
                if(mbr_result.min.x > mbr.min.x){
                    mbr_result.min.x = mbr.min.x;
                }
                if(mbr_result.min.y > mbr.min.y){
                    mbr_result.min.y = mbr.min.y;
                }
                if(mbr_result.max.x < mbr.max.x){
                    mbr_result.max.x = mbr.max.x;
                }
                if(mbr_result.max.y < mbr.max.y){
                    mbr_result.max.y = mbr.max.y;
                }
            }else{
                mbr_result = mbr;
                init_mbr = true;
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

            int64_t snap_id = snap_start;
            while(snap_id < snap_end){
                time_interval_left(snap_id, region_q, t_beg, t_end, processed_ids, rctIndex, r);
                t_beg = ++snap_id * rctIndex.period_snapshot;
            }

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

            int64_t snap_id = snap_start;
            while(snap_id < snap_end){
                time_interval_left_brute_force(snap_id, region_q, t_beg, t_end, processed_ids, rctIndex, r);
                t_beg = ++snap_id * rctIndex.period_snapshot;
            }

        }

        template<class RCTIndex>
        static void knn(const typename RCTIndex::size_type k, const util::geo::point& p_q,
                        const typename RCTIndex::size_type t_q, const RCTIndex &rctIndex,
                        std::vector<typename  RCTIndex::value_type> &r) {

            pq_knn_result_type pq_results, pq_candidates;
            pq_knn_movement_type pq_distances;
            auto snap_q = t_q / rctIndex.period_snapshot;
            const auto& snap =  rctIndex.snapshots[snap_q];
            //Init candidates
            pq_candidates.push(knn_element(snap.get_root(), p_q, util::geo::point{snap.get_lower0(), snap.get_lower1()}, util::geo::point{0, 0}));
            while(!pq_candidates.empty() && !(pq_results.size() >= k && pq_candidates.top().distance > pq_distances.top())){
                auto candidate = pq_candidates.top(); //copy is necessary
                pq_candidates.pop();
                if(candidate.is_point){
                    util::geo::point p;
                    if(search_object(candidate.id, t_q, rctIndex, p)){
                        knn_element knn_e(candidate.id, p_q, p);
                        insert_result(knn_e, pq_results, pq_distances, k);
                    }
                }else{
                    snap.enqueue_children(candidate, pq_candidates, p_q);
                }
            }
            double_t k_distance = 0;
            for (uint i = 0; i < k; i++) {
                r.emplace_back(pq_results.top().id);
                pq_results.pop();
            }


        }

        template<class RCTIndex>
        static bool MBR(const typename RCTIndex::size_type t_i, const typename RCTIndex::size_type t_j,
                        const typename RCTIndex::size_type oid,
                        const RCTIndex &rctIndex,
                        util::geo::region& mbr_result){

            auto traj_step = rctIndex.log_objects[oid].start_traj_step();
            typename RCTIndex::size_type t_start = traj_step.t;
            typename RCTIndex::size_type t_end = rctIndex.log_objects[oid].time_end();
            if (t_end < t_i || t_start > t_j) return false;

            bool mbr_init = false;
            auto aux_i = t_i;
            auto aux_j = t_j;
            if(aux_i <= t_start){
                aux_i = t_start+1;
                mbr_result = util::geo::region{traj_step.x, traj_step.y, traj_step.x, traj_step.y};
                mbr_init = true;
            };
            if(aux_j > t_end) aux_j = t_end;

            typename RCTIndex::size_type movement_i = 0, movement_j = 0, c_phrase_i = 0, c_phrase_j = 0,
                    ic_phrase_l = 0, ic_phrase_r = 0, delta_phrase_l = 0, delta_phrase_r = 0;

            rctIndex.log_objects[oid].time_to_movement(aux_i, aux_j, movement_i, movement_j);
            if (movement_i > 0 && movement_i <= movement_j) {
                //1. Phrases which involve the queried interval
                rctIndex.log_objects[oid].interval_phrases(movement_i, movement_j, c_phrase_i, c_phrase_j,
                                                           ic_phrase_l, delta_phrase_l,
                                                           ic_phrase_r, delta_phrase_r);

                //2. MBR of the contained phrases
                if(c_phrase_i <= c_phrase_j){
                    auto mbr = rctIndex.log_objects[oid].MBR(c_phrase_i, c_phrase_j);
                    update_MBR(mbr_result, mbr, mbr_init);
                }
                //3. MBR of the incomplete phrases
                if(delta_phrase_l || ic_phrase_l == ic_phrase_r){
                    //Check the first incomplete phrase
                    auto mbr_phrase = rctIndex.log_objects[oid].MBR(ic_phrase_l);
                    if(!mbr_init || !util::geo::contains(mbr_result, mbr_phrase)){
                        auto phrase_start = rctIndex.log_objects[oid].offsets[ic_phrase_l-1];
                        auto point_phrase = rctIndex.log_objects[oid].phrase_point(ic_phrase_l);
                        auto last_movement = std::min(rctIndex.log_objects[oid].last_movement(ic_phrase_l), movement_j);
                        auto move_ref_start =  phrase_start + delta_phrase_l+1; //count
                        auto move_ref_end = move_ref_start + (last_movement - movement_i); //count
                        auto mbr = rctIndex.log_reference.MBR(point_phrase.x, point_phrase.y, phrase_start, move_ref_start, move_ref_end);
                        update_MBR(mbr_result, mbr, mbr_init);
                    }

                }

                if(delta_phrase_r && ic_phrase_l < ic_phrase_r){
                    //Check the last incomplete phrase
                    auto mbr_phrase = rctIndex.log_objects[oid].MBR(ic_phrase_r);
                    if(!mbr_init || !util::geo::contains(mbr_result, mbr_phrase)){
                        auto phrase_start = rctIndex.log_objects[oid].offsets[ic_phrase_r-1];
                        auto point_phrase = rctIndex.log_objects[oid].phrase_point(ic_phrase_r);
                        auto start_movement = rctIndex.log_objects[oid].start_movement(ic_phrase_r);
                        auto move_ref_end = phrase_start+1 + (movement_j - start_movement); //count
                        auto mbr = rctIndex.log_reference.MBR(point_phrase.x, point_phrase.y, phrase_start, phrase_start+1, move_ref_end);
                        update_MBR(mbr_result, mbr, mbr_init);
                    }
                }
                return true;

            }
            //There is no movement, but if mbr_init=true, the interval of time considers the first position
            return mbr_init;

        }

        template<class RCTIndex>
        static void knn_trajectory(const typename RCTIndex::size_type k,
                                   const std::vector<uint32_t>& traj_x,
                                   const std::vector<uint32_t>& traj_y,
                                   const typename RCTIndex::size_type t_b,
                                   const typename RCTIndex::size_type t_e,
                                   const RCTIndex &rctIndex,
                                   std::vector<typename  RCTIndex::value_type> &r) {


            typedef std::unordered_map<typename RCTIndex::value_type, size_type> map_obj_heap_type;
            typedef struct {
                size_type t_b, t_e, ptr_bt, snap_id;
            } object_info_type;
            typedef util::heap_max_min<object_info_type, size_type> pq_object_type;

            //1. Building rmMq
            rmMq<uint32_t > rmMq_x(&traj_x);
            rmMq<uint32_t > rmMq_y(&traj_y);

            //Add the roots to the global heap and build the binary trees
            pq_knn_traj_element_type pq_global;
            std::vector<::util::geo::region> root_regions;
            size_type first = t_b / rctIndex.period_snapshot;
            size_type last = t_e / rctIndex.period_snapshot;
            for(size_type i = first; i <= last; ++i){
                //Building binary tree
                auto bt_begin = std::max(i*rctIndex.period_snapshot, t_b);
                auto bt_end = std::min((i+1)*rctIndex.period_snapshot-1, t_e);
                auto min_max_x = rmMq_x(bt_begin-t_b, bt_end-t_b);
                auto min_max_y = rmMq_y(bt_begin-t_b, bt_end-t_b);
                ::util::geo::region reg{min_max_x.min, min_max_y.min, min_max_x.max, min_max_y.max};
                root_regions.push_back(reg);
                //Adding the root
                const auto& snap =  rctIndex.snapshots[i];
                knn_traj_element root(snap.get_root(), reg, util::geo::point{snap.get_lower0(), snap.get_lower1()}, util::geo::point{0, 0}, i);
                snap.enqueue_children(root, pq_global, reg);
            }

            pq_knn_result_type pq_results;
            pq_knn_movement_type pq_distances;
            map_obj_heap_type object_to_heap;
            std::vector<pq_object_type> pq_objects;
            while(!pq_global.empty() && !(pq_results.size() >= k && pq_global.top().distance >= pq_distances.top())){
                auto candidate = pq_global.top();
                //std::cout << "New candidate: " << candidate.id << std::endl;
                pq_global.pop();
                //Header
                if(candidate.is_header){
                    //std::cout << "Search in pq_object" << std::endl;
                    //Search in the priority queue of the object
                    if(candidate.distance == candidate.max_distance){
                        knn_element knn_e(candidate.id, candidate.max_distance);
                        insert_result(knn_e, pq_results, pq_distances, k);
                    }else{
                        const auto it = object_to_heap.find(candidate.id);
                        auto &pq_object = pq_objects[it->second];
                        if (!pq_object.empty()) {
                            auto info = pq_object.top();
                            pq_object.pop();
                            if(info.t_e == info.t_b){
                                knn_element knn_e(candidate.id, candidate.max_distance);
                                insert_result(knn_e, pq_results, pq_distances, k);
                            }else if (info.t_e - info.t_b == 1) {
                                util::geo::point pb, pe;
                                if (search_object(candidate.id, info.t_b, rctIndex, pb)) {
                                    util::geo::point p{traj_x[info.t_b-t_b], traj_y[info.t_b-t_b]};
                                    auto dist = knn_support_helper_rtree::dmin(p, pb, pb);
                                    object_info_type info1{info.t_b, info.t_b, 0, info.snap_id};
                                    pq_object.push(info1, dist, dist);
                                };
                                if (search_object(candidate.id, info.t_e, rctIndex, pe)) {
                                    util::geo::point p{traj_x[info.t_e-t_b], traj_y[info.t_e-t_b]};
                                    auto dist = knn_support_helper_rtree::dmin(p, pe, pe);
                                    object_info_type info2{info.t_e, info.t_e, 0, info.snap_id};
                                    pq_object.push(info2, dist, dist);
                                };
                                auto min_max = pq_object.top_min_max();
                                pq_global.push(knn_traj_element(candidate.id, min_max.first, min_max.second));
                            } else {
                                //std::cout << ">1 element" << std::endl;
                                auto mid = (info.t_e - info.t_b) / 2 + info.t_b;
                                util::geo::region mbr_1, mbr_2;
                                distance_type dmin1, dmax1, dmin2, dmax2;
                                bool info1_ok = false, info2_ok = false;
                                object_info_type info1, info2;


                                //uto &bt = bts[info.snap_id - first];
                                if (MBR(info.t_b, mid, candidate.id, rctIndex, mbr_1)) {
                                    auto min_max_x = rmMq_x(info.t_b-t_b, mid-t_b);
                                    auto min_max_y = rmMq_y(info.t_b-t_b, mid-t_b);
                                    ::util::geo::region reg{min_max_x.min, min_max_y.min, min_max_x.max, min_max_y.max};
                                    dmin1 = knn_support_helper_rtree::dmin(reg, mbr_1.min, mbr_1.max);
                                    dmax1 = knn_support_helper_rtree::dmax(reg, mbr_1.min, mbr_1.max);
                                    info1 = object_info_type{info.t_b, mid, 0, info.snap_id};
                                    info1_ok = true;

                                }
                                if (MBR(mid + 1, info.t_e, candidate.id, rctIndex, mbr_2)) {
                                    auto min_max_x = rmMq_x(mid+1-t_b, info.t_e-t_b);
                                    auto min_max_y = rmMq_y(mid+1-t_b, info.t_e-t_b);
                                    ::util::geo::region reg{min_max_x.min, min_max_y.min, min_max_x.max, min_max_y.max};
                                    dmin2 = knn_support_helper_rtree::dmin(reg, mbr_2.min, mbr_2.max);
                                    dmax2 = knn_support_helper_rtree::dmax(reg, mbr_2.min, mbr_2.max);
                                    info2 = object_info_type{mid + 1, info.t_e, 0, info.snap_id};
                                    info2_ok = true;
                                }
                                //std::cout << "Check " << std::endl;
                                if (info1_ok && info2_ok) {
                                    if (dmin1 < dmax2) pq_object.push(info2, dmin2, dmax2);
                                    if (dmin2 < dmax1) pq_object.push(info1, dmin1, dmax1);
                                } else if (info1_ok) {
                                    pq_object.push(info1, dmin1, dmax1);
                                } else if (info2_ok) {
                                    pq_object.push(info2, dmin2, dmax2);
                                }
                                auto min_max = pq_object.top_min_max();
                                //std::cout << "Top min max" << std::endl;
                                pq_global.push(knn_traj_element(candidate.id, min_max.first, min_max.second));
                            }
                        }
                    }

                }else if(candidate.is_leaf){
                    const auto it = object_to_heap.find(candidate.id);
                    if(it == object_to_heap.end()){
                        pq_objects.push_back(pq_object_type());
                        auto snap_id = first;
                        while(snap_id <= last){
                            util::geo::region reg, mbr;
                            reg = root_regions[snap_id - first];
                            //bts[snap_id - first].get_region(0, reg);
                            auto b = std::max(snap_id*rctIndex.period_snapshot, t_b);
                            auto e = std::min((snap_id+1)*rctIndex.period_snapshot-1, t_e);
                            if(MBR(b, e, candidate.id, rctIndex, mbr)){
                                //std::cout << "MBR" << std::endl;
                                object_info_type info{b, e, 0, snap_id};
                                auto dmin = knn_support_helper_rtree::dmin(reg, mbr.min, mbr.max);
                                auto dmax = knn_support_helper_rtree::dmax(reg, mbr.min, mbr.max);
                                pq_objects[pq_objects.size()-1].push(info, dmin, dmax);
                            }
                            //std::cout << "Skip MBR" << std::endl;
                            ++snap_id;
                        }
                        if(!pq_objects[pq_objects.size()-1].empty()){
                            auto min_max = pq_objects[pq_objects.size()-1].top_min_max();
                            pq_global.push(knn_traj_element(candidate.id, min_max.first, min_max.second));
                        }
                        object_to_heap.insert({candidate.id, pq_objects.size()-1});
                    }

                }else{
                    util::geo::region reg = root_regions[candidate.snap_id - first];
                    rctIndex.snapshots[candidate.id].enqueue_children(candidate, pq_global, reg);
                }
            }

            for (uint i = 0; i < k; i++) {
                r.emplace_back(pq_results.top().id);
                pq_results.pop();
            }


        }
    };


}

#endif //RCT_SEARCH_OBJECT_SUPPORT_HPP
