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
// Created by Adrián on 09/03/2019.
//

#ifndef RCT_KNN_SUPPORT_HELPER_RTREE_HPP
#define RCT_KNN_SUPPORT_HELPER_RTREE_HPP

#include <geo_util.hpp>
#include <queue>
#include <SNode.h>

namespace rct {

    namespace knn_support_helper_rtree {



        typedef uint64_t size_type;
        typedef uint64_t distance_type;

        static distance_type dmin(const util::geo::point &_pq, const util::geo::point &_min, const util::geo::point &_max){
            distance_type x = 0, y = 0;
            if (_pq.x < _min.x) {
                x = _min.x - _pq.x;
            } else if (_pq.x > _max.x){
                x = _pq.x - _max.x;
            }
            if (_pq.y < _min.y) {
                y = _min.y - _pq.y;
            } else if(_pq.y > _max.y){
                y = _pq.y - _max.y;
            }
            return x * x + y * y;
        }

        static distance_type dmin(const util::geo::region &_rq, const util::geo::point &_min, const util::geo::point &_max){
            distance_type x = 0, y = 0;
            if (_rq.min.x > _max.x) {
                x = _rq.min.x - _max.x;
            } else if (_rq.max.x < _min.x){
                x = _min.x - _rq.max.x;
            }
            if (_rq.min.y > _max.y) {
                y = _rq.min.y - _max.y;
            } else if (_rq.max.y < _min.y){
                y = _min.y - _rq.max.y;
            }
            return x * x + y * y;
        }



        static distance_type dmax(const util::geo::point &_pq, const util::geo::point &_min, const util::geo::point &_max){
            distance_type x = 0, y = 0;
            if (_pq.x < _min.x + (_max.x - _min.x) / 2) {
                x = _max.x - _pq.x;
            } else {
                x = _pq.x - _min.x;
            }
            if(_pq.y < _min.y + (_max.y - _min.y)/2){
                y = _max.y - _pq.y;
            }else{
                y = _pq.y - _min.y;
            }
            return x * x + y * y;
        }

        static distance_type dmax(const util::geo::region &_rq, const util::geo::point &min, const util::geo::point &max){
            distance_type x = 0, y = 0;
            if (_rq.min.x > max.x) {
                x = _rq.max.x - min.x;
            } else if(_rq.max.x < min.x){
                x = max.x - _rq.min.x;
            }else{
                distance_type x1, x2;
                if(_rq.max.x > min.x){
                    x1 = _rq.max.x - min.x;
                }else{
                    x1 = min.x - _rq.max.x;
                }
                if(max.x > _rq.min.x){
                    x2 = max.x - _rq.min.x;
                }else{
                    x2 = _rq.min.x - max.x;
                }
                x = std::max(x1, x2);
            }
            if (_rq.min.y > max.y) {
                y = _rq.max.y - min.y;
            } else if(_rq.max.y < min.y){
                y = max.y - _rq.min.y;
            }else{
                distance_type y1, y2;
                if(_rq.max.y > min.y){
                    y1 = _rq.max.y - min.y;
                }else{
                    y1 = min.y - _rq.max.y;
                }
                if(max.y > _rq.min.y){
                    y2 = max.y - _rq.min.y;
                }else{
                    y2 = _rq.min.y - max.y;
                }
                y = std::max(y1, y2);
            }
            return x * x + y * y;
        }

        struct knn_element {

            bool is_point;
            size_type id;
            SRTree::SNode *ptr = nullptr;
            util::geo::point min;
            util::geo::point max;
            distance_type distance;
            distance_type max_distance;

            //size_type t;

            void common_knn_element(const util::geo::point &_pq, const util::geo::point &_min,
                                    const util::geo::point &_max){
                min = _min;
                max = _max;
                distance = dmin(_pq, min, max);
                max_distance = dmax(_pq, min, max);
            }

            knn_element(size_type _id, const util::geo::point &_p, distance_type _dist) {
                id = _id;
                min = _p;
                distance = _dist;
                // t = _t;
            }

            knn_element(size_type _id, distance_type _dist) {
                id = _id;
                distance = _dist;
                // t = _t;
            }

            /**** NOVOS METODOS ****/
            knn_element(const size_type _id, const util::geo::point &_pq, const util::geo::point &_p) {
                //std::cout << "KNN element: " << id << "p_q:" << _pq << std::endl;
                is_point = true;
                id = _id;
                //t = _t;
                distance_type x = 0, y = 0;
                if(_pq.x < _p.x){
                    x = _p.x - _pq.x;
                }else if(_pq.x > _p.x){
                    x =_pq.x - _p.x;
                }
                if(_pq.y < _p.y){
                    y = _p.y - _pq.y;
                }else if(_pq.y > _p.y){
                    y = _pq.y - _p.y;
                }
                distance = x*x + y*y;
            }

            knn_element(const size_type _id, const util::geo::point &_pq, const util::geo::point &_min,
                        const util::geo::point &_max) {
                is_point = true;
                ptr = nullptr;
                id = _id;
                common_knn_element(_pq, _min, _max);
            }

            knn_element(SRTree::SNode *_ptr, const util::geo::point &_pq, const util::geo::point &_min,
                        const util::geo::point &_max) {
                is_point = false;
                ptr = _ptr;
                id = 0;
                common_knn_element(_pq, _min, _max);
            }
        };

        struct knn_traj_element{
            bool is_header = false;
            bool is_leaf = false;
            size_type id;
            size_type snap_id;
            SRTree::SNode *ptr = nullptr;
            util::geo::point min;
            util::geo::point max;
            distance_type distance;
            distance_type max_distance;

            void common_knn_element(const util::geo::point &_pq, const util::geo::point &_min,
                                    const util::geo::point &_max){
                min = _min;
                max = _max;
                distance = dmin(_pq, _min, _max);
                max_distance = dmax(_pq, _min, _max);
            }
            void common_knn_element(const util::geo::region &_rq, const util::geo::point &_min,
                                    const util::geo::point &_max){
                min = _min;
                max = _max;
                distance = dmin(_rq, _min, _max);
                max_distance = dmax(_rq, _min, _max);
            }

            knn_traj_element(const size_type _id, const distance_type dmin,
                    const distance_type dmax){
                is_header = true;
                id = _id;
                distance = dmin;
                max_distance = dmax;

            }

            knn_traj_element(const size_type _id, const util::geo::point &_rq, const util::geo::point &_min,
                             const util::geo::point &_max, const size_type _snap = 0) {
                is_leaf = true;
                ptr = nullptr;
                id = _id;
                snap_id = _snap;
                common_knn_element(_rq, _min, _max);
            }

            knn_traj_element(const size_type _id, const util::geo::region &_rq, const util::geo::point &_min,
                             const util::geo::point &_max, const size_type _snap = 0) {
                is_leaf = true;
                ptr = nullptr;
                id = _id;
                snap_id = _snap;
                common_knn_element(_rq, _min, _max);
            }

            knn_traj_element(SRTree::SNode *_ptr, const util::geo::region &_rq,
                            const util::geo::point &_min, const util::geo::point &_max, const size_type _snap = 0) {
                is_leaf = false;
                ptr = _ptr;
                snap_id = _snap;
                common_knn_element(_rq, _min, _max);
            }

            knn_traj_element(SRTree::SNode *_ptr, const util::geo::point &_rq,
                             const util::geo::point &_min, const util::geo::point &_max, const size_type _snap = 0) {
                is_leaf = false;
                ptr = _ptr;
                snap_id = _snap;
                common_knn_element(_rq, _min, _max);
            }



        };


        class knn_element_comparator {
        public:
            bool operator()(const knn_element &lce, const knn_element &rce){
                if(lce.distance == rce.distance){
                    return lce.max_distance > rce.max_distance;
                }
                return lce.distance > rce.distance;
            }
        };

        class knn_result_comparator { //min
        public:
            bool operator()(const knn_element &lce, const knn_element &rce){
                if(lce.distance == rce.distance){
                    return lce.max_distance > rce.max_distance;
                }
                return lce.distance > rce.distance;
            }
        };

        class knn_traj_element_comparator { //min
        public:
            bool operator()(const knn_traj_element &lce, const knn_traj_element &rce){
                if(lce.distance == rce.distance){
                    return lce.max_distance > rce.max_distance;
                }
                return lce.distance >  rce.distance;
            }
        };

        typedef std::priority_queue<knn_element, std::vector<knn_element>, knn_result_comparator> pq_knn_result_type;
        typedef std::priority_queue<knn_traj_element, std::vector<knn_traj_element>, knn_traj_element_comparator> pq_knn_traj_element_type;
        typedef std::priority_queue<uint64_t , std::vector<uint64_t >, std::less<uint64_t >> pq_knn_movement_type;


    };

}

#endif //GRACT_KNN_SUPPORT_HELPER_HPP
