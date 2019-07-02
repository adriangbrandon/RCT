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

        struct knn_element {

            bool is_point;
            size_type id;
            SRTree::SNode *ptr = nullptr;
            util::geo::point min;
            util::geo::point max;
            size_type weight;
            distance_type distance;

            //size_type t;

            void common_knn_element(util::geo::point _pq, util::geo::point _min,
                                    util::geo::point _max){
                min = _min;
                max = _max;
                distance_type x = 0, y = 0;
                if (_pq.x < min.x) {
                    x = min.x - _pq.x;
                } else if (_pq.x > max.x){
                    x = _pq.x - max.x;
                }
                if (_pq.y < min.y) {
                    y = min.y - _pq.y;
                } else if(_pq.y > max.y){
                    y = _pq.y - max.y;
                }
                distance = x * x + y * y;
            }

            knn_element(size_type _id, util::geo::point _p, distance_type _dist) {
                id = _id;
                min = _p;
                distance = _dist;
                // t = _t;
                weight = (size_type) std::ceil(_dist);
            }

            /**** NOVOS METODOS ****/
            knn_element(size_type _id, util::geo::point _pq, util::geo::point _p) {
                //std::cout << "KNN element: " << id << "p_q:" << _pq << std::endl;
                is_point = true;
                id = _id;
                //t = _t;
                auto d_x = (distance_type) std::abs((int32_t) (_pq.x - _p.x));
                auto d_y = (distance_type) std::abs((int32_t) (_pq.y - _p.y));
                distance = d_x * d_x + d_y * d_y;
            }

            knn_element(size_type _id, util::geo::point _pq, util::geo::point _min,
                        util::geo::point _max) {
                is_point = true;
                ptr = nullptr;
                id = _id;
                common_knn_element(_pq, _min, _max);
            }

            knn_element(SRTree::SNode *_ptr, util::geo::point _pq, util::geo::point _min,
                        util::geo::point _max) {
                is_point = false;
                ptr = _ptr;
                id = 0;
                common_knn_element(_pq, _min, _max);
            }
        };

        class knn_element_comparator {
        public:
            bool operator()(const knn_element &lce, const knn_element &rce){
                if(lce.weight == rce.weight){
                    return lce.distance > rce.distance;
                }
                return lce.weight > rce.weight;
            }
        };

        class knn_result_comparator {
        public:
            bool operator()(const knn_element &lce, const knn_element &rce){
                return lce.distance > rce.distance;
            }
        };

        typedef std::priority_queue<knn_element, std::vector<knn_element>, knn_result_comparator> pq_knn_result_type;
        typedef std::priority_queue<knn_element, std::vector<knn_element>, knn_element_comparator> pq_knn_element_type;
        typedef std::priority_queue<uint64_t , std::vector<uint64_t >, std::less<uint64_t >> pq_knn_movement_type;


    };

}

#endif //GRACT_KNN_SUPPORT_HELPER_HPP
