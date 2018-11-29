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

namespace rct {


    namespace algorithm {

        /***
         * Everything Done.
            Error looking for: id=0 t=82
            Expected: 893, 273196
            Obtained: 892, 273195

         * @tparam RCTIndex
         * @param oid
         * @param t_q
         * @param rctIndex
         * @param r
         * @return
         */

        template<class RCTIndex>
        bool search_object(const typename RCTIndex::size_type oid, const typename RCTIndex::value_type t_q,
                           const RCTIndex &rctIndex, util::geo::point &r){

            auto traj_step = rctIndex.log_objects[oid].start_traj_step();
            if(traj_step.t > t_q) return false;
            if(traj_step.t == t_q){
                r = util::geo::point{traj_step.x, traj_step.y};
                return true;
            }
            typename RCTIndex::size_type idx_beg = 0, idx_end = 0;
            util::geo::movement movement_phrase;
            if(rctIndex.log_objects[oid].interval_ref(t_q, idx_beg, idx_end, movement_phrase)){
                util::geo::movement movement = rctIndex.log_reference.compute_movement(idx_beg, idx_end);
                r = util::geo::point{traj_step.x + movement_phrase.x + movement.x,
                                     traj_step.y + movement_phrase.y + movement.y};
                return true;
            };
            return false;

        }
    };
}

#endif //RCT_SEARCH_OBJECT_SUPPORT_HPP
