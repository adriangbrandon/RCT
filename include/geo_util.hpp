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

#ifndef GEO_UTIL_HPP
#define GEO_UTIL_HPP

#include <cstdint>

namespace util {

    namespace geo {

        struct point {
            uint32_t x;
            uint32_t y;
        };

        struct movement {
            int32_t x;
            int32_t y;
        };

        struct region {
            point min;
            point max;
        };

        struct traj_step {
            uint32_t t;
            uint32_t x;
            uint32_t y;
        };

        struct id_point {
            uint32_t id;
            uint32_t x;
            uint32_t y;
        };

        inline bool contains(const region &r_q, const point &p){
            return (r_q.min.x <= p.x && p.x <= r_q.max.x && r_q.min.y <= p.y && p.y <= r_q.max.y);
        }

        inline bool contains(const region &r_q, const region &r){
            return (r_q.min.x <= r.min.x && r.min.x <= r_q.max.x &&
                    r_q.min.y <= r.min.y && r.min.y <= r_q.max.y &&
                    r_q.min.x <= r.max.x && r.max.x <= r_q.max.x &&
                    r_q.min.y <= r.max.y && r.max.y <= r_q.max.y);
        }

        inline bool touches(const region &r_q, const region &r){
                return !(r.min.x > r_q.max.x || r.max.x < r_q.min.x
                        || r.min.y > r_q.max.y || r.max.y < r_q.min.y);
        }
    }
}

#endif //RCT_GEO_UTIL_HPP
