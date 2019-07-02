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
#include <iostream>
#include <cmath>

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

        struct mbr{
            movement min;
            movement max;
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


        std::ostream& operator<<(std::ostream& os, const point& pr)
        {
            return os << "(" << pr.x << ", " << pr.y << ")";
        }

        std::ostream& operator<<(std::ostream& os, const region& pr)
        {
            return os << "[" << pr.min << ", " << pr.max << "]";
        }

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

        inline region expand(const point &p,
                             const uint32_t max_speed, const uint32_t time,
                             const uint32_t max_x, const uint32_t max_y){

            auto distance = max_speed * time;
            region ret;
            if(p.x < distance){
                ret.min.x = 0;
            }else{
                ret.min.x = p.x - distance;
            }
            if(p.y < distance){
                ret.min.y = 0;
            }else{
                ret.min.y = p.y - distance;
            }

            if(p.x + distance > max_x){
                ret.max.x = max_x;
            }else{
                ret.max.x = p.x + distance;
            }
            if(p.y + distance > max_y){
                ret.max.y = max_y;
            }else{
                ret.max.y = p.y + distance;
            }
            return ret;
        }

        inline region expand(const region &r,
                             const uint32_t max_speed, const uint32_t time,
                             const uint32_t max_x, const uint32_t max_y){

            auto distance = max_speed * time;
            region ret;
            if(r.min.x < distance){
                ret.min.x = 0;
            }else{
                ret.min.x = r.min.x - distance;
            }
            if(r.min.y < distance){
                ret.min.y = 0;
            }else{
                ret.min.y = r.min.y - distance;
            }

            if(r.max.x + distance > max_x){
                ret.max.x = max_x;
            }else{
                ret.max.x = r.max.x + distance;
            }
            if(r.max.y + distance > max_y){
                ret.max.y = max_y;
            }else{
                ret.max.y = r.max.y + distance;
            }
            return ret;
        }

        inline double distance(const region &r, const point &p){
            double x = 0, y = 0;
            if(p.x < r.min.x){
                x = r.min.x - p.x;
            }else if(p.x > r.max.x){
                x =p.x - r.max.x;
            }
            if(p.y < r.min.y){
                y = r.min.y - p.y;
            }else if(p.y > r.max.y){
                y = p.y - r.max.y;
            }
            return std::sqrt(x * x + y * y);
        }

        inline double distance(const point &p, const point &p1){
            double x = 0, y = 0;
            x = std::abs((int32_t) (p1.x - p.x));
            y = std::abs((int32_t) (p1.y - p.y));
            return std::sqrt(x * x + y * y);
        }
    }
}

#endif //RCT_GEO_UTIL_HPP
