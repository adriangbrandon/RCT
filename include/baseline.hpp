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
// Created by Adrián on 12/1/21.
//

#ifndef RCT_BASELINE_HPP
#define RCT_BASELINE_HPP

#include <queue>
#include <vector>
#include <stdint.h>


namespace rct {

    class baseline {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef struct {
            value_type x;
            value_type y;
            value_type t;
        } tp_type;

    private:
        typedef struct {
            value_type oid;
            size_type distance;
        } oid_dist_type;

        class oid_dist_min { //min
        public:
            bool operator()(const oid_dist_type &lce, const oid_dist_type &rce){
                return lce.distance >  rce.distance;
            }
        };

        typedef std::priority_queue<oid_dist_type , std::vector<oid_dist_type>, oid_dist_min> pq_type;

        std::vector<std::vector<tp_type>> m_trajectories;

        static size_type dmin(const size_type &x1, const size_type &y1,
                              const size_type &x2, const size_type &y2){
            distance_type x = 0, y = 0;
            if (x1 < x2) {
                x = x2 - x1;
            } else if (x1 > x2){
                x = x1 - x2;
            }
            if (y1 < y2) {
                y = y2 - y1;
            } else if(y1 > y2){
                y = y1 - y2;
            }
            return x * x + y * y;
        }

    public:

        baseline(const std::string &file_name){
            std::ifstream input(file_name);
            value_type id, t, x, y;
            value_type old_id = -1;
            while (input) {
                input >> id >> t >> x >> y;
                if (input.eof()) break;
                if(id != old_id){
                    m_trajectories.emplace_back(std::vector<tp_type>());
                }

                m_trajectories[id].push_back(tp_type{x, y, t});
                old_id = id;
            }
            input.close();
        }

        void knn_traj(const size_type k, const std::vector<value_type>& traj_x,
                 const std::vector<value_type>& traj_y, const size_type t_b,
                 const size_type t_e, std::vector<value_type> &r){

            pq_type pq;
            for(value_type oid = 0; oid < m_trajectories.size(); ++oid){
                size_type max_distance = 0;
                bool add = false;
                for(size_type i = 0; i < m_trajectories[oid].size(); ++i){
                    const auto &tp = m_trajectories[oid][i];
                    if(tp.t > t_e) {
                        break;
                    }
                    if(tp.t >= t_b){
                        auto dis = dmin(tp.x, tp.y, traj_x[tp.t - t_b], traj_y[tp.t - t_b]);
                        if(dis > max_distance){
                            max_distance = dis;
                        }
                        add = true;
                    }

                }
                if(add){
                    oid_dist_type oid_dist{oid, max_distance};
                    pq.push(oid_dist);
                }

            }
            while(!pq.empty() && r.size() < k){
                r.push_back(pq.top().oid);
                pq.pop();
            }
        }

        void knn_int(const size_type k, const  value_type x,
                 const value_type y, const size_type t_b,
                 const size_type t_e, std::vector<value_type> &r){

            pq_type pq;
            for(value_type oid = 0; oid < m_trajectories.size(); ++oid){
                bool add = false;
                size_type min_distance = INT64_MAX;
                for(size_type i = 0; i < m_trajectories[oid].size(); ++i){
                    const auto &tp = m_trajectories[oid][i];
                    if(tp.t > t_e) {
                        break;
                    }
                    if(tp.t >= t_b){
                        auto dis = dmin(tp.x, tp.y, x, y);
                        if(dis < min_distance){
                            min_distance = dis;
                        }
                        add = true;
                    }
                }
                if(add){
                    oid_dist_type oid_dist{oid, min_distance};
                    pq.push(oid_dist);
                }

            }
            while(!pq.empty() && r.size() < k){
                r.push_back(pq.top().oid);
                pq.pop();
            }
        }


        size_type size_in_bytes(){
            size_type bytes = 0;
            for(size_type i = 0; i < m_trajectories.size(); ++i){
                bytes += sizeof(tp_type)*m_trajectories[i].size() + sizeof(size_type);
            }
            return bytes;
        }

    };



}

#endif //RCT_BASELINE_HPP
