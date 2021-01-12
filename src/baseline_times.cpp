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

#include <string>
#include <rct_index_rtree.hpp>
#include <rct_algorithm_rtree.hpp>
#include <baseline.hpp>

#define TIMES 1

using namespace std;
using namespace std::chrono;

int main(int argc, const char **argv) {


    size_t first_query_arg = 5;
    std::string dataset = argv[1];
    uint32_t t_queries = std::atoi(argv[2]);
    std::string queries = argv[3];

    rct::baseline m_baseline(dataset);


    std::vector<uint32_t> ids, ks;
    std::vector<uint32_t> t_starts;
    std::vector<uint32_t> t_ends;
    std::vector<util::geo::region> regions;
    std::vector<util::geo::point> points;
    int64_t type;
    uint32_t k, id, maxX, maxY, minX, minY, tstart, tend;

    if(t_queries == 0){
        ifstream finQ(queries);
        finQ >> type >> tstart >> tend >> minX >> minY >> k;
        t_starts.push_back(tstart);
        t_ends.push_back(tend);
        points.push_back(util::geo::point{minX, minY});
        ks.push_back(k);
        finQ >> type;
        while (finQ) {
            finQ >> type >> tstart >> tend >> minX >> minY >> k >> type;
            t_starts.push_back(tstart);
            t_ends.push_back(tend);
            points.push_back(util::geo::point{minX, minY});
            ks.push_back(k);
        }
        finQ.close();
        double_t t_find_knn_int = 0;
        for(uint64_t j = 0; j < TIMES; j++){
            auto start = util::time::user::now();
            for(uint64_t i = 0; i < t_starts.size(); i++){
                std::vector<uint32_t > res_knn;
                m_baseline.knn_int(ks[i], points[i].x, points[i].y, t_starts[i], t_ends[i], res_knn);
            }
            auto stop = util::time::user::now();
            auto milli = stop-start;
            double_t milli_query = milli/(double_t) t_starts.size();
            cout << "Time (µs) = " << milli << std::endl;
            cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
            t_find_knn_int += milli_query;
            sleep(1);
        }
        double_t avg_find_knn_int = t_find_knn_int/(double) TIMES;
        std::cout << "KNN interval (µs): " << t_find_knn_int << std::endl;
    }else{
        ifstream finQ(queries);
        finQ >> type >> tstart >> tend;
        auto length = tend - tstart + 1;
        std::vector<std::vector<uint32_t>> traj_x(1), traj_y(1);
        uint64_t c_query = 0;
        for (auto i = 0; i < length; ++i) {
            finQ >> minX >> minY;
            traj_x[c_query].push_back(minX);
            traj_y[c_query].push_back(minY);
        }
        finQ >> k >> type;
        t_starts.push_back(tstart);
        t_ends.push_back(tend);
        ks.push_back(k);
        while (finQ) {
            ++c_query;
            traj_x.emplace_back();
            traj_y.emplace_back();
            finQ >> type >> tstart >> tend;
            length = tend - tstart + 1;
            for (auto i = 0; i < length; ++i) {
                finQ >> minX >> minY;
                traj_x[c_query].push_back(minX);
                traj_y[c_query].push_back(minY);
            }
            finQ >> k >> type;
            t_starts.push_back(tstart);
            t_ends.push_back(tend);
            ks.push_back(k);
        }
        finQ.close();
        std::vector<double_t> knn_traj_times;
        //for(uint64_t steps_knn = 0; steps_knn <= 0; ++steps_knn){
        double_t t_find_knn_traj = 0;
        //uint64_t s = steps_knn;
        //if(s == 0) s = 1000;
        for(uint64_t j = 0; j < TIMES; j++){
            auto start = ::util::time::user::now();
            for(uint64_t i = 0; i < ks.size(); i++){
                std::vector<uint32_t > res_knn;
                m_baseline.knn_traj(ks[i], traj_x[i], traj_y[i], t_starts[i], t_ends[i], res_knn);
            }
            auto stop = ::util::time::user::now();
            auto milli = stop-start;
            double_t milli_query = milli/(double_t) t_starts.size();
            cout << "Time (µs) = " << milli << std::endl;
            cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
            t_find_knn_traj += milli_query;
            sleep(1);
        }
        double_t avg_find_knn_traj = t_find_knn_traj/(double) TIMES;
        std::cout << "KNN trajectory (µs): " << t_find_knn_traj << std::endl;
    }

    std::cout << "Everything is OK!" << std::endl;
//    dataset.clear();

}
