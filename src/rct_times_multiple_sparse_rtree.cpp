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

#define TIMES 1

using namespace std;
using namespace std::chrono;

int main(int argc, const char **argv) {

    std::string directory_queries = argv[4];

    if(!::util::file::end_slash(directory_queries)){
        directory_queries = directory_queries + "/";
    }

    std::vector<std::string> file_queries = {"so.txt",
                                             "traj.txt",
                                             "ts_s.txt",
                                             "ts_l.txt",
                                             "ti_s.txt",
                                             "ti_l.txt",
                                             "mbr.txt",
                                             "knn.txt",
                                             "knn_int.txt",
                                             "knn_traj.txt"};

    size_t first_query_arg = 4;
    std::string dataset_path = argv[1];
    double_t ratio = (double_t) atoi(argv[2])/(double_t) 100;
    uint32_t period = (uint32_t) atoi(argv[3]);
    std::string index_file = util::file::index_file("rct_index_multiple_sparse_rtree", argv, first_query_arg) + ".idx";
    std::cout << "Loading index: " << index_file << std::endl;
    rct::rct_index_rtree<rct::log_reference_sparse<>, rct::log_object_int_vector, rct::rlz_multiple_csa_bc_int64> m_rct_index;
    std::ifstream in(index_file);
    m_rct_index.load(in, dataset_path);
    in.close();
    std::cout << "Done" << std::endl;

    std::vector<uint32_t> ids, ks;
    std::vector<uint32_t> t_starts;
    std::vector<uint32_t> t_ends;
    std::vector<util::geo::region> regions;
    std::vector<util::geo::point> points;
    int64_t type;
    uint32_t id, maxX, maxY, minX, minY, tstart, tend;

    std::ifstream finQ(directory_queries + file_queries[0]);
    std::cout << directory_queries + file_queries[0] << std::endl;
    finQ >> type >> id >> tstart >> type;
    ids.push_back(id);
    t_starts.push_back(tstart);
    while (finQ) {
        finQ >> type >> id >> tstart >> type;
        ids.push_back(id);
        t_starts.push_back(tstart);
    }
    finQ.close();

    double_t t_object = 0;
    util::geo::point p;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            rct::algorithm::search_object(ids[i], t_starts[i], m_rct_index, p);
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_object += milli_query;
        sleep(1);
    }
    double_t avg_object = t_object/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.close();

    finQ.open(directory_queries + file_queries[1]);
    std::cout << directory_queries + file_queries[1] << std::endl;
    finQ >> type >> tstart >> tend >> id;
    ids.push_back(id);
    t_starts.push_back(tstart);
    t_ends.push_back(tstart + tend);
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> tend >> id >> type;
        ids.push_back(id);
        t_starts.push_back(tstart);
        t_ends.push_back(tstart+tend);
    }
    finQ.close();

    double_t t_traj = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<util::geo::traj_step> results;
            rct::algorithm::search_trajectory_fast(ids[i], t_starts[i], t_ends[i], m_rct_index, results);
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_traj += milli_query;
        sleep(1);
    }
    double_t avg_traj = t_traj/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();


    finQ.open(directory_queries + file_queries[2]);
    std::cout << directory_queries + file_queries[2] << std::endl;
    finQ >> type >> tstart >> minX >> maxX >> minY >> maxY;
    t_starts.push_back(tstart);
    regions.push_back(util::geo::region{util::geo::point{minX, minY}, util::geo::point{maxX, maxY}});
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> minX >> maxX >> minY >> maxY >> type;
        t_starts.push_back(tstart);
        regions.push_back(util::geo::region{util::geo::point{minX, minY}, util::geo::point{maxX, maxY}});
    }
    finQ.close();

    //std::ofstream log_ts("log_ts_RCT.txt");
    double_t t_slice_s = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<util::geo::id_point> res_t_s;
            rct::algorithm::time_slice(regions[i], t_starts[i], m_rct_index, res_t_s);
            //log_ts << " " << t_starts[i] << " " << regions[i] << "[results.size = " << res_t_s.size() << "]" << std::endl;
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_slice_s += milli_query;
        sleep(1);
    }
    double_t avg_slice_s = t_slice_s/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.open(directory_queries + file_queries[3]);
    std::cout << directory_queries + file_queries[3] << std::endl;
    finQ >> type >> tstart >> minX >> maxX >> minY >> maxY;
    t_starts.push_back(tstart);
    regions.push_back(util::geo::region{util::geo::point{minX, minY}, util::geo::point{maxX, maxY}});
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> minX >> maxX >> minY >> maxY >> type;
        t_starts.push_back(tstart);
        regions.push_back(util::geo::region{util::geo::point{minX, minY}, util::geo::point{maxX, maxY}});
    }
    finQ.close();

    double_t t_slice_l = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<util::geo::id_point> res_t_s;
            rct::algorithm::time_slice(regions[i], t_starts[i], m_rct_index, res_t_s);
            //log_ts << " " << t_starts[i] << " " << regions[i] << "[results.size = " << res_t_s.size() << "]" << std::endl;
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_slice_l += milli_query;
        sleep(1);
    }
    double_t avg_slice_l = t_slice_l/(double) TIMES;
    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();
    finQ.close();
    //log_ts.close();
    finQ.open(directory_queries + file_queries[4]);
    std::cout << directory_queries + file_queries[4] << std::endl;
    finQ >> type >> tstart >> tend >> minX >> maxX >> minY >> maxY;
    t_starts.push_back(tstart);
    t_ends.push_back(tend);
    regions.push_back(util::geo::region{minX, minY, maxX, maxY});
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> tend >> minX >> maxX >> minY >> maxY >> type;
        t_starts.push_back(tstart);
        t_ends.push_back(tend);
        regions.push_back(util::geo::region{minX, minY, maxX, maxY});
    }
    finQ.close();

    double_t t_interval_s = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_t_i;
            rct::algorithm::time_interval(regions[i], t_starts[i], t_ends[i], m_rct_index, res_t_i);
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_interval_s += milli_query;
        sleep(1);
    }
    double_t avg_interval_s = t_interval_s/(double) TIMES;


    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.open(directory_queries + file_queries[5]);
    std::cout << directory_queries + file_queries[5] << std::endl;
    finQ >> type >> tstart >> tend >> minX >> maxX >> minY >> maxY;
    t_starts.push_back(tstart);
    t_ends.push_back(tend);
    regions.push_back(util::geo::region{minX, minY, maxX, maxY});
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> tend >> minX >> maxX >> minY >> maxY >> type;
        t_starts.push_back(tstart);
        t_ends.push_back(tend);
        regions.push_back(util::geo::region{minX, minY, maxX, maxY});
    }
    finQ.close();

    double_t t_interval_l = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_t_i;
            rct::algorithm::time_interval(regions[i], t_starts[i], t_ends[i], m_rct_index, res_t_i);
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_interval_l += milli_query;
        sleep(1);
    }
    double_t avg_interval_l = t_interval_l/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.open(directory_queries + file_queries[6]);
    std::cout << directory_queries + file_queries[6] << std::endl;
    finQ >> type >> tstart >> tend >> id;
    ids.push_back(id);
    t_starts.push_back(tstart);
    t_ends.push_back(tend);
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> tend >> id >> type;
        ids.push_back(id);
        t_starts.push_back(tstart);
        t_ends.push_back(tend);
    }
    finQ.close();

    double_t t_mbr = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            util::geo::region mbr;
            rct::algorithm::MBR(t_starts[i], t_ends[i], ids[i], m_rct_index, mbr);
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_mbr += milli_query;
        sleep(1);
    }
    double_t avg_mbr = t_mbr/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.open(directory_queries + file_queries[7]);
    std::cout << directory_queries + file_queries[7] << std::endl;
    uint64_t k;
    finQ >> type >> tstart >> minX >> minY >> k;
    t_starts.push_back(tstart);
    points.push_back(util::geo::point{minX, minY});
    ks.push_back(k);
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> minX >> minY >> k >> type;
        t_starts.push_back(tstart);
        points.push_back(util::geo::point{minX, minY});
        ks.push_back(k);
    }
    finQ.close();
    double_t t_find_knn = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = util::time::user::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_knn;
            rct::algorithm::knn(ks[i], points[i], t_starts[i], m_rct_index, res_knn);
        }
        auto stop = util::time::user::now();
        auto milli = stop-start;
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (µs) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (µs) = " << milli_query << endl;
        t_find_knn += milli_query;
        sleep(1);
    }
    double_t avg_find_knn = t_find_knn/(double) TIMES;
    t_starts.clear();
    t_ends.clear();
    ks.clear();
    points.clear();

    finQ.open(directory_queries + file_queries[8]);
    std::cout << directory_queries + file_queries[8] << std::endl;
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
            rct::algorithm::knn_interval(ks[i], points[i], t_starts[i], t_ends[i], m_rct_index, res_knn);
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
    t_starts.clear();
    t_ends.clear();
    ks.clear();
    points.clear();
    ks.clear();

    finQ.open(directory_queries + file_queries[9]);
    std::cout << directory_queries + file_queries[9] << std::endl;
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
            rct::algorithm::knn_trajectory(ks[i], traj_x[i], traj_y[i],
                                               t_starts[i], t_ends[i], m_rct_index, res_knn);
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



    cout << "-----------------------------------------------------------------" << endl;
    cout << "Find_object (µs): " << avg_object <<  endl;
    cout << "Find_trajectory (µs): " << avg_traj <<  endl;
    cout << "Time Slice S (µs): " << avg_slice_s <<  endl;
    cout << "Time Slice L (µs): " << avg_slice_l <<  endl;
    cout << "Time Interval S (µs): " << avg_interval_s << endl;
    cout << "Time Interval L (µs): " << avg_interval_l <<  endl;
    cout << "MBR (µs): " << avg_mbr <<  endl;
    cout << "Knn (µs): " << avg_find_knn <<  endl;
    cout << "Knn int (µs): " << avg_find_knn_int <<  endl;
    cout << "Knn traj (µs): " << avg_find_knn_traj <<  endl;
    cout << "-----------------------------------------------------------------" << endl;

    std::string r = argv[2];
    std::ofstream out_res("RCT_multiple_sparse_rtree_" + r + "_" + util::file::remove_path(argv[1]) + ".res", std::ios_base::app);
    out_res <<  avg_object <<"," << avg_traj << "," << avg_slice_s << "," << avg_slice_l << ","
            << avg_interval_s << "," << avg_interval_l << "," << avg_mbr <<  "," << avg_find_knn <<  ","
            << avg_find_knn_int <<  "," << avg_find_knn_traj << endl;
    out_res.close();


    std::cout << "Everything is OK!" << std::endl;
    sdsl::write_structure<sdsl::HTML_FORMAT>(m_rct_index, util::file::index_file("rct_index_multiple_rtree", argv, 4) + ".html");
    sdsl::write_structure<sdsl::JSON_FORMAT>(m_rct_index, util::file::index_file("rct_index_multiple_rtree", argv, 4) + ".json");
    std::cout << "Done." << std::endl;

}
