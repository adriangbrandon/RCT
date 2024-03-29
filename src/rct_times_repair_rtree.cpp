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
#include <rct_index_grammar_rtree.hpp>
#include <rct_algorithm_rtree.hpp>

#define TIMES 1

using namespace std;
using namespace std::chrono;

int main(int argc, const char **argv) {


    size_t first_query_arg = 3;
    std::string dataset = argv[1];
    std::string index_file =  util::file::index_file("rct_index_repair_rtree", argv, first_query_arg)+ ".idx";
    rct::rct_index_grammar_rtree<rct::log_reference<>, rct::log_object_int_vector> m_rct_index;

    std::ifstream in(index_file);
    m_rct_index.load(in, dataset);
    in.close();
    std::cout << "Done" << std::endl;

    std::vector<uint32_t> ids;
    std::vector<uint32_t> t_starts;
    std::vector<uint32_t> t_ends;
    std::vector<util::geo::region> regions;
    std::vector<util::geo::point> points;
    int64_t type;
    uint32_t id, maxX, maxY, minX, minY, tstart, tend;

    std::ifstream finQ(argv[first_query_arg]);

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
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            rct::algorithm::search_object(ids[i], t_starts[i], m_rct_index, p);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_object += milli_query;
        sleep(10);
    }
    double_t avg_object = t_object/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.close();

    finQ.open(argv[first_query_arg+1]);
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
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<util::geo::traj_step> results;
            rct::algorithm::search_trajectory_fast(ids[i], t_starts[i], t_ends[i], m_rct_index, results);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_traj += milli_query;
        sleep(10);
    }
    double_t avg_traj = t_traj/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();


    finQ.open(argv[first_query_arg+2]);
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
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<util::geo::id_point> res_t_s;
            rct::algorithm::time_slice(regions[i], t_starts[i], m_rct_index, res_t_s);
            //log_ts << " " << t_starts[i] << " " << regions[i] << "[results.size = " << res_t_s.size() << "]" << std::endl;
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_slice_s += milli_query;
        sleep(10);
    }
    double_t avg_slice_s = t_slice_s/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.open(argv[first_query_arg+3]);
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
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<util::geo::id_point> res_t_s;
            rct::algorithm::time_slice(regions[i], t_starts[i], m_rct_index, res_t_s);
            //log_ts << " " << t_starts[i] << " " << regions[i] << "[results.size = " << res_t_s.size() << "]" << std::endl;
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_slice_l += milli_query;
        sleep(10);
    }
    double_t avg_slice_l = t_slice_l/(double) TIMES;
    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();
    finQ.close();
    //log_ts.close();
    finQ.open(argv[first_query_arg+4]);
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
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_t_i;
            rct::algorithm::time_interval(regions[i], t_starts[i], t_ends[i], m_rct_index, res_t_i);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_interval_s += milli_query;
        sleep(10);
    }
    double_t avg_interval_s = t_interval_s/(double) TIMES;

    /*double_t t_interval_brute_s = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_t_i;
            rct::algorithm::time_interval_brute_force(regions[i], t_starts[i], t_ends[i], m_rct_index, res_t_i);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_interval_brute_s += milli_query;
        sleep(10);
    }
    double_t avg_interval_brute_s = t_interval_brute_s/(double) TIMES;*/


    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.open(argv[first_query_arg+5]);
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
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_t_i;
            rct::algorithm::time_interval(regions[i], t_starts[i], t_ends[i], m_rct_index, res_t_i);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_interval_l += milli_query;
        sleep(10);
    }
    double_t avg_interval_l = t_interval_l/(double) TIMES;

    /*double_t t_interval_brute_l = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_t_i;
            rct::algorithm::time_interval(regions[i], t_starts[i], t_ends[i], m_rct_index, res_t_i);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_interval_brute_l += milli_query;
        sleep(10);
    }
    double_t avg_interval_brute_l = t_interval_brute_l/(double) TIMES;*/

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    finQ.open(argv[first_query_arg+6]);
    uint64_t k;
    std::vector<uint64_t> ks;
    finQ >> type >> tstart >> minX >> minY >> k;
    t_starts.push_back(tstart);
    points.push_back(util::geo::point{minX, minY});
    t_ends.push_back(k);
    finQ >> type;
    while (finQ) {
        finQ >> type >> tstart >> minX >> minY >> k >> type;
        t_starts.push_back(tstart);
        points.push_back(util::geo::point{minX, minY});
        t_ends.push_back(k);
    }
    finQ.close();
    double_t t_find_knn = 0;
    for(uint64_t j = 0; j < TIMES; j++){
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            std::vector<uint32_t > res_knn;
            rct::algorithm::knn(t_ends[i], points[i], t_starts[i], m_rct_index, res_knn);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_find_knn += milli_query;
        sleep(10);
    }
    double_t avg_find_knn = t_find_knn/(double) TIMES;
    t_starts.clear();
    t_ends.clear();
    points.clear();
    ks.clear();


    finQ.open(argv[first_query_arg+7]);
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
        auto start = high_resolution_clock::now();
        for(uint64_t i = 0; i < t_starts.size(); i++){
            util::geo::region mbr;
            rct::algorithm::MBR(t_starts[i], t_ends[i], ids[i], m_rct_index, mbr);
        }
        auto stop = high_resolution_clock::now();
        auto milli = duration_cast<milliseconds>(stop-start).count();
        double_t milli_query = milli/(double_t) t_starts.size();
        cout << "Time (ms) = " << milli << std::endl;
        cout << "Time per query (" << t_starts.size() << ")" << " (ms) = " << duration_cast<milliseconds>(stop-start).count()/((double)t_starts.size())<< endl;
        t_mbr += milli_query;
        sleep(10);
    }
    double_t avg_mbr = t_mbr/(double) TIMES;

    t_starts.clear();
    t_ends.clear();
    regions.clear();
    ids.clear();

    cout << "-----------------------------------------------------------------" << endl;
    cout << "Find_object (ms): " << avg_object <<  endl;
    cout << "Find_trajectory (ms): " << avg_traj <<  endl;
    cout << "Time Slice S (ms): " << avg_slice_s <<  endl;
    cout << "Time Slice L (ms): " << avg_slice_l <<  endl;
    cout << "Time Interval S (ms): " << avg_interval_s << endl;
    cout << "Time Interval L (ms): " << avg_interval_l <<  endl;
    cout << "Knn (ms): " << avg_find_knn <<  endl;
    cout << "MBR (ms): " << avg_mbr <<  endl;
    cout << "-----------------------------------------------------------------" << endl;

    struct decimal_comma : std::numpunct<char> {
        char do_decimal_point()   const { return ','; }  // separate with slash
    };
    std::cout.imbue(std::locale(std::cout.getloc(), new decimal_comma));
    cout << "TO EXCEL" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout <<  avg_object << endl;
    cout <<  avg_traj <<  endl;
    cout <<  avg_slice_s <<  endl;
    cout <<  avg_slice_l <<  endl;
    cout <<  avg_interval_s << endl;
    cout <<  avg_interval_l <<  endl;
    cout << avg_find_knn << endl;
    cout << avg_mbr << endl;
    cout << "-----------------------------------------------------------------" << endl;

    std::cout << "Everything is OK!" << std::endl;
    //    dataset.clear();

}
