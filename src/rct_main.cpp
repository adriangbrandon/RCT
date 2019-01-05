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


#include <cstdint>
#include <cstdlib>
#include <rct_index.hpp>
#include <rct_algorithm.hpp>
#include <vector>
#include <spiral_matrix_coder.hpp>

int main(int argc, const char* argv[]) {

    /*auto number = rct::spiral_matrix_coder::encode(404, -33497);
    auto pair = rct::spiral_matrix_coder::decode(number);
    std::cout << number << std::endl;
    std::cout << pair.first << ", " << pair.second << std::endl;
    exit(9);*/

    if(argc == 4){
        uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024*1024;
        uint32_t size_block_bytes = (uint32_t) atoi(argv[3]);
        std::string index_file = "rct_index_" + std::to_string(size_reference) + "_" + std::to_string(size_block_bytes) + ".idx";

        if(!util::file::file_exists(index_file)){
            std::cout << "Building index" << std::endl;
            auto t1 = util::time::user::now();
            rct::rct_index<2, rct::log_reference<>, rct::log_object_int_vector> m_rct_index(argv[1], size_reference, size_block_bytes, 120);
            auto t2 = util::time::user::now();
            std::cout << "User time: " << t2 - t1 << " µs" << std::endl;
            std::ofstream out("rct_index_" + std::to_string(size_reference) + "_" + std::to_string(size_block_bytes) + ".html");
            sdsl::write_structure<sdsl::HTML_FORMAT>(m_rct_index, out);
            sdsl::store_to_file(m_rct_index, index_file);
        }
        std::cout << "Loading index" << std::endl;
        rct::rct_index<2, rct::log_reference<>, rct::log_object_int_vector> m_rct_index;
        sdsl::load_from_file(m_rct_index, index_file);


        std::ifstream in(argv[1]);
        uint32_t id, t, x, y;
        util::geo::point r;
        /*
         * Error looking for: id=134 t=16678
         * Expected: 404, 248994
         * Obtained: 6950, 276802
         */



       /* rct::algorithm::search_object(134, 10972, m_rct_index, r);
        rct::algorithm::search_object(134, 16678, m_rct_index, r);
        while(in){
            in >> id >> t >> x >> y;
            if(in.eof()) break;
            rct::algorithm::search_object(id, t, m_rct_index, r);
            std::cout << "Obtained: " << r.x << ", " << r.y << std::endl;
            if(r.x != x || r.y != y){
                std::cout << "Error looking for: id=" << id << " t=" << t << std::endl;
                std::cout << "Expected: " << x << ", " << y << std::endl;
                std::cout << "Obtained: " << r.x << ", " << r.y << std::endl;
                exit(0);
            }
        }
        in.close();


        std::cout << "Search trajectory id=0 t_i=0 t_j=200: " << std::endl;
        std::vector<util::geo::traj_step> traj;
        rct::algorithm::search_trajectory(0, 0, 200, m_rct_index, traj);

        util::geo::region region{util::geo::point{890, 273100}, util::geo::point{900,273200}};
        //util::geo::region region{util::geo::point{10, 100}, util::geo::point{20, 200}};
        //rct::algorithm::time_interval()

        for(const auto &step : traj){
            std::cout << "t:" << step.t << " x:" << step.x << " y:" << step.y << std::endl;
        }

        std::vector<uint32_t >  results;
        rct::algorithm::time_interval(region, 100, 1000, m_rct_index, results);
        for(const auto &a : results){
            std::cout << "id: " << a << std::endl;
        }*/
    }

}

