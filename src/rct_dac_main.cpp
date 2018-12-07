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


int main(int argc, const char* argv[]) {

    if(argc == 4){
        uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024*1024;
        uint32_t size_block_bytes = (uint32_t) atoi(argv[3]);
        rct::rct_index<rct::log_reference<>, rct::log_object_int_vector> m_rct_index(argv[1], size_reference, size_block_bytes);
        std::ofstream out("rct_index_" + std::to_string(size_reference) + "_" + std::to_string(size_block_bytes) + ".html");
        sdsl::write_structure<sdsl::HTML_FORMAT>(m_rct_index, out);

        std::ifstream in(argv[1]);
        uint32_t id, t, x, y;
        util::geo::point r;
        while(in){
            in >> id >> t >> x >> y;
            if(in.eof() || id > 99) break;
            rct::algorithm::search_object(id, t, m_rct_index, r);
            std::cout << "Obtained: " << r.x << ", " << r.y << std::endl;
            if(r.x != x || r.y != y){
                std::cout << "Error looking for: id=" << id << " t=" << t << std::endl;
                std::cout << "Expected: " << x << ", " << y << std::endl;

                exit(0);
            }
        }
        out.close();
    }

}

