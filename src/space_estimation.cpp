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
// Created by Adrián on 1/2/21.
//


#include <rct_index.hpp>
#include <rmMq.hpp>

int main(int argc, const char* argv[]) {


    std::string index_file = argv[1];
    std::string dataset = argv[2];
    std::cout << "Loading index: " << index_file << std::endl;
    rct::rct_index<2, rct::log_reference_sparse<>, rct::log_object_int_vector> m_rct_index;
    sdsl::load_from_file(m_rct_index, index_file);
    auto c_size = sdsl::size_in_bytes(m_rct_index);
    sdsl::nullstream ns;
    auto new_size = m_rct_index.new_space(ns);

    std::ifstream in(dataset);
    int id, t, x, y;
    int prev_id = 0, prev_x, prev_t, prev_y;
    int i = 0;
    std::vector<uint64_t > x_values, y_values;
    while(in){
        in >> id >> t >> x >> y;
        if(in.eof()) break;
        if(i > 0){
            if(prev_id < id){
                //std::cout << v.m_id << " " << prev_t << " " << v.m_t << std::endl;
                rct::rmMq<> rmMq_x(&x_values);
                rct::rmMq<> rmMq_y(&y_values);
                new_size += sdsl::size_in_bytes(rmMq_x);
                new_size += sdsl::size_in_bytes(rmMq_y);
                x_values.clear();
                y_values.clear();
            }
        }
        x_values.push_back(x);
        y_values.push_back(y);
        prev_x = x;
        prev_y = y;
        prev_id = id;
        prev_t = t;
        ++i;
    }


    std::cout << "Current size: " << c_size << " new size: "
    << new_size << "(" << (new_size/(double) c_size)*100  << "%)" <<std::endl;


}