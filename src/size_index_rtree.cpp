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
// Created by Adrián on 07/03/2019.
//

#include <string>
#include <rct_index_rtree.hpp>
#include <rct_algorithm_rtree.hpp>

#define TIMES 1

using namespace std;
using namespace std::chrono;

int main(int argc, const char **argv) {


    size_t first_query_arg = 5;
    std::string dataset = argv[1];
    uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024*1024;
    uint32_t size_block_bytes = (uint32_t) atoi(argv[3]);
    uint32_t period = (uint32_t) atoi(argv[4]);
    std::string index_file =  util::file::index_file("rct_index_rtree", argv, first_query_arg)+ ".idx";
    std::cout << "Loading index: " << index_file << std::endl;
    rct::rct_index_rtree<rct::log_reference<>, rct::log_object_int_vector> m_rct_index;

    std::ifstream in(index_file);
    m_rct_index.load(in, dataset);
    in.close();
    std::cout << "Size: " << sdsl::size_in_bytes(m_rct_index) << std::endl;

}