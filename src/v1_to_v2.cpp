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
// Created by Adrián on 17/01/2019.
//

#include <rct_index.hpp>
#include <rct_index_rtree.hpp>

int main(int argc, const char* argv[]) {

    string dataset = argv[1];
    uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024*1024;
    uint32_t size_block_bytes = (uint32_t) atoi(argv[3]);
    uint32_t period = (uint32_t) atoi(argv[4]);
    std::string index_file = "rct_index_" + std::to_string(size_reference) + "_" + std::to_string(size_block_bytes)
                             + "_" + std::to_string(period) + ".idx";
    std::cout << "Loading index" << std::endl;
    rct::rct_index<2, rct::log_reference<>, rct::log_object_int_vector> m_rct_index;
    sdsl::load_from_file(m_rct_index, index_file);

    rct::rct_index_rtree<rct::log_reference<>, rct::log_object_int_vector> m_rct_index_rtree;
    m_rct_index_rtree.from_v1(m_rct_index, dataset);

    sdsl::util::clear(m_rct_index);

    std::cout << "Total objects: " << m_rct_index_rtree.total_objects << std::endl;

    std::string index_file_rtree = "rct_index_rtree_" + std::to_string(size_reference) + "_" + std::to_string(size_block_bytes)
                             + "_" + std::to_string(period) + ".idx";

    sdsl::store_to_file(m_rct_index_rtree, index_file_rtree);


}