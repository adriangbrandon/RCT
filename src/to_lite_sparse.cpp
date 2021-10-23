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
// Created by Adrián on 17/12/20.
//

#include <rct_index_rtree.hpp>


int main(int argc, const char* argv[]) {

    std::string dataset = argv[1];
    std::string index_file = util::file::index_file("rct_index_multiple_rtree", argv, argc) + ".idx";
    std::cout << index_file << std::endl;
    rct::rct_index_rtree<rct::log_reference<>, rct::log_object_int_vector, rct::rlz_multiple_csa_bc_int64> m_rct_index;
    std::ifstream in(index_file);
    m_rct_index.load(in, dataset);

    std::cout << "First index load" << std::endl;

    rct::rct_index_rtree<rct::log_reference_sparse<>, rct::log_object_lite_int_vector ,
                         rct::rlz_multiple_csa_bc_int64> m_rct_new_index;

    m_rct_new_index.to_lite(m_rct_index, dataset);
    std::string index_file_rtree =  util::file::index_file("rct_index_multiple_rtree_lite_sparse", argv, argc)+ ".idx";
    sdsl::store_to_file(m_rct_new_index, index_file_rtree);


}