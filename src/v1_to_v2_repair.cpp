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

#include <rct_index_grammar.hpp>
#include <rct_index_grammar_rtree.hpp>
#include <string>

int main(int argc, const char* argv[]) {

    string dataset = argv[1];
    std::string index_file =  util::file::index_file("rct_index_repair", argv, argc)+ ".idx";
    std::cout << "Loading index: " << index_file << std::endl;
    rct::rct_index_grammar<2, rct::log_reference<>, rct::log_object_int_vector> m_rct_index;
    sdsl::load_from_file(m_rct_index, index_file);

    rct::rct_index_grammar_rtree<rct::log_reference<>, rct::log_object_int_vector> m_rct_index_rtree;
    m_rct_index_rtree.from_v1(m_rct_index, dataset);

    sdsl::util::clear(m_rct_index);

    std::cout << "Total objects: " << m_rct_index_rtree.total_objects << std::endl;
    std::string index_file_rtree =  util::file::index_file("rct_index_repair_rtree", argv, argc)+ ".idx";
    sdsl::store_to_file(m_rct_index_rtree, index_file_rtree);


}