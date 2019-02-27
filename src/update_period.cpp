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
// Created by Adrián on 27/02/2019.
//

#include <string>
#include <rct_index.hpp>

int main(int argc, const char* argv[]) {

    std::string dataset_path = argv[1];
    double_t ratio = (double_t) atoi(argv[2])/(double_t) 100;
    uint32_t period = (uint32_t) atoi(argv[3]);
    std::string index_file = util::file::index_file("rct_index_multiple", argv, argc) + ".idx";
    std::cout << "Loading index: " << index_file << std::endl;
    rct::rct_index<2, rct::log_reference<>, rct::log_object_int_vector, rct::rlz_multiple_csa_bc_int64> m_rct_index;
    sdsl::load_from_file(m_rct_index, index_file);

    std::vector<uint32_t> new_periods = {240, 360, 720};

    for(auto p : new_periods){
        std::string dataset_name = util::file::remove_extension(util::file::remove_path(argv[1]));
        std::string new_index_file = "rct_index_multiple_" dataset_name + "_" + std::to_string(atoi(argv[2])) + "_" + std::to_string(p) + ".idx";
        rct::rct_index<2, rct::log_reference<>, rct::log_object_int_vector, rct::rlz_multiple_csa_bc_int64> new_rct_index(m_rct_index);
        new_rct_index.update_period_snapshot(p, dataset_path);
        sdsl::store_to_file(new_rct_index, new_index_file);
    }





}