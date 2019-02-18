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
// Created by Adrián on 12/12/2018.
//
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <spiral_matrix_coder.hpp>
#include <rlz_naive.hpp>
#include <geo_util.hpp>
#include <cmath>

int main(int argc, const char* argv[]) {

    if (argc == 2) {

        std::string dataset_name = argv[1];

        std::ifstream in(dataset_name);
        rct::rlz_csa_bc_int m_rlz;
        const auto size_ref = 10*1024*1024;
        const auto size_block = 1024;
        std::string ref = dataset_name + "_" + std::to_string(size_ref) + "_" + std::to_string(size_block) + ".ref";

        if(!util::file::file_exists(ref)){
            std::vector<uint32_t > input_reference;
            std::unordered_map<int, int> terminals_map;
            int max_diff_x = 0, max_diff_y = 0;
            int id, old_id = -1, t, x, old_x, y, old_y;
            int first_x, first_y;
            std::cout << "Array of movements: " << std::flush;
            while (in) {
                in >> id >> t >> x >> y;
                if (in.eof()) break;
                if (id == old_id) {
                    int diff_x = x - old_x;
                    int diff_y = y - old_y;
                    if (std::abs(x - first_x) > max_diff_x) {
                        max_diff_x = std::abs(x - first_x);
                    }
                    if (std::abs(y - first_y) > max_diff_y) {
                        max_diff_y = std::abs(y - first_y);
                    }
                    auto value = rct::spiral_matrix_coder::encode(diff_x, diff_y);
                    if (terminals_map.count(value) == 0) terminals_map[value] = 1;
                    input_reference.push_back(value);
                } else {
                    first_x = x;
                    first_y = y;
                }
                old_id = id;
                old_x = x;
                old_y = y;
            }
            rct::rlz_csa_bc_int rlz_aux(input_reference, std::vector<uint64_t >(),  size_ref, size_block, 0);
            m_rlz = rlz_aux;
            sdsl::store_to_file(m_rlz, ref);

        }else{
            sdsl::load_from_file(m_rlz, ref);
        }
        in.close();
        std::cout << "Done." << std::endl;

        in.open(dataset_name);
        int id = 0, old_id = -1, t, x, old_x, y, old_y;
        std::vector<uint32_t> movements;
        while (!in.eof() && id != -1) {
            in >> id >> t >> x >> y;
            if (in.eof()) id = (uint32_t) -1;
            //if (in.eof() || id > 99) id = (uint32_t) -1;
            if (id == old_id) {
                int32_t diff_x = x - old_x;
                int32_t diff_y = y - old_y;
                movements.push_back(rct::spiral_matrix_coder::encode(diff_x, diff_y));
            }
            if (id != old_id && old_id != -1) {
                std::cout << "Parsing: " << old_id << std::endl;
                std::vector<typename rct::rlz_csa_bc_int::factor_type> factors;
                m_rlz.init_factorization(&movements);
                while (m_rlz.has_next()) {
                    auto f = m_rlz.next();
                    //factors_log << "f.length: " << f.length << " f.offset: "  << f.offset << std::endl;
                    factors.push_back(f);
                }
                std::cout << "Factors: " << factors.size() << std::endl;
                std::cout << "Movements: " << movements.size() << std::endl;
                std::vector<uint32_t> result;
                std::cout << "Decompressing: " << std::flush;
                m_rlz.decompress(factors, result);
                if(result.size() != movements.size()){
                    std::cout << "Error size" << std::endl;
                    exit(0);
                }
                for(uint64_t i = 0; i < movements.size(); ++i){
                    if(movements[i] != result[i]){
                        std::cout << "Error at position: " << i << std::endl;
                        exit(0);
                    }
                }
                std::cout << "Ok!" << std::endl << std::endl;
                movements.clear();
                /*for(const auto &factor : factors){
                    std::cout << "offset: " << factor.offset << ", length: " << factor.length << std::endl;
                }
                std::cout << "Reference" << std::endl;
                m_rlz.print_reference(0, 256);
                std::cout << "Movements " << std::endl;
                for(const auto &move : movements){
                    std::cout << move << ", ";
                }
                std::cout << std::endl;*/

                //if(old_id == 0) exit(9);
            }
            old_id = id;
            old_x = x;
            old_y = y;
        }
    }

}