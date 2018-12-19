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
// Created by Adrián on 13/12/2018.
//

#include <cstdint>
#include <cstdlib>
#include <rct_index.hpp>
#include <rct_algorithm.hpp>
#include <runs_bitvector.hpp>

int main(int argc, const char* argv[]) {

    if (argc == 4) {
        uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024 * 1024;
        uint32_t size_block_bytes = (uint32_t) atoi(argv[3]);
        std::string index_file =
                "rct_index_" + std::to_string(size_reference) + "_" + std::to_string(size_block_bytes) + ".idx";
        rct::rct_index<2, rct::log_reference<>, rct::log_object_int_vector> m_rct_index;
        sdsl::load_from_file(m_rct_index, index_file);

        uint64_t byte1 = 0, byte2 = 0, byte3 = 0, byte4 = 0;
        double diff = 0;
        for( const auto &log : m_rct_index.log_objects){
            //log.stats(byte1, byte2, byte3, byte4, log.min_x_values);
            diff += log.diff_dac(log.min_x_values);
        }
        //std::cout << "Min X-> byte1: " << byte1 << " byte2: " << byte2 << " byte3: " << byte3 << " byte4: " << byte4 << std::endl;
        std::cout << "Min X: \t" << diff << std::endl;
        byte1 = 0, byte2 = 0, byte3 = 0, byte4 = 0;
        diff = 0;
        for( const auto &log : m_rct_index.log_objects){
            //log.stats(byte1, byte2, byte3, byte4, log.min_y_values);
            diff += log.diff_dac(log.min_y_values);
        }
        //std::cout << "Min Y-> byte1: " << byte1 << " byte2: " << byte2 << " byte3: " << byte3 << " byte4: " << byte4 << std::endl;
        std::cout << "Min Y: \t" << diff << std::endl;
        byte1 = 0, byte2 = 0, byte3 = 0, byte4 = 0;
        diff = 0;
        for( const auto &log : m_rct_index.log_objects){
            //log.stats(byte1, byte2, byte3, byte4, log.max_x_values);
            diff += log.diff_dac(log.max_x_values);
        }
       // std::cout << "Max X-> byte1: " << byte1 << " byte2: " << byte2 << " byte3: " << byte3 << " byte4: " << byte4 << std::endl;
        std::cout << "Max X: \t" << diff << std::endl;
        byte1 = 0, byte2 = 0, byte3 = 0, byte4 = 0;
        diff = 0;
        for( const auto &log : m_rct_index.log_objects){
            //log.stats(byte1, byte2, byte3, byte4, log.max_y_values);
            diff += log.diff_dac(log.max_y_values);
        }
        //std::cout << "Max Y-> byte1: " << byte1 << " byte2: " << byte2 << " byte3: " << byte3 << " byte4: " << byte4 << std::endl;
        std::cout << "Max Y: \t" << diff << std::endl;

        diff = 0;
        for( const auto &log : m_rct_index.log_objects){
            diff += log.diff_dac(log.x_values);
        }
        std::cout << "X: \t" << diff << std::endl;

        diff = 0;
        for( const auto &log : m_rct_index.log_objects){
            diff += log.diff_dac(log.y_values);
        }
        std::cout << "Y: \t" << diff << std::endl;

        diff = 0;
        for( const auto &log : m_rct_index.log_objects){
            diff += log.diff_dac(log.offsets);
        }
        std::cout << "Offsets: \t" << diff << std::endl;

        /*for(const auto &log : m_rct_index.log_objects){
            log.runs(log.disap);
        }*/


        /*diff = 0;
        for(const auto &log : m_rct_index.log_objects){
            diff += log.diff_oz_vector(log.info);
        }

        std::cout << "Disap OZ: " << diff << std::endl;

        diff = 0;
        for(const auto &log : m_rct_index.log_objects){
            diff += log.diff_runs_vector(log.info);
        }

        std::cout << "Disap Runs: " << diff << std::endl;

        diff = 0;
        for(const auto &log : m_rct_index.log_objects){
            diff += log.diff_sd_vector(log.info);
        }

        std::cout << "SD Runs: " << diff << std::endl;

        diff = 0;
        for(const auto &log : m_rct_index.log_objects){
            diff += log.diff_runs_bitvector(log.info);
        }
        std::cout << "Disap Runs Bitvector: " << diff << std::endl;

        diff = 0;
        for(const auto &log : m_rct_index.log_objects){
            diff += log.size_runs_bitvector(log.info);
        }
        std::cout << "Size Runs Bitvector: " << diff << std::endl;*/

       /* diff = 0;
        for(const auto &log : m_rct_index.log_objects){
            diff += log.diff_offset_vector(log.offsets);
        }
        std::cout << "Offsets Avg: " << diff << std::endl;

        for(uint64_t oid = 0; oid < m_rct_index.log_objects.size(); ++oid){
            rct::runs_bitvector<> m_runs(m_rct_index.log_objects[oid].disap);
            for(uint64_t i = 0; i < m_rct_index.log_objects[oid].disap.size(); ++i){
                if(m_runs[i] != m_rct_index.log_objects[oid].disap[i]){
                    std::cout << "Error at " << oid << ", " << i << std::endl;
                    exit(1);
                }
            }
        }
        std::cout << "Everythin ok!" << std::endl;*/
        //m_rct_index.log_objects[2].disap;

       /* sdsl::bit_vector bv = {1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1};

        rct::runs_bitvector<> m_runs(m_rct_index.log_objects[79].disap);
        rct::oz_vector<> m_oz(m_rct_index.log_objects[79].disap);
        //rct::runs_bitvector<> m_runs(bv);
        //rct::oz_vector<> m_oz(bv);
        std::cout << sdsl::size_in_bytes(m_runs) << std::endl;
        std::cout << sdsl::size_in_bytes(m_oz) << std::endl;*/

        diff = 0;
        for(const auto &log : m_rct_index.log_objects){
            diff += log.size_element(log.info);
        }
        std::cout << "Size Runs Bitvector: " << diff << std::endl;
    }

}
