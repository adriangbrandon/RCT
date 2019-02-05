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
// Created by Adrián on 04/02/2019.
//

#include <stdint.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <spiral_matrix_coder.hpp>
#include <map>
#include <util_mem.hpp>
#include <repair.hpp>
#include <split_repair_helper.hpp>
#include <rct_index_grammar.hpp>
#include <rct_algorithm.hpp>

using namespace rct;

uint32_t _translate_log(const std::vector<int64_t > &log, const std::vector<int64_t> &terminals,
                        int* &log_repair){

    std::map<int64_t, int32_t> map_terminals;
    for(int32_t i = 0; i < terminals.size(); i++){
        map_terminals.insert(std::pair<int64_t, int32_t>(terminals[i], i));
    }

    for(uint64_t i = 0; i < log.size(); i++){
        auto it_map_terminals = map_terminals.find(log[i]);
        if(it_map_terminals != map_terminals.end()){
            log_repair[i] = it_map_terminals->second;
        }else{
            std::cout << "error1" << std::endl;
            exit(20);
        }
    }

    uint32_t one_value = 0;
    auto it_map_terminals = map_terminals.find(1);
    if(it_map_terminals != map_terminals.end()){
        one_value = (uint32_t) it_map_terminals->second;
    }else{
        std::cout << "error2" << std::endl;
        exit(20);
    }
    return one_value;
}

void size_movement(int i, const repair &m_repair_algo, const uint64_t m_one_value, const uint64_t m_alpha, uint64_t &size) {
    //Skip values lower than 0
    if(i < m_one_value) return;
    while (i >= m_alpha) {
        size_movement(m_repair_algo.rules[i-m_alpha].left, m_repair_algo, m_one_value, m_alpha, size);
        i =  m_repair_algo.rules[i-m_alpha].right;
    }
    //Leaves
    if(i < m_alpha){
        //Translate and decode the positions
        ++size;
    }
}

int main(int argc, const char* argv[]) {


    std::ifstream in(argv[1]);
    std::vector<int64_t > container;

    uint32_t id, old_id = (uint32_t) -1, t, old_t = 0, x, old_x = 0, y, old_y = 0;
    //m_reap = std::vector<sdsl::bit_vector>(n_snapshots, sdsl::bit_vector(m_total_objects, 0));
    //m_disap = std::vector<sdsl::bit_vector>(n_snapshots, sdsl::bit_vector(m_total_objects, 0));
    std::cout << "Array of movements: " << std::flush;
    while (in) {
        in >> id >> t >> x >> y;
        if (in.eof()) break;
        if (id == old_id) {
            int32_t diff_x = x - old_x;
            int32_t diff_y = y - old_y;
            container.push_back(spiral_matrix_coder::encode(diff_x, diff_y));
        }else if(old_id != -1){
            container.push_back(-1);
        }

        old_id = id;
        old_x = x;
        old_y = y;
        old_t = t;
    }
    in.close();

    std::vector<int64_t > terminals;
    std::map<int64_t, u_char> map_terminals;
    split_repair_helper split_helper(&container, &terminals, &map_terminals);

    //1.Prepare log
    for(uint64_t i = 0; i < container.size(); ++i){
        if(container[i] > std::numeric_limits<int>::max() || container[i] < 0){
            split_helper.add_split_at(i);
        }else{
            if(!map_terminals.count(container[i])){
                terminals.push_back(container[i]);
                map_terminals.insert(std::pair<uint64_t , u_char>(container[i], 'a'));
            }
        }
    }
    sort(terminals.begin(), terminals.end());
    auto log_repair = new int[container.size()];
    uint64_t length_log = container.size();
    auto m_one_value = _translate_log(container, terminals, log_repair);

    //2. Running repair
    repair m_repair_algo;
    auto total_MB = util::memory::total_memory_megabytes() * 0.8;
    m_repair_algo.run(log_repair, length_log, terminals, (int) total_MB);

    //3. Decode log_repair
    auto m_alpha =  (uint32_t) m_repair_algo.alpha;
    auto m_translate_table = std::vector<uint64_t >(terminals.begin()+ m_one_value, terminals.end());

    std::string index_file = "rct_index_repair_120.idx";
    rct::rct_index_grammar<2, rct::log_reference<>, rct::log_object_int_vector> m_rct_index;
    std::ifstream index_stream(index_file);
    m_rct_index.load(index_stream);

    uint64_t phrases = 0;
    for(uint64_t i = 0; i < m_repair_algo.lenC; ++i){
        if(m_repair_algo.c[i] >= m_one_value){
            ++phrases;
        }
    }
    std::cout << "Phrases: " << phrases << std::endl;

    /*for(uint64_t i = 0; i < m_repair_algo.lenC; ++i){
        uint64_t size = 0;
        if(m_repair_algo.c[i] >= m_one_value){
            size_movement(m_repair_algo.c[i], m_repair_algo, m_one_value, m_alpha, size);
            std::cerr << size << std::endl;
        }else{
            std::cerr << "1" << std::endl;
        }

    }*/


}