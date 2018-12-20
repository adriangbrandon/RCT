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
// Created by Adrián on 18/12/2018.
//

#include <cstdint>
#include <cstdlib>
#include <runs_bitvector.hpp>

int main(int argc, const char* argv[]) {

    //sdsl::bit_vector bv = {0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0};
    sdsl::bit_vector bv = {1,1,1,1,0,0,0,0,
                           0,0,0,0,0,0,0,0,
                           0,0,0,1,1,1,1,0,
                           0,0,0,0,0,0,1,1,
                           1,1,0};
    //sdsl::bit_vector bv = {0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0};
    rct::runs_bitvector<> m_runs(bv);
    m_runs.print();
    for(uint64_t i = 0; i < bv.size(); ++i){
        if(m_runs[i] != bv[i]){
            std::cout << "Error at " << i << std::endl;
            exit(1);
        }
    }
    rct::succ_support_runs_bitvector_v2<1> m_succ;
    rct::succ_support_runs_bitvector_v2<0> m_succ_0;
    rct::rank_support_runs_bitvector<1> m_rank;
    sdsl::util::init_support(m_succ, &m_runs);
    sdsl::util::init_support(m_succ_0, &m_runs);
    sdsl::util::init_support(m_rank, &m_runs);
    for(uint64_t i = 0; i < m_runs.size(); ++i){
        std::cout << "i: " << i << " succ: " << m_succ(i) << std::endl;
    }
    m_succ_0(30);
    for(uint64_t i = 0; i < m_runs.size(); ++i){
        std::cout << "i: " << i << " succ0: " << m_succ_0(i) << std::endl;
    }
    /*for(uint64_t i = 0; i <= m_runs.size(); ++i){
        std::cout << "i: " <<i << " rank: " << m_rank(i) << std::endl;
    }*/
    /*uint64_t ones = 0, beg = 0, length = 0;
    m_rank(0, ones, beg, length);
    for(uint64_t i = 0; i < m_runs.size(); ++i){
        std::cout << "i: " <<i << " succ2: " << m_succ(i, ones, beg, length) << std::endl;
    }*/
   // m_succ_0(30);
    std::cout << "Everything OK!" << std::endl;
}