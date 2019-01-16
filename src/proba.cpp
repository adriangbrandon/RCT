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
// Created by Adrián on 16/01/2019.
//

#include <unistd.h>
#include <vector>
#include <iostream>

using namespace std;
using namespace std::chrono;


class elemento {

private:
    std::vector<uint64_t> data;

public:
    elemento(){
        data = std::vector<uint64_t>(100000, 10);
    }

    //! Copy constructor
    elemento(const elemento& o)
    {
        data = o.data;
    }

    //! Move constructor
    elemento(elemento&& o)
    {
        *this = std::move(o);
    }


    elemento &operator=(const elemento &o) {
        if (this != &o) {
           data = o.data;
        }
        return *this;
    }

    elemento &operator=(elemento &&o) {
        if (this != &o) {
            data = std::move(o.data);
            //m_succs_disap = std::move(o.m_succs_disap);
        }
        return *this;
    }

    uint64_t through(const uint64_t max) const{
        uint64_t a = 0;
        for(uint64_t i = 0; i < max; ++i){
            a += data[i];
            std::cerr << a;
        }
        return a;
    }
    void toSleep10(){
        sleep(10);
    }
};

int main(int argc, const char* argv[]) {

    std::vector<elemento> list;
    list.resize(50, elemento());
    uint64_t max = atoi(argv[1]);
    //t
    auto start1 = high_resolution_clock::now();
    for(uint64_t i = 0; i < max; ++i){
        elemento b = list[10];
        b.through(10);
        b.through(20);
        b.through(30);
        b.through(40);
        b.through(50);
    }
    auto end1 = high_resolution_clock::now();
    //t

    //t
    auto start2 = high_resolution_clock::now();
    for(uint64_t i = 0; i < max; ++i) {
        list[10].through(10);
        list[10].through(20);
        list[10].through(30);
        list[10].through(40);
        list[10].through(50);
    }
    auto end2 = high_resolution_clock::now();
    //t

    //t
    auto start3 = high_resolution_clock::now();
    for(uint64_t i = 0; i < max; ++i) {
        auto it = list.begin() + 10;
        it->through(10);
        it->through(20);
        it->through(30);
        it->through(40);
        it->through(50);
    }
    auto end3 = high_resolution_clock::now();
    //t


    auto start4 = high_resolution_clock::now();
    for(uint64_t i = 0; i < max; ++i) {
        for (auto e : list) {
            e.through(100);
        }
    }
    auto end4 = high_resolution_clock::now();

    auto start5 = high_resolution_clock::now();
    for(uint64_t i = 0; i < max; ++i) {
        for (auto &e : list) {
            e.through(100);
        }
    }
    auto end5 = high_resolution_clock::now();

    auto start6 = high_resolution_clock::now();
    for(uint64_t i = 0; i < max; ++i) {
        for (const auto &e : list) {
            e.through(100);
        }
    }
    auto end6 = high_resolution_clock::now();

    auto start7 = high_resolution_clock::now();
    for(uint64_t i = 0; i < max; ++i) {
        for (auto it = list.begin(); it != list.end(); ++it) {
            it->through(100);
        }
    }
    auto end7 = high_resolution_clock::now();



    std::cout << "Duration 1: "<< duration_cast<milliseconds>(end1-start1).count() << std::endl;
    std::cout << "Duration 2: " << duration_cast<milliseconds>(end2-start2).count() << std::endl;
    std::cout << "Duration 3: " << duration_cast<milliseconds>(end3-start3).count() << std::endl;
    std::cout << "Duration 4: "<< duration_cast<milliseconds>(end4-start4).count() << std::endl;
    std::cout << "Duration 5: " << duration_cast<milliseconds>(end5-start5).count() << std::endl;
    std::cout << "Duration 6: " << duration_cast<milliseconds>(end6-start6).count() << std::endl;
    std::cout << "Duration 7: " << duration_cast<milliseconds>(end7-start7).count() << std::endl;
    /********************
     * Duration 1: 428
     * Duration 2: 367
     * Duration 3: 372
     * Duration 4: 15261
     * Duration 5: 11585
     * Duration 6: 11712
     * Duration 7: 11081
     ********************/

}