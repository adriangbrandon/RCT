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
// Created by Adrián on 13/02/2019.
//

#include <iostream>
#include <fstream>
#include <sdsl/suffix_trees.hpp>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc > 1) {
        string file  = string(argv[1]);
        string ofile = file + ".sa";
        if (argc > 2) {
            ofile = argv[2];
        }
        int_vector<> sa;
        {
            int_vector<64> temp;
            load_vector_from_file(temp, file, 8);
            cout << "input elements = " << temp.size() << endl;
            if (temp.size()==0 or temp[temp.size()-1] != 0) {
                cout << "Add 0 to input `" << file << "`" << endl;
                temp.resize(temp.size()+1);
                temp[temp.size()-1] = 0;
                store_to_plain_array<uint64_t>(temp, file);
            }
        }
        qsufsort::construct_sa(sa, file.c_str(), 8);
        cout << "done" << endl;
        cout << "sa.size()="<<sa.size() << endl;
        cout << "sa.width()=" <<(int)sa.width() << endl;
        cout << "bit_compress..." << endl;
        util::bit_compress(sa);
        cout << "sa.width()=" <<(int)sa.width() << endl;
        store_to_file(sa, ofile);
    } else {
        cout << "Usage: " << argv[0] << " file [ofile]" << endl;
        cout << " Computes the SA from an array of 64-bit integers." << endl;
        cout << " Result is stored in `ofile`, or `file`.sa if `ofile`" << endl;
        cout << " is not specified." << endl;
    }
}