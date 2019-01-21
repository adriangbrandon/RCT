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


#include <vector>
#include <geo_util.hpp>
#include <rtree.hpp>

//
// Created by Adrián on 17/01/2019.
//
int main(int argc, char **argv) {

    std::vector <uint32_t> x_min = {0, 1, 3, 1, 2, 0, 0};
    std::vector <uint32_t> y_min = {1, 4, 5, 4, 2, 1, 1};
    std::vector <uint32_t> x_max = {6, 4, 8, 5, 5, 3, 6};
    std::vector <uint32_t> y_max = {4, 7, 6, 9, 5, 5, 3};
    std::vector <std::pair<uint32_t, util::geo::region>> data;

    for (size_t i = 0; i < x_min.size(); ++i) {
        util::geo::region r1{x_min[i], y_min[i], x_max[i], y_max[i]};
        data.push_back({i, r1});
    }

    /*size_t load_factor = 30;
    IStorageManager* memorymanager = StorageManager::createNewMemoryStorageManager();
    id_type indexIdentifier;
    ISpatialIndex* tree = RTree::createNewRTree(*memorymanager, 0.7, load_factor, load_factor, 2,
                                                SpatialIndex::RTree::RV_RSTAR, indexIdentifier);

    double plow[2], phigh[2];
    for(size_t i = 0; i < x_min.size(); ++i){
        plow[0]=x_min[i]; plow[1]=y_min[i];
        phigh[0]=x_max[i]; phigh[1]=y_max[i];
        region r1(x_min[i], y_min[i], x_max[i], y_max[i]);
        data.push_back({i, r1});
        Region r = Region(plow, phigh, 2);
        tree->insertData(0, 0, r, i);
    }
    size_t disc_factor = 1;
    size_t param = 1;
    int index_type = SRTREE_TYPE_REG;
    RTree::SRTree* tSRTree = new RTree::SRTree((RTree::RTree*)tree, disc_factor, index_type, param);

    delete tree;
    delete memorymanager;

    std::cout << "n: " << tSRTree->getNumObj() << std::endl;
    long* resultSRTree = new long[tSRTree->getNumObj()];
    int resultSRTreeSize = 0;
    tSRTree->range_query(2, 2, 2, 2, resultSRTree, &resultSRTreeSize);
    for(size_t i = 0; i < resultSRTreeSize; ++i){
        std::cout << resultSRTree[i] << ",";
    }
    std::cout << std::endl;
    delete tSRTree;*/

    rct::rtree ct_rtree(data);
    util::geo::region r_q{2, 2, 4, 4};
    std::vector <uint32_t> v;
    ct_rtree.intersection(r_q, v);
    for (uint64_t i = 0; i < v.size(); ++i) {
        std::cout << v[i] << ",";
    }
    std::cout << std::endl;
}
