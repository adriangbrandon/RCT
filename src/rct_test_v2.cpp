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

#include <string>
#include <rct_index_rtree.hpp>
#include <rct_algorithm_rtree.hpp>
#include <iostream>

int main(int argc, const char **argv) {

    std::string dataset = argv[1];
    uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024*1024;
    uint32_t size_block_bytes = (uint32_t) atoi(argv[3]);
    uint32_t period = (uint32_t) atoi(argv[4]);
    std::string index_file =  util::file::index_file("rct_index_rtree", argv, argc)+ ".idx";
    std::cout << "Loading index: " << index_file << std::endl;
    rct::rct_index_rtree<rct::log_reference<>, rct::log_object_int_vector> m_rct_index;

    std::ifstream in(index_file);
    m_rct_index.load(in, dataset);
    in.close();

   // Id: 162 TStart: 11572 TEnd: 13746
   /* std::vector<util::geo::traj_step> resultados_3;
    rct::algorithm::search_trajectory(162, 11572, 13746, m_rct_index, resultados_3);*/

    /* Consulta FAIL: tStart: 37982 tEnd:38782 minX:1383 maxX: 1703 minY:254629 maxY:254949
     * Valor esperado: id: 228
     * Valor obtido: id: 64
     *
     * */

    /*util::geo::region region{util::geo::point{1383, 254629}, util::geo::point{1703, 254949}};
    std::vector<uint32_t> resultados_2;
    rct::algorithm::time_interval(region, 37982, 38782, m_rct_index, resultados_2);
    for(const auto &r : resultados_2){
        std::cout << r << std::endl;
    }
    exit(10);*/

    struct ids_sort
    {
        inline bool operator() (const util::geo::id_point& p, const util::geo::id_point& q)
        {
            return (p.id < q.id);
        }
    };


    std::vector<std::string> queries_array;
    //queries_array.emplace_back("queries/traj.txt");
    queries_array.emplace_back("queries/ts_s.txt");
    queries_array.emplace_back("queries/ts_l.txt");
    queries_array.emplace_back("queries/ti_s.txt");
    queries_array.emplace_back("queries/ti_l.txt");

    /*Consulta MBR: oid: 313 tStart: 37566 tEnd:41408
It exists */

    //bool exist = m_index_ct.find_MBR(0, 2, 17, r);
    //bool exist = m_index_ct.find_MBR(0, 0, 17, r);
    //bool exist =  m_index_ct.find_MBR(0, 0, 1, r);
    //bool exist = m_index_ct.find_MBR(122, 17806, 17830, r);
    //bool exist =  m_index_ct.find_MBR(124, 28070, 28102, r);
    //bool exist =  m_index_ct.find_MBR(124, 28500, 28620, r);
    /*bool exist =  m_index_ct.find_MBR(313, 37566, 41408, r);
    std::cout << "exist: " << exist << std::endl;
    std::cout << "min_x: " << r.m_min_p.m_x << " min_y: " << r.m_min_p.m_y
    << " max_x: " << r.m_max_p.m_x << " max_y: " << r.m_max_p.m_y << std::endl;
    exit(20);*/

    int type, id;
    uint32_t t, minX, maxX, minY, maxY, x, y;
    // m_index1_ct.time_interval(279, 351, region(997, 276691, 1269, 277058));
    //exit(20);

    for(uint idx_q = 0; idx_q < queries_array.size(); idx_q++) {
        std::string query = queries_array[idx_q];
        std::cout << "Processing: " << query << std::endl;
        std::ifstream finQ(query);
        std::string results = query + "_r";
        std::ifstream finR(results);
        while (finQ) {
            finQ >> type;
            //Formato do arquivo imis1month
            if (!finQ.good()) continue; // skip newlines, etc.
            std::vector<util::geo::id_point> resultados;
            if (type == -1) {
                finQ >> t >> minX >> maxX >> minY >> maxY;
                std::cout << "T: " << t << " minX: " << minX << " maxX: " << maxX << " minY: " << minY << " maxY: " <<
                          maxY << std::endl;
                util::geo::region r{util::geo::point{minX, minY}, util::geo::point{maxX, maxY}};
                rct::algorithm::time_slice(r, t, m_rct_index, resultados);
                std::sort(resultados.begin(), resultados.end(), ids_sort());
                finR >> id;
                finR >> id;
                int index = 0;
                bool faltalog = false;
                while (id != -1) {
                    finR >> x >> y;
                    if (index < resultados.size() && (id != resultados[index].id
                        || x != resultados[index].x || y != resultados[index].y)) {
                        std::cout << "Consulta FAIL: t: " << t << " id:" << id << std::endl;
                        std::cout << "Valor esperado: id: " << id << "x: " << x << " y:" << y << std::endl;
                        std::cout << "Valor obtido: id: " << resultados[index].id << "x: " <<
                                  resultados[index].x << " y: " << resultados[index].y << std::endl;
                        exit(10);
                    }
                    finR >> id;
                    index++;
                }
                if (index == resultados.size()) {
                    std::cout << "Correcto" << std::endl;
                } else if (index > resultados.size()) {
                    std::cout << "Incorrecto. Faltan " << index - resultados.size() << " resultados" << std::endl;
                    exit(10);
                } else {
                    std::cout << "Incorrecto. Sobran " << resultados.size() - index << " resultados" << std::endl;
                    for(auto res : resultados){
                        std::cout << "id: " << res.id << " point: " << res.x << ", " << res.y << std::endl;
                    }
                    exit(10);
                }
                finQ >> type;
                resultados.clear();
                resultados.shrink_to_fit();
            } else if (type == -2) {
                uint tStart, tDiff, tResult;
                finQ >> tStart >> tDiff >> id;
                std::cout << "Id: " << id << " TStart: " << tStart << " TEnd: " << tStart+tDiff << std::endl;
                std::vector<util::geo::traj_step> resultados_3;
                rct::algorithm::search_trajectory_fast(id, tStart, tStart+ tDiff, m_rct_index, resultados_3);
                int index = 0;
                finR >> id;
                finR >> id;
                std::cout << "Resultados: " << resultados_3.size() << std::endl;
                while (id != -1) {
                    finR >> tResult >> x >> y;
                    if (index < resultados_3.size()
                        && (x != resultados_3[index].x || y != resultados_3[index].y
                            || tResult != resultados_3[index].t)) {
                        std::cout << "Index: " << index << std::endl;
                        std::cout << "Consulta FAIL: tStart: " << tStart  << "tEnd: " << tStart+tDiff << " id:" << id << std::endl;
                        std::cout << "Valor esperado: id: " << id << " t:" << tResult << " x: " << x << " y:" << y <<
                                  std::endl;
                        std::cout << "Valor obtido: id: " << id
                                  << " t: " << resultados_3[index].t
                                  << " x: " << resultados_3[index].x
                                  << " y: " << resultados_3[index].y << std::endl;
                        exit(10);
                    }
                    finR >> id;
                    index++;
                }
                std::cout << "Index: " << index << std::endl;
                if (index == resultados_3.size()) {
                    std::cout << "Correcto" << std::endl;
                } else if (index > resultados_3.size()) {
                    std::cout << "Incorrecto. Faltan " << index - resultados_3.size() << " resultados" << std::endl;
                    exit(10);
                } else {
                    std::cout << "Incorrecto. Sobran  resultados" << std::endl;
                    for (uint i = index; i < resultados_3.size(); i++) {
                        std::cout << "Resultado: " << resultados_3[i].t << " " << resultados_3[i].x << " " <<
                                  resultados_3[i].y << std::endl;
                    }
                    exit(10);
                }
                finQ >> type;
                resultados_3.clear();
                resultados_3.shrink_to_fit();
            } else if (type == -3) {
                uint tStart, tEnd, minX, maxX, minY, maxY;
                finQ >> tStart >> tEnd >> minX >> maxX >> minY >> maxY;
                util::geo::region region{util::geo::point{minX, minY}, util::geo::point{maxX, maxY}};
                // 1548 1820 281864 282231
                //2506 2686 1400 4123 279952 283629
                std::cout << "Consulta: tStart: " << tStart << " tEnd:" << tEnd << " minX:" << minX <<
                          " maxX: " << maxX << " minY:" << minY << " maxY:" << maxY << std::endl;

                std::vector<uint32_t> resultados_2;
                rct::algorithm::time_interval(region, tStart, tEnd, m_rct_index, resultados_2);
                std::sort(resultados_2.begin(), resultados_2.end());
                finR >> id;
                finR >> id;
                int index = 0;
                while (id != -1) {
                    if (index < resultados_2.size()
                        && id != resultados_2[index]) {
                        std::cout << "Consulta FAIL: tStart: " << tStart << " tEnd:" << tEnd << " minX:" << minX <<
                                  " maxX: " << maxX << " minY:" << minY << " maxY:" << maxY << std::endl;
                        std::cout << "Valor esperado: id: " << id << std::endl;
                        std::cout << "Valor obtido: id: " << resultados_2[index] << std::endl;
                        exit(10);
                    }
                    finR >> id;
                    index++;
                }
                if (index == resultados_2.size()) {
                    std::cout << "Correcto" << std::endl;
                } else if (index > resultados_2.size()) {
                    std::cout << "Incorrecto. Faltan " << index - resultados_2.size() << " resultados" << std::endl;
                    exit(10);
                }
                finQ >> type;
                resultados_2.clear();
                resultados_2.shrink_to_fit();
            }
        }
    }

    std::cout << "Everything is OK!" << std::endl;

//    dataset.clear();

}
