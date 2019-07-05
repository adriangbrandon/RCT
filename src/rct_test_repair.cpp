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
#include <rct_index_grammar.hpp>
#include <rct_algorithm.hpp>
#include <iostream>

int main(int argc, const char **argv) {

    std::string dataset = argv[1];
    std::string path_queries =  argv[2];
    std::string index_file =  util::file::index_file("rct_index_repair", argv, 2)+ ".idx";
    std::cout << "Loading index: " << index_file << std::endl;
    rct::rct_index_grammar<2, rct::log_reference<>, rct::log_object_int_vector> m_rct_index;
    sdsl::load_from_file(m_rct_index, index_file);

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
    queries_array.emplace_back("queries/traj.txt");
    queries_array.emplace_back("queries/ts_s.txt");
    queries_array.emplace_back("queries/ts_l.txt");
    queries_array.emplace_back("queries/ti_s.txt");
    queries_array.emplace_back("queries/ti_l.txt");
    queries_array.emplace_back("queries/knn.txt");
    queries_array.emplace_back("queries/mbr.txt");

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

                std::vector<uint32_t> resultados_2, resultados_brute_force;

                rct::algorithm::time_interval(region, tStart, tEnd, m_rct_index, resultados_2);
                //rct::algorithm::time_interval_brute_force(region, tStart, tEnd, m_rct_index, resultados_brute_force);
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
                /*std::sort(resultados_brute_force.begin(), resultados_brute_force.end());
                if(resultados_2.size() != resultados_brute_force.size()){
                    std::cout << "Brute force different size" << std::endl;
                    exit(10);
                }
                for(uint64_t j = 0; j < resultados_2.size(); ++j){
                    if(resultados_2[j] != resultados_brute_force[j]){
                        std::cout << "Brute force different values" << std::endl;
                        exit(10);
                    }
                }*/
                finQ >> type;
                resultados_2.clear();
                resultados_2.shrink_to_fit();
                resultados_brute_force.clear();
                resultados_brute_force.shrink_to_fit();
            }else if (type == -6){
                uint32_t t, x, y, K, xr, yr;
                finQ >> t >> x >> y >> K;
                util::geo::point p_q{x, y};
                // 1548 1820 281864 282231
                //2506 2686 1400 4123 279952 283629

                std::cout << "Executando: x=" << x << " y=" << y << " t=" << t << " k=" << K << std::endl;
                std::vector<uint32_t> resultados_2;
                //TODO: cambiar esto
                rct::algorithm::knn(K, p_q, t, m_rct_index, resultados_2);
                finR >> id;
                finR >> id;
                int index = 0;
                std::vector<uint> valoresEsperados;
                std::vector<uint64_t> x_knn, y_knn;
                std::vector<double_t> distances;
                std::map<double_t, std::vector<uint64_t>> map_distances_objects;
                while (id != -1) {
                    finR >> xr >> yr;
                    valoresEsperados.push_back(id);
                    double_t distance = util::geo::distance(util::geo::point{xr, yr}, p_q);
                    std::map<double_t, std::vector<uint64_t >>::iterator it = map_distances_objects.find(distance);
                    if(it == map_distances_objects.end()){
                        distances.push_back(distance);
                        std::vector<uint64_t> vector;
                        vector.push_back(id);
                        map_distances_objects.insert(std::pair<double_t , std::vector<uint64_t>>(distance, vector));
                    }else{
                        (*it).second.push_back(id);
                    }
                    finR >> id;
                    index++;
                }
                uint i = 0;
                for(uint64_t j = 0; j < distances.size(); j++){
                    std::map<double_t , std::vector<uint64_t >>::iterator it = map_distances_objects.find(distances[j]);
                    std::cout << std::endl;
                    std::cout << "distance: " << distances[j] << std::endl;
                    std::vector<uint64_t> vector = (*it).second;
                    if(j == distances.size()-1){
                        while(i < resultados_2.size()) {
                            std::cout << "i: " << i << std::endl;
                            std::cout << "j: " << j << std::endl;
                            std::cout << "size: " << vector.size() << std::endl;
                            std::cout << "id: " << resultados_2[i] << std::endl;
                            std::vector<uint64_t>::iterator it2 = std::find(vector.begin(), vector.end(),
                                                                            resultados_2[i]);
                            if (it2 == vector.end()) {
                                std::cout << "Error" << std::endl;
                                for(uint64_t ii = 0; ii < vector.size(); ii++){
                                    std::cout << "vector[" << ii <<  "]=" <<vector[ii] << std::endl;
                                }

                                std::cout << "valores esperados: ";
                                for(uint64_t ii = 0; ii < valoresEsperados.size(); ii++){
                                    std::cout << valoresEsperados[ii] << ", ";
                                }
                                std::cout << std::endl;

                                std::cout << "valores obtidos: ";
                                for(uint64_t ii = 0; ii < resultados_2.size(); ii++){
                                    std::cout << resultados_2[ii] << ", ";
                                }
                                std::cout << std::endl;
                                exit(11);
                            }
                            i++;
                        }
                    }else{
                        while(!vector.empty()){
                            std::cout << "i: " << i << std::endl;
                            std::cout << "j: " << j << std::endl;
                            std::cout << "size: " << vector.size() << std::endl;
                            std::cout << "id: " << resultados_2[i] << std::endl;
                            std::vector<uint64_t >::iterator it2 = std::find(vector.begin(), vector.end(), resultados_2[i]);
                            if(it2 == vector.end()){
                                std::cout << "Error " << std::endl;
                                for(uint64_t ii = 0; ii < vector.size(); ii++){
                                    std::cout << "vector[" << ii <<  "]=" <<vector[ii] << std::endl;
                                }
                                exit(10);
                            }else{
                                vector.erase(it2);
                            }
                            i++;
                        }
                    }
                }


                if (K == resultados_2.size()) {
                    std::cout << "Correcto" << std::endl;
                } else {
                    std::cout << "Consulta FAIL: x=" << x << " y=" << y << " t=" << t << " k=" << K << std::endl;
                    std::cout << "Esperados: " << valoresEsperados.size() << " obtidos: " << resultados_2.size() << std::endl;
                    exit(10);
                }
                finQ >> type;
            }else if (type == -7) {
                uint t_start, t_end, q_id;
                uint min_x, min_y, max_x, max_y;
                int wildcard;
                finQ >> t_start >> t_end >> q_id;
                std::cout << "MBR: " << q_id << " TStart: " << t_start << " TEnd: " << t_end << std::endl;
                finR >> wildcard >> min_x >> min_y >> max_x >> max_y >> wildcard;
                util::geo::region mbr_result;
                rct::algorithm::MBR(t_start, t_end, q_id, m_rct_index, mbr_result);
                if(mbr_result.min.x == min_x && mbr_result.min.y == min_y
                   && mbr_result.max.x == max_x && mbr_result.max.y == max_y){
                    std::cout << "Correcto" << std::endl;
                }else {
                    std::cout << "Obtained: " << mbr_result << std::endl;
                    util::geo::region expected{min_x, min_y, max_x, max_y};
                    std::cout << "Expected:" << expected << std::endl;
                    exit(10);
                }
                finQ >> type;
            }else if (type == -4) {
                uint q_id, t_q;
                uint x, y;
                int wildcard;
                finQ >> q_id >> t_q;
                std::cout << "SO: " << q_id << " T: " << t_q << std::endl;
                finR >> wildcard >> wildcard >> x >> y >> wildcard;
                util::geo::point p_result;
                rct::algorithm::search_object(q_id, t_q, m_rct_index, p_result);
                if(p_result.x == x && p_result.y == y){
                    std::cout << "Correcto" << std::endl;
                }else {
                    std::cout << "Obtained: " << p_result << std::endl;
                    util::geo::point expected{x, y};
                    std::cout << "Expected:" << expected << std::endl;
                    exit(10);
                }
                finQ >> type;
            }
        }
    }

    std::cout << "Everything is OK!" << std::endl;

//    dataset.clear();

}
