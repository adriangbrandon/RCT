//
// Created by adrian on 15/11/18.
//

#include <time_util.hpp>
#include <string>
#include <fstream>
#include <vector>
#include <spiral_matrix_coder.hpp>
#include <rlz_naive.hpp>

int main(int argc, const char* argv[]) {


    if(argc == 4){

        using t_factor = rct::rlz_csa_bc_int::factor_type;
        uint64_t size_movements = 0, size_rmqs = 0, size_disappeared = 0,
                size_length_phrase = 0, size_offset_phrase = 0, size_value_phrase = 0, size_rmq_phrase = 0;
        std::string dataset_name = argv[1];
        uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024*1024;
        uint32_t size_block_bytes = (uint32_t) atoi(argv[3]);
        std::ifstream in(dataset_name);
        std::ofstream factors_log("factors.log");
        std::vector<uint32_t > input_reference;
        uint32_t id, old_id=-1, t, x, old_x,  y, old_y;
        std::cout << "Array of movements: " << std::flush;
        while(in){
            in >> id >> t >> x >> y;
            if(in.eof() || id > 99) break;
            if(id == old_id){
                int32_t diff_x = x - old_x;
                int32_t diff_y = y - old_y;
                input_reference.push_back(rct::spiral_matrix_coder::encode(diff_x, diff_y));
            }
            old_id = id;
            old_x = x;
            old_y = y;
        }
        in.close();
        std::cout << "Done." << std::endl;
        std::cout << "RLZ: " << std::flush;
        rct::rlz_csa_bc_int rlz(input_reference, std::vector<uint64_t>(), size_reference, size_block_bytes, 0);
        std::cout << "Done." << std::endl;

        //Computing size of storing the reference
        x = 0, old_x = 0, y = 0, old_y = 0;
        bool increase_x = true, increase_y = true;
        uint64_t ones_x = 1, ones_y = 1;
        for(uint64_t i = 0; i < rlz.reference.size(); ++i){
            auto pair = rct::spiral_matrix_coder::decode(rlz.reference[i]);
            //Movements
            size_movements += 4 + std::abs(pair.first) + std::abs(pair.second);
            //Adding two bits per change
            if(increase_x && pair.first < 0){
                size_rmqs += 2;
                increase_x = false;
                ++ones_x;
            }else if (!increase_x && pair.first > 0){
                size_rmqs += 2;
                increase_x = true;
                ++ones_x;
            }
            if(increase_y && pair.second < 0){
                size_rmqs += 2;
                increase_y = false;
                ++ones_y;
            }else if (!increase_y && pair.second > 0){
                size_rmqs += 2;
                increase_y = true;
                ++ones_y;
            }
        }
        //Adding the bitmaps which mark the changes
        size_rmqs += std::ceil(ones_x*(2 + log2(rlz.reference.size()/(double_t) ones_x)));
        size_rmqs += std::ceil(ones_y*(2 + log2(rlz.reference.size()/(double_t) ones_y)));

        uint64_t total_length = input_reference.size();
        old_id=-1;
        std::vector<uint32_t > trajectory;
        uint64_t total_factors = 0;
        in.open(dataset_name);
        auto start = util::time::user::now();
        while(in){
            in >> id >> t >> x >> y;
            //if(in.eof()) break;
            if(in.eof() || id > 99) break;
            if(id == old_id){
                int32_t diff_x = x - old_x;
                int32_t diff_y = y - old_y;
                trajectory.push_back(rct::spiral_matrix_coder::encode(diff_x, diff_y));
            }
            if(id > old_id && old_id != -1){
                std::cout << "Parsing: " << id << std::endl;
                std::vector<t_factor> factors;
                rlz.init_factorization(&trajectory);
                while(rlz.has_next()){
                    auto f = rlz.next();
                    //factors_log << "f.length: " << f.length << " f.offset: "  << f.offset << std::endl;
                    factors.push_back(f);
                }
                std::vector<uint32_t > result;
                std::cout << "Decompressing: " << std::flush;
                rlz.decompress(factors, result);
                std::cout << "Done." << std::endl;
                std::cout << "Checking results: " << std::flush;
                if(result.size() == trajectory.size()){
                    for(uint64_t i = 0; i < result.size(); ++i){
                        if(result[i] != trajectory[i]){
                            std::cout << "Error - at position: " << i << std::endl;
                            std::cout << "Expected: " << trajectory[i] << std::endl;
                            std::cout << "Obtained: " << result[i] << std::endl;
                            exit(1);
                        }
                    }
                }else{
                    std::cout << "Error - different size: " << result.size() << " vs. " << trajectory.size() << std::endl;
                    exit(1);
                }
                std::cout << "Done." << std::endl;

                //Computing size of storing each trajectory
                if(!trajectory.empty()){
                    size_disappeared += trajectory.size();
                    size_length_phrase += std::ceil(factors.size()*(2 + log2(trajectory.size()/(double_t) factors.size())));
                    size_offset_phrase += 32 * factors.size();
                    size_value_phrase += 64 * factors.size();
                    size_rmq_phrase += 2 * factors.size();
                }else{
                    size_value_phrase += 64;
                }
                size_value_phrase += 32; //time_start
                trajectory.clear();
                total_factors += factors.size();
                std::cout << std::endl;
            }
            old_id = id;
            old_x = x;
            old_y = y;
        }
        auto end = util::time::user::now();
        in.close();
        factors_log.close();
        std::cout << "Parsing: done in " << end-start << " (µs)" << std::endl;



        std::cout << "---------------------STATS------------------------" << std::endl;

        std::cout << "Total factors: " << total_factors << std::endl;
        std::cout << "Total length: " << total_length << std::endl;
        std::cout << "                     Size                         " << std::endl;
        std::cout << "- Reference: " << size_movements + size_rmqs << " bits." << std::endl;
        std::cout << "    - MOVE : " << size_movements << " bits." << std::endl;
        std::cout << "    - RMQS : " << size_rmqs << " bits." << std::endl;

        std::cout << "- Trajectories: " << size_disappeared + size_length_phrase + size_offset_phrase + size_value_phrase + size_rmq_phrase << " bits." << std::endl;
        std::cout << "    - DISAP : " << size_disappeared << " bits." << std::endl;
        std::cout << "    - LENGTH: " << size_length_phrase << " bits." << std::endl;
        std::cout << "    - OFFSET: " << size_offset_phrase << " bits." << std::endl;
        std::cout << "    - VALUES: " << size_value_phrase << " bits." << std::endl;
        std::cout << "    - RMQS: " << size_rmq_phrase << " bits." << std::endl;

        std::cout << "--------------------------------------------------" << std::endl;

        std::string stats_file = dataset_name + "_" + std::to_string(size_block_bytes) + "_" + std::to_string(size_reference) + ".txt";
        std::ofstream stats_output(stats_file);
        stats_output << "[" << size_movements << ", " << size_rmqs << "," << size_disappeared << ", " << size_length_phrase
        << ", " << size_offset_phrase << ", " << size_value_phrase << ", " << size_rmq_phrase << "]" << std::endl;
        stats_output.close();

    }

}