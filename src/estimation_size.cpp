//
// Created by adrian on 15/11/18.
//

#include <string>
#include <fstream>
#include <vector>
#include <spiral_matrix_coder.hpp>
#include <rlz_naive.hpp>

int main(int argc, const char* argv[]) {


    if(argc == 3){

        std::string dataset_name = argv[1];
        uint32_t size_reference = (uint32_t) atoi(argv[2]) * 1024*1024;
        std::ifstream in(dataset_name);
        std::ofstream factors_log("factors.log");
        std::vector<uint32_t > input_reference;
        uint32_t id, old_id=-1, t, x, old_x,  y, old_y;
        std::cout << "Array of movements: " << std::flush;
        while(in){
            in >> id >> t >> x >> y;
            if(in.eof() || id > 100) break;
            if(id == old_id){
                int32_t diff_x = x - old_x;
                int32_t diff_y = y - old_y;
                input_reference.push_back(rct::spiral_matrix_coder::encode(diff_x, diff_y));
            }
            old_id = id;
            old_x = x;
            old_y = y;
        }
        std::cout << "Done." << std::endl;
        std::cout << "RLZ: " << std::flush;
        rct::rlz_naive<> rlz(input_reference, size_reference);
        std::cout << "Done." << std::endl;

        std::cout << "Parsing: " << std::flush;
        using t_factor = rct::rlz_naive<>::factor_type;
        std::vector<t_factor> factors;
        while(rlz.has_next()){
            auto f = rlz.next(input_reference);
            factors_log << "f.length: " << f.length << " f.offset" << std::endl;
            factors.push_back(f);
        }
        std::cout << "Done." << std::endl;

        std::cout << "Total factors: " << factors.size() << std::endl;

        std::vector<uint32_t > result;
        std::cout << "Decompressing: " << std::flush;
        rlz.decompress(factors, result);
        std::cout << "Done." << std::endl;

        std::cout << "Checking results: " << std::flush;
        if(result.size() == input_reference.size()){
            for(uint64_t i = 0; i < result.size(); ++i){
                if(result[i] != input_reference[i]){
                    std::cout << "Error - at position: " << i << std::endl;
                    exit(1);
                }
            }
        }else{
            std::cout << "Error - different size: " << result.size() << " vs. " << input_reference.size() << std::endl;
            exit(1);
        }
        std::cout << "Done." << std::endl;

    }

}