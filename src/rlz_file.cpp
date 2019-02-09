//
// Created by adrian on 14/11/18.
//

#include <vector>
#include <random>
#include <functional>
#include <rlz_naive.hpp>
#include <spiral_matrix_coder.hpp>

int main(int argc, const char* argv[]) {

    if(argc == 2){

        std::string dataset_name = argv[1];
        std::ifstream in(dataset_name);
        std::vector<uint64_t > movements;
        int max_diff_x = 0, max_diff_y = 0;
        int id, old_id=-1, t, x, old_x,  y, old_y;
        int first_x, first_y;
        std::cout << "Array of movements: " << std::flush;
        while(in){
            in >> id >> t >> x >> y;
            if(in.eof()) break;
            if(id == old_id){
                int diff_x = x - old_x;
                int diff_y = y - old_y;
                auto value = rct::spiral_matrix_coder::encode(diff_x, diff_y);
                movements.push_back(value);

            }else{
                first_x = x;
                first_y = y;
            }
            old_id = id;
            old_x = x;
            old_y = y;
        }
        in.close();


        std::ofstream out("movements.bin");
        out.write((char*) movements.data(), movements.size()*sizeof(uint64_t));
        out.close();

    }


}