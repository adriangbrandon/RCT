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
        std::vector<int > input_reference;
        std::unordered_map<int, int> terminals_map;
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
                if(std::abs(x - first_x) > max_diff_x){
                    max_diff_x = std::abs(x - first_x);
                }
                if(std::abs(y - first_y) > max_diff_y){
                    max_diff_y = std::abs(y - first_y);
                }
                int value = (int) rct::spiral_matrix_coder::encode(diff_x, diff_y);
                if(terminals_map.count(value) == 0) terminals_map[value] = 1;
                input_reference.push_back(value);
            }else{
                first_x = x;
                first_y = y;
            }
            old_id = id;
            old_x = x;
            old_y = y;
        }
        in.close();

        std::vector<int> terminals;
        for(auto it = terminals_map.begin(); it != terminals_map.end(); ++it){
            terminals.push_back(it->first);
        }
        std::sort(terminals.begin(), terminals.end());
        terminals_map.clear();


        for(uint32_t i = 0; i < terminals.size(); ++i){
            terminals_map[terminals[i]] = i;
        }


        int* array = new int[input_reference.size()];
        for(uint32_t i = 0; i < input_reference.size(); ++i){
            array[i] = terminals_map[input_reference[i]];
        }

        std::ofstream out("movements.bin");
        out.write((char*) array, input_reference.size()*sizeof(int));
        out.close();

        std::cout << "Max diff x: " << max_diff_x << std::endl;
        std::cout << "Max diff y: " << max_diff_y << std::endl;
    }


}