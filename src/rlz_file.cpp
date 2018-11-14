//
// Created by adrian on 14/11/18.
//

#include <vector>
#include <random>
#include <functional>
#include <rlz_naive.hpp>

int main(int argc, const char* argv[]) {

    if(argc == 2){

        std::string input_file = argv[1];



    }

    std::mt19937 rng;
    std::vector<uint32_t> trajectory(10);
    {
        std::uniform_int_distribution<uint64_t> distribution(2,20);
        auto dice = std::bind(distribution, rng);
        for (size_t j=0; j<trajectory.size(); ++j) {
            trajectory[j]=dice();
        }
    }

    for(const auto &v : trajectory){
        if(v == 0) exit(10);
    }

    std::cout << "Trajectory Done." << std::endl;
    rct::rlz_naive<> rlz(trajectory, 10);
    std::cout << "RLZ Done." << std::endl;

    while(rlz.has_next()){
        auto factor = rlz.next();
        std::cout << "{" << factor.offset << ", " << factor.length << "}" << std::endl;
    }







}