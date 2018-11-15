//
// Created by adrian on 13/11/18.
//


#include <vector>
#include <random>
#include <functional>
#include <rlz_naive.hpp>
#include <file_util.hpp>

int main(int argc, const char* argv[]) {

    std::mt19937 rng;
    std::vector<uint32_t> trajectory(100000000);
    {
        std::uniform_int_distribution<uint64_t> distribution(1,10);
        auto dice = std::bind(distribution, rng);
        for (size_t j=0; j<trajectory.size(); ++j) {
            trajectory[j]=dice();
        }
    }

    for(const auto &v : trajectory){
        if(v == 0) exit(10);
    }
    /*
    for(uint64_t i = 0; i < trajectory.size(); ++i){
        std::cout << "traj[" << i << "]=" << trajectory[i] << std::endl;
    }*/
    std::cout << "Trajectory Done." << std::endl;
    rct::rlz_naive<> rlz(trajectory, 10000);
    std::cout << "RLZ Done." << std::endl;

    std::vector<rct::rlz_naive<>::factor_type> factors;
    //uint64_t j = 0;
    while(rlz.has_next()){
        auto factor = rlz.next(trajectory);
        factors.push_back(factor);
    }
    std::cout << "Parsing Done." << std::endl;



    std::vector<uint32_t > result;
    rlz.decompress(factors, result);
    std::cout << "Decompress Done." << std::endl;

    for(uint64_t i = 0; i < trajectory.size(); ++i){
       // std::cout << "Checking position: " << i << std::endl;
        if(trajectory[i] != result[i]){
            std::cout << "Error at position: " << i << std::endl;
            std::cout << "Expeced: " << trajectory[i] << std::endl;
            std::cout << "Obtained: " << result[i] << std::endl;
            exit(1);
        }
    }

    std::cout << "OK!" << std::endl;


    util::file::write_to_file("proba.bin", trajectory);
    std::vector<uint32_t>  read_trajectory;
    util::file::read_from_file("proba.bin", read_trajectory);

    std::cout << trajectory.size() << " - " << read_trajectory.size() << std::endl;
    for(uint64_t i = 0; i < trajectory.size(); ++i){
        if(trajectory[i] != read_trajectory[i]){
            std::cout << "Error" << std::endl;
        }
    }




}