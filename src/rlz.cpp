//
// Created by adrian on 14/11/18.
//

#include <vector>
#include <random>
#include <functional>
#include <rlz_naive.hpp>
#include <spiral_matrix_coder.hpp>

int main(int argc, const char* argv[]) {

    if(argc == 4){

        std::string input = argv[1];
        uint64_t size_ref = (uint64_t) atoi(argv[2]) * 1024 * 1024;
        uint64_t size_block = (uint64_t) atoi(argv[3]);
        auto size = util::file::file_size(input) / sizeof(uint64_t);
        std::vector<uint64_t > movements(size);
        std::ifstream in(input);
        in.read((char*) movements.data(), movements.size()*sizeof(uint64_t));
        in.close();

        rct::rlz_csa_bc_int64 m_rlz(movements, std::vector<uint64_t >(),  size_ref, size_block, 0);
        m_rlz.init_factorization(&movements);

        std::vector<typename rct::rlz_csa_bc_int64::factor_type> factors;
        while (m_rlz.has_next()) {
            auto f = m_rlz.next();
            //factors_log << "f.length: " << f.length << " f.offset: "  << f.offset << std::endl;
            factors.push_back(f);
        }
        std::vector<uint64_t> result;
        std::cout << "Decompressing: " << std::flush;
        m_rlz.decompress(factors, result);
        if(result.size() != movements.size()){
            std::cout << "Error size" << std::endl;
            exit(1);
        }
        for(uint64_t i = 0; i < result.size(); ++i){
            if(result[i] != movements[i]){
                std::cout << "Error movement: " << i << std::endl;
                exit(1);
            }
        }
        std::cout << "Phrases: " << factors.size() << std::endl;

    }


}