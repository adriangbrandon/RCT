//
// Created by adrian on 13/11/18.
//

#include <rlz_int_vector.hpp>
#include <sdsl/int_vector.hpp>

int main(int argc, const char* argv[]) {


    sdsl::int_vector<> reference = {8,8,8,8,8,8,1, 2, 3, 4, 5, 6, 7};

    rct::rlz_int_vector<> m_rlz(reference);

    std::vector<uint32_t > container = {8,8,8,1,8,8,1,1,1,2,3,4,5,6,7,8,8,8,8,8,8,8};
    m_rlz.init_factorization(&container);
    std::vector<typename rct::rlz_int_vector<>::factor_type> factors;
    while (m_rlz.has_next()) {
        auto f = m_rlz.next();
        std::cout << "f.length: " << f.length << " f.offset: "  << f.offset << std::endl;
        factors.push_back(f);
    }
    std::cout << "Factors: " << factors.size() << std::endl;
    std::cout << "Movements: " << container.size() << std::endl;
    std::vector<uint32_t> result;
    std::cout << "Decompressing: " << std::flush;
    m_rlz.decompress(factors, result);
    if(result.size() != container.size()){
        std::cout << "Error size" << std::endl;
        exit(0);
    }
    for(uint64_t i = 0; i < container.size(); ++i){
        if(container[i] != result[i]){
            std::cout << "Error at position: " << i << std::endl;
            exit(0);
        }
    }


}