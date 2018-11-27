//
// Created by adrian on 15/11/18.
//

#ifndef RCT_RCT_INDEX_HPP
#define RCT_RCT_INDEX_HPP

#include <vector>
#include <log_object.hpp>
#include <geo_util.hpp>
#include "spiral_matrix_coder.hpp"

namespace rct {

    class rct_index {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;

    private:
        std::vector<log_object<>> m_log_objects;

    public:

        rct_index() = default;

        rct_index(const std::string &dataset_file, const size_type size_reference, const size_type size_block) {

            using t_factor = rct::rlz_naive<>::factor_type;
            std::ifstream in(dataset_file);
            std::ofstream factors_log("factors.log");
            std::vector<uint32_t> input_reference;
            uint32_t id, old_id = (uint32_t) -1, t, x, old_x = 0, y, old_y = 0;
            std::cout << "Array of movements: " << std::flush;
            while (in) {
                in >> id >> t >> x >> y;
                if (in.eof()) break;
                if (id == old_id) {
                    int32_t diff_x = x - old_x;
                    int32_t diff_y = y - old_y;
                    input_reference.push_back(spiral_matrix_coder::encode(diff_x, diff_y));
                }
                old_id = id;
                old_x = x;
                old_y = y;
            }
            in.close();
            std::cout << "Done." << std::endl;
            std::cout << "RLZ: " << std::flush;
            rct::rlz_naive<> rlz(input_reference, size_reference, size_block);
            std::cout << "Done." << std::endl;
            input_reference.clear();
            input_reference.shrink_to_fit();


            id = 0, old_id = (uint32_t) -1, old_x = 0, old_y = 0;
            in.open(dataset_file);
            std::vector<uint32_t> movements;
            std::vector<util::geo::traj_step> trajectory;
            while (!in.eof()) {
                in >> id >> t >> x >> y;
                //if (in.eof()) id = (uint32_t) -1;
                if(in.eof()) id = (uint32_t) -1;
                if (id == old_id) {
                    int32_t diff_x = x - old_x;
                    int32_t diff_y = y - old_y;
                    movements.push_back(rct::spiral_matrix_coder::encode(diff_x, diff_y));
                }
                if (id != old_id && old_id != -1) {
                    std::cout << "Parsing: " << old_id << std::endl;
                    std::vector<t_factor> factors;
                    rlz.init_factorization(&movements);
                    while (rlz.has_next()) {
                        auto f = rlz.next();
                        //factors_log << "f.length: " << f.length << " f.offset: "  << f.offset << std::endl;
                        factors.push_back(f);
                    }

                    m_log_objects.emplace_back(log_object<>(trajectory, factors));
                    movements.clear();
                    trajectory.clear();
                }
                old_id = id;
                old_x = x;
                old_y = y;
                trajectory.emplace_back(util::geo::traj_step{t, x, y});
            }
            std::cout << "Everything Done." << std::endl;
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += sdsl::serialize_vector(m_log_objects, out, child, "log_objects");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }


    };

}

#endif //RCT_RCT_INDEX_HPP
