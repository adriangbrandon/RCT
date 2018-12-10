//
// Created by adrian on 15/11/18.
//

#ifndef RCT_RCT_INDEX_HPP
#define RCT_RCT_INDEX_HPP

#include <vector>
#include <log_object.hpp>
#include <log_reference.hpp>
#include <geo_util.hpp>
#include "spiral_matrix_coder.hpp"
#include "time_util.hpp"

namespace rct {

    template <class t_log_reference = log_reference<>, class t_log_object = log_object<>, class t_rlz = rlz_csa_sada_int >
    class rct_index {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef t_log_reference log_reference_type;
        typedef t_log_object log_object_type;
        typedef t_rlz rlz_type;
        typedef typename rlz_type::factor_type factor_type;

    private:
        log_reference_type m_log_reference;
        std::vector<log_object_type> m_log_objects;

        void copy(const rct_index &o){
            m_log_reference = o.m_log_reference;
            m_log_objects = o.m_log_objects;
        }


    public:

        const log_reference_type &log_reference = m_log_reference;
        const std::vector<log_object_type> &log_objects = m_log_objects;

        rct_index() = default;

        rct_index(const std::string &dataset_file, const size_type size_reference, const size_type size_block) {

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
            rlz_type rlz(input_reference, size_reference, size_block);
            std::cout << "Done." << std::endl;
            input_reference.clear();
            input_reference.shrink_to_fit();

            std::vector<util::geo::movement> ref_movements;
            for(size_type i = 0; i < rlz.reference.size(); ++i){
                auto pair = rct::spiral_matrix_coder::decode(rlz.reference[i]);
                ref_movements.emplace_back(util::geo::movement{pair.first, pair.second});
            }
            m_log_reference = log_reference_type (ref_movements);
            ref_movements.clear();

            id = 0, old_id = (uint32_t) -1, old_x = 0, old_y = 0;
            in.open(dataset_file);
            std::vector<uint32_t> movements;
            std::vector<util::geo::traj_step> trajectory;
            auto start = util::time::user::now();
            while (!in.eof() && id != -1) {
                in >> id >> t >> x >> y;
                //if (in.eof()) id = (uint32_t) -1;
                if (in.eof() || id > 1) id = (uint32_t) -1;
                if (id == old_id) {
                    int32_t diff_x = x - old_x;
                    int32_t diff_y = y - old_y;
                    movements.push_back(rct::spiral_matrix_coder::encode(diff_x, diff_y));
                }
                if (id != old_id && old_id != -1) {
                    std::cout << "Parsing: " << old_id << std::endl;
                    std::vector<factor_type> factors;
                    rlz.init_factorization(&movements);
                    while (rlz.has_next()) {
                        auto f = rlz.next();
                        //factors_log << "f.length: " << f.length << " f.offset: "  << f.offset << std::endl;
                        factors.push_back(f);
                    }

                    m_log_objects.emplace_back(log_object_type(trajectory, factors));
                    movements.clear();
                    trajectory.clear();
                }
                old_id = id;
                old_x = x;
                old_y = y;
                trajectory.emplace_back(util::geo::traj_step{t, x, y});
            }
            auto end = util::time::user::now();
            std::cout << "Parsing in: " << end - start << " µs" << std::endl;
            std::cout << "Everything Done." << std::endl;

            //m_log_objects[0].print();
        }

        //! Copy constructor
        rct_index(const rct_index& o)
        {
            copy(o);
        }

        //! Move constructor
        rct_index(rct_index&& o)
        {
            *this = std::move(o);
        }


        rct_index &operator=(const rct_index &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        rct_index &operator=(rct_index &&o) {
            if (this != &o) {
                m_log_objects = o.m_log_objects;
                m_log_reference = o.m_log_reference;
            }
            return *this;
        }

        void swap(rct_index &o) {
            std::swap(m_log_objects, o.m_log_objects);
            m_log_reference.swap(o.m_log_reference);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            sdsl::write_member(m_log_objects.size(), out, child, "log_size");
            written_bytes += sdsl::serialize_vector(m_log_objects, out, child, "log_objects");
            written_bytes += m_log_reference.serialize(out, child, "log_reference");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            size_type log_size;
            sdsl::read_member(log_size, in);
            m_log_objects.resize(log_size);
            sdsl::load_vector(m_log_objects, in);
            m_log_reference.load(in);
        }


    };

}

#endif //RCT_RCT_INDEX_HPP
