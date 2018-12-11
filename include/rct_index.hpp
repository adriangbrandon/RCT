//
// Created by adrian on 15/11/18.
//

#ifndef RCT_RCT_INDEX_HPP
#define RCT_RCT_INDEX_HPP

#include <vector>
#include <snapshot.hpp>
#include <log_object.hpp>
#include <log_reference.hpp>
#include <geo_util.hpp>
#include <string>
#include "spiral_matrix_coder.hpp"
#include "time_util.hpp"

#define PRINT_STATS 1

namespace rct {

    template <uint64_t k = 2, class t_log_reference = log_reference<>, class t_log_object = log_object<>,
            class t_rlz = rlz_csa_sada_int >
    class rct_index {

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef t_log_reference log_reference_type;
        typedef t_log_object log_object_type;
        typedef t_rlz rlz_type;
        typedef typename rlz_type::factor_type factor_type;
        typedef snapshot<k> snapshot_type;

    private:
        size_type m_total_objects;
        size_type m_speed_max;
        size_type m_t_max;
        size_type m_x_max;
        size_type m_y_max;
        size_type m_level_max;
        size_type m_period_snapshot;
        log_reference_type m_log_reference;
        std::vector<log_object_type> m_log_objects;
        std::vector<snapshot_type> m_snapshots;

        void copy(const rct_index &o){
            m_total_objects = o.m_total_objects;
            m_speed_max = o.m_speed_max;
            m_t_max = o.m_t_max;
            m_x_max = o.m_x_max;
            m_y_max = o.m_y_max;
            m_level_max = o.m_level_max;
            m_period_snapshot = o.m_period_snapshot;
            m_log_reference = o.m_log_reference;
            m_log_objects = o.m_log_objects;
            m_snapshots = o.m_snapshots;
        }

        void _get_stats(const std::string &infile){
            m_total_objects = 0, m_speed_max = 0, m_t_max= 0, m_x_max = 0, m_y_max = 0;
            auto prev_id = (value_type) -1;
            int64_t prev_x, prev_y, prev_t, x, y, t, id, i = 0;
            std::ifstream in(infile);
            while(in){
                in >> id >> t >> x >> y;
                if(in.eof()) break;
                if(prev_id != id) ++m_total_objects;
                if(i > 0){
                    if(prev_id == id){
                        //std::cout << v.m_id << " " << prev_t << " " << v.m_t << std::endl;



                        uint64_t x_diff = std::abs((x - prev_x)) / (t - prev_t);
                        uint64_t y_diff = std::abs((y - prev_y)) / (t - prev_t);
                        if(x_diff > m_speed_max) {
                            m_speed_max = x_diff;
                        }
                        if(y_diff > m_speed_max) {
                            m_speed_max = y_diff;
                        }
                    }
                }
                if(x > m_x_max) m_x_max = x;
                if(y > m_y_max) m_y_max = y;
                if(t > m_t_max) m_t_max = t;
                prev_x = x;
                prev_y = y;
                prev_id = id;
                prev_t = t;
                ++i;
            }
            in.close();
            m_x_max++; m_y_max++;
            size_type max_xy = std::max(m_x_max, m_y_max);
            m_level_max = (size_type) std::ceil(log2(max_xy)/log2(k));
#if PRINT_STATS
            std::cout << "Total objects: " << m_total_objects << std::endl;
            std::cout << "Max speed: " << m_speed_max << std::endl;
            std::cout << "Max x: " << m_x_max << std::endl;
            std::cout << "Max y: " << m_y_max << std::endl;
            std::cout << "Max t: " << m_t_max << std::endl;
            std::cout << "Level max: " << m_level_max << std::endl;
#endif
        }


    public:

        const log_reference_type &log_reference = m_log_reference;
        const std::vector<log_object_type> &log_objects = m_log_objects;
        const std::vector<snapshot_type> &snapshots = m_snapshots;

        rct_index() = default;

        rct_index(const std::string &dataset_file, const size_type size_reference, const size_type size_block,
                  const size_type period_snapshot) {

            m_period_snapshot = period_snapshot;
            _get_stats(dataset_file);

            std::ifstream in(dataset_file);
            std::ofstream factors_log("factors.log");
            std::vector<uint32_t> input_reference;
            uint32_t id, old_id = (uint32_t) -1, t, x, old_x = 0, y, old_y = 0;
            size_type n_snapshots = util::math::ceil_div(m_t_max, m_period_snapshot);
            std::vector<k2_tree_representation_lite<k> > trees(n_snapshots, k2_tree_representation_lite<k>(m_level_max));
            std::cout << "Array of movements: " << std::flush;
            while (in) {
                in >> id >> t >> x >> y;
                if (in.eof()) break;
                if (id == old_id) {
                    int32_t diff_x = x - old_x;
                    int32_t diff_y = y - old_y;
                    input_reference.push_back(spiral_matrix_coder::encode(diff_x, diff_y));
                }
                if(t % m_period_snapshot == 0){
                    trees[t/period_snapshot].insertObject(x, y, id);
                }
                old_id = id;
                old_x = x;
                old_y = y;
            }
            in.close();
            std::cout << "Compressing snapshots. " << std::endl;
            m_snapshots = std::vector<snapshot<k>>(n_snapshots);
            for(size_type i = 0; i < n_snapshots; i++) {
                m_snapshots[i] = snapshot<k>(trees[i], m_total_objects);
            }
            trees.clear();
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
                if (in.eof()) id = (uint32_t) -1;
                //if (in.eof() || id > 99) id = (uint32_t) -1;
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
                m_total_objects = o.m_total_objects;
                m_speed_max = o.m_speed_max;
                m_t_max = o.m_t_max;
                m_x_max = o.m_x_max;
                m_y_max = o.m_y_max;
                m_level_max = o.m_level_max;
                m_period_snapshot = o.m_period_snapshot;
                m_log_reference = std::move(o.m_log_reference);
                m_log_objects = std::move(o.m_log_objects);
                m_snapshots = std::move(o.m_snapshots);
            }
            return *this;
        }

        void swap(rct_index &o) {
            std::swap(m_total_objects, o.m_total_objects);
            std::swap(m_speed_max, o.m_speed_max);
            std::swap(m_t_max, o.m_t_max);
            std::swap(m_x_max, o.m_x_max);
            std::swap(m_y_max, o.m_y_max);
            std::swap(m_level_max, o.m_level_max);
            std::swap(m_period_snapshot, o.m_period_snapshot);
            std::swap(m_log_objects, o.m_log_objects);
            m_log_reference.swap(o.m_log_reference);
            std::swap(m_snapshots, o.m_snapshots);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            written_bytes += sdsl::write_member(m_total_objects, out, child, "total_objects");
            written_bytes += sdsl::write_member(m_speed_max, out, child, "speed_max");
            written_bytes += sdsl::write_member(m_t_max, out, child, "t_max");
            written_bytes += sdsl::write_member(m_x_max, out, child, "x_max");
            written_bytes += sdsl::write_member(m_y_max, out, child, "y_max");
            written_bytes += sdsl::write_member(m_level_max, out, child, "level_max");
            written_bytes += sdsl::write_member(m_period_snapshot, out, child, "period_snapshot");
            //TODO: delete next line (log_size)
            sdsl::write_member(m_log_objects.size(), out, child, "log_size");
            written_bytes += sdsl::serialize_vector(m_log_objects, out, child, "log_objects");
            written_bytes += m_log_reference.serialize(out, child, "log_reference");
            written_bytes += sdsl::serialize_vector(m_snapshots, out, child, "snapshots");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            sdsl::read_member(m_total_objects, in);
            sdsl::read_member(m_speed_max, in);
            sdsl::read_member(m_t_max, in);
            sdsl::read_member(m_x_max, in);
            sdsl::read_member(m_y_max, in);
            sdsl::read_member(m_level_max, in);
            sdsl::read_member(m_period_snapshot, in);
            //TODO: delete these two lines and run m_log_objects.resize(m_total_objects);
            size_type log_size;
            sdsl::read_member(log_size, in);
            m_log_objects.resize(log_size);
            sdsl::load_vector(m_log_objects, in);
            m_log_reference.load(in);
            m_snapshots.resize(util::math::ceil_div(m_t_max, m_period_snapshot));
            sdsl::load_vector(m_snapshots, in);
        }


    };

}

#endif //RCT_RCT_INDEX_HPP
