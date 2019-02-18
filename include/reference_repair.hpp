//
// Created by adrian on 13/11/18.
//

#ifndef RCT_REFERENCE_REPAIR_HPP
#define RCT_REFERENCE_REPAIR_HPP

#include <vector>
#include <unordered_map>
#include <spiral_matrix_coder.hpp>
#include <split_repair_helper.hpp>
#include <sdsl/int_vector.hpp>
#include <repair.hpp>
#include <util_mem.hpp>
#include <geo_util.hpp>

namespace rct {


    template<class t_value = int64_t, class t_factor=  std::pair<uint64_t , uint64_t>>
    class reference_repair {

    public:
        typedef uint64_t size_type;
        typedef t_value value_type;
        typedef t_factor rlz_factor_type;
        typedef typename std::vector<value_type>::iterator iterator;
    private:
        std::vector<value_type> m_reference;
        std::unordered_map<value_type, rlz_factor_type> m_value_factor;
        std::vector<size_type> m_lengths;
        std::vector<util::geo::mbr> m_mbrs;

        std::vector<value_type> m_sequence;
        std::vector<value_type> m_translate_table;
        std::vector<std::pair<value_type, value_type>> m_rules;
        size_type m_alpha;
        size_type m_one_value;

        void copy(const reference_repair &n){
            m_reference = n.m_reference;
            m_value_factor = n.m_value_factor;
            m_lengths = n.m_lengths;
            m_mbrs = n.m_mbrs;
            m_alpha = n.m_alpha;
            m_one_value = n.m_one_value;
            m_translate_table = n.m_translate_table;
            m_sequence = n.m_sequence;
            m_rules = n.m_rules;
        }

        uint32_t _translate_log(const std::vector<value_type> &log, const std::vector<value_type> &terminals,
                                int* &log_repair){

            std::map<int64_t, int32_t> map_terminals;
            for(int32_t i = 0; i < terminals.size(); i++){
                map_terminals.insert(std::pair<value_type, int32_t>(terminals[i], i));
            }

            for(size_type i = 0; i < log.size(); i++){
                auto it_map_terminals = map_terminals.find(log[i]);
                if(it_map_terminals != map_terminals.end()){
                    log_repair[i] = it_map_terminals->second;
                }else{
                    std::cout << "error1" << std::endl;
                    exit(20);
                }
            }

            uint32_t one_value = 0;
            auto it_map_terminals = map_terminals.find(1);
            if(it_map_terminals != map_terminals.end()){
                one_value = (uint32_t) it_map_terminals->second;
            }else{
                std::cout << "error2" << std::endl;
                exit(20);
            }
            return one_value;
        }

        void decode_movements(value_type i, uint32_t &size, util::geo::movement &movement, util::geo::mbr &mbr) {
            //Skip values lower than 0
            if(i < m_one_value) return;
            while (i >= m_alpha) {
                decode_movements(m_rules[i-m_alpha].first, size, movement, mbr);
                i =  m_rules[i-m_alpha].second;
            }
            //Leaves
            if(i < m_alpha){
                //Translate and decode the positions
                auto value = m_translate_table[i-m_one_value];
                m_reference.emplace_back(value);
                auto pair = rct::spiral_matrix_coder::decode(value);
                movement.x += pair.first;
                movement.y += pair.second;
                if(mbr.min.x > movement.x){
                    mbr.min.x = movement.x;
                }
                if(mbr.min.y > movement.y){
                    mbr.min.y = movement.y;
                }
                if(mbr.max.x < movement.x){
                    mbr.max.x = movement.x;
                }
                if(mbr.max.y < movement.y){
                    mbr.max.y = movement.y;
                }
                ++size;
            }
        }

        void extend_factors(value_type i, uint32_t &start_pos, uint32_t &size) {
           // std::cout << "start_pos: " << start_pos << std::endl;
            //Skip values lower than 0
            if(i < m_one_value) return;
            if(i >= m_alpha) {
                uint32_t size_first = 0, pos_first = start_pos;
                extend_factors(m_rules[i-m_alpha].first, start_pos, size_first);
                if(m_rules[i-m_alpha].first >= m_alpha && m_value_factor.count(m_rules[i-m_alpha].first) == 0){
                    m_value_factor[m_rules[i-m_alpha].first] = rlz_factor_type{pos_first, size_first};
                }
                uint32_t size_second = 0, pos_second = start_pos;
                extend_factors(m_rules[i-m_alpha].second, start_pos, size_second);
                if(m_rules[i-m_alpha].second >= m_alpha && m_value_factor.count(m_rules[i-m_alpha].second) == 0){
                    m_value_factor[m_rules[i-m_alpha].second] = rlz_factor_type{pos_second, size_second};
                }
                size = size_first + size_second;
            }else{
                //Translate and decode the positions
                if(m_value_factor.count(i) == 0){
                    m_value_factor[i] = rlz_factor_type{start_pos, 1};
                }
                size=1;
                ++start_pos;
            }
        }

        void size_movements_rec(value_type i, size_type &size) {
            //Skip values lower than 0
            if(i < m_one_value) return;
            while (i >= m_alpha) {
                size_movements_rec(m_rules[i-m_alpha].first, size);
                i =  m_rules[i-m_alpha].second;
            }
            //Leaves
            if(i < m_alpha){
                //Translate and decode the positions
                ++size;
            }
        }


        void decode_all(){
            uint32_t size;
            uint32_t pos = 0;
            //Rules
            auto last_value = (int64_t) (m_rules.size() + m_alpha - 1);
            for(int64_t i = last_value; i >= 0; --i){
                if(m_value_factor.count(i) == 0){
                    size = 0;
                    util::geo::mbr mbr{0,0,0,0};
                    util::geo::movement movement{0,0};
                    decode_movements(i, size, movement, mbr);
                    m_value_factor[i] = rlz_factor_type{pos, size};
                    m_mbrs.push_back(mbr);
                    m_lengths.push_back(size);
                    size = 0;
                    extend_factors(i, pos, size);
                    //pos += size;
                }
            }

        }



        /*void decode_size(const uint64_t size){
            std::vector<value_type> sequence_sorted(m_sequence);
            std::sort(sequence_sorted.begin(), sequence_sorted.end(), std::greater<value_type >());
            uint32_t size_factor, pos = 0;
            for(uint32_t i = 0; i < sequence_sorted.size() && sizeof(value_type)*m_reference.size() < size; ++i){
                if(m_value_factor.count(sequence_sorted[i]) == 0){
                    size_factor = 0;
                    util::geo::mbr mbr{0,0,0,0};
                    util::geo::movement movement{0,0};
                    decode_movements(m_sequence[i], size_factor, movement, mbr);
                    m_value_factor[sequence_sorted[i]] = rlz_factor_type{pos, size_factor};
                    m_lengths.push_back(size_factor);
                    pos += size_factor;
                }
            }
        }*/

    public:


        const size_type &one_value = m_one_value;
        const std::vector<value_type> &sequence = m_sequence;
        const std::vector<size_type> &lengths = m_lengths;
        const std::vector<util::geo::mbr> &mbrs = m_mbrs;

        reference_repair() = default;

        reference_repair(const std::vector<value_type> &container, const std::vector<size_type> &lengths,
                         const size_type reference_size=0 , const size_type block_size = 0, const double_t ratio = 0){

            std::vector<value_type > terminals;
            std::map<value_type, u_char> map_terminals;
            split_repair_helper split_helper(&container, &terminals, &map_terminals);

            std::unordered_map<value_type, char> voc_map_container;
            std::vector<value_type> vocab;
            for(auto &value : container){
                if(value > 0 && voc_map_container.count(value) == 0){
                    voc_map_container[value] = 'a';
                    vocab.push_back(value);
                }
            }
            voc_map_container.clear();


            //1.Prepare log
            for(size_type i = 0; i < container.size(); ++i){
                if(container[i] > std::numeric_limits<int>::max() || container[i] < 0){
                  split_helper.add_split_at(i);
                }else{
                    if(!map_terminals.count(container[i])){
                        terminals.push_back(container[i]);
                        map_terminals.insert(std::pair<value_type , u_char>(container[i], 'a'));
                    }
                }
            }
            sort(terminals.begin(), terminals.end());
            auto log_repair = new int[container.size()];
            size_type length_log = container.size();
            m_one_value = _translate_log(container, terminals, log_repair);

            //2. Running repair
            repair m_repair_algo;
            auto total_MB = util::memory::total_memory_megabytes() * 0.8;
            m_repair_algo.run(log_repair, length_log, terminals, (int) total_MB);

            //3. Decode log_repair
            m_alpha =  (uint32_t) m_repair_algo.alpha;
            m_translate_table = std::vector<value_type>(terminals.begin()+ m_one_value, terminals.end());
            std::vector<value_type> reference;
            std::cout << "Decoding movements" << std::endl;

            m_sequence.resize(m_repair_algo.lenC);
            for(size_type i = 0; i < m_repair_algo.lenC; ++i){
                m_sequence[i] = m_repair_algo.c[i];
            }
            std::cout << "Sequence done" << std::endl;

            m_rules.resize(m_repair_algo.lenR);
            for(size_type i = 0; i < m_repair_algo.lenR; ++i){
                m_rules[i] = {m_repair_algo.rules[i].left, m_repair_algo.rules[i].right};
            };
            std::cout << "Rules done" << std::endl;
            m_repair_algo.clear();

            if(reference_size){
               // decode_size(reference_size);
            }else{
                decode_all();
            }
            //4. Adding values not included in the reference
            std::cout << "Adding extra values to the reference" << std::endl;

            for(auto &value : vocab){
                //std::cout << "value=" << value << std::endl;
                if(value > 0 && m_value_factor.count(value) == 0){
                    m_value_factor[value] = rlz_factor_type{(uint32_t) m_reference.size(),1};
                    m_reference.push_back(value);
                }
            }
            std::cout << "Reference Done. Size: " << m_reference.size() << std::endl;
        }

        size_type size_movements(value_type i){
            size_type size = 0;
            size_movements_rec(i, size);
            return size;
        }

        inline value_type operator[](const size_type& idx){
            return m_reference[idx];
        }

        inline value_type operator[](const size_type& idx) const{
            return m_reference[idx];
        }

        void* data(const size_type& idx){
            return &m_reference[idx];
        }

        rlz_factor_type get_factor(const size_type& v){
            return m_value_factor[v];
        }

        inline size_type size() const{
            return m_reference.size();
        }

        const iterator begin()
        {
            return m_reference.begin();
        }

        //! Iterator that points to the element after the last element of int_vector.
        /*! Time complexity guaranty is O(1).
         */
        const iterator end()
        {
            return m_reference.end();
        }

        //! Copy constructor
        reference_repair(const reference_repair& o)
        {
            copy(o);
        }

        //! Move constructor
        reference_repair(reference_repair&& o)
        {
            *this = std::move(o);
        }

        void swap(reference_repair& r)
        {
            if (this != &r) {
                std::swap(m_reference, r.m_reference);
                std::swap(m_mbrs, r.m_mbrs);
                std::swap(m_lengths, r.m_lengths);
                std::swap(m_value_factor, r.m_value_factor);
                std::swap(m_alpha, r.m_alpha);
                std::swap(m_one_value, r.m_one_value);
                std::swap(m_translate_table ,r.m_translate_table);
                std::swap(m_sequence , r.m_sequence);
                std::swap(m_rules, r.m_rules);
            }
        }

        reference_repair& operator=(const reference_repair& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        reference_repair& operator=(reference_repair&& v)
        {
            if (this != &v) {
                m_reference = std::move(v.m_reference);
                m_lengths = std::move(v.m_lengths);
                m_mbrs = std::move(v.m_mbrs);
                m_value_factor = std::move(v.m_value_factor);
                m_alpha = v.m_alpha;
                m_one_value = v.m_one_value;
                m_translate_table = std::move(v.m_translate_table);
                m_sequence = std::move(v.m_sequence);
                m_rules = std::move(v.m_rules);
            }
            return *this;
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            sdsl::write_member(m_reference.size(), out, child, "size_reference");
            written_bytes += sdsl::serialize_vector(m_reference, out, child, "reference");
            written_bytes += sdsl::write_member(m_alpha, out, child, "alpha");
            written_bytes += sdsl::write_member(m_one_value, out, child, "one_value");
            sdsl::write_member(m_translate_table.size(), out, child, "size_translate_table");
            written_bytes += sdsl::serialize_vector(m_translate_table, out, child, "translate_table");
            sdsl::write_member(m_sequence.size(), out, child, "size_sequence");
            written_bytes += sdsl::serialize_vector(m_sequence, out, child, "sequence");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in){
            size_type size_reference, size_translate_table, size_sequence;
            sdsl::read_member(size_reference, in);
            m_reference.resize(size_reference);
            sdsl::load_vector(m_reference, in);
            sdsl::read_member(m_alpha, in);
            sdsl::read_member(m_one_value, in);
            sdsl::read_member(size_translate_table, in);
            m_translate_table.resize(size_translate_table);
            sdsl::load_vector(m_translate_table, in);
            sdsl::read_member(size_sequence, in);
            m_sequence.resize(size_sequence);
            sdsl::load_vector(m_sequence, in);

        }

    };

}

#endif //RCT_REFERENCE_UNIFORM_SAMPLE_HPP
