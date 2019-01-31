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


    template<class t_value = int64_t>
    class reference_repair {

    public:
        typedef uint64_t size_type;
        typedef t_value value_type;
        typedef typename std::vector<value_type>::iterator iterator;
    private:
        std::vector<value_type> m_reference;
        size_type m_alpha;
        size_type m_one_value;
        std::vector<value_type> m_translate_table;

        void copy(const reference_repair &n){
            m_reference = n.m_reference;
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

        void decode_movements(int i, const repair &m_repair_algo) {
            //Skip values lower than 0
            if(i < m_one_value) return;
            while (i >= m_alpha) {
                decode_movements(m_repair_algo.rules[i-m_alpha].left, m_repair_algo);
                i =  m_repair_algo.rules[i-m_alpha].right;
            }
            //Leaves
            if(i < m_alpha){
                //Translate and decode the positions
                m_reference.emplace_back(m_translate_table[i-m_one_value]);
            }
        }


    public:

        reference_repair() = default;

        reference_repair(std::vector<value_type> &container){

            std::vector<value_type > terminals;
            std::map<value_type, u_char> map_terminals;
            split_repair_helper split_helper(&container, &terminals, &map_terminals);

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
            int total_MB = util::memory::total_memory_megabytes() * 0.8;
            m_repair_algo.run(log_repair, length_log, terminals, total_MB);

            //3. Decode log_repair
            m_alpha =  (uint32_t) m_repair_algo.alpha;
            m_translate_table = std::vector<value_type>(terminals.begin()+ m_one_value, terminals.end());
            std::vector<value_type> reference;
            std::cout << "Decoding movemnts" << std::endl;
            std::unordered_map<int, char> decoded;
            for(size_type i = 0; i < m_repair_algo.lenC; ++i){
                if(decoded.count(m_repair_algo.c[i]) == 0){
                    decode_movements(m_repair_algo.c[i], m_repair_algo);
                }
                decoded[m_repair_algo.c[i]] = 'a';
            }
            //4. Adding values not included in the reference
            std::cout << "Adding extra values to the reference" << std::endl;
            std::unordered_map<value_type, char> voc_reference;
            for(auto &value : m_reference){
                if(voc_reference.count(value) == 0) voc_reference[value] = 1;
            }
            for(auto &value : container){
                if(voc_reference.count(value) == 0){
                    voc_reference[value] = 1;
                    m_reference.push_back(value);
                }
            }
            std::cout << "Reference Done. Size: " << m_reference.size() << std::endl;
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
            }
            return *this;
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            sdsl::write_member(m_reference.size(), out, child, "size");
            written_bytes += sdsl::serialize_vector(m_reference, out, child, "reference");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in){
            size_type size;
            sdsl::read_member(size, in);
            m_reference.resize(size);
            sdsl::load_vector(m_reference, in);
        }

    };

}

#endif //RCT_REFERENCE_UNIFORM_SAMPLE_HPP
