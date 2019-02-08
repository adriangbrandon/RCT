//
// Created by adrian on 13/11/18.
//

#ifndef RCT_RLZ_REPAIR_V2_HPP
#define RCT_RLZ_REPAIR_V2_HPP

#include <file_util.hpp>
#include <unordered_map>
#include <reference_repair.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#define INITIAL_CAPACITY 1024

namespace rct {



    typedef struct rlz_factor {
        typedef uint32_t offset_type;
        typedef uint32_t length_type;
        offset_type offset;
        length_type length;
    } rlz_factor_type;


    template <class t_value = uint32_t, class t_reference = reference_repair<t_value, rlz_factor>,
              class t_csa = sdsl::csa_wt_int<> >
    class rlz_repair_v2 {

    public:

        typedef t_reference reference_type;
        typedef rlz_factor_type rlz_factor_type;
        typedef t_value value_type;
        typedef t_csa csa_type;
        typedef typename csa_type::size_type size_type;

    private:

        reference_type m_reference;
        size_type m_repair_pos = 0;
        size_type m_input_size = 0;
        size_type m_input_pos = 0;
        const std::vector<value_type>* m_input;

        void copy(const rlz_repair_v2 &n){
            m_reference = n.m_reference;
            m_repair_pos = n.m_repair_pos;
            m_input_size = n.m_input_size;
            m_input_pos = n.m_input_pos;
        }


    public:
        const reference_type &reference = m_reference;

        rlz_repair_v2() = default;

        rlz_repair_v2(std::vector<value_type> &container, const size_type reference_size, const size_type block_size = 0){

            m_reference = reference_type(container, reference_size);
            std::cout << "size reference: " << m_reference.size() << std::endl;
        }


        void init_factorization(const std::vector<value_type> *container){
            m_input = container;
            m_input_size = container->size();
            m_input_pos = 0;
        }


        rlz_factor_type next(){

            rlz_factor_type factor;
            auto repair_value = m_reference.sequence[m_repair_pos];
            if(repair_value < m_reference.one_value){
                //Value of input, it is a terminal
                factor =  m_reference.get_factor(m_input->at(m_input_pos));
                ++m_input_pos;
            }else{
                factor = m_reference.get_factor(repair_value);
                m_input_pos += factor.length;
            }
            ++m_repair_pos;
            if(m_input_pos == m_input_size) ++m_repair_pos;
            return factor;
        }


        inline bool has_next() {
            if(m_input_size == 0) ++m_repair_pos;
            return m_input_pos < m_input_size;
        }

        //! Copy constructor
        rlz_repair_v2(const rlz_repair_v2& o)
        {
            copy(o);
        }

        //! Move constructor
        rlz_repair_v2(rlz_repair_v2&& o)
        {
            *this = std::move(o);
        }


        void swap(rlz_repair_v2& r)
        {
            if (this != &r) {
                m_reference.swap(r.m_reference);
                std::swap(m_input_size, r.m_input_size);
                std::swap(m_input_pos, r.m_input_pos);
                std::swap(m_repair_pos, r.m_repair_pos);
            }
        }

        rlz_repair_v2& operator=(const rlz_repair_v2& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        rlz_repair_v2& operator=(rlz_repair_v2&& v)
        {
            if (this != &v) {
                m_reference = std::move(v.m_reference);
                m_input_size = std::move(v.m_input_size);
                m_input_pos = v.m_input_pos;
                m_repair_pos = v.m_repair_pos;
            }
            return *this;
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            m_reference.serialize(out, child, "reference");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in){
            m_reference.load(in);
        }

        void print_reference(const size_type i, const size_type j){
            for(size_type k = i; k < j; ++k){
                std::cout << m_reference[k] << ", ";
            }
            std::cout << std::endl;
        }
    };

    using rlz_repair_v2_csa_bc_int64 = rlz_repair_v2<int64_t , reference_repair<int64_t, rlz_factor_type>, sdsl::csa_bitcompressed<sdsl::int_alphabet<>>>;

}

#endif //RCT_RLZ_NAIVE_HPP
