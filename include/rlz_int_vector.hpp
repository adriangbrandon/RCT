//
// Created by adrian on 13/11/18.
//

#ifndef RCT_RLZ_INT_VECTOR_HPP
#define RCT_RLZ_INT_VECTOR_HPP

#include <file_util.hpp>
#include <unordered_map>
#include <reference_uniform_sample.hpp>
#include <reference_repair.hpp>
#include <reference_multiple_sequence.hpp>
#include <sdsl/suffix_arrays.hpp>

#define INITIAL_CAPACITY 1024

namespace rct {

    template <class t_value = uint32_t, class t_csa = sdsl::csa_wt_int<> >
    class rlz_int_vector {

    public:


        typedef uint32_t offset_type;
        typedef uint32_t length_type;
        typedef struct rlz_factor {
            offset_type offset;
            length_type length;
        } factor_type;
        typedef sdsl::int_vector<> reference_type;
        typedef t_value value_type;
        typedef t_csa csa_type;
        typedef typename csa_type::size_type size_type;
        typedef std::unordered_map<size_type, size_type> char2comp_type;

    private:

        reference_type m_reference;
        csa_type  m_csa;
        char2comp_type m_char2comp;

        size_type m_input_size = 0;
        size_type m_input_pos = 0;
        const std::vector<value_type>* m_input;

        void copy(const rlz_int_vector &n){
            m_reference = n.m_reference;
            m_csa = n.m_csa;
            m_input_size = n.m_input_size;
            m_input_pos = n.m_input_pos;
        }


    public:
        const reference_type &reference = m_reference;

        rlz_int_vector() = default;

        rlz_int_vector(const reference_type &reference){
            m_reference = reference;
            std::string file_rev_reference =  std::to_string(getpid()) + ".rev_ref";
            {
                std::map<size_type, size_type> D;
                // count occurrences of each symbol
                for(auto it = m_reference.begin(); it != m_reference.end(); ++it) {
                    auto value = *it;
                    D[value]++;
                }

                size_type index = 1;
                for(auto it = D.begin(); it != D.end(); ++it){
                    m_char2comp[it->first] = index;
                    ++index;
                }
                m_char2comp[0] = 0;
                sdsl::int_vector<> rev_reference;
                rev_reference.resize(m_reference.size());
                for(size_type i = 0; i < rev_reference.size() ; ++i){
                    auto value = *(m_reference.begin() + m_reference.size()-1-i);
                    rev_reference[i] = m_char2comp[value];
                }
                sdsl::store_to_plain_array<value_type>(rev_reference, file_rev_reference);
            }
            std::cout << "building SA " << std::endl;
            sdsl::construct(m_csa, file_rev_reference.c_str(), sizeof(value_type));
            std::cout << "done" << std::endl;
        }


        void init_factorization(const std::vector<value_type> *container){
            m_input = container;
            m_input_size = container->size();
            m_input_pos = 0;
        }


        rlz_factor next(){
            size_type start = 0;
            size_type end = m_csa.size()-1;
            size_type start_input = m_input_pos;
            while(m_input_pos < m_input_size){

                auto sym = m_input->at(m_input_pos);
                auto sym_comp = m_char2comp[sym];
                size_type res_start, res_end;
                if(start == 0 && end == m_csa.size()-1){
                    res_start = m_csa.C[sym_comp];
                    res_end = m_csa.C[sym_comp+1]-1;
                }else{
                    sdsl::backward_search(m_csa, start, end, sym_comp, res_start, res_end);
                }
                if(res_end < res_start){
                    length_type length = m_input_pos - start_input;
                    if(length == 0) {
                        ++m_input_pos;
                        offset_type offset = m_reference.size() - m_csa[res_start]-1;
                        return rlz_factor{offset, 1};
                    }
                    offset_type offset = m_reference.size() - (m_csa[start] + length);
                    return rlz_factor{offset, length};

                }else{
                    start = res_start;
                    end = res_end;
                    ++m_input_pos;
                }
            }
            length_type length = m_input_pos - start_input;
            offset_type offset = m_reference.size() - (m_csa[start] + length);
            return rlz_factor{offset, length};

        }


        inline bool has_next() const {
            return m_input_pos < m_input_size;
        }


         void decompress(const std::vector<rlz_factor> &factors, std::vector<value_type> &result){
            size_type pos = 0;
            result.resize(m_input_size);
            for(const auto &f : factors){
                /*std::cout << "pos: " << pos << "result.size(): " << result.size() << std::endl;
                std::cout << "f: " << f.offset << ", " << f.length << " reference.size(): " << m_reference.size() << std::endl;*/
                for(size_type i = f.offset; i < f.offset + f.length; ++i){
                    result[pos++] = m_reference[i];
                }
            }
        }

        //! Copy constructor
        rlz_int_vector(const rlz_int_vector& o)
        {
            copy(o);
        }

        //! Move constructor
        rlz_int_vector(rlz_int_vector&& o)
        {
            *this = std::move(o);
        }


        void swap(rlz_int_vector& r)
        {
            if (this != &r) {
                m_reference.swap(r.m_reference);
                m_csa.swap(r.m_csa);
                std::swap(m_input_size, r.m_input_size);
                std::swap(m_input_pos, r.m_input_pos);
            }
        }

        rlz_int_vector& operator=(const rlz_int_vector& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        rlz_int_vector& operator=(rlz_int_vector&& v)
        {
            if (this != &v) {
                m_reference = std::move(v.m_reference);
                m_csa = std::move(v.m_csa);
                m_input_size = std::move(v.m_input_size);
                m_input_pos = v.m_input_pos;
            }
            return *this;
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            m_reference.serialize(out, child, "reference");
            m_csa.serialize(out, child, "csa");
            sdsl::write_member(m_input_pos, out, child, "input_pos");
            sdsl::write_member(m_input_size, out, child, "input_size");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in){
            size_type size;
            m_reference.load(in);
            m_csa.load(in);
            sdsl::read_member(m_input_pos, in);
            sdsl::read_member(m_input_size, in);
        }

        void print_reference(const size_type i, const size_type j){
            for(size_type k = i; k < j; ++k){
                std::cout << m_reference[k] << ", ";
            }
            std::cout << std::endl;
        }
    };

    //using rlz_csa_int64 = rlz_naive<uint64_t , reference_uniform_sample<uint64_t>, sdsl::int_vector<>>;

}

#endif //RCT_RLZ_NAIVE_HPP
