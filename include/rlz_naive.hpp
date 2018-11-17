//
// Created by adrian on 13/11/18.
//

#ifndef RCT_RLZ_NAIVE_HPP
#define RCT_RLZ_NAIVE_HPP

#include <file_util.hpp>
#include <reference_uniform_sample.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

#define INITIAL_CAPACITY 1024

namespace rct {

    template <class t_value = uint32_t, class t_reference = reference_uniform_sample<t_value>,
              class t_fm_index = sdsl::csa_wt_int<> >
    class rlz_naive {

    public:


        typedef uint32_t offset_type;
        typedef uint32_t length_type;
        typedef struct rlz_factor {
            offset_type offset;
            length_type length;
        } factor_type;
        typedef t_reference reference_type;
        typedef t_value value_type;
        typedef t_fm_index fm_index_type;
        typedef typename fm_index_type::size_type size_type;

    private:

        reference_type m_reference;
        fm_index_type  m_fm_index;

        size_type m_input_size = 0;
        size_type m_input_pos = 0;
        const std::vector<value_type>* m_input;

        void copy(const rlz_naive &n){
            m_reference = n.m_reference;
            m_fm_index = n.m_fm_index;
            m_input_size = n.m_input_size;
            m_input_pos = n.m_input_pos;
        }


    public:
        const reference_type &reference = m_reference;

        rlz_naive(std::vector<value_type> &container, const size_type reference_size){

            m_reference = reference_type(container, reference_size);
            sdsl::int_vector<> rev_reference;
            //std::vector<uint32_t> rev_reference;
           // rev_reference.width((uint8_t) (sizeof(value_type)*8));
           /* for(size_type i = 0; i < m_reference.size(); ++i){
                std::cout << "ref[" << i << "]=" << m_reference[i] << std::endl;
            }*/
            rev_reference.resize(m_reference.size());
            std::reverse_copy(m_reference.begin(), m_reference.end(), rev_reference.begin());
           // util::file::write_to_file("rev_ref.txt", rev_reference);
            sdsl::construct_im(m_fm_index, rev_reference);
        }


        void init_factorization(const std::vector<value_type> *container){
            m_input = container;
            m_input_size = container->size();
            m_input_pos = 0;
        }


        rlz_factor next(){
            size_type start = 0;
            size_type end = m_fm_index.size()-1;
            size_type start_input = m_input_pos;
            while(m_input_pos < m_input_size){

                auto sym = m_input->at(m_input_pos);
                size_type res_start, res_end;
                if(start == 0 && end == m_fm_index.size()-1){
                    auto sym_comp = m_fm_index.char2comp[sym];
                    res_start = m_fm_index.C[sym_comp];
                    res_end = m_fm_index.C[sym_comp+1]-1;
                }else{
                    sdsl::backward_search(m_fm_index, start, end, sym, res_start, res_end);
                }
                if(res_end <= res_start){
                    length_type length = m_input_pos - start_input;
                    if(length == 0) {
                        ++m_input_pos;
                        offset_type offset = m_reference.size() - m_fm_index[res_start]-1;
                        return rlz_factor{offset, 1};
                    }
                    offset_type offset = m_reference.size() - (m_fm_index[start] + length);
                    return rlz_factor{offset, length};

                }else{
                    start = res_start;
                    end = res_end;
                    ++m_input_pos;
                }
            }
            length_type length = m_input_pos - start_input;
            offset_type offset = m_reference.size() - (m_fm_index[start] + length);
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
                std::memcpy(&result[pos], m_reference.data(f.offset), f.length * sizeof(value_type));
                pos += f.length;
            }
        }


        void swap(rlz_naive& r)
        {
            if (this != &r) {
                m_reference.swap(r.m_reference);
                m_fm_index.swap(r.m_fm_index);
                std::swap(m_input_size, r.m_input_size);
                std::swap(m_input_pos, r.m_input_pos);
            }
        }

        rlz_naive& operator=(const rlz_naive& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        rlz_naive& operator=(rlz_naive&& v)
        {
            if (this != &v) {
                m_reference = std::move(v.m_reference);
                m_fm_index = std::move(v.m_fm_index);
                m_input_size = std::move(v.m_input_size);
                m_input_pos = v.m_input_pos;
            }
            return *this;
        }
    };
}

#endif //RCT_RLZ_NAIVE_HPP
