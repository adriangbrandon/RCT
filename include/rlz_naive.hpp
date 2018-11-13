//
// Created by adrian on 13/11/18.
//

#ifndef RCT_RLZ_NAIVE_HPP
#define RCT_RLZ_NAIVE_HPP

#include <reference_uniform_sample.hpp>
#include <sdsl/suffix_arrays.hpp>

namespace rct {

    template <class t_reference = reference_uniform_sample, class t_fm_index = sdsl::csa_wt<>>
    class rlz_naive {

    public:

        typedef uint32_t size_type;
        typedef uint32_t offset_type;
        typedef uint32_t length_type;
        typedef struct rlz_factor {
            offset_type offset;
            length_type length;
        } factor_type;
        typedef t_reference reference_type;
        typedef typename reference_type::value_type value_type;
        typedef t_fm_index fm_index_type;

    private:

        reference_type m_reference;
        fm_index_type  m_fm_index;
        std::vector<value_type> m_input;
        size_type m_input_pos = 0;

        void copy(const rlz_naive &n){
            m_reference = n.m_reference;
            m_fm_index = n.m_fm_index;
            m_input = n.m_input;
            m_input_pos = n.m_input_pos;
        }


    public:

        rlz_naive(std::vector<value_type> &container, const size_type reference_size){

            m_reference = reference_type(container, reference_size);
            std::vector<value_type > rev_reference;
            rev_reference.resize(m_reference.size());
            std::reverse_copy(m_reference.begin(), m_reference.end(), rev_reference.begin());
            rev_reference[rev_reference.size()-1] = 0;
            sdsl::construct_im(m_fm_index, rev_reference, (uint8_t) sizeof(value_type));
            m_input = std::move(container);
        }


        rlz_factor next(){
            size_type start = 0;
            size_type end = m_fm_index.size()-1;
            size_type start_input = m_input_pos;
            while(m_input_pos != end){
                auto sym = m_fm_index.char2comp[m_input[m_input_pos]];
                size_type res_start, res_end;
                if(start == 0 && end == m_fm_index.size()-1){
                    res_start = m_fm_index.C[sym];
                    res_end = m_fm_index.C[sym+1]-1;
                }else{
                    sdsl::backward_search(m_fm_index, start, end, sym, res_start, res_end);
                }
                if(res_end <= res_start){
                    length_type length = m_input_pos - start_input;
                    if(length == 0) ++m_input_pos;
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
            return m_input_pos < m_input.size();
        }


        void swap(rlz_naive& r)
        {
            if (this != &r) {
                m_reference.swap(r.m_reference);
                m_fm_index.swap(r.m_fm_index);
                std::swap(m_input, r.m_input);
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
                m_input = std::move(v.m_input);
                m_input_pos = v.m_input_pos;
            }
            return *this;
        }
    };
}

#endif //RCT_RLZ_NAIVE_HPP
