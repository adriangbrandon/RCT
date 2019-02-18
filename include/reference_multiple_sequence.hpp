//
// Created by adrian on 13/11/18.
//

#ifndef RCT_REFERENCE_MULTIPLE_SEQUENCE_HPP
#define RCT_REFERENCE_MULTIPLE_SEQUENCE_HPP

#include <vector>
#include <unordered_map>
#include <sdsl/int_vector.hpp>
#include <sdsl/csa_bitcompressed.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <file_util.hpp>

namespace rct {


    template<class t_value = uint64_t, class t_sa = sdsl::csa_bitcompressed<sdsl::int_alphabet<>>>
    class reference_multiple_sequence {

    public:
        typedef uint64_t size_type;
        typedef t_value value_type;
        typedef t_sa sa_type;
        typedef typename std::vector<value_type>::iterator iterator;
        typedef std::unordered_map<size_type, size_type> char2comp_type;
    private:
        std::vector<value_type> m_reference;

        void copy(const reference_multiple_sequence &n){
            m_reference = n.m_reference;
        }

        void prepare_compute_factors(sa_type &csa_bit, char2comp_type &char2comp, const std::vector<value_type> &reference){
            sdsl::util::clear(csa_bit);
            sdsl::util::clear(char2comp);
            std::string file_rev_reference = "rev_reference_"+std::to_string(getpid())+".sdsl";
            {
                std::map<size_type, size_type> D;
                // count occurrences of each symbol
                for(auto it = reference.begin(); it != reference.end(); ++it) {
                    auto value = *it;
                    D[value]++;
                    if(value ==  0){
                        std::cout << "Contains zero symbol (reference)" << std::endl;
                        exit(0);
                    }

                }

                size_type index = 1;
                for(auto it = D.begin(); it != D.end(); ++it){
                    char2comp[it->first] = index;
                    ++index;
                }
                char2comp[0] = 0;
                sdsl::int_vector<> rev_reference;
                rev_reference.resize(reference.size());
                for(size_type i = 0; i < rev_reference.size() ; ++i){
                    auto value = *(reference.begin() + reference.size()-1-i);
                    rev_reference[i] = char2comp[value];
                }
                /*for(auto v : rev_reference){
                    if(v == 0){
                        std::cout << "Contains zero symbol (rev_reference)" << std::endl;
                        exit(0);
                    }
                }*/
                sdsl::store_to_plain_array<value_type>(rev_reference, file_rev_reference);
            }
            //std::cout << "building SA " << std::endl;
            sdsl::construct(csa_bit, file_rev_reference.c_str(), sizeof(value_type));
            //std::cout << "done" << std::endl;
            util::file::remove_file(file_rev_reference);
        }

        bool count_factors(iterator &it, const iterator &end_input, sa_type &csa_bit, char2comp_type &char2comp, size_type &counter,
                           std::vector<value_type> &reference){

            size_type start = 0;
            size_type end = csa_bit.size()-1;
            size_type pos = 0;
            const iterator start_input = it;
            while(it != end_input){
                auto sym = *it;
                //std::cout << "Symbol: " << sym << std::endl;
                if(char2comp.count(sym) == 0){
                    if(start_input != it) ++counter;
                    reference.push_back(sym);
                    prepare_compute_factors(csa_bit, char2comp, reference);
                    ++counter;
                    return true;
                }
                auto sym_comp = char2comp[sym];
                size_type res_start, res_end;
                if(start == 0 && end == csa_bit.size()-1){
                    res_start = csa_bit.C[sym_comp];
                    res_end = csa_bit.C[sym_comp+1]-1;
                }else{
                    /*std::cout << "pos: " << pos << std::endl;
                    std::cout << "start: " << start << " end: " << end << std::endl;*/
                    sdsl::backward_search(csa_bit, start, end, sym_comp, res_start, res_end);
                    //std::cout << "res_start: " << res_start << " res_end: " << res_end << std::endl;
                }
                if(res_end < res_start){
                    ++counter;
                    return true;
                }else{
                    start = res_start;
                    end = res_end;
                    ++it; ++pos;
                }
            }
            ++counter;
            return false;
        }

        size_type compute_factors(const iterator &beg,  const iterator &end, std::vector<value_type> &reference){
            sa_type csa_bit;
            char2comp_type char2comp;
            prepare_compute_factors(csa_bit, char2comp, reference);
            size_type factors = 0;
            auto it = beg;
            while(count_factors(it, end, csa_bit,char2comp, factors, reference)){
                //nothing to do
            }
            return factors;
        }

    public:

        reference_multiple_sequence() = default;

        reference_multiple_sequence(std::vector<value_type> &container, const std::vector<size_type> &lengths,
                                    const size_type reference_size=0 , const size_type block_size = 0, const double_t ratio = 0.08){


            assert(!lengths.empty());
            size_type traj_done = 1;
            iterator beg = container.begin();
            iterator end = container.begin() + lengths[0];
            m_reference.resize(lengths[0]);
            std::copy(beg, end, m_reference.begin());
            std::cout << "size: " << lengths.size() << std::endl;
            while(traj_done < lengths.size()){
                if(lengths[traj_done] > 0){
                    beg = end;
                    end = beg + lengths[traj_done];
                    for(auto ita = beg; ita != end; ++ita){
                        if(*ita == 0){
                            std::cout << "Contains zero." << std::endl;
                            exit(1);
                        }
                    }
                    std::cout << "length: " << lengths[traj_done] << std::endl;
                    std::vector<value_type> reference(m_reference);
                    auto factors = compute_factors(beg, end, reference);
                    double_t r = factors / (double_t) lengths[traj_done];
                    std::cout << "R: " << r << std::endl;
                    if(r > ratio){
                        std::cout << "new reference" << std::endl;
                        auto size_ref = std::distance(m_reference.begin(), m_reference.end());
                        m_reference.resize(m_reference.size() + lengths[traj_done]);
                        std::copy(beg, end, m_reference.begin()+size_ref);
                        auto pos = 0;
                        for(auto v : m_reference){
                            if(v == 0){
                                std::cout << "Reference contains 0 after update" << std::endl;
                                std::cout << "Position: " << pos << std::endl;
                                exit(0);
                            }
                            ++pos;
                        }
                    }else{
                        m_reference = reference;
                    }
                    std::cout << "size reference: " << m_reference.size() << std::endl;
                }
                ++traj_done;

            }
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
        reference_multiple_sequence(const reference_multiple_sequence& o)
        {
            copy(o);
        }

        //! Move constructor
        reference_multiple_sequence(reference_multiple_sequence&& o)
        {
            *this = std::move(o);
        }

        void swap(reference_multiple_sequence& r)
        {
            if (this != &r) {
                std::swap(m_reference, r.m_reference);
            }
        }

        reference_multiple_sequence& operator=(const reference_multiple_sequence& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        reference_multiple_sequence& operator=(reference_multiple_sequence&& v)
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

#endif //RCT_REFERENCE_MULTIPLE_SEQUENCE_HPP
