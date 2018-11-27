//
// Created by adrian on 11/04/17.
//

#ifndef RCT_RMQ_SUCCINCT_CT_HPP
#define RCT_RMQ_SUCCINCT_CT_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/bp_support.hpp>

namespace rct {

    template <bool t_min = true>
    class rmq_succinct_ct {

    public:
        typedef typename sdsl::bit_vector::size_type size_type;
        typedef typename sdsl::bit_vector::size_type value_type;

    private:
        sdsl::bit_vector m_bp;
        sdsl::bp_support_sada<> m_bp_support;
        size_type m_size;

    private:
        template<class t_rac>
        void build_cartesian_tree(const t_rac& vec){

            typedef typename t_rac::size_type size_type;
            m_bp.resize(2*(vec.size()+1));      // resize bit vector for balanced pare ntheses to 2(n+1) bits
            sdsl::util::set_to_value(m_bp, 0);
            std::stack<typename t_rac::value_type> vec_stack;

            size_type k=1;
            m_bp[0]=1; //extra root
            for (size_type i=0; i < vec.size(); ++i) {
                sdsl::int_vector<>::value_type l = vec[i];
                if (t_min) {
                    while (vec_stack.size() > 0 and l < vec_stack.top()) {
                        vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
                    }
                } else {
                    while (vec_stack.size() > 0 and l > vec_stack.top()) {
                        vec_stack.pop(); ++k; /*bp[k++] = 0; bp is already initialized to zero*/ // writing a closing parenthesis
                    }

                }
                vec_stack.push(l);
                m_bp[k++] = 1;
            }

        }

    public:

        rmq_succinct_ct(){};

        template<class t_rac>
        rmq_succinct_ct(const t_rac vec){
            build_cartesian_tree(vec);
            sdsl::util::init_support(m_bp_support, &m_bp);
        }

        size_type size(){
            return m_bp.size();
        }

        rmq_succinct_ct& operator=(const rmq_succinct_ct& rm) {
            if (this != &rm) {
                m_bp = rm.m_bp;
                m_bp_support = rm.m_bp_support;
                m_bp_support.set_vector(&m_bp);
            }
            return *this;
        }

        rmq_succinct_ct& operator=(rmq_succinct_ct&& rm) {
            if (this != &rm) {
                m_bp = std::move(rm.m_bp);
                m_bp_support = std::move(rm.m_bp_support);
                m_bp_support.set_vector(&m_bp);
            }
            return *this;
        }

        //! Swap method
        /*! Swaps the content of the two data structure.
         *  You have to use set_vector to adjust the supported bit_vector.
         *  \param bp_support Object which is swapped.
         */
        void swap(rmq_succinct_ct& p) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            m_bp.swap(p.m_bp);
            sdsl::util::swap_support(m_bp_support, p.m_bp_support, &m_bp, &p.m_bp);
        }


        size_type operator()(const size_type l, const size_type r) const{
            size_type i = m_bp_support.select(l+2)-1;
            size_type j = m_bp_support.select(r+2);
            size_type rmq = m_bp_support.rmq(i, j);
            return m_bp_support.rank(rmq)-1;
        }

        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_bp.serialize(out, child, "bp");
            written_bytes += m_bp_support.serialize(out, child, "bp_support");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in){
            m_bp.load(in);
            m_bp_support.load(in, &m_bp);
        }


    };
}

#endif //RLZ_REFERENCE_RMQ_SUCCINCT_CT_HPP
