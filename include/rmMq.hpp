//
// Created by adrian on 11/04/17.
//

#ifndef CONTACT_RMMQ_SUCCINCT_CT_HPP
#define CONTACT_RMMQ_SUCCINCT_CT_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/bp_support.hpp>
#include <rmq_succinct_ct.hpp>
#include <sdsl/sd_vector.hpp>

namespace rct {

    template<class t_value = uint64_t, class t_sampling = sdsl::sd_vector<>>
    class rmMq {

    public:
        typedef typename sdsl::bit_vector::size_type size_type;
        typedef t_value value_type;
        typedef t_sampling sampling_type;
        typedef typename sampling_type::rank_1_type rank_sampling_type;
        typedef typename sampling_type::select_1_type select_sampling_type;
        typedef struct {
            value_type min;
            value_type max;
        } result_type;
        typedef struct {
            size_type min;
            size_type max;
            bool min_ok = false;
            bool max_ok = false;
        } pos_result_type;

    private:
        rmq_succinct_ct<true> m_rmq;
        rmq_succinct_ct<false> m_rMq;

        bool m_start_increasing;
        sampling_type m_sampling;
        rank_sampling_type m_rank;
        select_sampling_type m_select;
        const std::vector<value_type>* m_v;

    public:
        void set_vector(const std::vector<value_type>* vec){
            m_v = vec;
        }

    private:
        void construction(){

            if(m_v->size() < 1) return;
            bool increase;
            std::vector<size_type> ones;
            std::vector<value_type> max, min;
            for(size_type i = 0; i < m_v->size()-1; ++i){
                int64_t diff = static_cast<int64_t>(m_v->at(i+1)) - static_cast<int64_t>(m_v->at(i));
                if(i == 0){
                    increase = diff >= 0;
                    m_start_increasing = increase;
                }else{
                    if(increase && diff < 0){
                        increase = false;
                        ones.push_back(i);
                        max.push_back(m_v->at(i));
                    }else if(!increase && diff > 0){
                        increase = true;
                        ones.push_back(i);
                        min.push_back(m_v->at(i));
                    }
                }
            }
            ones.push_back(m_v->size());
            m_sampling = sampling_type(ones.begin(), ones.end());
            sdsl::util::init_support(m_rank, &m_sampling);
            sdsl::util::init_support(m_select, &m_sampling);
            m_rmq = rmq_succinct_ct<true>(min);
            m_rMq = rmq_succinct_ct<false>(max);
        }

        inline size_type rank_to_max(const size_type rank) const{
            return (rank - m_start_increasing)/2 + m_start_increasing;

        }

        inline size_type rank_to_min(const size_type rank) const{
            return (rank - !m_start_increasing)/2 + !m_start_increasing;
        }

        inline size_type select_max(const size_type i) const{
            return m_select(2*i - m_start_increasing);
        }

        inline size_type select_min(const size_type i) const{
            return m_select(2*i - !m_start_increasing);
        }

    public:

        rmMq(){};

        rmMq(const std::vector<value_type> *vec){
            set_vector(vec);
            construction();
        }

        size_type size(){
            return m_v->size();
        }

        //! Copy constructor
        rmMq(const rmMq& o)
        {
            copy(o);
        }

        //! Move constructor
        rmMq(rmMq&& o)
        {
            *this = std::move(o);
        }

        rmMq& operator=(const rmMq& o) {
            if (this != &o) {
                m_rmq = o.m_rmq;
                m_rMq = o.m_rMq;
                m_start_increasing = o.m_start_increasing;
                m_sampling = o.m_sampling;
                m_rank = o.m_rank;
                m_rank.set_vector(&m_sampling);
                m_select = o.m_select;
                m_select.set_vector(&m_sampling);
            }
            return *this;
        }

        rmMq& operator=(rmMq&& o) {
            if (this != &o) {
                m_rmq = std::move(o.m_rmq);
                m_rMq = std::move(o.m_rMq);
                m_start_increasing = std::move(o.m_start_increasing);
                m_sampling = std::move(o.m_sampling);
                m_rank = std::move(o.m_rank);
                m_rank.set_vector(&m_sampling);
                m_select = std::move(o.m_select);
                m_select.set_vector(&m_sampling);
            }
            return *this;
        }

        //! Swap method
        /*! Swaps the content of the two data structure.
         *  You have to use set_vector to adjust the supported bit_vector.
         *  \param bp_support Object which is swapped.
         */
        void swap(rmMq& o) {
            m_rmq.swap(o.m_rmq);
            m_rMq.swap(o.m_rMq);
            std::swap(m_start_increasing, o.m_start_increasing);
            std::swap(m_sampling, o.m_sampling);
            sdsl::util::swap_support(m_rank, o.m_rank, &m_sampling, &o.m_sampling);
            sdsl::util::swap_support(m_select, o.m_select, &m_sampling, &o.m_sampling);
        }

        value_type min(const size_type l, const size_type r) const {
            size_type min_l = 0, min_r = 0;
            auto rank_l = m_rank(l+1);
            if(rank_l){
                min_l = rank_to_min(rank_l)+1;
            }
            auto rank_r = m_rank(r);
            if(rank_r){
                min_r = rank_to_min(rank_r);
            }
            value_type minimum = m_v->at(l);
            if(min_l > 0 && min_r > 0 && min_l <= min_r){
                size_type min_pos = m_rmq(min_l-1, min_r-1);
                auto pos = select_min(min_pos+1);
                if(m_v->at(pos) < minimum){
                    minimum = m_v->at(pos);
                }
            }
            if(m_v->at(r) < minimum){
                minimum = m_v->at(r);
            }
            return minimum;
        }

        value_type max(const size_type l, const size_type r) const {
            size_type max_l = 0, max_r = 0;
            auto rank_l = m_rank(l+1);
            if(rank_l){
                max_l = rank_to_min(rank_l)+1;
            }
            auto rank_r = m_rank(r);
            if(rank_r){
                max_r = rank_to_min(rank_r);
            }
            value_type maximum = m_v->at(l);
            if(max_l > 0 && max_r > 0 && max_l <= max_r){
                size_type max_pos = m_rMq(max_l-1, max_r-1);
                auto pos = select_max(max_pos+1);
                if(m_v->at(pos) > maximum){
                    maximum = m_v->at(pos);
                }
            }
            if(m_v->at(r) > maximum){
                maximum = m_v->at(r);
            }
            return maximum;
        }

        pos_result_type locals(const size_type l, const size_type r) const{
            size_type min_l = 0, max_l = 0, min_r = 0, max_r = 0;
            pos_result_type result;
            auto rank_l = m_rank(l+1);
            if(rank_l){
                min_l = rank_to_min(rank_l)+1;
                max_l = rank_to_max(rank_l)+1;
            }else{
                min_l = 1;
                max_l = 1;
            }
            auto rank_r = m_rank(r);
            if(rank_r){
                min_r = rank_to_min(rank_r);
                max_r = rank_to_max(rank_r);
            }
            if(min_l > 0 && min_r > 0 && min_l <= min_r){
                if(min_l < min_r) {
                    size_type min_pos = m_rmq(min_l - 1, min_r - 1);
                    result.min = select_min(min_pos + 1);
                }else{
                    result.min = select_min(min_l);//TODO: revisar
                }
                result.min_ok = true;
            }

            if(max_l > 0 && max_r > 0 && max_l <= max_r) {
                if (max_l < max_r) {
                    size_type max_pos = m_rMq(max_l - 1, max_r - 1);
                    result.max = select_max(max_pos + 1);
                } else {
                    result.max = select_max(max_l);//TODO: revisar
                }
                result.max_ok = true;
            }
            return result;
        }

        result_type operator()(const size_type l, const size_type r) const{
            size_type min_l = 0, max_l = 0, min_r = 0, max_r = 0;
            auto rank_l = m_rank(l+1);
            if(rank_l){
                min_l = rank_to_min(rank_l)+1;
                max_l = rank_to_max(rank_l)+1;
            }else{
                min_l = 1;
                max_l = 1;
            }
            auto rank_r = m_rank(r);
            if(rank_r){
                min_r = rank_to_min(rank_r);
                max_r = rank_to_max(rank_r);
            }
            value_type minimum = m_v->at(l);
            if(min_l > 0 && min_r > 0 && min_l <= min_r){
                if(min_l < min_r) {
                    size_type min_pos = m_rmq(min_l - 1, min_r - 1);
                    auto pos = select_min(min_pos + 1);
                    if (m_v->at(pos) < minimum) {
                        minimum = m_v->at(pos);
                    }
                }else{
                    auto pos = select_min(min_l);//TODO: revisar
                    if(m_v->at(pos) < minimum){
                        minimum = m_v->at(pos);
                    }
                }
            }
            if(m_v->at(r) < minimum){
                minimum = m_v->at(r);
            }


            value_type maximum = m_v->at(l);
            if(max_l > 0 && max_r > 0 && max_l <= max_r){
                if(max_l < max_r){
                    size_type max_pos = m_rMq(max_l-1, max_r-1);
                    auto pos = select_max(max_pos+1);
                    if(m_v->at(pos) > maximum){
                        maximum = m_v->at(pos);
                    }
                }else{
                    auto pos = select_max(max_l);//TODO: revisar
                    if(m_v->at(pos) > maximum){
                        maximum = m_v->at(pos);
                    }
                }

            }
            if(m_v->at(r) > maximum){
                maximum = m_v->at(r);
            }

            return result_type{minimum, maximum};
        }

        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
            sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_rmq.serialize(out, child, "rmq");
            written_bytes += m_rMq.serialize(out, child, "rMq");
            written_bytes += sdsl::write_member(m_start_increasing, out, child, "start_increasing");
            written_bytes += m_sampling.serialize(out, child, "sampling");
            written_bytes += m_rank.serialize(out, child, "rank");
            written_bytes += m_select.serialize(out, child, "select");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in){
            m_rmq.load(in);
            m_rMq.load(in);
            sdsl::load(m_start_increasing, in);
            m_sampling.load(in);
            m_rank.load(in, &m_sampling);
            m_select.load(in, &m_sampling);
        }


    };
}

#endif //RLZ_REFERENCE_RMQ_SUCCINCT_CT_HPP
