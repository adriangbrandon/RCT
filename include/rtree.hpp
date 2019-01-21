//
// Created by adrian on 13/07/18.
//

#ifndef SUCCINCTCT_RTREE_HPP
#define SUCCINCTCT_RTREE_HPP

#include <vector>
#include <queue>
#include <geo_util.hpp>
#include <SRTree.h>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>

using namespace SpatialIndex;

namespace rct {



    class rtree {

        struct queue_element {
            util::geo::region r;
            uint64_t distance;
            SNode* ptr;
            uint32_t id;
            int32_t pl0, pl1;
        };

        class compare_queue_region {
        public:

            bool operator()(const queue_element& lhs, const queue_element& rhs) const
            {
                return lhs.distance > rhs.distance;
            }
        };

        class compare_queue_object {
        public:

            bool operator()(const queue_element& lhs, const queue_element& rhs) const
            {
                return lhs.distance < rhs.distance;
            }
        };

    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;
        typedef queue_element queue_element_type;
        typedef std::priority_queue<queue_element, std::vector<queue_element>, compare_queue_region> queue_region_type;
        typedef std::priority_queue<queue_element, std::vector<queue_element>, compare_queue_object> queue_object_type;

    private:
        RTree::SRTree* m_tSRTree = nullptr;

    private:
        void copy(const rtree& p){
           m_tSRTree = p.m_tSRTree;
        }

        /*void add_to_region_queue(const std::vector<int> &l0, const std::vector<int> &l1, const std::vector<int> &h0,
                                   const std::vector<int> &h1, const std::vector<SNode*> &nodes,
                                   queue_region_type &queue_region, const util::geo::point &p_q) const{
            for(size_type i = 0; i < nodes.size(); ++i){
                queue_element qe;
                util::geo::region r{l0[i], l1[i], h0[i], h1[i]};
                qe.r = r;
                qe.ptr = nodes[i];
                qe.pl0 = l0[i];
                qe.pl1 = l1[i];
                qe.distance = util::distance_knn(r, p_q);
                queue_region.push(qe);
            }
        }

        void add_to_object_queue(const std::vector<int> &l0, const std::vector<int> &l1, const std::vector<int> &h0,
                                 const std::vector<int> &h1, const std::vector<int> &ids,
                                 queue_object_type &queue_object, const point &p_q, size_type &n_disap,
                                 const sdsl::bit_vector disap) const{
            for(size_type i = 0; i < ids.size(); ++i){
                queue_element qe;
                region r(l0[i], l1[i], h0[i], h1[i]);
                qe.r = r;
                qe.distance = util::distance_knn(r, p_q);
                qe.id = ids[i];
                if(disap[qe.id]) ++n_disap;
                queue_object.push(qe);
            }
        }*/

    public:

        rtree(){};
        rtree(const std::vector<std::pair<value_type, util::geo::region>> &data, const size_t capacity = 30){

            IStorageManager* memorymanager = StorageManager::createNewMemoryStorageManager();
            id_type indexIdentifier;
            ISpatialIndex* tree = RTree::createNewRTree(*memorymanager, 0.7, capacity, capacity, 2,
                                                        SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
            double plow[2], phigh[2];
            for(size_t i = 0; i < data.size(); ++i){
                auto region = data[i].second;
                plow[0]=region.min.x; plow[1]=region.min.y;
                phigh[0]=region.max.x; phigh[1]=region.max.y;
                Region r = Region(plow, phigh, 2);
                tree->insertData(0, 0, r, data[i].first);
            }
            size_t disc_factor = 1;
            size_t param = 1;
            int index_type = SRTREE_TYPE_COMP;
            m_tSRTree = new RTree::SRTree((RTree::RTree*)tree, disc_factor, index_type, param);
            delete tree;
            delete memorymanager;
        }

        void intersection(const util::geo::region &r_q, std::vector<value_type> &ret) const {
            int size = 0;
            long * results = new long[m_tSRTree->getNumObj()];
            m_tSRTree->range_query(r_q.min.x, r_q.min.y, r_q.max.x, r_q.max.y, results, &size);
            ret = std::vector<value_type>(results, results+size);
            delete [] results;
        }


       /* void knn_candidates(const size_type k, const point &p_q, util::pq_knn_candidate &candidates,
                            const sdsl::bit_vector &disap) const{

            int pl0 = m_tSRTree->lower0;
            int pl1 = m_tSRTree->lower1;

            //First nodes
            queue_region_type priority_region_queue;
            queue_object_type priority_object_queue;
            queue_element qe;
            qe.ptr = m_tSRTree->root;
            qe.pl0 = pl0;
            qe.pl1 = pl1;
            priority_region_queue.push(qe);
            size_type n_disap = 0;
            while(!priority_region_queue.empty() && (priority_object_queue.size() < (k + n_disap)
                  || priority_object_queue.top().distance >= priority_region_queue.top().distance)){

                auto element = priority_region_queue.top();
                priority_region_queue.pop();
                std::vector<int> l0, l1, h0, h1, ids;
                std::vector<SNode*> nodes;
                element.ptr->children_region(element.pl0, element.pl1, l0, l1, h0, h1, nodes, ids);
                if(ids.empty()){
                    add_to_region_queue(l0, l1, h0, h1, nodes, priority_region_queue, p_q);
                }else{
                    add_to_object_queue(l0, l1, h0, h1, ids, priority_object_queue, p_q, n_disap, disap);
                }
            }
            while(!priority_object_queue.empty()){
                queue_element_type c = priority_object_queue.top();
                util::pq_knn_candidate_element candidate(c.id, 0, 0, 0, 0, c.distance);
                candidates.push(candidate);
                priority_object_queue.pop();
            }
        }*/


        //! Assignment move operation
        rtree& operator=(rtree&& p) {
            if (this != &p) {
                m_tSRTree = p.m_tSRTree;
            }
            return *this;
        }

        //! Assignment operator
        rtree& operator=(const rtree& snapshot)
        {
            if (this != &snapshot) {
                copy(snapshot);
            }
            return *this;
        }

        //! Copy constructor
        rtree(const rtree& p)
        {
            copy(p);
        }

        //! Move constructor
        rtree(rtree&& snapshot)
        {
            *this = std::move(snapshot);
        }

        //! Swap method
        /*! Swaps the content of the two data structure.
         *  You have to use set_vector to adjust the supported bit_vector.
         *  \param bp_support Object which is swapped.
         */
        void swap(rtree& p)
        {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            std::swap(m_tSRTree, p.m_tSRTree);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_tSRTree->memory_usage();
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

    };
}


#endif //SUCCINCTCT_RTREE_HPP
