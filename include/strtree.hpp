//
// Created by adrian on 13/07/18.
//

#ifndef SUCCINCTCT_STRTREE_HPP
#define SUCCINCTCT_STRTREE_HPP

#include <vector>
#include <queue>
#include <geometry.hpp>
#include <SRTree.h>

using namespace SpatialIndex;

namespace rct {

    class strtree {

        struct queue_element {
            region r;
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

        class MyDataStream : public IDataStream
        {
        public:
            MyDataStream(std::vector<std::pair<uint32_t, region>> &data) : m_pNext(0)
            {
                processed = 0;
                m_data = data;

                if (!data.empty()){
                    id_type id;
                    double low[2], high[2];
                    id = data[processed].first;
                    low[0] = data[processed].second.m_min_p.m_x;
                    low[1] = data[processed].second.m_min_p.m_y;
                    high[0] = data[processed].second.m_max_p.m_x;
                    high[1] = data[processed].second.m_max_p.m_y;
                    Region r = Region(low, high, 2);
                    m_pNext = new RTree::Data(0, 0, r, id);
                }
            }

            virtual ~MyDataStream()
            {
                if (m_pNext != 0) delete m_pNext;
            }

            virtual IData* getNext()
            {

                if (m_pNext == 0) return 0;

                RTree::Data* ret = m_pNext;
                m_pNext = 0;

                if (m_data.empty()) {
                    throw Tools::IllegalArgumentException("Empty data.");
                }

                processed++;
                if(m_data.size() != processed){
                    id_type id;
                    double low[2], high[2];
                    id = m_data[processed].first;
                    low[0] = m_data[processed].second.m_min_p.m_x;
                    low[1] = m_data[processed].second.m_min_p.m_y;
                    high[0] = m_data[processed].second.m_max_p.m_x;
                    high[1] = m_data[processed].second.m_max_p.m_y;
                    Region r = Region(low, high, 2);
                    m_pNext = new RTree::Data(0, 0, r, id);
                }
                return ret;
            }

            virtual bool hasNext() throw (Tools::NotSupportedException)
            {
                return (m_data.size() > 0 && m_data.size() > processed);
            }

            virtual uint32_t size() throw (Tools::NotSupportedException)
            {
                return (uint32_t) m_data.size();
            }

            virtual void rewind() throw (Tools::NotSupportedException)
            {
                id_type id;
                double low[2], high[2];

                if (m_pNext != 0)
                {
                    delete m_pNext;
                    m_pNext = 0;
                }

                processed = 0;

                if (m_data.empty())
                    throw Tools::IllegalArgumentException("Empty data.");

                id = m_data[processed].first;
                low[0] = m_data[processed].second.m_min_p.m_x;
                low[1] = m_data[processed].second.m_min_p.m_y;
                high[0] = m_data[processed].second.m_max_p.m_x;
                high[1] = m_data[processed].second.m_max_p.m_y;
                Region r = Region(low, high, 2);
                m_pNext = new RTree::Data(0, 0, r, id);
            }

            std::vector<std::pair<uint32_t, region>> m_data;
            RTree::Data* m_pNext;
            uint32_t processed;
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
        void copy(const strtree& p){
           m_tSRTree = p.m_tSRTree;
        }

        void add_to_region_queue(const std::vector<int> &l0, const std::vector<int> &l1, const std::vector<int> &h0,
                                   const std::vector<int> &h1, const std::vector<SNode*> &nodes,
                                   queue_region_type &queue_region, const point &p_q){
            for(size_type i = 0; i < nodes.size(); ++i){
                queue_element qe;
                region r(l0[i], l1[i], h0[i], h1[i]);
                qe.r = r;
                qe.ptr = nodes[i];
                qe.pl0 = l0[i];
                qe.pl1 = l1[i];
                qe.distance = static_cast<uint64_t >(std::ceil(util::distance(r, p_q)));
                queue_region.push(qe);
            }
        }

        void add_to_object_queue(const std::vector<int> &l0, const std::vector<int> &l1, const std::vector<int> &h0,
                                 const std::vector<int> &h1, const std::vector<int> &ids,
                                 queue_object_type &queue_object, const point &p_q, size_type &n_disap,
                                 const sdsl::bit_vector disap){
            for(size_type i = 0; i < ids.size(); ++i){
                queue_element qe;
                region r(l0[i], l1[i], h0[i], h1[i]);
                qe.r = r;
                qe.distance = static_cast<uint64_t >(std::ceil(util::distance(r, p_q)));
                qe.id = ids[i];
                if(disap[qe.id]) ++n_disap;
                queue_object.push(qe);
            }
        }

    public:

        strtree(){};
        strtree(std::vector<std::pair<value_type, region>> &data, const size_t capacity = 30){

            IStorageManager* memorymanager = StorageManager::createNewMemoryStorageManager();
            id_type indexIdentifier;
            MyDataStream stream(data);
            ISpatialIndex* tree = RTree::createAndBulkLoadNewRTree(
                    RTree::BLM_STR, stream,
                    *memorymanager, 0.7, capacity, capacity, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
            size_t disc_factor = 1;
            size_t param = 1;
            int index_type = SRTREE_TYPE_COMP;
            m_tSRTree = new RTree::SRTree((RTree::RTree*)tree, disc_factor, index_type, param);
            std::cout << m_tSRTree->memory_usage() << std::endl;
            delete tree;
            delete memorymanager;
        }

        void intersection(const region &r_q, std::vector<value_type> &ret){
            int size = 0;
            long * results = new long[m_tSRTree->getNumObj()];
            m_tSRTree->range_query(r_q.m_min_p.m_x, r_q.m_min_p.m_y, r_q.m_max_p.m_x, r_q.m_max_p.m_y, results, &size);
            ret = std::vector<value_type>(results, results+size);
            delete [] results;
        }

        void knn_candidates(const size_type k, const point &p_q, std::vector<queue_element_type> &candidates,
                            const sdsl::bit_vector &disap){

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
            candidates.resize(priority_object_queue.size());
            for(size_type i = 1; i <= candidates.size(); ++i){
                candidates[candidates.size()-i] = priority_object_queue.top();
                priority_object_queue.pop();
            }
        }


        //! Assignment move operation
        strtree& operator=(strtree&& p) {
            if (this != &p) {
                std::swap(m_tSRTree, p.m_tSRTree);
            }
            return *this;
        }

        //! Assignment operator
        strtree& operator=(const strtree& snapshot)
        {
            if (this != &snapshot) {
                copy(snapshot);
            }
            return *this;
        }

        //! Copy constructor
        strtree(const strtree& p)
        {
            copy(p);
        }

        //! Move constructor
        strtree(strtree&& snapshot)
        {
            *this = std::move(snapshot);
        }

        //! Swap method
        /*! Swaps the content of the two data structure.
         *  You have to use set_vector to adjust the supported bit_vector.
         *  \param bp_support Object which is swapped.
         */
        void swap(strtree& p)
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
