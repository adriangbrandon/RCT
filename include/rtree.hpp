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


    public:
        typedef uint64_t size_type;
        typedef uint32_t value_type;

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


        inline SRTree::SNode* get_root() const {
            return m_tSRTree->root;
        }

        inline value_type get_lower0() const {
            return m_tSRTree->lower0;
        }

        inline value_type get_lower1() const {
            return m_tSRTree->lower1;
        }

        template<class Candidate, class Candidates>
        void enqueue_children(const Candidate &candidate, Candidates &candidates, const util::geo::point &p_q) const {
            std::vector<int32_t > lx, ly, hx, hy, ids;
            std::vector<SNode*> nodes;
            candidate.ptr->children_region(candidate.min.x, candidate.min.y, lx, ly, hx, hy, nodes, ids);
            //std::cout << "Children of: [(" << candidate.min.x << ", " << candidate.min.y << "), (" << candidate.max.x << ", " << candidate.max.y << ")]" << std::endl;
            if(ids.empty()){
                //Enqueue Regions
                for(size_type i = 0; i < nodes.size(); ++i){
                    //std::cout << "Enqueue Region: " << " [(" << lx[i] << ", " << ly[i] << "), (" << hx[i] << ", " << hy[i] << ")]" << std::endl;
                    candidates.push(Candidate(nodes[i], p_q, util::geo::point{(value_type) lx[i], (value_type) ly[i]}, util::geo::point{(value_type) hx[i], (value_type) hy[i]}));
                }
            }else{
                //Enqueue Objects
                for(size_type i = 0; i < ids.size(); ++i){
                    //std::cout << "Enqueue Object: " << ids[i] << " [(" << lx[i] << ", " << ly[i] << "), (" << hx[i] << ", " << hy[i] << ")]" << std::endl;
                    candidates.push(Candidate(ids[i], p_q, util::geo::point{(value_type) lx[i], (value_type) ly[i]}, util::geo::point{(value_type) hx[i], (value_type) hy[i]}));
                }
            }
        }


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
