#pragma once
#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <climits>
#include <array>
#include <unordered_map>
#include <cstddef>

#include <omp.h>
#include <cstdlib>
#include <unistd.h>

#include "utils.hpp"

class Graph
{
    public:
        Graph(): nv_(-1), ne_(-1),
                 edge_indices_(nullptr), edge_list_(nullptr),
                 edge_active_(nullptr), mate_(nullptr), D_(nullptr), M_(nullptr) 
        {}
                
        Graph(GraphElem nv): 
            nv_(nv), ne_(-1), 
            edge_list_(nullptr), edge_active_(nullptr)
        {
            edge_indices_   = new GraphElem[nv_+1];
            M_              = new EdgeTuple[nv_];
            D_              = new GraphElem[nv_*2];
            std::fill(D_, D_ + nv_*2, -1);
#ifdef USE_OMP_OFFLOAD
#pragma omp target enter data map(to:this[:1])
#pragma omp target enter data map(alloc:edge_indices_[0:nv_+1])
#pragma omp target enter data map(alloc:D_[0:nv_*2])
#pragma omp target enter data map(alloc:M_[0:nv_])
#endif
        }

        Graph(GraphElem nv, GraphElem ne): 
            nv_(nv), ne_(ne) 
        {
            edge_indices_   = new GraphElem[nv_+1];
            edge_list_      = new Edge[ne_];
            edge_active_    = new EdgeActive[ne_];
            M_              = new EdgeTuple[nv_];
            D_              = new GraphElem[nv_*2];
#ifdef USE_OMP_OFFLOAD
#pragma omp target enter data map(to:this[:1])
#pragma omp target enter data map(alloc:edge_indices_[0:nv_+1])
#pragma omp target enter data map(alloc:edge_list_[0:ne_])
#pragma omp target enter data map(alloc:edge_active_[0:ne_])
#pragma omp target enter data map(alloc:D_[0:nv_*2])
#pragma omp target enter data map(alloc:M_[0:nv_])
#endif
        }

        ~Graph() 
        {
            delete [] edge_indices_;
            delete [] edge_list_;
            delete [] edge_active_;
            delete [] D_;
            delete [] M_;
        }
       
        /* 
        #ifdef USE_OMP_OFFLOAD                  
        Graph(const Graph &g) 
        {
            nv_ = g.nv_;
            ne_ = g.ne_;
            memcpy(edge_list_, g.edge_list_, sizeof(Edge)*ne_); 
            memcpy(edge_indices_, g.edge_indices_, sizeof(GraphElem)*(nv_+1)); 
            memcpy(edge_active_, g.edge_active_, sizeof(EdgeActive)*ne_); 
        }
        #endif
        */

        Graph(const Graph &other) = delete;
        Graph& operator=(const Graph& d) = delete;
 
        void set_edge_index(GraphElem const vertex, GraphElem const e0)
        {
#if defined(DEBUG_BUILD)
            assert((vertex >= 0) && (vertex <= nv_));
            assert((e0 >= 0) && (e0 <= ne_));
            edge_indices_.at(vertex) = e0;
#else
            edge_indices_[vertex] = e0;
#endif
        } 
        
        void edge_range(GraphElem const vertex, GraphElem& e0, 
                GraphElem& e1) const
        {
            e0 = edge_indices_[vertex];
            e1 = edge_indices_[vertex+1];
        } 

        void set_nedges(GraphElem ne) 
        { 
            ne_ = ne;
            edge_list_      = new Edge[ne_];
            edge_active_    = new EdgeActive[ne_];
#ifdef USE_OMP_OFFLOAD
#pragma omp target enter data map(alloc:edge_list_[0:ne_])
#pragma omp target enter data map(alloc:edge_active_[0:ne_])
#endif
        }

        GraphElem get_nv() const { return nv_; }
        GraphElem get_ne() const { return ne_; }
       
        // return edge and active info
        // ----------------------------
       
        Edge const& get_edge(GraphElem const index) const
        { return edge_list_[index]; }
         
        Edge& set_edge(GraphElem const index)
        { return edge_list_[index]; }       
                
        EdgeActive const& get_active_edge(GraphElem const index) const
        { return edge_active_[index]; }
         
        EdgeActive& get_active_edge(GraphElem const index)
        { return edge_active_[index]; }
 
        // print edge list (with weights)
        void print(bool print_weight = true) const
        {
            if (ne_ < MAX_PRINT_NEDGE)
            {
                for (GraphElem i = 0; i < nv_; i++)
                {
                    GraphElem e0, e1;
                    edge_range(i, e0, e1);
                    if (print_weight) { // print weights (default)
                        for (GraphElem e = e0; e < e1; e++)
                        {
                            Edge const& edge = get_edge(e);
                            std::cout << i << " " << edge.tail_ << " " << edge.weight_ << std::endl;
                        }
                    }
                    else { // don't print weights
                        for (GraphElem e = e0; e < e1; e++)
                        {
                            Edge const& edge = get_edge(e);
                            std::cout << i << " " << edge.tail_ << std::endl;
                        }
                    }
                }
            }
            else
            {
                std::cout << "Graph size is {" << nv_ << ", " << ne_ << 
                    "}, which will overwhelm STDOUT." << std::endl;
            }
        }

        void print_M() const
        {
            std::cout << "Matched vertices: " << std::endl;
            for (GraphElem i = 0; i < nv_; i++)
                std::cout << M_[i].ij_[0] << " ---- " << M_[i].ij_[1] << std::endl;
        }
         
        // if mate[mate[v]] == v then
        // we're good
        void check_results()
        {
            bool success = true;
            for (GraphElem i = 0; i < nv_; i++)
            {                
                if ((mate_[mate_[M_[i].ij_[0]]] != M_[i].ij_[0])
                        || (mate_[mate_[M_[i].ij_[1]]] != M_[i].ij_[1]))
                {
                    std::cout << "Validation FAILED." << std::endl; 
                    std::cout << "mate_[mate_[" << M_[i].ij_[0] << "]] != " << M_[i].ij_[0] << std::endl;
                    std::cout << "mate_[mate_[" << M_[i].ij_[1] << "]] != " << M_[i].ij_[1] << std::endl;
                    success = false;

                }
            }
            if (success)
                std::cout << "Validation SUCCESS." << std::endl;
        }
        
        // print statistics about edge distribution
        void print_stats()
        {
          std::vector<GraphElem> pdeg(nv_, 0);
          for (GraphElem v = 0; v < nv_; v++)
          {
            GraphElem e0, e1;
            edge_range(v, e0, e1);
            for (GraphElem e = e0; e < e1; e++)
              pdeg[v] += 1;
          }

          std::sort(pdeg.begin(), pdeg.end());
          GraphWeight loc = (GraphWeight)(nv_ + 1)/2.0;
          GraphElem median;
          if (fmod(loc, 1) != 0)
            median = pdeg[(GraphElem)loc]; 
          else
            median = (pdeg[(GraphElem)floor(loc)] + pdeg[((GraphElem)floor(loc)+1)]) / 2;
          GraphElem spdeg = std::accumulate(pdeg.begin(), pdeg.end(), 0);
          GraphElem mpdeg = *(std::max_element(pdeg.begin(), pdeg.end()));
          std::transform(pdeg.cbegin(), pdeg.cend(), pdeg.cbegin(),
              pdeg.begin(), std::multiplies<GraphElem>{});

          GraphElem psum_sq = std::accumulate(pdeg.begin(), pdeg.end(), 0);

          GraphWeight paverage = (GraphWeight) spdeg / nv_;
          GraphWeight pavg_sq  = (GraphWeight) psum_sq / nv_;
          GraphWeight pvar     = std::abs(pavg_sq - (paverage*paverage));
          GraphWeight pstddev  = sqrt(pvar);

          std::cout << std::endl;
          std::cout << "--------------------------------------" << std::endl;
          std::cout << "Graph characteristics" << std::endl;
          std::cout << "--------------------------------------" << std::endl;
          std::cout << "Number of vertices: " << nv_ << std::endl;
          std::cout << "Number of edges: " << ne_ << std::endl;
          std::cout << "Maximum number of edges: " << mpdeg << std::endl;
          std::cout << "Median number of edges: " << median << std::endl;
          std::cout << "Expected value of X^2: " << pavg_sq << std::endl;
          std::cout << "Variance: " << pvar << std::endl;
          std::cout << "Standard deviation: " << pstddev << std::endl;
          std::cout << "--------------------------------------" << std::endl;
        }
        
        // deactivate edge x -- y
        inline void deactivate_edge(GraphElem x, GraphElem y)
        {
            GraphElem e0, e1;
            edge_range(x, e0, e1);
            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = get_active_edge(e);
                if (edge.edge_->tail_ == y && edge.active_)
                {
                    edge.active_ = false;
                    break;
                }
            }
        }

        bool is_matched(const GraphElem v)
        {
            auto found = std::find_if(M_, M_ + nv_, 
                    [&](EdgeTuple const& et) 
                    { return ((et.ij_[0] == v) || (et.ij_[1] == v)); });
            if (found == M_ + nv_)
                return false;
            return true;
        }
        
        // maximal edge matching using OpenMP
        void maxematch()
        {
            Edge* max_edges = new Edge[nv_];
#ifdef USE_OMP_OFFLOAD
#pragma omp target enter data map(alloc:max_edges[0:nv_])
#pragma omp target teams distribute parallel for \
            map(always, tofrom:mate_[0:nv_])
#else
            #pragma omp parallel for default(shared) schedule(static)
#endif
            for (GraphElem i = 0; i < nv_; i++)
            {
                GraphElem e0, e1;
                edge_range(i, e0, e1);
                for (GraphElem e = e0; e < e1; e++)
                {
                    EdgeActive& edge = get_active_edge(e);
                    if (edge.active_)
                    {
                        if (edge.edge_->weight_ > max_edges[i].weight_)
                            max_edges[i] = *edge.edge_;
                        // break tie using vertex index
                        if (edge.edge_->weight_ == max_edges[i].weight_)
                        {
                            if (edge.edge_->tail_ > max_edges[i].tail_)
                                max_edges[i] = *edge.edge_;
                        }
                    }
                }
                mate_[i] = max_edges[i].tail_;
                const GraphElem y = mate_[i];
                // initiate matching request
                if (y != -1)
                {
                  GraphElem mate_y = mate_[y]; 
                  if (mate_y == i)
                  {
                    D_[i*2    ] = i;
                    D_[i*2 + 1] = y;
                    EdgeTuple et(i, y, max_edges[i].weight_); 
                    M_[i] = et;
                    // mark y<->i inactive, because its matched
                    deactivate_edge(y, i);
                    deactivate_edge(i, y);
                  }
                }
                else // invalidate all neigboring vertices 
                {
                    for (GraphElem e = e0; e < e1; e++)
                    {
                        EdgeActive& edge = get_active_edge(e);
                        if (edge.active_)
                        {
                          const GraphElem z = edge.edge_->tail_;
                          edge.active_ = false;
                          deactivate_edge(z, i); // invalidate x -- z
                        }
                    }
                }
            }
#ifdef USE_OMP_OFFLOAD
#pragma omp target exit data map(from:mate_[0:nv_])
#pragma omp target exit data map(from:D_[0:nv_*2])
#pragma omp target exit data map(from:M_[0:nv_])
#endif
        } 

        EdgeActive *edge_active_;
        GraphElem *edge_indices_;
        Edge *edge_list_;
   
        GraphElem get_num_vertices() { return nv_;};
        GraphElem get_num_edges() {return ne_;};
        GraphElem* get_index_ranges() {return edge_indices_;};
        void* get_edge_list() {return edge_list_;};

    private:
        GraphElem nv_, ne_;
        GraphElem* mate_;
        GraphElem* D_;
        EdgeTuple* M_;
};

// read in binary edge list files using POSIX I/O
class BinaryEdgeList
{
    public:
        BinaryEdgeList() : 
            M_(-1), N_(-1)
        {}
        
        // read a file and return a graph
        Graph* read(std::string binfile)
        {
            std::ifstream file;

            file.open(binfile.c_str(), std::ios::in | std::ios::binary); 

            if (!file.is_open()) 
            {
                std::cout << " Error opening file! " << std::endl;
                std::abort();
            }

            // read the dimensions 
            file.read(reinterpret_cast<char*>(&M_), sizeof(GraphElem));
            file.read(reinterpret_cast<char*>(&N_), sizeof(GraphElem));
#ifdef EDGE_AS_VERTEX_PAIR
            GraphElem weighted;
            file.read(reinterpret_cast<char*>(&weighted), sizeof(GraphElem));
            N_ *= 2;
#endif
            // create local graph
            Graph *g = new Graph(M_, N_);

            uint64_t tot_bytes=(M_+1)*sizeof(GraphElem);
            ptrdiff_t offset = 2*sizeof(GraphElem);

            if (tot_bytes < INT_MAX)
                file.read(reinterpret_cast<char*>(&g->edge_indices_[0]), tot_bytes);
            else 
            {
                int chunk_bytes=INT_MAX;
                uint8_t *curr_pointer = (uint8_t*) &g->edge_indices_[0];
                uint64_t transf_bytes = 0;

                while (transf_bytes < tot_bytes)
                {
                    file.read(reinterpret_cast<char*>(&curr_pointer[offset]), chunk_bytes);
                    transf_bytes += chunk_bytes;
                    offset += chunk_bytes;
                    curr_pointer += chunk_bytes;

                    if ((tot_bytes - transf_bytes) < INT_MAX)
                        chunk_bytes = tot_bytes - transf_bytes;
                } 
            }    

            N_ = g->edge_indices_[M_] - g->edge_indices_[0];
            g->set_nedges(N_);
            tot_bytes = N_*(sizeof(Edge));
            offset = 2*sizeof(GraphElem) + (M_+1)*sizeof(GraphElem) + g->edge_indices_[0]*(sizeof(Edge));

#if defined(GRAPH_FT_LOAD)
            ptrdiff_t currpos = file.tellg();
            ptrdiff_t idx = 0;
            GraphElem* vidx = (GraphElem*)malloc(M_ * sizeof(GraphElem));

            const int num_sockets = (GRAPH_FT_LOAD == 0) ? 1 : GRAPH_FT_LOAD;
            printf("Read file from %d sockets\n", num_sockets);
            int n_blocks = num_sockets;

            GraphElem NV_blk_sz = M_ / n_blocks;
            GraphElem tid_blk_sz = omp_get_num_threads() / n_blocks;

            #pragma omp parallel
            {
                for (int b=0; b<n_blocks; b++) 
                {

                    long NV_beg = b * NV_blk_sz;
                    long NV_end = std::min(M_, ((b+1) * NV_blk_sz) );
                    int tid_doit = b * tid_blk_sz;

                    if (omp_get_thread_num() == tid_doit) 
                    {
                        // for each vertex within block
                        for (GraphElem i = NV_beg; i < NV_end ; i++) 
                        {
                            // ensure first-touch allocation
                            // read and initialize using your code
                            vidx[i] = idx;
                            const GraphElem vcount = g->edge_indices_[i+1] - g->edge_indices_[i];
                            idx += vcount;
                            file.seekg(currpos + vidx[i] * sizeof(Edge), std::ios::beg);
                            file.read(reinterpret_cast<char*>(&g->edge_list_[vidx[i]]), sizeof(Edge) * (vcount));
                        }
                    }
                }
            }
            free(vidx);
#else
            if (tot_bytes < INT_MAX)
                file.read(reinterpret_cast<char*>(&g->edge_list_[0]), tot_bytes);
            else 
            {
                int chunk_bytes=INT_MAX;
                uint8_t *curr_pointer = (uint8_t*)&g->edge_list_[0];
                uint64_t transf_bytes = 0;

                while (transf_bytes < tot_bytes)
                {
                    file.read(reinterpret_cast<char*>(&curr_pointer[offset]), tot_bytes);
                    transf_bytes += chunk_bytes;
                    offset += chunk_bytes;
                    curr_pointer += chunk_bytes;

                    if ((tot_bytes - transf_bytes) < INT_MAX)
                        chunk_bytes = (tot_bytes - transf_bytes);
                } 
            }   
#endif

            file.close();

            for(GraphElem i=1;  i < M_+1; i++)
                g->edge_indices_[i] -= g->edge_indices_[0];   
            g->edge_indices_[0] = 0;
            
            // store reference to edge and active info
            // and ensure weights are positive
            for (GraphElem i=0; i < N_; i++)
                g->edge_active_[i].edge_ = &(g->edge_list_[i]);
            
            return g;
        }
    private:
        GraphElem M_, N_;
};

// RGG graph
class GenerateRGG
{
    public:
        GenerateRGG(GraphElem nv):
            nv_(nv), rn_(0)
        {
            // calculate r(n)
            GraphWeight rc = sqrt((GraphWeight)log(nv_)/(GraphWeight)(PI*nv_));
            GraphWeight rt = sqrt((GraphWeight)2.0736/(GraphWeight)nv_);
            rn_ = (rc + rt)/(GraphWeight)2.0;
            
            assert(((GraphWeight)1.0) > rn_);
        }

        Graph* generate(bool isLCG, bool unitEdgeWeight = true, GraphWeight randomEdgePercent = 0.0)
        {
            std::vector<GraphWeight> X, Y;

            X.resize(nv_);
            Y.resize(nv_);

            // create graph, edge list to be populated later
            Graph *g = new Graph(nv_);
            
            // measure the time to generate random numbers
            double st = omp_get_wtime();

            if (!isLCG) {
                // set seed (declared an extern in utils)
                seed = (unsigned)reseeder(1);

#if defined(PRINT_RANDOM_XY_COORD)
                #pragma omp parallel for
                for (GraphElem i = 0; i < nv_; i++) {
                    X[i] = genRandom<GraphWeight>(0.0, 1.0);
                    Y[i] = genRandom<GraphWeight>(0.0, 1.0);
                    std::cout << "X, Y: " << X[i] << ", " << Y[i] << std::endl;
                }
#else
                #pragma omp parallel for
                for (GraphElem i = 0; i < nv_; i++) {
                    X[i] = genRandom<GraphWeight>(0.0, 1.0);
                    Y[i] = genRandom<GraphWeight>(0.0, 1.0);
                }
#endif
            }
            else { // LCG
                // X | Y
                // e.g seeds: 1741, 3821
                // create LCG object
                // seed to generate x0
                LCG xr(/*seed*/1, X.data(), nv_); 
                
                // generate random numbers between 0-1
                xr.generate();

                // rescale xr further between lo-hi
                // and put the numbers in Y taking
                // from X[n]
                xr.rescale(Y.data(), nv_, 0);

#if defined(PRINT_RANDOM_XY_COORD)
                        for (GraphElem i = 0; i < nv_; i++) {
                            std::cout << "X, Y: " << X[i] << ", " << Y[i] << std::endl;
                        }
#endif
            }
                 
            double et = omp_get_wtime();
            double tt = et - st;
                
            std::cout << "Average time to generate " << nv_ 
                << " random numbers using LCG (in s): " 
                << tt << std::endl;

            // edges
            std::vector<EdgeTuple> edgeList;

#if defined(CHECK_NUM_EDGES)
            GraphElem numEdges = 0;
#endif
            for (GraphElem i = 0; i < nv_; i++) {
                for (GraphElem j = i + 1; j < nv_; j++) {
                    // euclidean distance:
                    // 2D: sqrt((px-qx)^2 + (py-qy)^2)
                    GraphWeight dx = X[i] - X[j];
                    GraphWeight dy = Y[i] - Y[j];
                    GraphWeight ed = sqrt(dx*dx + dy*dy);
                    // are the two vertices within the range?
                    if (ed <= rn_) {
                       if (!unitEdgeWeight) {
                            edgeList.emplace_back(i, j, ed);
                            edgeList.emplace_back(j, i, ed);
                        }
                        else {
                            edgeList.emplace_back(i, j);
                            edgeList.emplace_back(j, i);
                        }
#if defined(CHECK_NUM_EDGES)
                        numEdges += 2;
#endif

                        g->edge_indices_[i+1]++;
                        g->edge_indices_[j+1]++;
                    }
                }
            }

            // add random edges based on 
            // randomEdgePercent 
            if (randomEdgePercent > 0.0) {
                const GraphElem pnedges = (edgeList.size()/2);
                // extra #edges
                const GraphElem nrande = ((GraphElem)(randomEdgePercent * (GraphWeight)pnedges)/100);

#if defined(PRINT_EXTRA_NEDGES)
                int extraEdges = 0;
#endif

                unsigned rande_seed = (unsigned)(time(0)^getpid());
                GraphWeight weight = 1.0;
                std::hash<GraphElem> reh;

                // cannot use genRandom if it's already been seeded
                std::default_random_engine re(rande_seed); 
                std::uniform_int_distribution<GraphElem> IR, JR; 
                std::uniform_real_distribution<GraphWeight> IJW; 

                for (GraphElem k = 0; k < nrande; k++) {

                    // randomly pick start/end vertex and target from my list
                    const GraphElem i = (GraphElem)IR(re, std::uniform_int_distribution<GraphElem>::param_type{0, (nv_- 1)});
                    const GraphElem j = (GraphElem)JR(re, std::uniform_int_distribution<GraphElem>::param_type{0, (nv_- 1)});

                    if (i == j) 
                        continue;

                    // check for duplicates prior to edgeList insertion
                    auto found = std::find_if(edgeList.begin(), edgeList.end(), 
                            [&](EdgeTuple const& et) 
                            { return ((et.ij_[0] == i) && (et.ij_[1] == j)); });

                    // OK to insert, not in list
                    if (found == std::end(edgeList)) { 

                        // calculate weight
                        if (!unitEdgeWeight) {
                            GraphWeight dx = X[i] - X[j];
                            GraphWeight dy = Y[i] - Y[j];
                            weight = sqrt(dx*dx + dy*dy);
                        }

#if defined(PRINT_EXTRA_NEDGES)
                        extraEdges += 2;
#endif
#if defined(CHECK_NUM_EDGES)
                        numEdges += 2;
#endif                       
                        edgeList.emplace_back(i, j, weight);
                        edgeList.emplace_back(j, i, weight);
                        g->edge_indices_[i+1]++;
                        g->edge_indices_[j+1]++;
                    }
                }

#if defined(PRINT_EXTRA_NEDGES)
                std::cout << "Adding extra " << (extraEdges/2) << " edges while trying to incorporate " 
                    << randomEdgePercent << "%" << " extra edges globally." << std::endl;
#endif
            } // end of (conditional) random edges addition

            // set graph edge indices
            std::partial_sum(g->edge_indices_, g->edge_indices_ + (nv_+1), g->edge_indices_);
             
            for(GraphElem i = 1; i < nv_+1; i++)
                g->edge_indices_[i] -= g->edge_indices_[0];   
            g->edge_indices_[0] = 0;

            g->set_edge_index(0, 0);
            for (GraphElem i = 0; i < nv_; i++)
                g->set_edge_index(i+1, g->edge_indices_[i+1]);
            
            const GraphElem nedges = g->edge_indices_[nv_] - g->edge_indices_[0];
            g->set_nedges(nedges);
            
            // set graph edge list
            // sort edge list
            auto ecmp = [] (EdgeTuple const& e0, EdgeTuple const& e1)
            { return ((e0.ij_[0] < e1.ij_[0]) || ((e0.ij_[0] == e1.ij_[0]) && (e0.ij_[1] < e1.ij_[1]))); };

            if (!std::is_sorted(edgeList.begin(), edgeList.end(), ecmp)) {
#if defined(DEBUG_PRINTF)
                std::cout << "Edge list is not sorted." << std::endl;
#endif
                std::sort(edgeList.begin(), edgeList.end(), ecmp);
            }
#if defined(DEBUG_PRINTF)
            else
                std::cout << "Edge list is sorted!" << std::endl;
#endif
  
            GraphElem ePos = 0;
            for (GraphElem i = 0; i < nv_; i++) {
                GraphElem e0, e1;

                g->edge_range(i, e0, e1);
#if defined(DEBUG_PRINTF)
                if ((i % 100000) == 0)
                    std::cout << "Processing edges for vertex: " << i << ", range(" << e0 << ", " << e1 <<
                        ")" << std::endl;
#endif
                for (GraphElem j = e0; j < e1; j++) {
                    Edge &edge = g->set_edge(j);
                    g->edge_active_[j].edge_ = &edge;

                    assert(ePos == j);
                    assert(i == edgeList[ePos].ij_[0]);
                    
                    edge.tail_ = edgeList[ePos].ij_[1];
                    edge.weight_ = edgeList[ePos].w_;

                    ePos++;
                }
            }
                        
#if defined(CHECK_NUM_EDGES)
            const GraphElem ne = g->get_ne();
            assert(ne == numEdges);
#endif
            edgeList.clear();
            
            X.clear();
            Y.clear();
            
            return g;
        }

        GraphWeight get_d() const { return rn_; }
        GraphElem get_nv() const { return nv_; }

    private:
        GraphElem nv_;
        GraphWeight rn_;
};
#endif
