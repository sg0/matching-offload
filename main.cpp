#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "graph.hpp"

static std::string inputFileName;
static GraphElem nvRGG = 0;
static int generateGraph = 0;

static GraphWeight randomEdgePercent = 0.0;
static bool randomNumberLCG = false;

static void parseCommandLine(const int argc, char * const argv[]);

int main(int argc, char *argv[])
{
    double t0, t1, td, td0, td1;
    int max_threads, req_threads;
    
    max_threads = omp_get_max_threads();
    req_threads = omp_get_num_threads();

    std::cout << "Maximum threads: " << max_threads << ", requested threads: " << req_threads << std::endl;
    
    // command line options
    parseCommandLine(argc, argv);
 
    Graph* g = nullptr;
    
    td0 = omp_get_wtime();

    // generate graph only supports RGG as of now
    if (generateGraph) 
    { 
        GenerateRGG gr(nvRGG);
        g = gr.generate(randomNumberLCG, isUnitEdgeWeight, randomEdgePercent);

        std::cout << "Generated Random Geometric Graph with d: " << gr.get_d() << std::endl;
        const GraphElem nv = g->get_nv();
        const GraphElem ne = g->get_ne();
        std::cout << "Number of vertices: " << nv << std::endl;
        std::cout << "Number of edges: " << ne << std::endl;
    }
    else // read input graph
    {
        BinaryEdgeList rm;
        if (readBalanced == true)
            g = rm.read_balanced(me, nprocs, ranksPerNode, inputFileName);
        else
            g = rm.read(me, nprocs, ranksPerNode, inputFileName);
    }
        
    g->print_dist_stats();
    assert(g != nullptr);

    td1 = omp_get_wtime();
    td = td1 - td0;

    if (!generateGraph)
        std::cout << "Time to read input file and create distributed graph (in s): " 
            << td << std::endl;
    else
        std::cout << "Time to generate distributed graph of " 
            << nvRGG << " vertices (in s): " << td << std::endl;

    // invoke matching
    t0 = omp_get_wtime();
    g->maxematch();
    t1 = omp_get_wtime();
    double p_tot = t1 - t0;

    std::cout << "Execution time (in s) for maximal edge matching: " 
        << t_tot << std::endl;

#if defined(CHECK_RESULTS)    
    g->check_results();
#endif
#if defined(PRINT_RESULTS)    
    g->print_M();
#endif
 
    return 0;
}

void parseCommandLine(int argc, char** const argv)
{
  int ret;
  optind = 1;
  bool help_text = false;

  if (argc == 1)
  {
      nvRGG = DEFAULT_NV;
      generateGraph = (nvRGG > 0)? true : false; 
  }
  else
  {
      while ((ret = getopt(argc, argv, "f:n:lp:h")) != -1) 
      {
          switch (ret) {
              case 'f':
                  inputFileName.assign(optarg);
                  break;
              case 'n':
                  nvRGG = atol(optarg);
                  if (nvRGG > 0)
                      generateGraph = true; 
                  break;
              case 'l':
                  randomNumberLCG = true;
                  break;
              case 'p':
                  randomEdgePercent = atof(optarg);
                  break;
              case 'h':
                  std::cout << "Set OMP_NUM_THREADS (max threads reported: " << omp_get_max_threads() << ") and affinity." << std::endl;
                  std::cout << "Usage [1] (use real-world file): ./neve_threads [-l] [-f /path/to/binary/file.bin] (see README)" << std::endl;
                  std::cout << "Usage [2] (use synthetic graph): ./neve_threads [-n <#vertices>] [-l] [-p <\% extra edges>]" << std::endl;
                  help_text = true;
                  break;
              default:
                  std::cout << "Please check the passed options." << std::endl;
                  break;
          }
      }
  }

  if (help_text)
      std::exit(EXIT_SUCCESS);

  if (!generateGraph && inputFileName.empty()) 
  {
      std::cerr << "Must specify a binary file name with -f or provide parameters for generating a graph." << std::endl;
      std::abort();
  }
   
  if (!generateGraph && randomNumberLCG) 
  {
      std::cerr << "Must specify -n <#vertices> for graph generation using LCG." << std::endl;
      std::abort();
  } 
   
  if (!generateGraph && (randomEdgePercent > 0.0)) 
  {
      std::cerr << "Must specify -n <#vertices> for graph generation first to add random edges to it." << std::endl;
      std::abort();
  } 
  
  if (generateGraph && ((randomEdgePercent < 0.0) || (randomEdgePercent >= 100.0))) 
  {
      std::cerr << "Invalid random edge percentage for generated graph!" << std::endl;
      std::abort();
  }
} // parseCommandLine
