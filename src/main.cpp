#define HEAD_INFO

#include "sfmt/SFMT.h"
#include "head.h"

class Argument{
public:
    int k;
    string dataset;
    double epsilon;
    string model;
    double T;
    int threshold;
    int rrNum;
};

#include "graph.h"
#include "infgraph.h"
#include "RRset.h"

/*
 * 1. run (imm.cpp): process parameters, call run_with_parameter
 * 2. run_with_parameter (imm.cpp) : call Imm::InfluenceMaximize
 * 3. Imm::InfluenceMaximize (imm.h): call infgraph :: g.build_hyper_graph_r and g.build_seed_set
 * 4. infgraph (infgraph.h) (derived class of graph.h) : call build hypernodes
 * 5. Build RRset (continuous_rrset.h and discrete_rrset.h):
 * #ifdef CONTINUOUS (CONTINUOUS is defined in Makefile)
            BuildHypergraphNode(random_number[i], i, arg );
            Implement in continuous_rrset.h
    #endif
    #ifdef DISCRETE (DISCRETE is defined in Makefile)
            BuildHypergraphNode(random_number[i], i );
            Implement in discrete_rrset.h
    #endif
 */


void run_with_parameter(InfGraph &g, const Argument & arg)
{
    cout << "--------------------------------------------------------------------------------" << endl;
    /*
     * R is the number of RRsets
     * threshold: only output to files if the node's influence is greater than the threshold
     */
    cout << arg.dataset << " k=" << arg.k << " epsilon=" << arg.epsilon <<   " " << arg.model <<" "<<arg.threshold <<endl;

    RRset::InfluenceMaximize(g, arg);
//    RRset::getRRsets(g, arg.rrNum, arg, arg.threshold);
//    RRset::getRRsets(g, arg, arg.threshold);

    //INFO(g.seedSet);
    //INFO(g.InfluenceHyperGraph());
    Timer::show();
}
void Run(int argn, char **argv)
{
    Argument arg;


    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-help") || argv[i] == string("--help") || argn == 1)
        {
            cout << "./tim -dataset *** -epsilon *** -k ***  -model IC|LT|TR|CONT " << endl;
            return ;
        }
        if (argv[i] == string("-dataset"))
            arg.dataset = argv[i + 1];
        if (argv[i] == string("-epsilon"))
            arg.epsilon = atof(argv[i + 1]);
        if (argv[i] == string("-T"))
            arg.T = atof(argv[i + 1]);
        if (argv[i] == string("-k"))
            arg.k = atoi(argv[i + 1]);
        if (argv[i] == string("-model"))
            arg.model = argv[i + 1];
        if (argv[i] == string("-R"))
            arg.rrNum = atoi(argv[i + 1]);
        if (argv[i] == string("-threshold"))
            arg.threshold = atoi(argv[i + 1]);
    }
    ASSERT(arg.dataset != "");
    ASSERT(arg.model == "IC" || arg.model == "LT" || arg.model == "TR" || arg.model=="CONT");

    string graph_file;
    if (arg.model == "IC")
        graph_file = arg.dataset + "graph_ic.inf";
    else if (arg.model == "LT")
        graph_file = arg.dataset + "graph_lt.inf";
    else if (arg.model == "TR")
        graph_file = arg.dataset + "graph_tr.inf";
    else if (arg.model == "CONT")
        graph_file = arg.dataset + "graph_cont.inf";
    else
    ASSERT(false);

    InfGraph g(arg.dataset, graph_file);


    if (arg.model == "IC")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "LT")
        g.setInfuModel(InfGraph::LT);
    else if (arg.model == "TR")
        g.setInfuModel(InfGraph::IC);
    else if (arg.model == "CONT")
        g.setInfuModel(InfGraph::CONT);
    else
    ASSERT(false);

    INFO(arg.T);

    run_with_parameter(g, arg);
}


int main(int argn, char **argv)
{
    __head_version = "v1";
    OutputInfo info(argn, argv);


    Run( argn, argv );
}


