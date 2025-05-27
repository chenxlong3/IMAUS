#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"
#include "preprocessing.h"

void init_random_seed()
{
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
}

void init_random_seed(double rand_seed)
{
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(rand_seed);
}

int main(int argc, char* argv[])
{

    TArgument Arg(argc, argv);

    if (Arg._probDist == PROB_DIST_ERROR)
    {
        LogInfo("The input probability distribution is not supported:", Arg._probDistStr);
        LogInfo("The supported probability distribution: weights, wc, uniform, skewed");
        return 0;
    }

    if (Arg._func == FUNC_ERROR)
    {
        LogInfo("The input func is not supported: ", Arg._funcStr);
        LogInfo("The supported func: format, im");
    }

    init_random_seed(Arg._rand_seed);

    const std::string infilename = Arg._dir + "/" + Arg._graphname;
    if (Arg._func == FORMAT)
    {
        // Format the graph
        preprocess(infilename, Arg);
        return 0;
    }
    Timer mainTimer("main");
    // Load the reverse graph
    Graph graph;
    TResult tRes;
    Graph candEdges;
    
    GraphBase::LoadGraph(graph, infilename);
    int seedSize = Arg._seedsize;
    int edgeSize = Arg._edgesize;
    auto delta = Arg._delta;
    if (delta < 0) delta = 1.0 / graph.size();

    std::cout << "edgeSize k=" << edgeSize << std::endl;
    int probDist = GraphBase::LoadGraphProbDist(infilename);
    // if func == EVAL
    Arg.build_cand_probseeds_filenames();
    if (Arg._func == PREP_CAND) {
        Arg.build_cand_probseeds_filenames();
        prepare_cand_edges(graph, Arg._probseeds_filename, Arg);
        return 0;
    }
    // if (Arg._seedMode == "IM")
    // {
    //     candedges_filename = Arg._dir + "/candEdges_" + Arg._graphname;
    //     probseeds_filename = Arg._dir + "/probseeds_" + Arg._graphname;
    // }
    // else {
    //     candedges_filename = Arg._dir + "/candEdges_" + Arg._seedMode + "_" + Arg._graphname;
    //     probseeds_filename = Arg._dir + "/probseeds_" + Arg._seedMode + "_" + Arg._graphname;
    // }
    if (Arg._func == EVAL)
    {
        Arg._resultFolder = "./result/evaluation";
        Arg.build_outfilename(edgeSize, (ProbDist)probDist, graph, Arg._func);

        std::vector<UVWEdge> edges;
        TIO::read_UVWEdges(Arg._outFileName, edges, "./result");
        LogInfo("Number of edges", edges.size());
        AddEdges(graph, edges);
        TAlg tAlg(graph, tRes, candEdges);
        tAlg.set_vanilla_sample(1);
        tAlg.set_prob_dist((ProbDist)probDist); // Set propagation model
        tAlg.set_probseed(Arg._probseeds_filename);
        
        double res = tAlg.estimateInfStoppingRules(0.01, delta);
        TIO::WriteInfluence(Arg._outFileName, res, Arg._resultFolder);
        return 0;
    }
    
    
    TIO::LoadGraphStruct(Arg._candedges_filename, candEdges, false);
    // Initialize a result object to record the results
    TAlg tAlg(graph, tRes, candEdges);
    tAlg.set_vanilla_sample(Arg._vanilla);
    tAlg.set_prob_dist((ProbDist)probDist); // Set propagation model
    tAlg.set_probseed(Arg._probseeds_filename);
    
    Arg.build_outfilename(edgeSize, (ProbDist)probDist, graph, Arg._func);
    std::cout << "--- The Begin of " << Arg._outFileName << " ---\n";
    if (Arg._func == SANDWICH)
    {
        tAlg.sandwich(edgeSize, Arg._eps, delta, Arg._gamma);
    }
    // else if (Arg._func == SANDWICHP)
    // {
    //     tAlg.sandwichP(edgeSize, Arg._eps, delta, Arg._gamma);
    // }
    else if (Arg._func == RAND)
    {
        tAlg.selectRandom(edgeSize);
    }
    else if (Arg._func == OUTDEG) {
        tAlg.selectOutDeg(edgeSize);
    }
    else if (Arg._func == PROB) {
        tAlg.selectProb(edgeSize);
    }
    else if (Arg._func == GREEDY) {
        tAlg.selectGreedy(edgeSize, Arg._num_samples);
        // tAlg.selectGreedyNode(edgeSize, Arg._num_samples);
    }
    else if (Arg._func == AIS) {
        tAlg.selectAIS(edgeSize, Arg._gamma);
    }
    else if (Arg._func == SUBSIM) {
        tAlg.selectIM(edgeSize, Arg._eps, delta);
    }
    // else if (Arg._func == IM) {
    //     if (!Arg._hist)
    //     {
    //         tAlg.subsimOnly(seedSize, Arg._eps, delta);
    //     }
    //     else
    //     {
    //         std::cout <<"HIST is invoked." <<std::endl;
    //         if (seedSize < 10)
    //         {
    //             tAlg.subsimWithTrunc(seedSize, Arg._eps, delta);
    //         }
    //         else
    //         {
    //             tAlg.subsimWithHIST(seedSize, Arg._eps, delta);
    //         }
    //     }
    // }

    TIO::WriteResult(Arg._outFileName, tRes, Arg._resultFolder);
    TIO::write_UVWEdges(Arg._outFileName, Arg._resultFolder, tRes);
    std::cout << "---The End of " << Arg._outFileName << "---\n";
    tAlg.RefreshHypergraph();
    return 0;
}