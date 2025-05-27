#pragma once

// Preprocess function
static void prepare_seed_edges_IM(std::string infilename, TArgument Arg)
{
    Graph graph;
    GraphBase::LoadGraph(graph, infilename);
    int probDist = GraphBase::LoadGraphProbDist(infilename);
    TResult tRes;
    Graph placeHold;
    TAlg tAlg(graph, tRes, placeHold);
    tAlg.set_vanilla_sample(Arg._vanilla);
    tAlg.set_prob_dist((ProbDist)probDist); // Set propagation model
    auto delta = Arg._delta;
    if (delta < 0) delta = 1.0 / graph.size();
    int seedSize = Arg._seedsize;
    // Generating seeds
    LogInfo("Generating seeds");
    if (!Arg._hist)
    {
        tAlg.subsimOnly(seedSize, Arg._eps, delta);
    }
    else
    {
        std::cout <<"HIST is invoked." <<std::endl;
        if (seedSize < 10)
        {
            tAlg.subsimWithTrunc(seedSize, Arg._eps, delta);
        }
        else
        {
            tAlg.subsimWithHIST(seedSize, Arg._eps, delta);
        }
    }
    tAlg.RefreshHypergraph();

    std::string outfilename = Arg._graphname + "_num" + std::to_string(seedSize) + "_IM";
    TIO::WriteOrderSeeds(outfilename, tRes, Arg._dir);

    auto vecSeeds = tRes.get_seed_vec();
    
    // candidate edges
    std::vector<double> vecWeightedIndegree(graph.size(), 0.0);
    std::vector<double> vecWeightedOutdegree(graph.size(), 0.0);
    std::vector<uint32_t> vecOutDegree(graph.size(), 0);
    // compute weighted in-degree and out-degree
    for (size_t i = 0; i < graph.size(); i++) {
        for (auto& nbr : graph[i]) {
            auto u = nbr.first;
            auto w = nbr.second;
            vecWeightedIndegree[i] += w;    // since the graph is reversed
            vecWeightedOutdegree[u] += w;
            vecOutDegree[u]++;
        }
    }
    Graph candEdges((uint32_t)graph.size());
    std::vector<bool> vecboolSeed(graph.size(), false);
    // init vecboolSeed
    for (auto& seed : vecSeeds) {
        vecboolSeed[seed] = true;
    }
    uint32_t edge_cnt = 0;
    for (size_t i = 0; i < graph.size(); ++i)
    {
        std::unordered_map<uint32_t, bool> vecboolVst4Seed;
        for (auto& nbr : graph[i]) {
            auto u = nbr.first;
            // if u is a seed node
            if (vecboolSeed[u]) {
                vecboolVst4Seed[u] = true;
            }
        }
        // for every seed
        for (auto& seed : vecSeeds) {
            if (vecboolVst4Seed[seed] || seed == i) {
                continue;
            }
            uint32_t targ_indeg = graph[i].size(), st_outdeg = vecOutDegree[seed];
            double avg_w_in, avg_w_out;
            
            if (targ_indeg == 0)
                avg_w_in = dsfmt_gv_genrand_open_close();
            if (st_outdeg == 0) avg_w_out = dsfmt_gv_genrand_open_close();
            if (targ_indeg != 0 && st_outdeg != 0) {
                avg_w_in = vecWeightedIndegree[i] / targ_indeg;
                avg_w_out = vecWeightedOutdegree[seed] / st_outdeg;
            }
            double w = (avg_w_in + avg_w_out) / 2;
            // w = 0.6;
            if (w > 1.0) {
                ExitMessage("probability over 1, please check");
            }
            candEdges[i].push_back(std::make_pair(seed, w));
            edge_cnt++;
        }
        // sort candEdges[i]
        std::sort(candEdges[i].begin(), candEdges[i].end(), smaller_first);
    }
    LogInfo("Total number of candidate edges", edge_cnt);
    LogInfo("--- Saving the candidate edges ---");
    std::string candedges_filepath = Arg._dir + "/candEdges_" + outfilename;
    TIO::SaveGraphStruct(candedges_filepath, candEdges, false);
    return;
}

// Preprocess to generate random seeds and prepare cand edges
static void prepare_seed_edges_RAND(std::string infilename, TArgument Arg) {
    Graph graph;
    GraphBase::LoadGraph(graph, infilename);
    int probDist = GraphBase::LoadGraphProbDist(infilename);
    TResult tRes;
    Graph placeHold;
    auto delta = Arg._delta;
    if (delta < 0) delta = 1.0 / graph.size();
    int seedSize = Arg._seedsize;
    // Generating seeds randomly
    std::vector<bool> vecboolSeed(graph.size(), false);
    std::vector<uint32_t> vecSeeds;
    while (vecSeeds.size() < seedSize) {
        uint32_t seed = dsfmt_gv_genrand_uint32() % graph.size();
        if (!vecboolSeed[seed]) {
            vecboolSeed[seed] = true;
            vecSeeds.push_back(seed);
        }
    }

    std::string outfilename = Arg._graphname + "_num" + std::to_string(seedSize) + "_RAND";
    std::vector<double> vecProb(seedSize, 0.0);
    // Assign the probability of each seed randomly
    for (int i = 0; i < seedSize; i++) {
        vecProb[i] = dsfmt_gv_genrand_open_close();
    }
    tRes.set_seed_vec(vecSeeds);
    TIO::WriteOrderSeeds(outfilename, tRes, Arg._dir);
    
    // candidate edges
    std::vector<double> vecWeightedIndegree(graph.size(), 0.0);
    std::vector<double> vecWeightedOutdegree(graph.size(), 0.0);
    std::vector<uint32_t> vecOutDegree(graph.size(), 0);
    // compute weighted in-degree and out-degree
    for (size_t i = 0; i < graph.size(); i++) {
        for (auto& nbr : graph[i]) {
            auto u = nbr.first;
            auto w = nbr.second;
            vecWeightedIndegree[i] += w;    // since the graph is reversed
            vecWeightedOutdegree[u] += w;
            vecOutDegree[u]++;
        }
    }
    Graph candEdges((uint32_t)graph.size());

    uint32_t edge_cnt = 0;
    for (size_t i = 0; i < graph.size(); ++i)
    {
        std::unordered_map<uint32_t, bool> vecboolVst4Seed;
        for (auto& nbr : graph[i]) {
            auto u = nbr.first;
            // if u is a seed node
            if (vecboolSeed[u]) {
                vecboolVst4Seed[u] = true;
            }
        }
        // for every seed
        for (auto& seed : vecSeeds) {
            if (vecboolVst4Seed[seed] || seed == i) {
                continue;
            }
            uint32_t targ_indeg = graph[i].size(), st_outdeg = vecOutDegree[seed];
            double avg_w_in, avg_w_out;
            
            if (targ_indeg == 0)
                avg_w_in = dsfmt_gv_genrand_open_close();
            if (st_outdeg == 0) avg_w_out = dsfmt_gv_genrand_open_close();
            if (targ_indeg != 0 && st_outdeg != 0) {
                avg_w_in = vecWeightedIndegree[i] / targ_indeg;
                avg_w_out = vecWeightedOutdegree[seed] / st_outdeg;
            }
            double w = (avg_w_in + avg_w_out) / 2;
            // w = 0.6;
            if (w > 1.0) {
                ExitMessage("probability over 1, please check");
            }
            candEdges[i].push_back(std::make_pair(seed, w));
            edge_cnt++;
        }
        // sort candEdges[i]
        std::sort(candEdges[i].begin(), candEdges[i].end(), smaller_first);
    }
    LogInfo("Total number of candidate edges", edge_cnt);
    LogInfo("--- Saving the candidate edges ---");
    std::string candedges_filepath = Arg._dir + "/candEdges_" + outfilename;
    TIO::SaveGraphStruct(candedges_filepath, candEdges, false);
    return;
}

static void preprocess(std::string infilename, TArgument Arg) {
    GraphBase::FormatGraph(infilename, Arg._probDist, Arg._wcVar, Arg._probEdge, Arg._skewType);
    prepare_seed_edges_IM(infilename, Arg);
    prepare_seed_edges_RAND(infilename, Arg);
}

// prepare candidate edges for different seeds
static void prepare_cand_edges(const Graph& g, const std::string seed_filepath, TArgument Arg) {
    std::vector<double> vecWeightedIndegree(g.size(), 0.0);
    std::vector<double> vecWeightedOutdegree(g.size(), 0.0);
    std::vector<uint32_t> vecOutDegree(g.size(), 0);
    std::string outfilename = Arg._graphname + "_num" + std::to_string(Arg._seedsize) + "_" + Arg._seedMode;
    Nodelist seeds;
    std::vector<double> vecProb;
    TIO::ReadProbSeeds(seed_filepath, seeds, vecProb);
    // compute weighted in-degree and out-degree
    for (size_t i = 0; i < g.size(); i++) {
        for (auto& nbr : g[i]) {
            auto u = nbr.first;
            auto w = nbr.second;
            vecWeightedIndegree[i] += w;    // since the graph is reversed
            vecWeightedOutdegree[u] += w;
            vecOutDegree[u]++;
        }
    }
    Graph candEdges((uint32_t)g.size());
    std::vector<bool> vecboolSeed(g.size(), false);
    // init vecboolSeed
    for (auto& seed : seeds) {
        vecboolSeed[seed] = true;
    }
    uint32_t edge_cnt = 0;
    for (size_t i = 0; i < g.size(); ++i)
    {
        std::unordered_map<uint32_t, bool> vecboolVst4Seed;
        for (auto& nbr : g[i]) {
            auto u = nbr.first;
            // if u is a seed node
            if (vecboolSeed[u]) {
                vecboolVst4Seed[u] = true;
            }
        }
        // for every seed
        for (auto& seed : seeds) {
            if (vecboolVst4Seed[seed] || seed == i) {
                continue;
            }
            uint32_t targ_indeg = g[i].size(), st_outdeg = vecOutDegree[seed];
            double avg_w_in, avg_w_out;
            
            if (targ_indeg == 0)
                avg_w_in = dsfmt_gv_genrand_open_close();
            if (st_outdeg == 0) avg_w_out = dsfmt_gv_genrand_open_close();
            if (targ_indeg != 0 && st_outdeg != 0) {
                avg_w_in = vecWeightedIndegree[i] / targ_indeg;
                avg_w_out = vecWeightedOutdegree[seed] / st_outdeg;
            }
            double w = (avg_w_in + avg_w_out) / 2;
            if (w > 1.0) {
                ExitMessage("probability over 1, please check");
            }
            candEdges[i].push_back(std::make_pair(seed, w));
            edge_cnt++;
        }
        std::sort(candEdges[i].begin(), candEdges[i].end(), smaller_first);
    }
    LogInfo("Total number of candidate edges", edge_cnt);
    LogInfo("--- Saving the candidate edges ---");
    std::string candedges_filepath = Arg._dir + "/candEdges_" + outfilename;
    TIO::SaveGraphStruct(candedges_filepath, candEdges, false);
    return;
}