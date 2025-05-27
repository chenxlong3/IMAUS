#pragma once

class Alg
{
private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE;
    /// _numRRsets: number of RR sets.
    size_t _numCandEdges = 0;
    size_t _numRRsets = 0;
    /// Upper bound in the last round for __mode=1.
    double _boundLast = DBL_MAX;
    /// The minimum upper bound among all rounds for __model=2.
    double _boundMin = DBL_MAX;
    /// Two hyper-graphs, one is used for selecting seeds and the other is used for validating influence.
    THyperGraph _hyperGraph, _hyperGraphVldt;
    /// Result object.
    TResult& _res;
    /// Seed set.
    Nodelist _vecSeed;
    std::vector<double> _vecProb;
    std::unordered_map<uint32_t, double> _mapSeedProb;
    std::vector<bool> _vecBoolSeed;
    std::vector<UVWEdge> _vecResEdges;
    std::vector<UVWEdge> _vecResEdgesLB;
    std::vector<UVWEdge> _vecResEdgesUB;

    std::vector<std::pair<uint32_t, uint32_t>> _vecIdx4Edges;
    std::vector<size_t> _vecStIdx4Edges;

    // Two-hop Inf for seeds
    std::unordered_map<uint32_t, double> _mapSeed2HopInf;

    Graph& _candEdges;
    ProbDist _probDist = WC;

    double _baseNumRRsets = 0.0;

    std::vector<uint32_t> _vecOutDegree;
    std::vector<uint32_t> _vecVldtInf;

    /// Maximum coverage by lazy updating.
    double MaxCoverVanilla(const int targetSize);
    double MaxCoverOutDegPrority(const int targetSize);
    double MaxCoverIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize);
    double MaxCoverSentinelSet(const int targetSize, const int totalTargetSize);

    /// Maximum coverage by maintaining the top-k marginal coverage.
    double MaxCoverTopK(const int targetSize);
    /// Maximum coverage.
    double MaxCover(const int targetSize);

    // Edge selection PMC
    double EdgeSelectionPMC(const int targetSize, std::string mode);
    double EdgeSelectionLB(const int targetSize);

public:
    Alg(Graph& graph, TResult& tRes, Graph& candEdges) : _hyperGraph(graph), _hyperGraphVldt(graph), _res(tRes), _candEdges(candEdges)
    {
        _numV = _hyperGraph.get_nodes();
        _numE = _hyperGraph.get_edges();
        _vecOutDegree = std::vector<uint32_t>(_numV);
        _vecBoolSeed = std::vector<bool>(_numV, false);
        _vecStIdx4Edges = std::vector<size_t>(_numV);
        for (auto &nbrs : graph)
        {
            for (auto &node : nbrs)
            {
                _vecOutDegree[node.first]++;
            }
        }
        if (_candEdges.size() > 0) {
            for (int i=0; i<this->_numV; i++) {
                this->_vecStIdx4Edges[i] = this->_numCandEdges;
                _numCandEdges += this->_candEdges[i].size();
                for (int j=0; j<this->_candEdges[i].size(); j++) {
                    this->_vecIdx4Edges.push_back(std::make_pair(i, j));
                }
            }
        }
        
    }
    ~Alg()
    {
    }
    
    /// Set cascade model.
    void set_prob_dist(const ProbDist weight);
    void set_vanilla_sample(const bool isVanilla);
    void set_probseed(std::string seedfile);
    void set_cand_edges_prob(std::string candEdgesProb);

    void RefreshHypergraph()
    {
        _hyperGraph.RefreshHypergraph();
        _hyperGraphVldt.RefreshHypergraph();
    }
    /// Evaluate influence spread for the seed set constructed
    double EfficInfVldtAlg();
    /// Evaluate influence spread for a given seed set
    double EfficInfVldtAlg(const Nodelist vecSeed);


    double estimateRRSize();

    double subsimOnly(const int targetSize, const double epsilon, const double delta);
    double subsimWithTrunc(const int targetSize, const double epsilon, const double delta);
    double IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr);

    double FindFixSub(const int targetSize, const int totalTargetSize, const double epsilon, const double delta);
    double FindRemSet(const int targetSize, const double epsilon, const double targeEpsilon, const double delta);
    double FindDynamSub(const int totalTargetSize, const double epsilon, const double delta);

    double subsimWithHIST(const int targetSize, const double epsilon, const double delta);

    // sandwich
    double sandwich(const int targetSize, const double epsilon, const double delta, const double gamma);

    double boundMaximize(const int targetSize, const double epsilon, const double delta, const std::string mode);
    void InitIdx4Edges();
    // sandwich+
    void compTwoHopInf();
    double sandwichP(const int targetSize, const double epsilon, const double delta, const double gamma);
    
    // Estimate influence increment
    double estimateInfIncforSand(std::vector<UVWEdge>& edges, const double epsilon, const double delta);
    // EVAL
    double estimateInfStoppingRules(double eps, double delta);
    // RAND
    void selectRandom(const int targetSize);
    // OUTDEG
    void selectOutDeg(const int targetSize);
    // PROB
    void selectProb(const int targetSize);
    // Greedy
    double inf_estimation(const size_t num_samples);
    void selectGreedy(const int targetSize, const size_t num_samples);
    void selectGreedyNode(const int targetSize, const size_t num_samples);
    // AIS
    void selectAIS(const int targetSize, const double gamma);
    // OPIM
    void selectIM(const int targetSize, double epsilon, double delta);
};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;
