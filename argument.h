#pragma once


class Argument
{
public:
    // Function parameter.
    // format: format graph
    // im: influence maximization
    std::string _funcStr = "im";
    FuncType _func = IM;

    // The number of nodes to be selected. Default is 50.
    int _seedsize = 50;
    int _edgesize = 100;

    // For the uniform setting, every edge has the same diffusion probability.
    float _probEdge = float(0.1);

    // Error threshold 1-1/e-epsilon.
    double _eps = 0.05;

    // Failure probability delta. Default is 1/#nodes.
    double _delta = -1.0;

    // Estimation error
    double _gamma = 0.05;

    // Graph name. Default is "facebook".
    std::string _graphname = "facebook";

    std::string _candEdgesProb = "orig";

    // Probability distribution
    // weights: graph data with weights
    // wc: wc setting
    // uniform: uniform setting
    // skewed: skewed distribution
    std::string _probDistStr = "wc";
    ProbDist _probDist = WC;
    std::string _skewType = "exp";

    std::string _seedDist = "uniRAND";

    // seed mode
    std::string _seedMode = "IM";

    // Directory
    std::string _dir = "graphInfo";

    // Result folder
    std::string _resultFolder = "result";

    // File name of the result
    std::string _outFileName;

    std::string _candedges_filename, _probseeds_filename;

    // wc variant
    double _wcVar = 1.0;

    // sample RR set with the vanilla method
    bool _vanilla = false;

    // use hist algorithm
    bool _hist = false;
    //Random seed
    int _rand_seed = 0;
    // Number of samples for greedy
    size_t _num_samples = 1000;
    // For evaluation
    std::string _method = "";

    Argument(int argc, char* argv[])
    {
        std::string param, value;

        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') break;

            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');

            if (!param.compare("-func")) _funcStr = value;
            else if (!param.compare("-seedsize")) _seedsize = stoi(value);
            else if (!param.compare("-edgesize")) _edgesize = stoi(value);
            else if (!param.compare("-num_samples")) _num_samples = stoi(value);
            else if (!param.compare("-eps")) _eps = stod(value);
            else if (!param.compare("-delta")) _delta = stod(value);
            else if (!param.compare("-gamma")) _gamma = stod(value);
            else if (!param.compare("-gname")) _graphname = value;
            else if (!param.compare("-dir")) _dir = value;
            else if (!param.compare("-outpath")) _resultFolder = value;
            else if (!param.compare("-pdist")) _probDistStr = value;
            else if (!param.compare("-pedge")) _probEdge = stof(value);
            else if (!param.compare("-wcvariant")) _wcVar = stod(value);
            else if (!param.compare("-skew")) _skewType = value;
            else if (!param.compare("-vanilla")) _vanilla = (value == "1");
            else if (!param.compare("-hist")) _hist = (value == "1");
            else if (!param.compare("-rand_seed")) _rand_seed = stoi(value);
            else if (!param.compare("-method")) _method = value;
            else if (!param.compare("-seed_dist")) _seedDist = value;
            else if (!param.compare("-seed_mode")) _seedMode = value;
            else if (!param.compare("-candEdges_prob")) _candEdgesProb = value;
        }

        if (_wcVar <= 0)
        {
            //wrong input
            _wcVar = 1.0;
        }

        decode_func_type();
        decode_prob_dist();
    }



    void build_outfilename(int seedSize, ProbDist dist, Graph& graph, FuncType func=IM, std::string method="SANDWICH")
    {
        std::string distStr; 

        if (dist == WEIGHTS)
        {
            _probDistStr = "weights";
        }
        else if (dist == WC)
        {
            _probDistStr = "wc";
        }
        else if (dist == UNIFORM)
        {
            _probDistStr = "uniform";

            for (int i = 0; i < graph.size(); i++)
            {
                if (graph[i].size() > 0 )
                {
                    _probEdge = graph[i][0].second;
                    break;
                }
            }
        }
        else
        {
            _probDistStr = "skewed";
        }
        if (_func == IM) {
            _outFileName = TIO::BuildOutFileName(_graphname, "subsim", seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps);
        }
        else if (func == SANDWICH) {
            _outFileName = TIO::BuildOutFileName(_graphname, "SANDWICH", seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps, this->_gamma, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == SANDWICHP) {
            _outFileName = TIO::BuildOutFileName(_graphname, "SANDWICHP", seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps, this->_gamma, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == RAND) {
            _outFileName = TIO::BuildOutFileName(_graphname, "RAND", seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == OUTDEG) {
            _outFileName = TIO::BuildOutFileName(_graphname, "OUTDEG", seedSize, _probDistStr, _probEdge, 0, 0.0, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == PROB) {
            _outFileName = TIO::BuildOutFileName(_graphname, "PROB", seedSize, _probDistStr, _probEdge, 0, 0.0, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == GREEDY) {
            _outFileName = TIO::BuildOutFileName(_graphname, "GREEDY", seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == SUBSIM)
        {
            _outFileName = TIO::BuildOutFileName(_graphname, "SUBSIM", seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == AIS)
        {
            _outFileName = TIO::BuildOutFileName(_graphname, "AIS", seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0, _gamma, _seedsize, _seedDist, _seedMode, _candEdgesProb);
        }
        else if (func == EVAL)
        {
            if(this->_method == "subsim" || this->_method == "SUBSIM") {
                _outFileName = TIO::BuildOutFileName(_graphname, this->_method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
            }
            else if (this->_method == "RAND" || this->_method == "GREEDY") {
                _outFileName = TIO::BuildOutFileName(_graphname, this->_method, seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
            }
            else if (this->_method == "OUTDEG" || this->_method == "PROB")
            {
                _outFileName = TIO::BuildOutFileName(_graphname, this->_method, seedSize, _probDistStr, _probEdge, 0, 0.0, 0.0, _seedsize, _seedDist, _seedMode, _candEdgesProb);
            }
            else if (this->_method == "AIS") {
                _outFileName = TIO::BuildOutFileName(_graphname, this->_method, seedSize, _probDistStr, _probEdge, this->_rand_seed, 0.0, this->_gamma, _seedsize, _seedDist, _seedMode, _candEdgesProb);
            }
            else if (this->_method == "SANDWICH" || this->_method == "SANDWICHP")
            {
                _outFileName = TIO::BuildOutFileName(_graphname, this->_method, seedSize, _probDistStr, _probEdge, this->_rand_seed, this->_eps, this->_gamma, _seedsize, _seedDist, _seedMode, _candEdgesProb);
            }
        }

        return ;
    }

    // Fill candedges and probseeds filenames from Arg
    void build_cand_probseeds_filenames() {
        this->_candedges_filename = this->_dir + "/" + "candEdges_" + this->_graphname + "_num" + std::to_string(this->_seedsize) + "_" + this->_seedMode;
        this->_probseeds_filename = this->_dir + "/" + "seed_" + this->_graphname + "_num" + std::to_string(this->_seedsize) + "_" + this->_seedMode + "_" + this->_seedDist;
    }

    void decode_prob_dist()
    {
        if (_probDistStr == "wc")
        {
            _probDist = WC;
        }
        else if (_probDistStr == "uniform")
        {
            _probDist = UNIFORM;
        }
        else if (_probDistStr == "skewed")
        {
            _probDist = SKEWED;
        }
        else if (_probDistStr == "weights")
        {
            _probDist = WEIGHTS;
        }
        else 
        {
            _probDist = PROB_DIST_ERROR;
        }
    }

    void decode_func_type()
    {
        if (_funcStr == "format")
        {
            _func = FORMAT;
        }
        else if (_funcStr == "im")
        {
            _func = IM;
        }
        else if (_funcStr == "SUBSIM")
        {
            _func = SUBSIM;
        }
        else if (_funcStr == "OUTDEG")
        {
            _func = OUTDEG;
        }
        // PROB
        else if (_funcStr == "PROB")
        {
            _func = PROB;
        }
        // RAND
        else if (_funcStr == "RAND")
        {
            _func = RAND;
        }
        // AIS
        else if (_funcStr == "AIS")
        {
            _func = AIS;
        }
        // SANDWICH
        else if (_funcStr == "SANDWICH")
        {
            _func = SANDWICH;
        }
        else if (_funcStr == "SANDWICHP")
        {
            _func = SANDWICHP;
        }
        // GREEDY
        else if (_funcStr == "GREEDY")
        {
            _func = GREEDY;
        }
        else if (_funcStr == "EVAL")
        {
            _func = EVAL;
        }
        else if (_funcStr == "PREP_CAND")
        {
            _func = PREP_CAND;
        }
        else
        {
            _func = FUNC_ERROR;
        }
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;
