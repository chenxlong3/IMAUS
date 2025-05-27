#pragma once

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#endif

class IoController
{
public:
    static void mkdir_absence(const char* outFolder)
    {
#if defined(_WIN32)
        CreateDirectoryA(outFolder, nullptr); // can be used on Windows
#else
        mkdir(outFolder, 0733); // can be used on non-Windows
#endif
    }

    /// Save a serialized file
    template <class T>
    static void SaveFile(const std::string filename, const T& output)
    {
        std::ofstream outfile(filename, std::ios::binary);

        if (!outfile.eof() && !outfile.fail())
        {
            StreamType res;
            serialize(output, res);
            outfile.write(reinterpret_cast<char*>(&res[0]), res.size());
            outfile.close();
            res.clear();
            std::cout << "Save file successfully: " << filename << '\n';
        }
        else
        {
            std::cout << "Save file failed: " + filename << '\n';
            exit(1);
        }
    }

    /// Load a serialized file
    template <class T>
    static void LoadFile(const std::string filename, T& input)
    {
        std::ifstream infile(filename, std::ios::binary);

        if (!infile.eof() && !infile.fail())
        {
            infile.seekg(0, std::ios_base::end);
            const std::streampos fileSize = infile.tellg();
            infile.seekg(0, std::ios_base::beg);
            std::vector<uint8_t> res(fileSize);
            infile.read(reinterpret_cast<char*>(&res[0]), fileSize);
            infile.close();
            input.clear();
            auto it = res.cbegin();
            input = deserialize<T>(it, res.cend());
            res.clear();
        }
        else
        {
            std::cout << "Cannot open file: " + filename << '\n';
            exit(1);
        }
    }

    /// Save graph structure to a file
    static void SaveGraphStruct(const std::string graphName, const Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";

        if (isReverse) postfix = ".vec.rvs.graph";

        const std::string filename = graphName + postfix;
        SaveFile(filename, vecGraph);
    }

    /// Load graph structure from a file
    static void LoadGraphStruct(const std::string graphName, Graph& vecGraph, const bool isReverse)
    {
        std::string postfix = ".vec.graph";

        if (isReverse) postfix = ".vec.rvs.graph";

        const std::string filename = graphName + postfix;
        LoadFile(filename, vecGraph);
    }

    static void SaveGraphProbDist(const std::string graphName, int dist)
    {
        std::ofstream outFile(graphName + ".probdist");
        outFile<< dist;
    }


    static int LoadGraphProbDist(const std::string graphName)
    {
        std::string filename = graphName + ".probdist";
        std::ifstream infile(filename);
        int probDist = WEIGHTS;

        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return probDist;
        }

        infile >> probDist;
        infile.close();
        std::cout << "probability distribution: " << probDist << std::endl;
        return probDist;
    }

    /// Get out-file name
    static std::string BuildOutFileName(const std::string graphName, const std::string algName, const int edgesize,
                                        const std::string probDist, const float probEdge, const int rand_seed = 0, const double eps = 0.1,
                                        const double gamma = 0.05,
                                        const int seedsize=50, const std::string seedDist="uni", const std::string seedMode="IM",
                                        const std::string cand_edges_prob="orig")
    {
        if (probDist == "uniform")
        {
            return graphName + "_" + std::to_string(rand_seed) + "_" + algName + "_k" + std::to_string(edgesize) + "_" + std::to_string(eps) + "_gamma" + std::to_string(gamma) + "_" + probDist + "_" + std::
                   to_string(probEdge) + "_seeds" + seedMode + std::to_string(seedsize) + "_" + seedDist + "_candEdges" + cand_edges_prob;
        }

        return graphName + "_" + std::to_string(rand_seed) + "_" + algName + "_k" + std::to_string(edgesize) + "_" + std::to_string(eps) + "_gamma" + std::to_string(gamma) + "_" + probDist + "_seeds" + seedMode + std::to_string(seedsize) + "_" + seedDist + "_candEdges" + cand_edges_prob;
    }

    /// Print the results
    static void WriteResult(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
    {
        const auto approx = resultObj.get_approximation();
        const auto approx_LB = resultObj.get_approximation_LB(); // approx_LB
        const auto approx_UB = resultObj.get_approximation_UB(); // approx_UB
        const auto runTime = resultObj.get_running_time();
        const auto samplingTime = resultObj.get_sampling_time();
        const auto selectionTime = resultObj.get_selection_time();
        const auto influence = resultObj.get_influence();
        const auto influenceOriginal = resultObj.get_influence_original();
        const auto seedSize = resultObj.get_seed_size();
        const auto RRsetsSize = resultObj.get_RRsets_size();
        const auto influence_LB = resultObj.get_influence_LB();
        const auto influence_UB = resultObj.get_influence_UB();
        const auto influence_original_LB = resultObj.get_influence_original_LB();
        const auto influence_original_UB = resultObj.get_influence_original_UB();
        const auto LB_value = resultObj.get_LB_value();
        const auto UB_value = resultObj.get_UB_value();
        const auto estimate_time = resultObj.get_estimate_time();
        const auto bound_max_time = resultObj.get_bound_max_time();

        std::cout << "   --------------------" << std::endl;
        std::cout << "  |Approx.: " << approx << std::endl;
        std::cout << "  |Time (sec): " << runTime << std::endl;
        std::cout << "  |Influence: " << influence << std::endl;
        std::cout << "  |Self-estimated influence: " << influenceOriginal << std::endl;
        std::cout << "  |#Seeds: " << seedSize << std::endl;
        std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
        std::cout << "   --------------------" << std::endl;
        mkdir_absence((outFolder + "/info/").c_str());
        std::ofstream outFileNew(outFolder + "/info/" + outFileName);

        if (outFileNew.is_open())
        {
            outFileNew << "Approx.: " << approx << std::endl;
            outFileNew << "Approx. LB: " << approx_LB << std::endl;
            outFileNew << "Approx. UB: " << approx_UB << std::endl;
            outFileNew << "Time (sec): " << runTime << std::endl;
            outFileNew << "Sampling Time (sec): " << samplingTime << std::endl;
            outFileNew << "Selection Time (sec): " << selectionTime << std::endl;
            outFileNew << "Estimate Time (sec): " << estimate_time << std::endl;
            outFileNew << "Bound Max Time (sec): " << bound_max_time << std::endl;
            outFileNew << "Influence_LB: " << influence_LB << std::endl;
            outFileNew << "Influence_UB: " << influence_UB << std::endl;
            outFileNew << "LB value: " << LB_value << std::endl;
            outFileNew << "UB value: " << UB_value << std::endl;
            outFileNew << "Self-estimated LB: " << influence_original_LB << std::endl;
            outFileNew << "Self-estimated UB: " << influence_original_UB << std::endl;
            outFileNew << "Self-estimated influence: " << influenceOriginal << std::endl;
            outFileNew << "#Seeds: " << seedSize << std::endl;
            outFileNew << "#RR sets: " << RRsetsSize << std::endl;
            outFileNew.close();
        }
        LogInfo("Finish writing the result to", outFolder + "/info/" + outFileName);
    }

    /// Print the seeds
    static void WriteOrderSeeds(const std::string& outFileName, const TResult& resultObj, const std::string& outFolder)
    {
        auto vecSeed = resultObj.get_seed_vec();
        mkdir_absence(outFolder.c_str());
        const auto outpath = outFolder + "/seed";
        mkdir_absence(outpath.c_str());
        std::ofstream outFile(outpath + "/seed_" + outFileName);

        for (auto i = 0; i < vecSeed.size(); i++)
        {
            outFile << vecSeed[i] << '\n';
        }

        outFile.close();
    }
    
    // Write Influence
    static void WriteInfluence(const std::string& outFileName, const double inf, const std::string& outFolder)
    {
        mkdir_absence(outFolder.c_str());
        std::ofstream outFile(outFolder + "/" + outFileName);
        outFile << inf << '\n';
        outFile.close();
        LogInfo("Finish writing the result to:", outFolder + "/" + outFileName);
    }

    /// Write Order Probabilistic Seeds
    static void WriteOrderProbSeeds(const std::string& outFileName, const std::string& outFolder, const TResult& resultObj, const std::vector<double>& vecProb)
    {
        auto vecSeed = resultObj.get_seed_vec();
        mkdir_absence(outFolder.c_str());
        std::ofstream outFile(outFolder + "/" + outFileName);

        for (auto i = 0; i < vecSeed.size(); i++) {
            outFile << vecSeed[i] << " " << vecProb[i] << '\n';
        }

        outFile.close();
        return;
    }

    // Read probabilistic seeds from file
    static void ReadProbSeeds(const std::string& outFilePath, std::vector<uint32_t>& vecSeed, std::vector<double>& vecProb)
    {
        uint32_t seed;
        double prob;
        std::ifstream inFile(outFilePath);

        while (inFile >> seed >> prob) {
            vecSeed.push_back(seed);
            vecProb.push_back(prob);
        }

        inFile.close();
        return;
    }

    static void write_UVWEdges(const std::string& outFileName, const std::string& outFolder, const TResult& resultObj) {
        const auto res = resultObj.get_res_edges();
        const auto outpath = outFolder + "/edges";
        const auto filename = outpath + "/edges_" + outFileName;
        mkdir_absence(outpath.c_str());
        std::ofstream outfile(filename);
        for (auto& tuple : res) {
            outfile << tuple.first << '\t' << tuple.second.first << '\t' << tuple.second.second << '\n';
        }
        outfile.close();
        LogInfo("Finish saving the edges to " + filename);
        return;
    }

    static void read_UVWEdges(const std::string& outFileName, std::vector<UVWEdge>& res, const std::string& outFolder) {
        uint32_t u, v;
        double w;
        const auto outpath = outFolder + "/edges";
        std::ifstream inFile(outpath + "/edges_" + outFileName);
        LogInfo("edges path", outpath + "/edges_" + outFileName);
        while (inFile >> u >> v >> w) {
            res.push_back(make_pair(u, Edge(v, w)));
        }
        inFile.close();
        return;
    }
    
};

using TIO = IoController;
using PIO = std::shared_ptr<IoController>;
