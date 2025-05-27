#include "stdafx.h"

// code from OPIM
double Alg::MaxCoverVanilla(const int targetSize)
{
    // optimization with minimum upper bound among all rounds [Default].
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    // degMap: map degree to the nodes with this degree
    RRsets degMap(maxDeg + 1);

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                // Top-k influential nodes constructed
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// In the case of identical marginal，the node with largest out degree is chosen
double Alg::MaxCoverOutDegPrority(const int targetSize)
{
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

// max cover used in the IM-sentinel Phase
double Alg::MaxCoverIMSentinel(std::vector<uint32_t> &seedSet, const int targetSize)
{
    // seedSet: the sentinel set obtained in the Sentinel Set Selection Phase
    std::unordered_set<uint32_t> subSeedSet(seedSet.begin(), seedSet.end());
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto node : seedSet)
    {
        sumInf += coverage[node];
        _vecSeed.push_back(node);
        coverage[node] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[node])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                coverage[nodeIdx]--;
            }
        }
    }

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = targetSize;
                auto degBound = deg;
                FRset vecBound(targetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
                // std::cout << ">>>[greedy-lazy] influence: " << finalInf << ", min-bound: " << _boundMin <<
                //           ", last-bound: " << _boundLast << '\n';
                return finalInf;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverSentinelSet(const int targetSize, const int totalTargetSize)
{
    //targetSize: the size of the sentinel set
    //totalTargetSize: the total number of the seed set
    _boundLast = DBL_MAX, _boundMin = DBL_MAX;
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        if (coverage[i] == 0) continue;

        degMap[coverage[i]].push_back(i);
    }

    size_t sumInf = 0;
    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    _vecSeed.clear();

    for (auto deg = maxDeg; deg > 0; deg--) // Enusre deg > 0
    {
        auto& origVecNode = degMap[deg];
        std::vector<std::pair<uint32_t, uint32_t>> vecPair;

        for (auto idx = origVecNode.size(); idx--;)
        {
            auto node = origVecNode[idx];
            auto nodeCoverage = coverage[node];

            if (deg > nodeCoverage)
            {
                degMap[nodeCoverage].push_back(node);
                continue;
            }

            vecPair.push_back(std::make_pair(_vecOutDegree[node], node));
        }

        degMap.pop_back();

        if (vecPair.size() == 0)
        {
            continue;
        }

        /* sort nodes by their out-degre in ascending order */
        sort(vecPair.begin(), vecPair.end());
        std::vector<uint32_t> newVecNode;

        for (auto &nodePair : vecPair)
        {
            newVecNode.push_back(nodePair.second);
        }

        degMap.push_back(newVecNode);
        auto &vecNode = degMap[deg];

        for (auto idx = vecNode.size(); idx--;)
        {
            auto argmaxIdx = vecNode[idx];
            const auto currDeg = coverage[argmaxIdx];

            if (deg > currDeg)
            {
                degMap[currDeg].push_back(argmaxIdx);
                continue;
            }

            if (true)
            {
                // Find upper bound
                auto topk = totalTargetSize;
                auto degBound = deg;
                FRset vecBound(totalTargetSize);
                // Initialize vecBound
                auto idxBound = idx + 1;

                while (topk && idxBound--)
                {
                    vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                }

                while (topk && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (topk && idxBound--)
                    {
                        vecBound[--topk] = coverage[degMap[degBound][idxBound]];
                    }
                }

                MakeMinHeap(vecBound);
                // Find the top-k marginal coverage
                auto flag = topk == 0;

                while (flag && idxBound--)
                {
                    const auto currDegBound = coverage[degMap[degBound][idxBound]];

                    if (vecBound[0] >= degBound)
                    {
                        flag = false;
                    }
                    else if (vecBound[0] < currDegBound)
                    {
                        MinHeapReplaceMinValue(vecBound, currDegBound);
                    }
                }

                while (flag && --degBound)
                {
                    idxBound = degMap[degBound].size();

                    while (flag && idxBound--)
                    {
                        const auto currDegBound = coverage[degMap[degBound][idxBound]];

                        if (vecBound[0] >= degBound)
                        {
                            flag = false;
                        }
                        else if (vecBound[0] < currDegBound)
                        {
                            MinHeapReplaceMinValue(vecBound, currDegBound);
                        }
                    }
                }

                _boundLast = double(accumulate(vecBound.begin(), vecBound.end(), size_t(0)) + sumInf) * _numV / _numRRsets;

                if (_boundMin > _boundLast) _boundMin = _boundLast;
            }

            if (_vecSeed.size() >= targetSize)
            {
                goto afterGreedy;
            }

            sumInf += currDeg;
            _vecSeed.push_back(argmaxIdx);
            _vecVldtInf.push_back(sumInf);
            coverage[argmaxIdx] = 0;

            for (auto edgeIdx : _hyperGraph._FRsets[argmaxIdx])
            {
                if (edgeMark[edgeIdx]) continue;

                edgeMark[edgeIdx] = true;

                for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
                {
                    if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                    coverage[nodeIdx]--;
                }
            }
        }

        degMap.pop_back();
    }

afterGreedy:
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    // std::cout << "  >>>[greedy-lazy] influence: " << finalInf << ", seed set: " << _vecSeed.size() << ", min-bound: " << _boundMin <<
    //           ", last-bound: " << _boundLast << std::endl;

    if (_vecSeed.size() == targetSize)
    {
        // if the sample size is sufficiently large, the result is reliable
        if (_numRRsets > 1000)
        {
            return finalInf;
        }
    }

    // if covering all the RR-sets, the sentinel set may include some nodes which cover only a small number of RR-sets.
    // such nodes should not be included.
    uint32_t threshold = 0.9 * sumInf;
    for (int i = _vecSeed.size() - 1; i > 0; i--)
    {
        if (_vecVldtInf[i - 1] >= threshold)
        {
            degMap[0].push_back(_vecSeed[i]);
            _vecSeed.pop_back();
        }
        else
        {
            break;
        }
    }

    // std::cout << "seedset size reaching 0.9 coverage: " << _vecSeed.size() << std::endl;
    // the following code is to select the nodes with large out-degree. 
    int seedSetSize = (degMap[0].size() > targetSize) ? targetSize : degMap[0].size();
    std::vector<std::pair<uint32_t, uint32_t>> vecHeap;

    for (int i = 0; i < seedSetSize; i++)
    {
        auto &node = degMap[0][i];
        vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
    }

    std::make_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    const auto nodeNum = degMap[0].size();


    for (int i = seedSetSize; i < nodeNum; i++)
    {
        uint32_t node = degMap[0][i];
        uint32_t currDeg = _vecOutDegree[node];

        if (currDeg > vecHeap[0].first)
        {
            std::pop_heap(vecHeap.begin(), vecHeap.end());
            vecHeap.pop_back();
            vecHeap.push_back(std::make_pair(_vecOutDegree[node], node));
            std::push_heap(vecHeap.begin(), vecHeap.end());
        }
    }

    std::sort_heap(vecHeap.begin(), vecHeap.end(), GreaterPair);
    std::unordered_set<uint32_t> seedHashSet(_vecSeed.begin(), _vecSeed.end());

    for (auto &node : vecHeap)
    {
        if (_vecSeed.size() >= targetSize)
        {
            break;
        }

        if (seedHashSet.find(node.second) != seedHashSet.end())
        {
            std::cout << "node exist" << std::endl;
            continue;
        }

        _vecSeed.push_back(node.second);
        seedHashSet.insert(node.second);

        if (_vecVldtInf.size() < _vecSeed.size())
        {
            _vecVldtInf.push_back(sumInf);
        }
    }

    return 1.0 * _numV; // All RR sets are covered.
}

double Alg::MaxCoverTopK(const int targetSize)
{
    FRset coverage(_numV, 0);
    size_t maxDeg = 0;

    for (auto i = _numV; i--;)
    {
        const auto deg = _hyperGraph._FRsets[i].size();
        coverage[i] = deg;

        if (deg > maxDeg) maxDeg = deg;
    }

    RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree

    for (auto i = _numV; i--;)
    {
        //if (coverage[i] == 0) continue;
        degMap[coverage[i]].push_back(i);
    }

    Nodelist sortedNode(_numV); // sortedNode: record the sorted nodes in ascending order of degree
    Nodelist nodePosition(_numV); // nodePosition: record the position of each node in the sortedNode
    Nodelist degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
    uint32_t idxSort = 0;
    size_t idxDegree = 0;

    for (auto& nodes : degMap)
    {
        degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
        idxDegree++;

        for (auto& node : nodes)
        {
            nodePosition[node] = idxSort;
            sortedNode[idxSort++] = node;
        }
    }

    // check if an edge is removed
    std::vector<bool> edgeMark(_numRRsets, false);
    // record the total of top-k marginal gains
    size_t sumTopk = 0;

    for (auto deg = maxDeg + 1; deg--;)
    {
        if (degreePosition[deg] <= _numV - targetSize)
        {
            sumTopk += deg * (degreePosition[deg + 1] - (_numV - targetSize));
            break;
        }

        sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
    }

    _boundMin = 1.0 * sumTopk;
    _vecSeed.clear();
    size_t sumInf = 0;

    /*
    * sortedNode: position -> node
    * nodePosition: node -> position
    * degreePosition: degree -> position (start position of this degree)
    * coverage: node -> degree
    * e.g., swap the position of a node with the start position of its degree
    * swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
    */
    for (auto k = targetSize; k--;)
    {
        const auto seed = sortedNode.back();
        sortedNode.pop_back();
        const auto newNumV = sortedNode.size();
        sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
        sumInf += coverage[seed];
        _vecSeed.push_back(seed);
        coverage[seed] = 0;

        for (auto edgeIdx : _hyperGraph._FRsets[seed])
        {
            if (edgeMark[edgeIdx]) continue;

            edgeMark[edgeIdx] = true;

            for (auto nodeIdx : _hyperGraph._RRsets[edgeIdx])
            {
                if (coverage[nodeIdx] == 0) continue; // This node is seed, skip

                const auto currPos = nodePosition[nodeIdx]; // The current position
                const auto currDeg = coverage[nodeIdx]; // The current degree
                const auto startPos = degreePosition[currDeg]; // The start position of this degree
                const auto startNode = sortedNode[startPos]; // The node with the start position
                // Swap this node to the start position with the same degree, and update their positions in nodePosition
                std::swap(sortedNode[currPos], sortedNode[startPos]);
                nodePosition[nodeIdx] = startPos;
                nodePosition[startNode] = currPos;
                // Increase the start position of this degree by 1, and decrease the degree of this node by 1
                degreePosition[currDeg]++;
                coverage[nodeIdx]--;

                // If the start position of this degree is in top-k, reduce topk by 1
                if (startPos >= newNumV - targetSize) sumTopk--;
            }
        }

        _boundLast = 1.0 * (sumInf + sumTopk);

        if (_boundMin > _boundLast) _boundMin = _boundLast;
    }

    _boundMin *= 1.0 * _numV / _numRRsets;
    _boundLast *= 1.0 * _numV / _numRRsets;
    const auto finalInf = 1.0 * sumInf * _numV / _numRRsets;
    std::cout << "  >>>[greedy-topk] influence: " << finalInf << ", min-bound: " << _boundMin <<
              ", last-bound: " << _boundLast << '\n';
    return finalInf;
}

double Alg::MaxCover(const int targetSize)
{
    if (targetSize >= 1000) return MaxCoverTopK(targetSize);

    return MaxCoverVanilla(targetSize);
}

void Alg::set_prob_dist(const ProbDist dist)
{
    _probDist = dist;
    _hyperGraph.set_prob_dist(dist);
    _hyperGraphVldt.set_prob_dist(dist);
}

void Alg::set_vanilla_sample(const bool isVanilla)
{
    if (isVanilla)
    {
        std::cout << "Vanilla sampling method is used" << std::endl;
    }

    _hyperGraph.set_vanilla_sample(isVanilla);
    _hyperGraphVldt.set_vanilla_sample(isVanilla);
}

void Alg::set_probseed(const std::string seedfile) {
    // check existence of seedfile
    TIO::ReadProbSeeds(seedfile, this->_vecSeed, this->_vecProb);
    assert(this->_vecSeed.size() > 0);
    for (int i = 0; i < this->_vecSeed.size(); ++i) {
        this->_mapSeedProb[this->_vecSeed[i]] = this->_vecProb[i];
        // this->_vecBoolSeed[this->_vecSeed[i]] = true;
    }
    return;
}

void Alg::set_cand_edges_prob(const std::string cand_edges_prob="orig") {
    if (cand_edges_prob == "orig") {
        return;
    }
    if (cand_edges_prob == "wc")
    {
        // Traverse the candidate edges
        for (uint32_t i=0; i<this->_candEdges.size(); ++i) {
            // Traverse the nodes in the candidate edge
            for (uint32_t j=0; j<this->_candEdges[i].size(); ++j) {
                this->_candEdges[i][j].second = 1.0 / (double)this->_hyperGraph._graph[i].size();
            }
        }
        LogInfo("Finish setting the candidate edges probability with WC");
    }
    
}

double Alg::EfficInfVldtAlg()
{
    return EfficInfVldtAlg(_vecSeed);
}

double Alg::EfficInfVldtAlg(const Nodelist vecSeed)
{
    Timer EvalTimer("Inf. Eval.");
    std::cout << "  >>>Evaluating influence in [0.99,1.01]*EPT with prob. 99.9% ...\n";
    const auto inf = _hyperGraph.EfficInfVldtAlg(vecSeed);
    //const auto inf = _hyperGraphVldt.EfficInfVldtAlg(vecSeed);
    std::cout << "  >>>Done! influence: " << inf << ", time used (sec): " << EvalTimer.get_total_time() << '\n';
    return inf;
}

double Alg::estimateRRSize()
{
    const int sampleNum = 100;
    _hyperGraph.BuildRRsets(sampleNum);
    double avg= _hyperGraph.HyperedgeAvg();
    _hyperGraph.RefreshHypergraph();
    return avg;
}

double Alg::subsimOnly(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR); // R1
        _hyperGraphVldt.BuildRRsets(numR); // R2
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCover(targetSize);
        time2 += timerSubsim.get_operation_time();
        const auto infVldt = _hyperGraphVldt.CalculateInf(_vecSeed);

        const auto degVldt = infVldt * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double avgSize = _hyperGraph.HyperedgeAvg();

        if (currApprox >= approx - epsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Influence via R2: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
            return 0;
        }
    }

    return 0.0;
}

int decideMultiple(int ratio, int numRRsets)
{
    int multiple = 1;

    if (numRRsets < 100)
    {
        return multiple;
    }

    if (ratio >= 32)
    {
        multiple = 8;
    }
    else if (ratio >= 16)
    {
        multiple = 4;
    }
    else if (ratio >= 4)
    {
        multiple = 2;
    }
    else
    {
        multiple = 1;
    }

    return multiple;
}

double Alg::subsimWithTrunc(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
    const auto numRbase = size_t(3 * log(1 / delta));
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    double infVldt = 0.0;
    int multiple = 1;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        const auto infSelf = MaxCoverOutDegPrority(targetSize);
        time2 += timerSubsim.get_operation_time();
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        infVldt = _hyperGraphVldt.EvalSeedSetInf(connSet, _numRRsets * multiple);
        time4 += timerSubsim.get_operation_time();
        const auto degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();
        // if truncRRset is more efficient, increase the size of R2 in next iteration
        int ratio = fullRRSize / truncRRSize;
        multiple = decideMultiple(ratio, _numRRsets);

        if (currApprox >= approx - epsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Time for full RR sets: " << time1  << std::endl;
            std::cout << "==>Time for truncated RR set: " << time4 << std::endl;
            std::cout << "==>Time for greedy: " << time2 << std::endl;
            return 0;
        }
    }

    return 0.0;
}

double Alg::IncreaseR2(std::unordered_set<uint32_t> &connSet, double a, double upperOPT, double targetAppr)
{
    size_t vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
    size_t R1RRsets = _hyperGraph.get_RR_sets_size();
    int multiple = 4;
    double estimateAppr = 0.0;
    double lowerSelect = 0;
    int maxMultiple = 3;
    double infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
    double degVldt = infVldt * vldtRRsets / _numV;
    lowerSelect = (pow2(sqrt(degVldt * multiple + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0) ;
    estimateAppr = (lowerSelect / (multiple * vldtRRsets)) / (upperOPT / R1RRsets);

    if (estimateAppr < targetAppr)
    {
        return 0.0;
    }

    _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, vldtRRsets * multiple);
    infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
    vldtRRsets = _hyperGraphVldt.get_RR_sets_size();
    degVldt = infVldt *  vldtRRsets / _numV;
    lowerSelect = (pow2(sqrt(degVldt + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0);
    double newAppr = (lowerSelect / vldtRRsets) / (upperOPT / R1RRsets);
    return (newAppr > targetAppr) ? newAppr : 0.0;
}

double Alg::FindRemSet(const int targetSize, const double epsilon, const double targetEpsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    size_t subSeedSetSize = _vecSeed.size();
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    // delta for upper bound on the number of RR sets
    const double delta_upper = delta / 3.0;
    const double alpha = sqrt(log(3.0 / delta_upper));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize - subSeedSetSize) + log(3.0 / delta_upper)));
    //const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
    const auto numRbase = size_t(_baseNumRRsets);
    const auto maxNumR = size_t(2.0 * _numV * pow2(alpha + beta) / targetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 3.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    double infVldt = 0.0;
    double currApprox = 0.0;
    double infSelf = 0.0;
    int multiple = 1;
    std::unordered_set<uint32_t> subSeedSet(_vecSeed.begin(), _vecSeed.end());
    std::vector<uint32_t> vecSubSeed(_vecSeed.begin(), _vecSeed.end());

    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        _hyperGraph.BuildRRsetsEarlyStop(subSeedSet, numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        infSelf = MaxCoverIMSentinel(vecSubSeed, targetSize);
        time2 += timerSubsim.get_operation_time();
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        infVldt = _hyperGraphVldt.EvalSeedSetInf(connSet, _numRRsets * multiple);
        time4 += timerSubsim.get_operation_time();
        const auto degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        const auto upperDegOPT = upperBound * _numRRsets / _numV;
        const auto lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;
        const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
        const auto currApprox = lowerSelect / upperOPT;

        std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.EvalHyperedgeAvg();

        // if truncRRset is more efficient, increase the size of R2 in next iteration
        int ratio = fullRRSize / truncRRSize;
        multiple = decideMultiple(ratio, _numRRsets);

        if (currApprox >= approx - targetEpsilon)
        {
            _res.set_approximation(currApprox);
            _res.set_running_time(timerSubsim.get_total_time());
            _res.set_influence(infVldt);
            _res.set_influence_original(infSelf);
            _res.set_seed_vec(_vecSeed);
            _res.set_RR_sets_size(_numRRsets * 2);
            std::cout << "==>Time for full RR in IM-Sentinel phase: " << time1  << std::endl;
            std::cout << "==>Time for truncated RR in IM-Sentinel phase: " << time4 << std::endl;
            std::cout << "==>Time for greedy in IM-Sentinel phase: " << time2 << std::endl;
            std::cout << "==>Influence via R2 in IM-Sentinel phase: " << infVldt << ", time: " << _res.get_running_time() << '\n';
            return 0;
        }
    }

    return 0.0;
}

double Alg::FindDynamSub(const int totalTargetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    const double e = exp(1);
    const double x = (1.0 - 1.0 / totalTargetSize);
    const int minSubSize = ceil(log(1 - epsilon) / log(x));
    const double alpha = sqrt(log(6.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, totalTargetSize) + log(6.0 / delta)));
    const auto numRbasePrevious = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize);
    const auto numRbase = size_t(_baseNumRRsets);
    
    // the successful probability of at least 1-delta/3
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / totalTargetSize / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 3.0 / delta);
    const double a2 = log(numIter * 6.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;
    double time4 = 0.0;
    int multiple = 1;
    double infVldt = 0.0;
    bool firstRound = true;


    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerSubsim.get_operation_time();
        //build R1
        _hyperGraph.BuildRRsets(numR); // R1
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerSubsim.get_operation_time();
        _vecVldtInf.clear();
        int targetSize = firstRound ? (totalTargetSize / 4) : (totalTargetSize / 8);
        firstRound = false;
        const auto infSelf = MaxCoverSentinelSet(targetSize, totalTargetSize);
        time2 += timerSubsim.get_operation_time();
        std::vector<double> vecAppro(_vecSeed.size());
        int lastPos = 0;
        bool found = false;

        if (_vecSeed.size() < targetSize)
        {
            //low influence
            continue;
        }

        double calcAppr = 0.0;

        for (int i = _vecSeed.size() - 1; i >= 0; i--)
        {
            infVldt = _vecVldtInf[i];
            double lowerDeg = infVldt;
            double upperDeg = _boundMin;
            double a = log(numIter * 6.0  / delta);
            double lower = pow2(sqrt(lowerDeg + a * 2.0 / 9.0) - sqrt(a / 2.0)) - a / 18.0;
            double upper = pow2(sqrt(upperDeg + a1 / 2.0) + sqrt(a1 / 2.0));
            upper = (upper > numR) ? numR : upper;
            vecAppro[i] = lower / upper;
            calcAppr = (1 - pow(x, i + 1) - epsilon) * 1.2;

            if (vecAppro[i] > calcAppr)
            {
                found = true;
                lastPos = i;
                break;
            }
        }


        if (!found)
        {
            lastPos = 0;
        }

        size_t setSize = (lastPos + 1);
        setSize = setSize > 10 ? setSize : 10;
        setSize = (setSize > targetSize) ? targetSize : setSize;
        setSize = (setSize < minSubSize) ? minSubSize : setSize;
        std::vector<uint32_t> dynSeedSet(_vecSeed.begin(), _vecSeed.begin() + setSize);
        _vecSeed.clear();
        _vecSeed.assign(dynSeedSet.begin(), dynSeedSet.end());
        std::unordered_set<uint32_t> connSet(_vecSeed.begin(), _vecSeed.end());
        _hyperGraphVldt.RefreshHypergraph();
        _hyperGraphVldt.BuildRRsetsEarlyStop(connSet, _numRRsets * multiple);
        infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
        time4 += timerSubsim.get_operation_time();
        double degVldt = infVldt * multiple * _numRRsets / _numV;
        auto upperBound = _boundMin;

        double upperDegOPT = upperBound * _numRRsets / _numV;
        double lowerSelect = (pow2(sqrt(degVldt + a2 * 2.0 / 9.0) - sqrt(a2 / 2.0)) - a2 / 18.0) / multiple;

        if (lowerSelect < 0)
        {
            lowerSelect = 1.0 * _vecSeed.size() / _numV * _numRRsets * multiple;
        }

        double upperOPT = pow2(sqrt(upperDegOPT + a1 / 2.0) + sqrt(a1 / 2.0));
        upperOPT = (upperOPT > _numRRsets) ? _numRRsets : upperOPT;
        const auto currApprox = lowerSelect / upperOPT;
        std::cout << "lower bound: " << (lowerSelect * _numV / (_numRRsets)) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
        std::cout << "-->SUBSIM (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                  " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
        const double approx = 1 - pow(x, _vecSeed.size());
        double targetAppr = approx - epsilon;

        if (currApprox >= targetAppr)
        {
            goto succ;
        }

        if (_numRRsets < 100)
        {
            continue;
        }

        double fullRRSize = _hyperGraph.HyperedgeAvg();
        double truncRRSize = _hyperGraphVldt.HyperedgeAvg();

        if (fullRRSize / truncRRSize < 2)
        {
            continue;
        }

        double lowerThreshold = (upperOPT * _numV / _numRRsets) * targetAppr;

        if ((1.0 * infVldt / multiple) > lowerThreshold && lowerThreshold > 0)
        {
            double newAppr = IncreaseR2(connSet, a2, upperOPT, targetAppr);
            time4 += timerSubsim.get_operation_time();

            if (newAppr > targetAppr)
            {
                std::cout << "increase R2 successfully" << std::endl;
                infVldt = _hyperGraphVldt.CalculateInfEarlyStop();
                goto succ;
            }
        }
    }

succ:
    std::cout << "==>Time for full RR in SentinelSet phase: " << time1  << std::endl;
    std::cout << "==>Time for truncated RR in SentinelSet phase: " << time4 << std::endl;
    std::cout << "==>Time for greedy in SentinelSet phase: " << time2 << std::endl;
    std::cout << "==>size of sentinel set: " << _vecSeed.size() << ", inf: " << infVldt << std::endl;
    std::cout << "==>total time for SentinelSet phase: " << timerSubsim.get_total_time() << std::endl;
    return 0.0;
}

double Alg::subsimWithHIST(const int targetSize, const double epsilon, const double delta)
{
    Timer timerSubsim("SUBSIM");
    _baseNumRRsets = 3 * log(1 / delta);
    
    std::cout << std::endl;
    std::cout << "Sentinel Set Selection Phase" << std::endl;
    FindDynamSub(targetSize, epsilon / 2, delta / 2);
    _hyperGraph.RefreshHypergraph();
    _hyperGraphVldt.RefreshHypergraph();

    std::cout << std::endl;
    std::cout << "IM-Sentinel Phase" << std::endl;
    FindRemSet(targetSize, epsilon / 2, epsilon, delta / 2);
    _res.set_running_time(timerSubsim.get_total_time());
    return 0.0;
}

// ---------- IMAPS ----------
double Alg::EdgeSelectionPMC(const int targetSize, std::string mode="UB") {
    // Use a heap to store the <node, prob coverage> pair
    std::priority_queue<std::pair<uint32_t, double>, std::vector<std::pair<uint32_t, double>>, CompareBySecondDouble> heap;
    std::vector<double> coverage(this->_numV, 0.0);
    std::vector<double> RR_utility(this->_hyperGraph._vecUtility);

    for(uint32_t i=0; i<this->_numV; i++) {
        // store coverage
        // if the node is not covered 
        if (this->_candEdges[i].empty())
            continue;
        double prob_for_i = this->_candEdges[i].back().second;
        assert(prob_for_i > 0);
        double node_i_cov = prob_for_i * this->_hyperGraph._vecNodeGain[i];
        heap.push(std::make_pair(i, node_i_cov));
        coverage[i] = node_i_cov;
    }
    std::vector<UVWEdge> vec_edges;
    
    uint32_t max_idx;
    double cov_num = 0.0;
    double update_time = 0.0;
    // LogInfo("Start selection");
    while (vec_edges.size() < targetSize)
    {
        std::pair<uint32_t, double> top = heap.top();
        heap.pop();

        // Lazy Update
        if (top.second > coverage[top.first]) {
            // Update coverage of top
            top.second = coverage[top.first];
            heap.push(top);
            continue;
        }
        max_idx = top.first;
        cov_num += coverage[max_idx];
        // LogInfo("node", max_idx);
        // LogInfo("coverage", coverage[max_idx]);
        
        double origin_cov = coverage[max_idx];
        // selection
        size_t top_seed = this->_candEdges[max_idx].back().first;   // the seed of the target edge
        double top_prob = this->_candEdges[max_idx].back().second;      // the probability of the target edge

        this->_candEdges[max_idx].pop_back();       // pop the edge for node max_idx
        double origin_prob = top_prob;
        if (mode == "UB")
        {
            vec_edges.push_back(std::make_pair(top_seed, Edge(max_idx, origin_prob)));            
        }
        else if (mode == "LB") {
            origin_prob /= _mapSeedProb[top_seed];
            vec_edges.push_back(std::make_pair(top_seed, Edge(max_idx, origin_prob)));
        }

        std::vector<size_t>& RRsets_cov_by_max_node = this->_hyperGraph._FRsets[max_idx];
        // LogInfo("Updating");
        for (uint32_t j=0; j<RRsets_cov_by_max_node.size(); j++) {
            size_t RRset_idx = RRsets_cov_by_max_node[j];
            Nodelist& node_list = this->_hyperGraph._RRsets[RRset_idx];
            // timer.get_operation_time();
            for (uint32_t& node: node_list){
                if (node != max_idx && !this->_candEdges[node].empty()) {
                    double cur_node_prob = this->_candEdges[node].back().second;
                    coverage[node] -= (RR_utility[RRset_idx] * top_prob * cur_node_prob);
                }
            }
            // update_time += timer.get_operation_time();
            // update RR set utility
            RR_utility[RRset_idx] *= (1 - top_prob);
        }

        //update edges pointing to the same node
        coverage[max_idx] = 0.0;
        if (!this->_candEdges[max_idx].empty()) {
            double new_prob = this->_candEdges[max_idx].back().second;
            // double test = 0;
            // for (uint32_t j=0; j<RRsets_cov_by_max_node.size(); j++) {
            //     size_t RRset_idx = RRsets_cov_by_max_node[j];
            //     // coverage[max_idx] += live_prob_per_RRset[RRset_idx];
            //     test += Fnew_prob * live_prob_per_RRset[RRset_idx];
            // }
            coverage[max_idx] = origin_cov / top_prob * new_prob * (1-top_prob);
            heap.push(std::make_pair(max_idx, coverage[max_idx]));
        }
        double extra_val = 0.0;
        uint32_t topk = 0;
        std::vector<std::pair<uint32_t, double>> vec_pop;
        // LogInfo("Fetching Topk");
        // while (topk < targetSize)
        // {
        //     std::pair<uint32_t, double> top = heap.top();
            
        //     vec_pop.push_back(top);
        //     heap.pop();
        //     if (top.second > coverage[top.first])
        //     {
        //         heap.push(std::make_pair(top.first, coverage[top.first]));
        //         continue;
        //     }
        //     extra_val += top.second;
        //     topk++;
        // }
        // _boundLast = (extra_val + cov_num) * _numV / _numRRsets;
        // if (_boundMin > _boundLast) _boundMin = _boundLast;

        // for (auto &node : vec_pop)
        // {
        //     heap.push(node);
        // }
    }
    if (mode == "UB")
    {
        this->_vecResEdgesUB = vec_edges;
    }
    else if (mode == "LB") {
        this->_vecResEdgesLB = vec_edges;
    }
    // LogInfo("Post process: adding cand edges back");
    // Post-process
    for (int i = vec_edges.size()-1; i >= 0; i--)
    {
        // LogInfo(i);
        auto& uvw = vec_edges[i];
        uint32_t u = uvw.first, v = uvw.second.first;
        double w = uvw.second.second;

        if (mode == "LB") {
            w *= _mapSeedProb[u];
        }

        // std::cout << u << " " << v << " " << w << std::endl;
        this->_candEdges[v].push_back(std::make_pair(u,w));
    }

    LogInfo("cov num", cov_num);
    return cov_num  * this->_numV / this->_numRRsets;
}

double Alg::EdgeSelectionLB(const int targetSize) { 
    std::priority_queue<std::pair<uint32_t, double>, std::vector<std::pair<uint32_t, double>>, CompareBySecondDouble> heap;
    std::vector<double> coverage(this->_numCandEdges, 0.0);
    std::vector<double> RR_utility(this->_hyperGraph._vecUtility);

    // initialize coverage for each edge
    for (uint32_t i=0; i<this->_numCandEdges; i++) {
        auto& ij = this->_vecIdx4Edges[i];
        uint32_t v = ij.first, u = this->_candEdges[v][ij.second].first;
        double puv = this->_candEdges[v][ij.second].second;
        double pu = this->_mapSeedProb[u];
        coverage[i] = puv * pu * this->_hyperGraph._vecEdgeGain[i];
        heap.push(std::make_pair(i, coverage[i]));
    }
    std::vector<UVWEdge> vec_edges;
    
    uint32_t max_idx;
    double cov_num = 0.0;
    double update_time = 0.0;
    // LogInfo("Start selection");
    while (vec_edges.size() < targetSize)
    {
        std::pair<uint32_t, double> top = heap.top();
        heap.pop();

        // Lazy Update
        if (top.second > coverage[top.first]) {
            // Update coverage of top
            top.second = coverage[top.first];
            heap.push(top);
            continue;
        }
        max_idx = top.first;
        cov_num += coverage[max_idx];
        // LogInfo("edge idx", max_idx);
        // LogInfo("coverage", coverage[max_idx]);
        
        double origin_cov = coverage[max_idx];
        // selection
        auto& ij = this->_vecIdx4Edges[max_idx];
        uint32_t top_v = ij.first;
        uint32_t top_seed = this->_candEdges[top_v][ij.second].first;   // the seed of the target edge
        double top_prob = this->_candEdges[top_v][ij.second].second;      // the probability of the target edge

        vec_edges.push_back(std::make_pair(top_seed, Edge(top_v, top_prob)));
        // std::cout << top_seed << " " << top_v << " " << top_prob << std::endl;
        std::vector<size_t>& RRsets_cov_by_max_node = this->_hyperGraph._FRsets4Edges[max_idx];
        // LogInfo("Updating");
        // LogInfo("RR sets to update", RRsets_cov_by_max_node.size());
        for (uint32_t j=0; j<RRsets_cov_by_max_node.size(); j++) {
            size_t RRset_idx = RRsets_cov_by_max_node[j];
            Nodelist& edge_list = this->_hyperGraph._RRsets4Edges[RRset_idx];
            // timer.get_operation_time();
            for (uint32_t& e_idx: edge_list){
                if (e_idx == max_idx) continue;
                auto v = this->_vecIdx4Edges[e_idx].first;
                auto& uw = this->_candEdges[v][this->_vecIdx4Edges[e_idx].second];
                auto u = uw.first; auto w = uw.second;
                if (u != top_seed) {
                    double a = RR_utility[RRset_idx] / ((1-_mapSeedProb[u]*this->_hyperGraph._acprob_by_seed[RRset_idx][u])*(1-_mapSeedProb[top_seed]*this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]));
                    double b = _mapSeedProb[top_seed] * top_prob * (1-this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]);
                    double c = _mapSeedProb[u] * w * (1-this->_hyperGraph._acprob_by_seed[RRset_idx][u]);
                    
                    coverage[e_idx] -= (a*b*c);
                    // if (a < 0.0 || a>1.0) {
                    //     std::cout << "a: " << a << std::endl;
                    //     std::cout << "RR_utility[RRset_idx]: " << RR_utility[RRset_idx] << std::endl;
                    //     std::cout << "this->_hyperGraph._acprob_by_seed[RRset_idx][u]: " << this->_hyperGraph._acprob_by_seed[RRset_idx][u] << std::endl;
                    //     std::cout << "this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]: " << this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed] << std::endl;
                    // }
                    // if (a < 0.0 || a>1.0) {
                    //     std::cout << "a: " << a << std::endl;
                    //     std::cout << "RR_utility[RRset_idx]: " << RR_utility[RRset_idx] << std::endl;
                    //     std::cout << "this->_hyperGraph._acprob_by_seed[RRset_idx][u]: " << this->_hyperGraph._acprob_by_seed[RRset_idx][u] << std::endl;
                    //     std::cout << "this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]: " << this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed] << std::endl;
                    // }
                    assert(a >= 0.0 && a <= 1.01);
                    assert(b >= 0.0 && b <= 1.01);
                    assert(c >= 0.0 && c <= 1.01);
                } else {
                    double a = RR_utility[RRset_idx] / (1-_mapSeedProb[top_seed]*this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]);
                    double b = _mapSeedProb[top_seed] * w * top_prob * (1-this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]);
                    coverage[e_idx] -= (a*b);
                    assert(a >= 0.0 && a <= 1.01);
                    assert(b >= 0.0 && b <= 1.01);
                }
            }
            // update_time += timer.get_operation_time();
            double prev_u = RR_utility[RRset_idx];
                
            // update RR set utility
            RR_utility[RRset_idx] = RR_utility[RRset_idx] / (1-_mapSeedProb[top_seed]*this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]);
            this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed] += (1-this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed])*top_prob;
            assert(this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed] >= 0.0 && this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed] <= 1.0);
            RR_utility[RRset_idx] *= (1-_mapSeedProb[top_seed]*this->_hyperGraph._acprob_by_seed[RRset_idx][top_seed]);
            // if (prev_u > 0.5) {
            //     LogInfo("RRset_idx", RRset_idx);
            //     LogInfo("before update", prev_u);
            //     LogInfo("after update", RR_utility[RRset_idx]);
            // }
            assert(RR_utility[RRset_idx] >= 0.0 && RR_utility[RRset_idx] <= 1.0);
        }
        // std::cout << "max_idx: " << max_idx << std::endl;
        coverage[max_idx] = 0.0;
    }

    this->_vecResEdgesLB = vec_edges;
    LogInfo("Post process: clear acprob");
    for (uint32_t i=0; i<this->_numRRsets; i++) {
        for (auto& kv: this->_hyperGraph._acprob_by_seed[i]) {
            if (kv.second < -0.5) continue;
            kv.second = 0.0;
        }
    }
    // Post-process
    // for (int i = vec_edges.size()-1; i >= 0; i--)
    // {
    //     // LogInfo(i);
    //     auto& uvw = vec_edges[i];
    //     uint32_t u = uvw.first, v = uvw.second.first;
    //     double w = uvw.second.second;

    //     // std::cout << u << " " << v << " " << w << std::endl;
    //     this->_candEdges[v].push_back(std::make_pair(u,w));
    // }

    LogInfo("cov num", cov_num);
    return cov_num  * this->_numV / this->_numRRsets;
}

void Alg::InitIdx4Edges() {
    for (uint32_t i=0; i<this->_numV; i++) {
        for (uint32_t j=0; j<this->_candEdges[i].size(); j++) {
            this->_vecIdx4Edges.push_back(std::make_pair(i, j));
        }
    }
    return;
}

double Alg::boundMaximize(const int targetSize, const double epsilon, const double delta, const std::string mode="UB") {
    Timer timerboundMax("boundMaximize");
    const double e = exp(1);
    const double approx = 1 - 1.0 / e;
    const double alpha = sqrt(log(12.0 / delta));
    const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(12.0 / delta)));
    const auto numRbase = std::max(size_t(2.0 * pow2((1 - 1 / e) * alpha + beta)), this->_numRRsets);
    const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / 10 / pow2(epsilon)) + 1;
    const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
    const double a1 = log(numIter * 6.0 / delta);
    const double a2 = log(numIter * 6.0 / delta);
    double time1 = 0.0, time2 = 0.0, time3 = 0.0;

    std::cout << std::endl;
    for (auto idx = 1; idx <= numIter; idx++)
    {
        const auto numR = numRbase << (idx-1);
        if (numR < this->_numRRsets) continue;
        std::cout << "Iteration: " << idx << " RR set: " << numR << std::endl;
        timerboundMax.get_operation_time();
        // LogInfo("Building RR sets");
        _hyperGraph.BuildRRsets(numR); // R1
        _hyperGraphVldt.BuildRRsets(numR); // R2
        // LogInfo("Updating Utilities");
        _hyperGraph.updateUtility(_mapSeedProb);
        _hyperGraphVldt.updateUtility(_mapSeedProb);
        _numRRsets = _hyperGraph.get_RR_sets_size();
        time1 += timerboundMax.get_operation_time();
        // LogInfo("Edge Selection");
        
        time2 += timerboundMax.get_operation_time();
        double infVldt = 0.0;
        if (mode == "UB") {
            const auto infSelf = EdgeSelectionPMC(targetSize, mode);
            time2 += timerboundMax.get_operation_time();
            infVldt = _hyperGraphVldt.ComputeInfUpperBound(_vecResEdgesUB, _mapSeedProb);
            const auto degVldt = infVldt * _numRRsets / _numV;
            // auto upperBound = _boundMin;
            auto upperBound = infSelf;
            const auto upperDegOPT = upperBound * _numRRsets / _numV / approx;
            const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
            const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
            const auto currApprox = lowerSelect / upperOPT;
            std::cout << "InfSelf: " << infSelf << ", infVldt: " << infVldt << std::endl;
            std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
            std::cout << "-->IMAPS (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                    " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
            double avgSize = _hyperGraph.HyperedgeAvg();

            if (currApprox >= approx - epsilon)
            {
                _res.set_mode(mode);
                _res.set_approximation_UB(currApprox);
                _res.set_running_time_UB(timerboundMax.get_total_time());
                _res.set_influence_original(infSelf);
                _res.set_edges_bounds(_vecResEdgesUB);
                _res.set_UB_value(infVldt);
                _res.set_sampling_time_UB(time1);
                _res.set_selection_time_UB(time2);
                _res.set_RR_sets_size_UB(_numRRsets * 2);
                std::cout << "==>Influence via R2: " << infVldt << ", time: " << _res.get_running_time() << '\n';
                std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
                return 0.0;
            }
        }
        else if (mode == "LB") {
            LogInfo("--- update coverage graph for edges ---");
            this->_hyperGraph.updateCoverage4Edges(this->_candEdges, this->_vecStIdx4Edges, this->_numCandEdges);
            time2 += timerboundMax.get_operation_time();
            LogInfo("--- update coverage graph for edges complete ---");
            const auto infSelf = EdgeSelectionLB(targetSize);
            time2 += timerboundMax.get_operation_time();
            infVldt = _hyperGraphVldt.ComputeInfLowerBound(_vecResEdgesLB, _mapSeedProb);
            const auto degVldt = infVldt * _numRRsets / _numV;
            // auto upperBound = _boundMin;
            auto upperBound = infSelf;
            const auto upperDegOPT = upperBound * _numRRsets / _numV / approx;
            const auto lowerSelect = pow2(sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
            const auto upperOPT = pow2(sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0));
            const auto currApprox = lowerSelect / upperOPT;
            std::cout << "InfSelf: " << infSelf << ", infVldt: " << infVldt << std::endl;
            std::cout << "lower bound: " << (lowerSelect * _numV / _numRRsets) << ", upperBound: " << (upperOPT * _numV / _numRRsets) << std::endl;
            std::cout << "-->IMAPS (" << idx << "/" << numIter << ") approx. (max-cover): " << currApprox <<
                    " (" << infSelf / upperBound << "), #RR sets: " << _numRRsets << '\n';
            double avgSize = _hyperGraph.HyperedgeAvg();

            if (currApprox >= approx - epsilon)
            {
                _res.set_mode(mode);
                _res.set_approximation_LB(currApprox);
                _res.set_running_time_LB(timerboundMax.get_total_time());
                _res.set_influence_original(infSelf);
                _res.set_edges_bounds(_vecResEdgesLB);
                _res.set_LB_value(infVldt);
                _res.set_sampling_time_LB(time1);
                _res.set_selection_time_LB(time2);
                _res.set_RR_sets_size_LB(_numRRsets * 2);
                std::cout << "==>Influence via R2: " << infVldt << ", time: " << _res.get_running_time() << '\n';
                std::cout << "==>Time for RR sets and greedy: " << time1 << ", " << time2 << '\n';
                return 0.0;
            }
        }
    }

    return 0.0;
}

// double Alg::lowerBoundMaximize(const int targetSize, const double epsilon, const double delta) {
//     Timer timerSubsim("lowerBoundMaximize");
//     const double e = exp(1);
//     const double approx = 1 - 1.0 / e;
//     const double alpha = sqrt(log(6.0 / delta));
//     const double beta = sqrt((1 - 1 / e) * (logcnk(_numV, targetSize) + log(6.0 / delta)));
//     const auto numRbase = _numRRsets;
//     const auto maxNumR = size_t(2.0 * _numV * pow2((1 - 1 / e) * alpha + beta) / targetSize / pow2(epsilon)) + 1;
//     const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
//     const double a1 = log(numIter * 3.0 / delta);
//     const double a2 = log(numIter * 3.0 / delta);
//     double time1 = 0.0, time2 = 0.0, time3 = 0.0;
//     return 0.0;
// }

// Influence Increment Estimation for Sandwich
// double Alg::estimateInfIncforSand(std::vector<UVWEdge>& edges, const double eps, const double delta) {
//     uint32_t k_seed = this->_vecSeed.size();
//     assert(k_seed > 0);
//     assert(k_seed == this->_vecProb.size());
//     double Gamma = 2*(1 + eps)*(1 + eps/3)*log(2.0 / delta) / (eps*eps);
//     size_t theta=this->_numRRsets;
//     double Sigma=0.0;
// }


// Influence Estimation of probabilistic seeds
double Alg::estimateInfStoppingRules(double eps, double delta) {
    this->_hyperGraphVldt.RefreshHypergraph();
    uint32_t k_seed = this->_vecSeed.size();
    uint32_t currProgress = 0;
    Timer EvalTimer("evaluation");
    assert(k_seed > 0);
    double Gamma = 2*(1 + eps)*(1 + eps/3)*log(2.0 / delta) / (eps*eps);

    size_t theta=0;
    double Sigma=0.0;
    LogInfo("Gamma", Gamma);
    while (Sigma < Gamma) {
        this->_hyperGraphVldt.BuildRRsets(this->_hyperGraphVldt.get_RR_sets_size()+1);
        RRset& rrset = this->_hyperGraphVldt._RRsets.back();
        double not_ac_prob = 1.0;
        for (auto node : rrset) {
            if (_mapSeedProb.find(node) != _mapSeedProb.end()) {
                not_ac_prob *= (1 - _mapSeedProb[node]);
            }
        }
        Sigma += (1.0 - not_ac_prob);
        theta++;
        if (Sigma * 100 > Gamma * currProgress) {
            const auto evalTime = EvalTimer.get_operation_time();
            if (evalTime > 100)
                std::cout << "\tEstimation-Progress at: " << currProgress << "%, " << "time used: " << evalTime << std::endl;
            currProgress += 20;
        }
    }
    double est_inf = this->_numV * Gamma / theta;
    LogInfo("Total Number of RR sets", theta);
    LogInfo("Estimated Influence", est_inf);
    return est_inf;
}

// Sandwich
double Alg::sandwich(const int targetSize, const double epsilon, const double delta, const double gamma) {
    double bound_max_time=0.0, estimate_time=0.0;
    Timer timerSandwich("Sandwich");
    
    LogInfo("--- Maximize UB function ---");
    this->boundMaximize(targetSize, epsilon, delta, "UB");
    bound_max_time += timerSandwich.get_operation_time();
    timerSandwich.get_operation_time();
    // print first ten elements of edge LB and edge UB
    // for(int i=0; i<10; i++) {
    //     std::cout << "i: " << i << std::endl;
    //     std::cout << "UB: " << this->_vecResEdgesUB[i].first << " " << this->_vecResEdgesUB[i].second.first << " " << this->_vecResEdgesUB[i].second.second << std::endl;
    // }
    LogInfo("--- Maximize LB function ---");
    this->boundMaximize(targetSize, epsilon, delta, "LB");
    bound_max_time += timerSandwich.get_operation_time();


    this->_hyperGraphVldt.set_vanilla_sample(true);
    AddEdges(_hyperGraphVldt._graph, this->_vecResEdgesUB);
    double infUB = this->estimateInfStoppingRules(gamma, delta);
    PopEdges(_hyperGraphVldt._graph, this->_vecResEdgesUB);

    AddEdges(_hyperGraphVldt._graph, this->_vecResEdgesLB);
    double infLB = this->estimateInfStoppingRules(gamma, delta);
    PopEdges(_hyperGraphVldt._graph, this->_vecResEdgesLB);
    estimate_time = timerSandwich.get_operation_time();
    
    LogInfo("Influence of UB", infUB);
    LogInfo("Influence of LB", infLB);

    _res.set_mode("total");
    _res.set_running_time(timerSandwich.get_total_time());
    _res.set_influence_LB(infLB);
    _res.set_influence_UB(infUB);
    _res.set_estimate_time(estimate_time);
    _res.set_bound_max_time(bound_max_time);
    
    if (infUB < infLB) {
        _res.set_res_edges(_vecResEdgesLB);
    } else {
        _res.set_res_edges(_vecResEdgesUB);
    }
    _res.set_combined_res_two_bounds();
    return 0.0;
}

// RAND
void Alg::selectRandom(const int targetSize) {
    while(this->_vecResEdges.size() < targetSize) {
        int rand_node = rand() % this->_numV;
        if (!this->_candEdges[rand_node].empty()) {
            this->_vecResEdges.push_back(std::make_pair(this->_candEdges[rand_node].back().first, Edge(rand_node, this->_candEdges[rand_node].back().second)));
        }
    }
    _res.set_res_edges(_vecResEdges);
    return;
}

// OUTDEG
void Alg::selectOutDeg(const int targetSize) {
    Timer timerOutdeg("OUTDEG");
    std::vector<std::vector<uint32_t>> outdegMap(this->_numV);
    uint32_t maxOutDeg = 0;
    for (uint32_t i=0; i<this->_numV; i++) {
        outdegMap[this->_vecOutDegree[i]].push_back(i);
        if (maxOutDeg < this->_vecOutDegree[i]) {
            maxOutDeg = this->_vecOutDegree[i];
        }
    }
    while (this->_vecResEdges.size() < targetSize) {
        for (int i=maxOutDeg; i>=0; i--) {
            if (outdegMap[i].empty()) {
                continue;
            }
            for (auto node : outdegMap[i]) {
                if (!this->_candEdges[node].empty()) {
                    this->_vecResEdges.push_back(std::make_pair(this->_candEdges[node].back().first, Edge(node, this->_candEdges[node].back().second)));
                }
                if (this->_vecResEdges.size() >= targetSize) {
                    _res.set_mode("total");
                    _res.set_approximation(-1);
                    _res.set_running_time(timerOutdeg.get_total_time());
                    _res.set_influence(0);
                    _res.set_influence_original(0);
                    _res.set_RR_sets_size(_numRRsets * 2);
                    _res.set_res_edges(_vecResEdges);
                    return;
                }
            }
        }
    }
    _res.set_mode("total");
    _res.set_approximation(-1);
    _res.set_running_time(timerOutdeg.get_total_time());
    _res.set_influence(0);
    _res.set_influence_original(0);
    _res.set_RR_sets_size(_numRRsets * 2);
    _res.set_res_edges(_vecResEdges);
    return;
}

void Alg::selectGreedy(const int targetSize, const size_t num_samples) {
    std::priority_queue<std::pair<size_t, double>, std::vector<std::pair<size_t, double>>, CompareBySecondDoubleSizet> heap;
    std::vector<double> edge_influence;
    std::vector<std::pair<uint32_t, uint32_t>> idx_ij;
    // Traverse candidate edges and initialize edge influence
    size_t edge_cnt = 0;
    Timer timerGreedy("Greedy");

    double time_sampling = 0.0;
    double cur_sigma = this->_hyperGraph.BuildRRsetsEstimation(num_samples, this->_mapSeedProb);
    LogInfo("--- Initialize influence ---");
    for (uint32_t i=0; i<this->_numV; i++) {
        for (uint32_t j=0; j<this->_candEdges[i].size(); j++) {
            auto& edge = this->_candEdges[i][j];
            uint32_t u = edge.first;
            double w = edge.second;
            UVWEdge uvw_edge = std::make_pair(u, Edge(i, w));

            AddEdge(this->_hyperGraph._graph, uvw_edge);
            timerGreedy.get_operation_time();
            double inf = this->_hyperGraph.BuildRRsetsEstimation(num_samples, this->_mapSeedProb);
            time_sampling += timerGreedy.get_operation_time();
            PopEdge(this->_hyperGraph._graph, uvw_edge);
            heap.push(std::make_pair(edge_cnt++, inf - cur_sigma));
            // edge_influence.push_back(inf);
            idx_ij.push_back(std::make_pair(i, j));
        }
        assert(heap.size() == idx_ij.size());
        uint16_t currProgress = 0;
        if (i % 1000 == 0)
        {
            LogInfo("1000 passed");
            LogInfo("Edges count", edge_cnt);
            LogInfo("Time sampling", time_sampling);
        }
        
        // if (this->_numV * 100 > i * currProgress) {
        //     const auto evalTime = timerGreedy.get_operation_time();
        //     if (evalTime > 100)
        //         std::cout << "\tEstimation-Progress at: " << currProgress << "%, " << "time used: " << evalTime << std::endl;
        //     currProgress += 20;
        // }
        
    }

    LogInfo("--- Start selection ---");
    while (this->_vecResEdges.size() < targetSize)
    {
        std::pair<size_t, double> top = heap.top();
        heap.pop();

        std::pair<uint32_t, uint32_t> ij = idx_ij[top.first];
        uint32_t v = ij.first, u = this->_candEdges[v][ij.second].first;
        double w = this->_candEdges[v][ij.second].second;
        UVWEdge uvw_edge = std::make_pair(u, Edge(v, w));
        this->_vecResEdges.push_back(uvw_edge);
        AddEdge(this->_hyperGraph._graph, uvw_edge);
        cur_sigma += top.second;
        std::cout << u << " " << v << std::endl;
        LogInfo("Marginal Influence", top.second);

        size_t prev_top = heap.top().first;
        uint32_t cnt = 0;
        while (cnt == 0 || heap.top().first != prev_top)
        {
            top = heap.top();
            heap.pop();
            prev_top = top.first;
            ij = idx_ij[top.first];
            v = ij.first; u = this->_candEdges[v][ij.second].first;
            w = this->_candEdges[v][ij.second].second;
            
            UVWEdge uvw_edge = std::make_pair(u, Edge(v, w));
            AddEdge(this->_hyperGraph._graph, uvw_edge);
            timerGreedy.get_operation_time();
            double inf = this->_hyperGraph.BuildRRsetsEstimation(num_samples, this->_mapSeedProb);
            time_sampling += timerGreedy.get_operation_time();
            PopEdge(this->_hyperGraph._graph, uvw_edge);
            heap.push(std::make_pair(prev_top, inf-cur_sigma));
            ++cnt;
        }
    }
    
    _res.set_mode("total");
    _res.set_approximation(-1);
    _res.set_running_time(timerGreedy.get_total_time());
    // _res.set_influence(heap.top().second);
    _res.set_influence_original(cur_sigma);
    _res.set_RR_sets_size(_numRRsets);
    _res.set_res_edges(_vecResEdges);
    return;
}

void Alg::selectGreedyNode(const int targetSize, const size_t num_samples) {
    std::priority_queue<std::pair<size_t, double>, std::vector<std::pair<size_t, double>>, CompareBySecondDoubleSizet> heap;
    std::vector<double> edge_influence;
    // std::vector<std::pair<uint32_t, uint32_t>> idx_ij;
    // Traverse candidate edges and initialize edge influence
    size_t edge_cnt = 0;
    Timer timerGreedy("Greedy");

    double time_sampling = 0.0;
    double cur_sigma = this->_hyperGraph.BuildRRsetsEstimation(num_samples, this->_mapSeedProb);
    LogInfo("--- Initialize influence ---");
    for (uint32_t i=0; i<this->_numV; i++) {
        auto& edge = this->_candEdges[i].back();
        uint32_t u = edge.first;
        double w = edge.second;
        UVWEdge uvw_edge = std::make_pair(u, Edge(i, w));

        AddEdge(this->_hyperGraph._graph, uvw_edge);
        timerGreedy.get_operation_time();
        double inf = this->_hyperGraph.BuildRRsetsEstimation(num_samples, this->_mapSeedProb);
        time_sampling += timerGreedy.get_operation_time();
        PopEdge(this->_hyperGraph._graph, uvw_edge);
        heap.push(std::make_pair(i, inf-cur_sigma));

        uint16_t currProgress = 0;
        if (i % 1000 == 0) {
            LogInfo("1000 passed");
            LogInfo("Time sampling", time_sampling);
        }
        
        // if (this->_numV * 100 > i * currProgress) {
        //     const auto evalTime = timerGreedy.get_operation_time();
        //     if (evalTime > 100)
        //         std::cout << "\tEstimation-Progress at: " << currProgress << "%, " << "time used: " << evalTime << std::endl;
        //     currProgress += 20;
        // }
        
    }

    LogInfo("--- Start selection ---");
    while (this->_vecResEdges.size() < targetSize)
    {
        std::pair<size_t, double> top = heap.top();
        heap.pop();

        // std::pair<uint32_t, uint32_t> ij = idx_ij[top.first];
        uint32_t v = top.first, u = this->_candEdges[v].back().first;
        double w = this->_candEdges[v].back().second;

        UVWEdge uvw_edge = std::make_pair(u, Edge(v, w));
        this->_vecResEdges.push_back(uvw_edge);
        AddEdge(this->_hyperGraph._graph, uvw_edge);
        cur_sigma += top.second;
        
        std::cout << u << " " << v << std::endl;
        LogInfo("Influence", top.second);

        this->_candEdges[v].pop_back();
        size_t prev_top = v;
        if (!this->_candEdges[v].empty())
        {
            u = this->_candEdges[v].back().first;
            w = this->_candEdges[v].back().second;
            auto new_uvw_edge = std::make_pair(u, Edge(v, w));
            AddEdge(this->_hyperGraph._graph, new_uvw_edge);
            heap.push(std::make_pair(v, this->_hyperGraph.BuildRRsetsEstimation(num_samples, this->_mapSeedProb) - cur_sigma));
            PopEdge(this->_hyperGraph._graph, new_uvw_edge);
            prev_top = v;
        }
        
        while (heap.top().first != prev_top)
        {
            top = heap.top();
            heap.pop();
            prev_top = top.first;

            v = top.first; u = this->_candEdges[v].back().first;
            w = this->_candEdges[v].back().second;
            
            UVWEdge uvw_edge = std::make_pair(u, Edge(v, w));
            AddEdge(this->_hyperGraph._graph, uvw_edge);
            timerGreedy.get_operation_time();
            double inf = this->_hyperGraph.BuildRRsetsEstimation(num_samples, this->_mapSeedProb);
            time_sampling += timerGreedy.get_operation_time();
            PopEdge(this->_hyperGraph._graph, uvw_edge);
            heap.push(std::make_pair(prev_top, inf - cur_sigma));
        }
    }
    
    _res.set_mode("total");
    _res.set_approximation(-1);
    _res.set_running_time(timerGreedy.get_total_time());
    // _res.set_influence(heap.top().second);
    _res.set_influence_original(cur_sigma);
    _res.set_RR_sets_size(_numRRsets);
    _res.set_res_edges(_vecResEdges);
    return;
}

void Alg::selectAIS(const int targetSize, const double gamma) {
    Timer timerAIS("AIS");
    double init_inf = this->estimateInfStoppingRules(gamma, 1.0 / (double)this->_numV);
    std::vector<bool> vec_bool_seed (this->_numV, false);
    std::vector<bool> vec_bool_cov_by_S (this->_hyperGraphVldt._RRsets.size(), false);
    std::vector<uint32_t> marginal_cov (this->_numV, 0);
    std::priority_queue<std::pair<uint32_t, double>, std::vector<std::pair<uint32_t, double>>, CompareBySecondDouble> heap;
    for (auto seed : this->_vecSeed) {
        vec_bool_seed[seed] = true;
        for (auto RR_idx: this->_hyperGraphVldt._FRsets[seed]) {
            vec_bool_cov_by_S[RR_idx] = true;
        }
    }
    for (uint32_t i=0; i<this->_numV; i++) {
        if (vec_bool_seed[i]) {
            continue;
        }
        for (auto RR_idx: this->_hyperGraphVldt._FRsets[i]) {
            if (!vec_bool_cov_by_S[RR_idx]) {
                marginal_cov[i] += 1;
            }
        }
        if (!this->_candEdges[i].empty()) {
            heap.push(std::make_pair(i, marginal_cov[i] * this->_candEdges[i].back().second));
        }
    }
    // For k iterations
    while (this->_vecResEdges.size() < targetSize)
    {
        auto top = heap.top();
        heap.pop();
        auto node = top.first;
        auto seed = this->_candEdges[node].back().first;
        auto p = this->_candEdges[node].back().second;
        if (top.second > marginal_cov[top.first]*p) {
            top.second = marginal_cov[top.first]*p;
            heap.push(top);
            continue;
        }
        
        this->_vecResEdges.push_back(std::make_pair(seed, Edge(node, p)));
        this->_candEdges[node].pop_back();

        // update the RR sets
        for (auto RR_idx: this->_hyperGraphVldt._FRsets[node]) {
            if (!vec_bool_cov_by_S[RR_idx]) {
                if (dsfmt_gv_genrand_open_close() < p) {
                    vec_bool_cov_by_S[RR_idx] = true;
                    for (auto node: this->_hyperGraphVldt._RRsets[RR_idx]) {
                        if (!vec_bool_seed[node]) {
                            marginal_cov[node] -= 1;
                        }
                    }
                }
            }
        }
        if (!this->_candEdges[node].empty()) {
            heap.push(std::make_pair(node, marginal_cov[node] * this->_candEdges[node].back().second));
        }
    }
    
    _res.set_res_edges(_vecResEdges);
    _res.set_running_time(timerAIS.get_total_time());
    return;
}

void Alg::selectProb(const int targetSize) {
    Timer timerProb("prob");
    // build a heap of <node, max_prob>
    std::priority_queue<std::pair<uint32_t, double>, std::vector<std::pair<uint32_t, double>>, CompareBySecondDouble> heap;
    // Traverse the nodes and push the back of candedges into heap
    for (uint32_t i=0; i<this->_numV; i++) {
        if (!this->_candEdges[i].empty()) {
            heap.push(std::make_pair(i, this->_candEdges[i].back().second));
        }
    }
    // for k rounds, select the node with the max prob
    for (int i=0; i<targetSize; i++) {
        std::pair<uint32_t, double> top = heap.top();
        heap.pop();
        uint32_t node = top.first;
        double prob = top.second;
        this->_vecResEdges.push_back(std::make_pair(this->_candEdges[node].back().first, Edge(node, prob)));
        // update the heap
        this->_candEdges[node].pop_back();
        if (!this->_candEdges[node].empty()) {
            heap.push(std::make_pair(node, this->_candEdges[node].back().second));
        }
    } 
    _res.set_mode("total");
    _res.set_approximation(-1);
    _res.set_running_time(timerProb.get_total_time());
    // _res.set_influence(heap.top().second);
    // _res.set_influence_original(heap.top().second);
    // _res.set_RR_sets_size(_numRRsets);
    _res.set_res_edges(_vecResEdges);
}

void Alg::selectIM(const int targetSize, double epsilon, double delta) {
    this->subsimOnly(targetSize, epsilon, delta);
    // traverse vecSeed
    for (auto seed : this->_vecSeed) {
        if (this->_candEdges[seed].empty()) {
            continue;
        }
        auto& edge = this->_candEdges[seed].back();
        this->_vecResEdges.push_back(std::make_pair(edge.first, Edge(seed, edge.second)));
    }
    _res.set_res_edges(_vecResEdges);
    return;
}

/// Improve by 2hop
void Alg::compTwoHopInf() {
    // one hop influence vector
    std::vector<double> oneHopInf(this->_numV, 0.0);
    std::vector<bool> bool_seed(this->_numV, false);

    // probability map from a node to its in-seed
    std::vector<std::unordered_map<uint32_t, double>> probMap(this->_numV);

    for (auto seed : this->_vecSeed) {
        bool_seed[seed] = true;
    }
    // traverse every node to compute the one hop influence
    for (uint32_t i=0; i<this->_numV; i++) {
        for (auto& nbr : this->_hyperGraph._graph[i]) {
            if (bool_seed[nbr.first]) {
                oneHopInf[i] = oneHopInf[i] + (1 - oneHopInf[i]) * nbr.second;
                probMap[i][nbr.first] = nbr.second;
            }
        }
        if (bool_seed[i]) {
            oneHopInf[i] = this->_mapSeedProb[i] + (1-this->_mapSeedProb[i]) * oneHopInf[i];
        }
    }

    // Traverse every seed to compute the two hop influence
    for (auto seed : this->_vecSeed) {
        if (this->_mapSeedProb[seed] == 1) {
            this->_mapSeed2HopInf[seed] = 1.0;
            continue;
        }

        double twoHopInf = 0.0;
        for (auto& nbr : this->_hyperGraph._graph[seed]) {
            double conditional_prob = oneHopInf[nbr.first];
            if (probMap[nbr.first].find(seed) != probMap[nbr.first].end()) {
                double inv_prob = probMap[nbr.first][seed];
                double self_seed_prob = 0.0;
                if (bool_seed[nbr.first]) {
                    self_seed_prob = this->_mapSeedProb[nbr.first];
                }
                double ac_by_the_seed = this->_mapSeedProb[seed] * inv_prob;
                double tmp_prod = (1 - (oneHopInf[nbr.first] - self_seed_prob) / (1 - self_seed_prob)) / (1 - ac_by_the_seed);
                conditional_prob = oneHopInf[nbr.first] - (1-self_seed_prob)*ac_by_the_seed*tmp_prod;
            }
            twoHopInf += (1-twoHopInf) * conditional_prob * nbr.second;
        }
        
        this->_mapSeed2HopInf[seed] = this->_mapSeedProb[seed] + (1-this->_mapSeedProb[seed]) * twoHopInf;
        assert(this->_mapSeed2HopInf[seed] >= this->_mapSeedProb[seed] && this->_mapSeed2HopInf[seed] <= 1.0);
    }
    LogInfo("--- Finish two hop computation ---");
    return;
}

double Alg::sandwichP(const int targetSize, const double epsilon, const double delta, const double gamma) {
    double bound_max_time=0.0, estimate_time=0.0, twohop_time=0.0;
    Timer timerSandwich("SandwichP");
    
    this->compTwoHopInf();
    twohop_time = timerSandwich.get_operation_time();
    
    this->boundMaximize(targetSize, epsilon, delta, "UB");
    bound_max_time += timerSandwich.get_operation_time();

    candEdgesforLB(this->_candEdges, this->_mapSeed2HopInf);

    timerSandwich.get_operation_time();
    this->boundMaximize(targetSize, epsilon, delta, "LB");
    bound_max_time += timerSandwich.get_operation_time();

    this->_hyperGraphVldt.set_vanilla_sample(true);
    AddEdges(_hyperGraphVldt._graph, this->_vecResEdgesUB);
    double infUB = this->estimateInfStoppingRules(gamma, delta);
    PopEdges(_hyperGraphVldt._graph, this->_vecResEdgesUB);

    AddEdges(_hyperGraphVldt._graph, this->_vecResEdgesLB);
    double infLB = this->estimateInfStoppingRules(gamma, delta);
    PopEdges(_hyperGraphVldt._graph, this->_vecResEdgesLB);
    estimate_time = timerSandwich.get_operation_time();
    
    LogInfo("Influence of UB", infUB);
    LogInfo("Influence of LB", infLB);

    _res.set_mode("total");
    _res.set_influence_LB(infLB);
    _res.set_influence_UB(infUB);
    _res.set_estimate_time(estimate_time);
    _res.set_bound_max_time(bound_max_time);
    
    if (infUB < infLB) {
        _res.set_res_edges(_vecResEdgesLB);
    } else {
        _res.set_res_edges(_vecResEdgesUB);
    }
    _res.set_combined_res_two_bounds();
    _res.set_running_time(_res.get_running_time() + twohop_time);
    return 0.0;
}