#pragma once

class ResultInfo
{
private:
    double __RunningTime = -1.0, __Influence = -1.0, __InfluenceOriginal = -1.0, __Approx = -1.0;
    double __RunningTime_UB = -1.0, __Influence_UB = -1.0, __InfluenceOriginal_UB = -1.0, __UB_value, __Approx_UB = -1.0, __bound_max_time = -1.0, __estimate_time = -1.0;
    double __RunningTime_LB = -1.0, __Influence_LB = -1.0, __InfluenceOriginal_LB = -1.0, __LB_value, __Approx_LB = -1.0;
    double __samplingTime_UB = -1.0, __selectionTime_UB = -1.0, __samplingTime_LB = -1.0, __selectionTime_LB = -1.0;
    double __samplingTime = -1.0, __selectionTime = -1.0;
    int __SeedSize = 0;
    size_t __RRsetsSize = 0, __RRsetsSize_LB = 0, __RRsetsSize_UB = 0;
    Nodelist __VecSeed;
    std::string __mode = "total";
    std::vector<UVWEdge> __VecResEdges, __VecResEdgesLB, __VecResEdgesUB;
public:
    ResultInfo()
    {
    }

    ~ResultInfo()
    {
    }

    /// Get running time.
    double get_running_time() const
    {
        return __RunningTime;
    }

    double get_sampling_time() const
    {
        return __samplingTime;
    }

    double get_selection_time() const
    {
        return __selectionTime;
    }

    /// Get influence spread.
    double get_influence() const
    {
        return __Influence;
    }

    /// Get influence of LB
    double get_influence_LB() const
    {
        return __Influence_LB;
    }

    /// Get Influence of UB
    double get_influence_UB() const
    {
        return __Influence_UB;
    }

    double get_influence_original_LB() const
    {
        return __InfluenceOriginal_LB;
    }

    /// Get Influence of UB
    double get_influence_original_UB() const
    {
        return __InfluenceOriginal_UB;
    }

    /// Get LB value
    double get_LB_value() const {
        return __LB_value;
    }
    double get_UB_value() const {
        return __UB_value;
    }

    double get_estimate_time() const
    {
        return __estimate_time;
    }

    double get_bound_max_time() const
    {
        return __bound_max_time;
    }

    /// Get self-estimated influence spread.
    double get_influence_original() const
    {
        return __InfluenceOriginal;
    }

    /// Get approximation guarantee.
    double get_approximation() const
    {
        return __Approx;
    }

    /// Get approximation LB
    double get_approximation_LB() const
    {
        return __Approx_LB;
    }

    /// Get approximation UB
    double get_approximation_UB() const
    {
        return __Approx_UB;
    }

    /// Get seed sets.
    Nodelist get_seed_vec() const
    {
        return __VecSeed;
    }

    /// Get seed size.
    int get_seed_size() const
    {
        return __SeedSize;
    }

    /// Get the number of RR sets.
    size_t get_RRsets_size() const
    {
        return __RRsetsSize;
    }

    /// Get edge set
    std::vector<UVWEdge> get_res_edges() const
    {
        return __VecResEdges;
    }

    /// Set running time.
    void set_running_time(const double value)
    {
        __RunningTime = value;
    }
    
    /// Set running time LB
    void set_running_time_LB(const double value)
    {
        __RunningTime_LB = value;
    }

    void set_running_time_UB(const double value)
    {
        __RunningTime_UB = value;
    }

    /// set sampling time
    void set_sampling_time(const double value)
    {
        __samplingTime = value;
    }

    /// set sampling time LB
    void set_sampling_time_LB(const double value)
    {
        __samplingTime_LB = value;
    }

    /// set sampling time UB
    void set_sampling_time_UB(const double value)
    {
        __samplingTime_UB = value;
    }
    /// set selection time
    void set_selection_time(const double value)
    {
        // LB or UB
        if (__mode == "UB") {
            __selectionTime_UB = value;
        } else if (__mode == "LB") {
            __selectionTime_LB = value;
        } else {
            __selectionTime = value;
        }
    }
    /// set selection time UB
    void set_selection_time_UB(const double value)
    {
        __selectionTime_UB = value;
    }

    /// set selection time LB
    void set_selection_time_LB(const double value)
    {
        __selectionTime_LB = value;
    }

    /// Set influence spread.
    void set_influence(const double value)
    {
        // LB or UB
        __Influence = value;
    }

    void set_LB_value(const double value)
    {
        __LB_value = value;
    }

    void set_UB_value(const double value)
    {
        __UB_value = value;
    }

    /// Set influence of LB
    void set_influence_LB(const double value)
    {
        __Influence_LB = value;
    }
    void set_influence_UB(const double value)
    {
        __Influence_UB = value;
    }

    /// Set self-estimated influence spread
    void set_influence_original(const double value)
    {
        // LB or UB
        if (__mode == "UB") {
            __InfluenceOriginal_UB = value;
        } else if (__mode == "LB") {
            __InfluenceOriginal_LB = value;
        } else {
            __InfluenceOriginal = value;
        }
    }

    /// Set approximation guarantee.
    void set_approximation(const double value)
    {
        __Approx = value;
    }

    /// Set approximation guarantee LB
    void set_approximation_LB(const double value)
    {
        __Approx_LB = value;
    }

    /// Set approximation guarantee UB
    void set_approximation_UB(const double value)
    {
        __Approx_UB = value;
    }

    /// Set seed sets.
    void set_seed_vec(Nodelist& vecSeed)
    {
        __VecSeed = vecSeed;
        set_seed_size((int)vecSeed.size());
    }

    /// Set seed size.
    void set_seed_size(const int value)
    {
        __SeedSize = value;
    }

    /// Set the number of RR sets.
    void set_RR_sets_size(const size_t value)
    {
        // LB or UB
        if (__mode == "UB") {
            __RRsetsSize_UB = value;
        } else if (__mode == "LB") {
            __RRsetsSize_LB = value;
        } else {
            __RRsetsSize = value;
        }
    }

    /// set number of RR sets UB
    void set_RR_sets_size_UB(const size_t value)
    {
        __RRsetsSize_UB = value;
    }
    /// set number of RR sets LB
    void set_RR_sets_size_LB(const size_t value)
    {
        __RRsetsSize_LB = value;
    }

    /// Set the mode
    void set_mode(const std::string mode)
    {
        __mode = mode;
    }

    /// Set the edges
    void set_res_edges(std::vector<UVWEdge>& vecResEdges)
    {
        __VecResEdges = vecResEdges;
    }
    /// Set the edges of Lower Bound
    void set_edges_bounds(std::vector<UVWEdge>& vecResEdges)
    {
        if (__mode == "UB") {
            __VecResEdgesUB = vecResEdges;
        } else if (__mode == "LB") {
            __VecResEdgesLB = vecResEdges;
        }
    }

    void set_estimate_time(const double value)
    {
        __estimate_time = value;
    }

    void set_bound_max_time(const double value)
    {
        __bound_max_time = value;
    }

    void set_combined_res_two_bounds() {
        this->__Approx = std::min(this->__Approx_LB, this->__Approx_UB);
        this->__samplingTime = this->__samplingTime_LB + this->__samplingTime_UB;
        this->__selectionTime = this->__selectionTime_LB + this->__selectionTime_UB;
        this->__RRsetsSize = std::max(this->__RRsetsSize_LB, this->__RRsetsSize_UB);
        this->__RunningTime = this->__samplingTime + this->__selectionTime + this->__estimate_time;
    }
};

typedef ResultInfo TResult;
typedef std::shared_ptr<ResultInfo> PResult;
