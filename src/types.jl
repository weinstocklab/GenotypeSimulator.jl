"""
    CoalescentNode

Simple coalescent tree structure representing a node in the genealogy.

# Fields
- `id::Int`: Unique identifier for the node
- `parent::Union{CoalescentNode, Nothing}`: Parent node (nothing for root)
- `children::Vector{CoalescentNode}`: Child nodes (empty for leaves)
- `time::Float64`: Time of the node (0.0 for present-day samples)
- `mutations::Vector{Int}`: Positions of mutations on the branch leading to this node
"""
mutable struct CoalescentNode
    id::Int
    parent::Union{CoalescentNode, Nothing}
    children::Vector{CoalescentNode}
    time::Float64
    mutations::Vector{Int}
    
    CoalescentNode(id::Int) = new(id, nothing, CoalescentNode[], 0.0, Int[])
    CoalescentNode(id::Int, parent, children, time, mutations) = 
        new(id, parent, children, time, mutations)
end

"""
    DemographicEpoch

Piecewise-constant effective population size on `[t_start, t_end)`, measured
backwards in generations.
"""
struct DemographicEpoch
    t_start::Float64
    t_end::Float64
    ne::Float64

    function DemographicEpoch(t_start, t_end, ne)
        t_start = Float64(t_start)
        t_end = Float64(t_end)
        ne = Float64(ne)
        @assert t_start >= 0 "Epoch start time must be non-negative"
        @assert t_end > t_start "Epoch end time must be greater than start time"
        @assert ne > 0 "Epoch effective population size must be positive"
        new(t_start, t_end, ne)
    end
end

"""
    DemographyModel

Ordered piecewise-constant demographic schedule for time-varying coalescent rates.
"""
struct DemographyModel
    epochs::Vector{DemographicEpoch}

    function DemographyModel(epochs::Vector{DemographicEpoch})
        @assert !isempty(epochs) "DemographyModel must contain at least one epoch"
        @assert epochs[1].t_start == 0.0 "DemographyModel must start at time 0"
        prev_end = epochs[1].t_start
        for epoch in epochs
            @assert epoch.t_start == prev_end "Demography epochs must be contiguous and sorted"
            prev_end = epoch.t_end
        end
        new(epochs)
    end
end

"""
    constant_demography(ne::Real) -> DemographyModel

Convenience constructor for a constant-size population.
"""
function constant_demography(ne::Real)
    return DemographyModel([DemographicEpoch(0.0, Inf, ne)])
end

"""
    recent_bottleneck_demography(ne_present::Real, ne_bottleneck::Real,
                                 t_start::Real, t_end::Real;
                                 ne_ancestral::Real=ne_present) -> DemographyModel

Three-epoch demography with a recent bottleneck on `[t_start, t_end)`.
"""
function recent_bottleneck_demography(ne_present::Real, ne_bottleneck::Real,
                                      t_start::Real, t_end::Real;
                                      ne_ancestral::Real=ne_present)
    t_start = Float64(t_start)
    t_end = Float64(t_end)
    @assert t_start > 0 "Bottleneck start time must be > 0"
    @assert t_end > t_start "Bottleneck end time must be greater than start time"

    return DemographyModel([
        DemographicEpoch(0.0, t_start, ne_present),
        DemographicEpoch(t_start, t_end, ne_bottleneck),
        DemographicEpoch(t_end, Inf, ne_ancestral),
    ])
end

"""
    PopulationParams

Parameters for human population genetic simulation.

# Fields
- `ne::Int`: Effective population size
- `mutation_rate::Float64`: Mutation rate per base pair per generation
- `recombination_rate::Float64`: Recombination rate per base pair per generation
- `sequence_length::Int`: Length of simulated sequence in base pairs
- `sample_size::Int`: Number of diploid individuals to sample
"""
struct PopulationParams
    ne::Int
    mutation_rate::Float64
    recombination_rate::Float64
    sequence_length::Int
    sample_size::Int
    
    function PopulationParams(ne, mutation_rate, recombination_rate, sequence_length, sample_size)
        @assert ne > 0 "Effective population size must be positive"
        @assert mutation_rate >= 0 "Mutation rate must be non-negative"
        @assert recombination_rate >= 0 "Recombination rate must be non-negative"
        @assert sequence_length > 0 "Sequence length must be positive"
        @assert sample_size > 0 "Sample size must be positive"
        
        new(ne, mutation_rate, recombination_rate, sequence_length, sample_size)
    end
end

"""
    HUMAN_PARAMS

Default parameters for human population genetic simulation of a 1Mb region.
- Effective population size: 10,000
- Mutation rate: 1.25e-8 per bp per generation
- Recombination rate: 1e-8 per bp per generation
- Sequence length: 1,000,000 bp (1Mb)
- Sample size: 100 individuals
"""
const HUMAN_PARAMS = PopulationParams(
    10_000,         # Ne = 10,000
    1.25e-8,        # mutation rate
    1e-8,           # recombination rate  
    1_000_000,      # 1Mb sequence
    100             # 100 individuals (200 haplotypes)
)
