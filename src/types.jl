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