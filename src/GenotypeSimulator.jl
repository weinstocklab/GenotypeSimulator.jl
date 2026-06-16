module GenotypeSimulator

using Random
using Distributions
using StatsBase
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
import Statistics: cor, mean, std, var

# Export main types and functions
export PopulationParams, HUMAN_PARAMS, DemographicEpoch, DemographyModel
export constant_demography, recent_bottleneck_demography
export CoalescentNode
export simulate_genotypes, haplotypes_to_diploid, extract_genotypes, extract_genotypes_arg
export calculate_stats, save_genotypes, tajimas_d, print_stats, calculate_spectral_stats, SpectralStats
export simulate_coalescent_times, build_coalescent_tree
export RecombinationMap, uniform_recombination_map, periodic_hotspot_recombination_map
export build_arg_tree, simulate_with_recombination, simulate_genotypes_marginal
export SelectionModel, NeutralSelection, DirectionalSelection, BalancingSelection, BackgroundSelection
export SelectionParameters, simulate_with_selection, allele_frequency_trajectory, effective_population_size_with_selection
export count_nodes, tree_height, count_mutations, get_all_mutation_positions
export add_mutations!, add_mutations_arg!, sample_recombination_position, allele_frequency_spectrum
export simulate_genotypes_memory_optimized, haplotypes_to_diploid_optimized, memory_benchmark
export SparseGenotypes, BitPackedGenotypes, choose_optimal_representation
export clear_mutations!

# Export memory optimization functions
export MemoryOptimizationConfig, get_memory_config, set_memory_config!, auto_configure_memory_optimization!
export simulate_genotypes_adaptive, simulate_with_recombination_adaptive, simulate_with_selection_adaptive
export haplotypes_to_diploid_adaptive, memory_efficient_simulation_pipeline
export estimate_memory_usage, print_memory_config
export CompactCoalescentNode, CompactTree, MemoryPool, get_memory_pool

# Include submodules
include("types.jl")
include("coalescent.jl")
include("recombination.jl")
include("selection.jl")
include("mutations.jl")
include("genotypes.jl")
include("statistics.jl")
include("io.jl")

# Memory optimization modules
include("memory_optimized_types.jl")
include("memory_optimized_genotypes.jl")
include("memory_optimization_config.jl")
include("adaptive_simulation.jl")

end # module GenotypeSimulator
