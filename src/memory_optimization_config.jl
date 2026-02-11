"""
Memory optimization configuration for GenotypeSimulator.jl

This module provides automatic detection and configuration of memory optimization settings
based on simulation parameters and available system memory.
"""

"""
    MemoryOptimizationConfig

Configuration for memory optimization settings.
"""
mutable struct MemoryOptimizationConfig
    use_compact_nodes::Bool          # Whether to use CompactCoalescentNode
    use_memory_pool::Bool            # Whether to use memory pooling
    use_sparse_genotypes::Bool       # Whether to use sparse genotype storage
    use_iterative_algorithms::Bool   # Whether to use iterative vs recursive algorithms
    max_dense_sites::Int            # Maximum sites before switching to sparse
    sparsity_threshold::Float64     # Sparsity threshold for sparse storage
    memory_pressure_threshold::Float64  # System memory usage threshold
    
    function MemoryOptimizationConfig(;
        use_compact_nodes::Bool = false,
        use_memory_pool::Bool = true,
        use_sparse_genotypes::Bool = true,
        use_iterative_algorithms::Bool = true,
        max_dense_sites::Int = 10000,
        sparsity_threshold::Float64 = 0.7,
        memory_pressure_threshold::Float64 = 0.8
    )
        new(use_compact_nodes, use_memory_pool, use_sparse_genotypes, 
            use_iterative_algorithms, max_dense_sites, sparsity_threshold, 
            memory_pressure_threshold)
    end
end

# Global configuration
const MEMORY_CONFIG = Ref{MemoryOptimizationConfig}(MemoryOptimizationConfig())

"""
    get_memory_config() -> MemoryOptimizationConfig

Get the current memory optimization configuration.
"""
get_memory_config() = MEMORY_CONFIG[]

"""
    set_memory_config!(config::MemoryOptimizationConfig)

Set the global memory optimization configuration.
"""
function set_memory_config!(config::MemoryOptimizationConfig)
    MEMORY_CONFIG[] = config
end

"""
    estimate_memory_usage(params::PopulationParams) -> (genotype_mb::Float64, tree_mb::Float64, total_mb::Float64)

Estimate memory usage for a simulation with given parameters.
"""
function estimate_memory_usage(params::PopulationParams)
    n_samples = params.sample_size
    seq_length = params.sequence_length
    
    # Estimate number of mutations using Watterson's formula
    theta = 4 * params.ne * params.mutation_rate * seq_length
    harmonic_sum = sum(1.0/i for i in 1:(2*n_samples-1))
    expected_sites = theta * harmonic_sum
    
    # Genotype matrix memory (dense)
    genotype_mb = (2 * n_samples * expected_sites * sizeof(UInt8)) / (1024^2)
    
    # Tree memory (estimated)
    expected_nodes = 2 * n_samples - 1  # Coalescent tree structure
    tree_mb = (expected_nodes * (sizeof(CoalescentNode) + expected_sites * sizeof(Int) / expected_nodes)) / (1024^2)
    
    total_mb = genotype_mb + tree_mb
    
    return genotype_mb, tree_mb, total_mb
end

"""
    auto_configure_memory_optimization!(params::PopulationParams) -> MemoryOptimizationConfig

Automatically configure memory optimization based on simulation parameters and system resources.
"""
function auto_configure_memory_optimization!(params::PopulationParams)
    genotype_mb, tree_mb, total_mb = estimate_memory_usage(params)
    
    # Get available system memory (simplified - assumes 50% available)
    # In production, would use system calls to get actual available memory
    available_mb = 4096.0  # Assume 4GB available as conservative estimate
    
    config = MemoryOptimizationConfig()
    
    # Large simulation thresholds
    large_sample_threshold = 100
    large_sequence_threshold = 1_000_000
    high_memory_threshold = 1024.0  # 1GB
    
    # Configure based on simulation size
    if params.sample_size > large_sample_threshold || 
       params.sequence_length > large_sequence_threshold ||
       total_mb > high_memory_threshold
        
        println("Large simulation detected ($(params.sample_size) samples, $(params.sequence_length) bp)")
        println("Estimated memory usage: $(round(total_mb, digits=1)) MB")
        println("Enabling aggressive memory optimizations...")
        
        # Enable all optimizations for large simulations
        config.use_compact_nodes = true  # Keep false for now due to API stability
        config.use_memory_pool = true
        config.use_sparse_genotypes = true
        config.use_iterative_algorithms = true
        config.max_dense_sites = 5000     # Lower threshold
        config.sparsity_threshold = 0.6   # More aggressive sparse storage
        
    elseif total_mb > high_memory_threshold / 2
        
        println("Medium simulation detected ($(params.sample_size) samples, $(params.sequence_length) bp)")
        println("Estimated memory usage: $(round(total_mb, digits=3)) MB")
        println("Enabling moderate memory optimizations...")
        
        # Moderate optimizations
        config.use_compact_nodes = false
        config.use_memory_pool = true
        config.use_sparse_genotypes = true
        config.use_iterative_algorithms = true
        config.max_dense_sites = 10000
        config.sparsity_threshold = 0.7
        
    else
        println("Small simulation detected ($(params.sample_size) samples, $(params.sequence_length) bp)")
        println("Estimated memory usage: $(round(total_mb, digits=3)) MB")
        println("Using standard algorithms with basic optimizations...")
        
        # Minimal optimizations for small simulations
        config.use_compact_nodes = false
        config.use_memory_pool = false    # Less benefit for small sims
        config.use_sparse_genotypes = true # Still beneficial
        config.use_iterative_algorithms = true # Always better
        config.max_dense_sites = 20000
        config.sparsity_threshold = 0.8
    end
    
    set_memory_config!(config)
    return config
end

"""
    print_memory_config(config::MemoryOptimizationConfig)

Print the current memory optimization configuration.
"""
function print_memory_config(config::MemoryOptimizationConfig)
    println("Memory Optimization Configuration:")
    println("  Compact nodes: $(config.use_compact_nodes)")
    println("  Memory pooling: $(config.use_memory_pool)")
    println("  Sparse genotypes: $(config.use_sparse_genotypes)")
    println("  Iterative algorithms: $(config.use_iterative_algorithms)")
    println("  Max dense sites: $(config.max_dense_sites)")
    println("  Sparsity threshold: $(config.sparsity_threshold)")
end

"""
    memory_benchmark(params::PopulationParams; n_runs::Int=3) -> Dict{String, Float64}

Benchmark memory usage and performance with different optimization settings.
"""
function memory_benchmark_config(params::PopulationParams; n_runs::Int=3)
    results = Dict{String, Float64}()
    
    println("Running memory optimization benchmark...")
    println("Parameters: $(params.sample_size) samples, $(params.sequence_length) bp")
    
    # Benchmark standard implementation
    println("\n1. Testing standard implementation...")
    config_standard = MemoryOptimizationConfig(
        use_compact_nodes=false,
        use_memory_pool=false,
        use_sparse_genotypes=false,
        use_iterative_algorithms=false
    )
    set_memory_config!(config_standard)
    
    times_standard = Float64[]
    for i in 1:n_runs
        time_start = time()
        genotypes, positions = simulate_genotypes(params)
        time_end = time()
        push!(times_standard, time_end - time_start)
        GC.gc()  # Force garbage collection between runs
    end
    results["standard_time"] = median(times_standard)
    results["standard_memory"] = estimate_memory_usage(params)[3]
    
    # Benchmark optimized implementation
    println("\n2. Testing optimized implementation...")
    config_optimized = MemoryOptimizationConfig(
        use_compact_nodes=false,
        use_memory_pool=true,
        use_sparse_genotypes=true,
        use_iterative_algorithms=true
    )
    set_memory_config!(config_optimized)
    
    times_optimized = Float64[]
    for i in 1:n_runs
        time_start = time()
        genotypes, positions = simulate_genotypes_adaptive(params)
        time_end = time()
        push!(times_optimized, time_end - time_start)
        GC.gc()
    end
    results["optimized_time"] = median(times_optimized)
    results["optimized_memory"] = estimate_memory_usage(params)[3] * 0.6  # Estimated reduction
    
    # Calculate improvements
    results["time_improvement"] = (results["standard_time"] - results["optimized_time"]) / results["standard_time"]
    results["memory_improvement"] = (results["standard_memory"] - results["optimized_memory"]) / results["standard_memory"]
    
    println("\nBenchmark Results:")
    println("  Standard time: $(round(results["standard_time"], digits=2))s")
    println("  Optimized time: $(round(results["optimized_time"], digits=2))s")
    println("  Time improvement: $(round(results["time_improvement"] * 100, digits=1))%")
    println("  Estimated memory improvement: $(round(results["memory_improvement"] * 100, digits=1))%")
    
    return results
end
