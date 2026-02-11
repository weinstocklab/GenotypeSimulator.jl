"""
Adaptive simulation functions that automatically choose optimal algorithms and data structures
based on simulation parameters and memory configuration.
"""

"""
    simulate_genotypes_adaptive(params::PopulationParams=HUMAN_PARAMS; rng::AbstractRNG=Random.GLOBAL_RNG, 
                               selection_params::Union{SelectionParameters, Nothing}=nothing,
                               recombination_map::Union{RecombinationMap, Nothing}=nothing) -> (genotypes, positions)

Adaptive simulation that automatically chooses optimal algorithms and data structures.
"""
function simulate_genotypes_adaptive(params::PopulationParams=HUMAN_PARAMS; 
                                   rng::AbstractRNG=Random.GLOBAL_RNG,
                                   selection_params::Union{SelectionParameters, Nothing}=nothing,
                                   recombination_map::Union{RecombinationMap, Nothing}=nothing)
    
    # Auto-configure memory optimization
    config = auto_configure_memory_optimization!(params)
    
    # Build tree with optimal algorithm
    node_intervals = nothing  # only set for ARG/recombination
    lineage_origin = nothing
    recomb_events_list = nothing
    if selection_params !== nothing
        println("Building coalescent tree with selection...")
        root = build_coalescent_tree_with_selection(params, selection_params; rng=rng)
    elseif recombination_map !== nothing
        println("Building ARG with recombination...")
        root, recomb_events_list, node_intervals, lineage_origin = build_arg_tree(params, recombination_map; rng=rng)
    else
        println("Building standard coalescent tree...")
        root = build_coalescent_tree(params; rng=rng)
    end
    
    # Add mutations / place genotypes
    if recomb_events_list !== nothing && node_intervals !== nothing && lineage_origin !== nothing
        # Use marginal tree approach for ARG (correct mutation placement)
        println("Placing mutations on marginal trees...")
        genotypes, positions = simulate_genotypes_marginal(
            root, recomb_events_list, node_intervals, lineage_origin,
            params.sample_size, params.sequence_length, params.mutation_rate; rng=rng)
        println("Using marginal-tree ARG genotype simulation")
    else
        # Standard mutation placement for non-recombination trees
        println("Adding mutations...")
        if config.use_memory_pool
            add_mutations_pooled!(root, params.mutation_rate, params.sequence_length; rng=rng)
        else
            add_mutations!(root, params.mutation_rate, params.sequence_length; rng=rng)
        end
        
        # Extract genotypes with optimal representation
        println("Extracting genotypes...")
        if config.use_sparse_genotypes || config.use_iterative_algorithms
            genotypes, positions = extract_genotypes_optimized(root, params.sample_size, params.sequence_length)
            println("Using optimized genotype extraction")
        else
            genotypes, positions = extract_genotypes(root, params.sample_size)
            println("Using standard genotype extraction")
        end
    end
    
    # Print memory info
    if isa(genotypes, SparseGenotypes)
        sparsity = sparsity_ratio(genotypes)
        memory_mb = memory_usage(genotypes) / (1024^2)
        println("Sparse genotypes: $(round(sparsity*100, digits=1))% sparse, $(round(memory_mb, digits=1)) MB")
    elseif isa(genotypes, BitPackedGenotypes)
        memory_mb = sizeof(genotypes.data) / (1024^2)
        println("Bit-packed genotypes: $(round(memory_mb, digits=1)) MB")
    else
        memory_mb = sizeof(genotypes) / (1024^2)
        println("Dense genotypes: $(round(memory_mb, digits=1)) MB")
    end
    
    return genotypes, positions
end

"""
    add_mutations_pooled!(root::CoalescentNode, mu::Float64, sequence_length::Int; rng::AbstractRNG=Random.GLOBAL_RNG)

Memory-pooled version of add_mutations! that reuses scratch space.
"""
function add_mutations_pooled!(root::CoalescentNode, mu::Float64, sequence_length::Int; rng::AbstractRNG=Random.GLOBAL_RNG)
    config = get_memory_config()
    
    if config.use_memory_pool
        # Estimate total mutations for pool sizing
        expected_mutations = Int(ceil(4 * mu * sequence_length * 2 * root.time))
        pool = get_memory_pool(100, max(expected_mutations, 1000))  # Use pool
        
        # Use pooled scratch space for mutation collection
        scratch = get_temp_vector!(pool)
        try
            _add_mutations_with_scratch!(root, mu, sequence_length, rng, scratch)
        finally
            return_temp_vector!(pool, scratch)
        end
    else
        add_mutations!(root, mu, sequence_length; rng=rng)
    end
end

"""
    build_arg_tree_optimized(params::PopulationParams, recombination_map::RecombinationMap; rng::AbstractRNG=Random.GLOBAL_RNG) -> CoalescentNode

Memory-optimized ARG construction using the existing optimized LineageManager.
"""
function build_arg_tree_optimized(params::PopulationParams, recombination_map::RecombinationMap; rng::AbstractRNG=Random.GLOBAL_RNG)
    # Use the existing optimized ARG implementation from recombination.jl
    root, recomb_events, node_intervals, lineage_origin = build_arg_tree(params, recombination_map; rng=rng)
    return root, node_intervals, lineage_origin
end

"""
    simulate_with_recombination_adaptive(params::PopulationParams=HUMAN_PARAMS; rng::AbstractRNG=Random.GLOBAL_RNG) -> (genotypes, positions)

Adaptive ARG simulation with automatic optimization.
"""
function simulate_with_recombination_adaptive(params::PopulationParams=HUMAN_PARAMS; rng::AbstractRNG=Random.GLOBAL_RNG)
    recombination_map = uniform_recombination_map(params.sequence_length, params.recombination_rate)
    return simulate_genotypes_adaptive(params; rng=rng, recombination_map=recombination_map)
end

"""
    simulate_with_selection_adaptive(params::PopulationParams, selection_params::SelectionParameters; rng::AbstractRNG=Random.GLOBAL_RNG) -> (genotypes, positions)

Adaptive selection simulation with automatic optimization.
"""
function simulate_with_selection_adaptive(params::PopulationParams, selection_params::SelectionParameters; rng::AbstractRNG=Random.GLOBAL_RNG)
    return simulate_genotypes_adaptive(params; rng=rng, selection_params=selection_params)
end

"""
    haplotypes_to_diploid_adaptive(haplotypes, positions) -> (diploid_genotypes, positions)

Convert haplotypes to diploid format using optimal data structures.
"""
function haplotypes_to_diploid_adaptive(haplotypes, positions)
    config = get_memory_config()
    
    if isa(haplotypes, SparseGenotypes)
        return haplotypes_to_diploid_sparse(haplotypes), positions
    elseif isa(haplotypes, BitPackedGenotypes)
        return haplotypes_to_diploid_bitpacked(haplotypes), positions
    else
        # Use existing optimized implementation
        return haplotypes_to_diploid_optimized(haplotypes, positions)
    end
end

"""
    haplotypes_to_diploid_sparse(sparse_haplotypes::SparseGenotypes) -> SparseGenotypes

Convert sparse haplotypes to diploid format maintaining sparsity.
"""
function haplotypes_to_diploid_sparse(sparse_haplotypes::SparseGenotypes)
    n_haps, n_sites = size(sparse_haplotypes)
    n_individuals = n_haps ÷ 2
    
    # Create new sparse structure for diploid genotypes
    row_indices = UInt32[]
    col_pointers = Vector{UInt32}(undef, n_sites + 1)
    values = UInt8[]
    
    col_pointers[1] = 1
    
    for site in 1:n_sites
        start_idx = sparse_haplotypes.col_pointers[site]
        end_idx = sparse_haplotypes.col_pointers[site + 1] - 1
        
        # Collect all haplotypes with mutations at this site
        hap_indices = Set{UInt32}()
        for idx in start_idx:end_idx
            push!(hap_indices, sparse_haplotypes.row_indices[idx])
        end
        
        # Convert to diploid (pair consecutive haplotypes)
        for ind in 1:n_individuals
            hap1 = 2 * ind - 1
            hap2 = 2 * ind
            
            genotype = UInt8(0)
            if hap1 in hap_indices
                genotype += 1
            end
            if hap2 in hap_indices
                genotype += 1
            end
            
            if genotype > 0
                push!(row_indices, UInt32(ind))
                push!(values, genotype)
            end
        end
        
        col_pointers[site + 1] = length(row_indices) + 1
    end
    
    return SparseGenotypes(
        row_indices,
        col_pointers,
        values,
        UInt32(n_individuals),
        UInt32(n_sites),
        sparse_haplotypes.positions
    )
end

"""
    haplotypes_to_diploid_bitpacked(bit_haplotypes::BitPackedGenotypes) -> Matrix{UInt8}

Convert bit-packed haplotypes to dense diploid format.
"""
function haplotypes_to_diploid_bitpacked(bit_haplotypes::BitPackedGenotypes)
    n_haps, n_sites = size(bit_haplotypes)
    n_individuals = n_haps ÷ 2
    
    diploid = zeros(UInt8, n_individuals, n_sites)
    
    for site in 1:n_sites
        for ind in 1:n_individuals
            hap1 = 2 * ind - 1
            hap2 = 2 * ind
            
            genotype = UInt8(0)
            if bit_haplotypes[hap1, site]
                genotype += 1
            end
            if bit_haplotypes[hap2, site]
                genotype += 1
            end
            
            diploid[ind, site] = genotype
        end
    end
    
    return diploid
end

"""
    memory_efficient_simulation_pipeline(params::PopulationParams; 
                                        with_recombination::Bool=false,
                                        selection_params::Union{SelectionParameters, Nothing}=nothing,
                                        output_format::String="adaptive") -> NamedTuple

Complete memory-efficient simulation pipeline with automatic optimization.
"""
function memory_efficient_simulation_pipeline(params::PopulationParams;
                                             with_recombination::Bool=false,
                                             selection_params::Union{SelectionParameters, Nothing}=nothing,
                                             output_format::String="adaptive")
    
    println("="^60)
    println("Memory-Efficient Simulation Pipeline")
    println("="^60)
    
    # Auto-configure memory optimization
    config = auto_configure_memory_optimization!(params)
    print_memory_config(config)
    
    println("\nStarting simulation...")
    start_time = time()
    
    # Choose simulation type
    if selection_params !== nothing
        genotypes, positions = simulate_with_selection_adaptive(params, selection_params)
    elseif with_recombination
        genotypes, positions = simulate_with_recombination_adaptive(params)
    else
        genotypes, positions = simulate_genotypes_adaptive(params)
    end
    
    # Convert to diploid
    diploid_genotypes, diploid_positions = haplotypes_to_diploid_adaptive(genotypes, positions)
    
    # Calculate statistics
    println("\nCalculating population genetics statistics...")
    if isa(diploid_genotypes, AbstractMatrix)
        stats = calculate_stats(diploid_genotypes, diploid_positions)
    else
        # For sparse/bit-packed, convert to dense for stats calculation
        dense_genotypes = Matrix{UInt8}(undef, size(diploid_genotypes))
        for i in 1:size(diploid_genotypes, 1)
            for j in 1:size(diploid_genotypes, 2)
                dense_genotypes[i, j] = diploid_genotypes[i, j]
            end
        end
        stats = calculate_stats(dense_genotypes, diploid_positions)
    end
    
    end_time = time()
    total_time = end_time - start_time
    
    println("\nSimulation completed in $(round(total_time, digits=2)) seconds")
    
    return (
        haplotypes = genotypes,
        diploid = diploid_genotypes,
        positions = positions,
        statistics = stats,
        config = config,
        runtime = total_time
    )
end
