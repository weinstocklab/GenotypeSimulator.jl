"""
Memory-optimized genotype extraction functions

This module provides memory-efficient alternatives to the standard genotype extraction
that significantly reduce RAM usage during simulation.
"""

# Import base functions to extend them
using ..GenotypeSimulator: haplotypes_to_diploid, calculate_stats
using Statistics: cor

"""
    merge_sorted_mutations!(dest::Vector{UInt32}, src1::Vector{UInt32}, src2::Vector{Int})

Efficiently merge two sorted mutation vectors into dest. src1 is assumed sorted,
src2 (node.mutations) may need sorting. Result is stored in dest (sorted, unique values).
"""
function merge_sorted_mutations!(dest::Vector{UInt32}, src1::Vector{UInt32}, src2::Vector{Int})
    empty!(dest)
    
    # Handle edge cases
    if isempty(src1) && isempty(src2)
        return dest
    elseif isempty(src1)
        for val in src2
            push!(dest, UInt32(val))
        end
        sort!(dest)
        # Remove duplicates
        if length(dest) > 1
            j = 1
            for i in 2:length(dest)
                if dest[i] != dest[j]
                    j += 1
                    dest[j] = dest[i]
                end
            end
            resize!(dest, j)
        end
        return dest
    elseif isempty(src2)
        append!(dest, src1)
        return dest
    end
    
    # Pre-sort src2 if needed (node.mutations should already be sorted from add_mutations!)
    src2_sorted = issorted(src2)
    if !src2_sorted
        src2 = sort(src2)
    end
    
    # Two-pointer merge
    i, j = 1, 1
    n1, n2 = length(src1), length(src2)
    
    while i <= n1 && j <= n2
        val1 = src1[i]
        val2 = UInt32(src2[j])
        
        if val1 < val2
            push!(dest, val1)
            i += 1
        elseif val1 > val2
            push!(dest, val2)
            j += 1
        else  # Equal - skip duplicate
            push!(dest, val1)
            i += 1
            j += 1
        end
    end
    
    # Append remaining elements
    while i <= n1
        push!(dest, src1[i])
        i += 1
    end
    
    while j <= n2
        push!(dest, UInt32(src2[j]))
        j += 1
    end
    
    return dest
end

# Method overload for cor to handle SparseGenotypes
"""
    cor(sg::SparseGenotypes) -> Matrix{Float64}

Compute correlation matrix for SparseGenotypes by converting to dense matrix.
"""
function cor(sg::SparseGenotypes)
    # Convert to dense matrix for correlation calculation
    dense_matrix = Matrix{Float64}(undef, Int(sg.n_rows), Int(sg.n_cols))
    
    # Fill with zeros first
    fill!(dense_matrix, 0.0)
    
    # Fill in the non-zero values
    for j in 1:Int(sg.n_cols)
        start_idx = sg.col_pointers[j]
        end_idx = sg.col_pointers[j + 1] - 1
        
        for idx in start_idx:end_idx
            row = sg.row_indices[idx]
            dense_matrix[row, j] = Float64(sg.values[idx])
        end
    end
    
    return cor(dense_matrix)
end

"""
    extract_genotypes_optimized(root::CoalescentNode, n_samples::Int, sequence_length::Int) -> (Union{Matrix{UInt8}, SparseGenotypes, BitPackedGenotypes}, Vector{UInt32})

Memory-optimized version of extract_genotypes that:
1. Uses smaller data types (UInt32 instead of Int)
2. Reuses temporary objects from memory pool
3. Uses iterative instead of recursive traversal
4. Automatically chooses optimal storage format
"""
function extract_genotypes_optimized(root::CoalescentNode, n_samples::Int, sequence_length::Int)
    # Get all mutation positions
    positions = get_all_mutation_positions_optimized(root)
    n_sites = length(positions)
    
    if n_sites == 0
        return zeros(UInt8, 2 * n_samples, 0), UInt32[]
    end
    
    # Get memory pool
    pool = get_memory_pool(n_samples, n_sites)
    
    # Use pre-allocated buffer if it fits, otherwise allocate new
    if size(pool.genotype_buffer, 1) >= 2 * n_samples && size(pool.genotype_buffer, 2) >= n_sites
        genotypes = view(pool.genotype_buffer, 1:(2 * n_samples), 1:n_sites)
        fill!(genotypes, 0)
    else
        genotypes = zeros(UInt8, 2 * n_samples, n_sites)
    end
    
    # Create position lookup for faster access
    pos_to_idx = Dict{UInt32, Int}()
    for (idx, pos) in enumerate(positions)
        pos_to_idx[UInt32(pos)] = idx
    end
    
    # Iterative traversal to avoid stack overflow and reduce allocations
    extract_genotypes_iterative!(root, genotypes, pos_to_idx, n_samples, pool)
    
    # Convert to optimal representation
    if genotypes isa SubArray
        genotypes_copy = Matrix{UInt8}(genotypes)
        optimal_genotypes, optimal_positions = choose_optimal_representation(genotypes_copy, Int.(positions))
    else
        optimal_genotypes, optimal_positions = choose_optimal_representation(genotypes, Int.(positions))
    end
    
    return optimal_genotypes, UInt32.(optimal_positions)
end

"""
    extract_genotypes_iterative!(root::CoalescentNode, genotypes::AbstractMatrix{UInt8}, 
                                pos_to_idx::Dict{UInt32, Int}, n_samples::Int, pool::MemoryPool)

Iterative genotype extraction using sorted vectors instead of Sets for better performance.
Avoids expensive Set hashing operations by using two-pointer merge on sorted vectors.
"""
function extract_genotypes_iterative!(root::CoalescentNode, genotypes::AbstractMatrix{UInt8},
                                     pos_to_idx::Dict{UInt32, Int}, n_samples::Int, pool::MemoryPool)
    # Stack for iterative traversal: (node, inherited_mutations_vector)
    # Use preallocated vectors to avoid allocations
    stack = Tuple{CoalescentNode, Vector{UInt32}}[]
    mutation_buffers = Vector{UInt32}[UInt32[] for _ in 1:10]  # Small pool of buffers
    buffer_idx = 1
    
    # Start with root and empty mutation vector
    root_mutations = UInt32[]
    push!(stack, (root, root_mutations))
    
    while !isempty(stack)
        node, inherited_mutations = pop!(stack)
        
        # Get current mutations for this node (use merge with node.mutations)
        current_mutations = buffer_idx <= length(mutation_buffers) ? mutation_buffers[buffer_idx] : UInt32[]
        buffer_idx += 1
        
        merge_sorted_mutations!(current_mutations, inherited_mutations, node.mutations)
        
        if isempty(node.children)  # Leaf node
            hap_id = node.id
            if hap_id <= 2 * n_samples
                # Set genotypes for this haplotype
                for mut_pos in current_mutations
                    site_idx = get(pos_to_idx, mut_pos, 0)
                    if site_idx > 0
                        genotypes[hap_id, site_idx] = 1
                    end
                end
            end
            
            # Clear buffer for reuse
            empty!(current_mutations)
            buffer_idx -= 1
        else
            # Internal node - add children to stack with copies of current mutations
            for (i, child) in enumerate(node.children)
                if i == 1
                    # First child reuses current_mutations
                    push!(stack, (child, current_mutations))
                else
                    # Subsequent children need copies
                    child_mutations = buffer_idx <= length(mutation_buffers) ? mutation_buffers[buffer_idx] : UInt32[]
                    buffer_idx += 1
                    empty!(child_mutations)
                    append!(child_mutations, current_mutations)
                    push!(stack, (child, child_mutations))
                end
            end
        end
    end
end

"""
    get_all_mutation_positions_optimized(root::CoalescentNode) -> Vector{UInt32}

Memory-optimized version that uses UInt32 and avoids unnecessary allocations.
"""
function get_all_mutation_positions_optimized(root::CoalescentNode)
    positions = Set{UInt32}()
    
    # Iterative traversal
    stack = CoalescentNode[root]
    
    while !isempty(stack)
        node = pop!(stack)
        
        # Add mutations from this node
        for pos in node.mutations
            push!(positions, UInt32(pos))
        end
        
        # Add children to stack
        for child in node.children
            push!(stack, child)
        end
    end
    
    return sort!(collect(positions))
end

"""
    haplotypes_to_diploid_optimized(haplotypes::Union{Matrix{UInt8}, SparseGenotypes, BitPackedGenotypes}) -> Union{Matrix{UInt8}, SparseGenotypes}

Memory-optimized diploid conversion that handles different input formats efficiently.
"""
function haplotypes_to_diploid_optimized(haplotypes::Matrix{UInt8})
    n_haps, n_sites = size(haplotypes)
    n_individuals = n_haps ÷ 2
    
    diploid_genotypes = zeros(UInt8, n_individuals, n_sites)
    
    @inbounds for j in 1:n_sites
        for i in 1:n_individuals
            hap1_idx = 2 * i - 1
            hap2_idx = 2 * i
            diploid_genotypes[i, j] = haplotypes[hap1_idx, j] + haplotypes[hap2_idx, j]
        end
    end
    
    return diploid_genotypes
end

function haplotypes_to_diploid_optimized(haplotypes::SparseGenotypes)
    n_haps = Int(haplotypes.n_rows)
    n_sites = Int(haplotypes.n_cols)
    n_individuals = n_haps ÷ 2
    
    # For sparse haplotypes, result might also be sparse
    diploid_values = UInt8[]
    diploid_rows = UInt32[]
    diploid_cols = UInt32[]
    
    for j in 1:n_sites
        for i in 1:n_individuals
            hap1_idx = 2 * i - 1
            hap2_idx = 2 * i
            
            val1 = haplotypes[hap1_idx, j]
            val2 = haplotypes[hap2_idx, j]
            diploid_val = val1 + val2
            
            if diploid_val > 0
                push!(diploid_values, diploid_val)
                push!(diploid_rows, UInt32(i))
                push!(diploid_cols, UInt32(j))
            end
        end
    end
    
    # Convert to dense if not very sparse
    sparsity = 1.0 - length(diploid_values) / (n_individuals * n_sites)
    
    if sparsity < 0.7
        # Convert to dense
        diploid_dense = zeros(UInt8, n_individuals, n_sites)
        for (idx, (i, j)) in enumerate(zip(diploid_rows, diploid_cols))
            diploid_dense[i, j] = diploid_values[idx]
        end
        return diploid_dense
    else
        # Keep as sparse - create new SparseGenotypes directly
        # Convert to CSC format manually
        row_indices = Vector{UInt32}()
        col_pointers = Vector{UInt32}(undef, n_sites + 1)
        values_csc = Vector{UInt8}()
        
        col_pointers[1] = 1
        for j in 1:n_sites
            for (idx, (row, col)) in enumerate(zip(diploid_rows, diploid_cols))
                if col == j
                    push!(row_indices, row)
                    push!(values_csc, diploid_values[idx])
                end
            end
            col_pointers[j + 1] = length(row_indices) + 1
        end
        
        return SparseGenotypes(
            row_indices,
            col_pointers,
            values_csc,
            UInt32(n_individuals),
            UInt32(n_sites),
            UInt32.(collect(1:n_sites))
        )
    end
end

function haplotypes_to_diploid_optimized(haplotypes::BitPackedGenotypes)
    n_haps = Int(haplotypes.n_rows)
    n_sites = Int(haplotypes.n_cols)
    n_individuals = n_haps ÷ 2
    
    diploid_genotypes = zeros(UInt8, n_individuals, n_sites)
    
    @inbounds for j in 1:n_sites
        for i in 1:n_individuals
            hap1_idx = 2 * i - 1
            hap2_idx = 2 * i
            
            val1 = haplotypes[hap1_idx, j] ? 1 : 0
            val2 = haplotypes[hap2_idx, j] ? 1 : 0
            diploid_genotypes[i, j] = val1 + val2
        end
    end
    
    return diploid_genotypes
end

"""
    simulate_genotypes_memory_optimized(params::PopulationParams = HUMAN_PARAMS) -> (Union{Matrix{UInt8}, SparseGenotypes, BitPackedGenotypes}, Vector{UInt32})

Memory-optimized version of simulate_genotypes that uses all available optimizations.
"""
function simulate_genotypes_memory_optimized(params::PopulationParams = HUMAN_PARAMS)
    println("Simulating genotypes (memory-optimized) for $(params.sample_size) individuals...")
    println("Sequence length: $(params.sequence_length) bp")
    println("Mutation rate: $(params.mutation_rate)")
    println("Effective population size: $(params.ne)")
    
    # Build coalescent tree
    println("Building coalescent tree...")
    root = build_coalescent_tree(params.sample_size, params.ne)
    
    # Add mutations
    println("Adding mutations...")
    add_mutations!(root, params)
    
    # Extract genotypes with memory optimization
    println("Extracting genotypes (optimized)...")
    genotypes, positions = extract_genotypes_optimized(root, params.sample_size, params.sequence_length)
    
    println("Simulation complete!")
    println("Generated $(length(positions)) variant sites")
    
    # Report memory usage
    if genotypes isa SparseGenotypes
        println("Using sparse representation ($(round(sparsity_ratio(genotypes) * 100, digits=1))% sparse)")
        println("Memory usage: $(memory_usage(genotypes)) bytes")
    elseif genotypes isa BitPackedGenotypes
        println("Using bit-packed representation")
        println("Memory usage: $(sizeof(genotypes.data) + sizeof(genotypes.positions)) bytes")
    else
        println("Using dense representation")
        println("Genotype matrix: $(size(genotypes))")
    end
    
    return genotypes, positions
end

"""
    memory_benchmark(params::PopulationParams)

Compare memory usage between standard and optimized implementations.
"""
function memory_benchmark(params::PopulationParams)
    println("Memory Benchmark for $(params.sample_size) individuals, $(params.sequence_length) bp")
    println("=" * 70)
    
    # Benchmark standard implementation
    println("Standard implementation:")
    GC.gc()  # Force garbage collection
    memory_before = Base.gc_bytes()
    
    genotypes_std, positions_std = simulate_genotypes(params)
    genotypes_diploid_std = haplotypes_to_diploid(genotypes_std)
    
    GC.gc()
    memory_after = Base.gc_bytes()
    memory_std = memory_after - memory_before
    
    println("  Memory used: $(round(memory_std / 1024^2, digits=2)) MB")
    println("  Variants: $(length(positions_std))")
    println("  Matrix size: $(size(genotypes_diploid_std))")
    
    # Benchmark optimized implementation
    println("\nOptimized implementation:")
    GC.gc()
    memory_before = Base.gc_bytes()
    
    genotypes_opt, positions_opt = simulate_genotypes_memory_optimized(params)
    genotypes_diploid_opt = haplotypes_to_diploid_optimized(genotypes_opt)
    
    GC.gc()
    memory_after = Base.gc_bytes()
    memory_opt = memory_after - memory_before
    
    println("  Memory used: $(round(memory_opt / 1024^2, digits=2)) MB")
    println("  Variants: $(length(positions_opt))")
    
    # Calculate improvement
    if memory_std > 0
        improvement = (memory_std - memory_opt) / memory_std * 100
        println("\nMemory reduction: $(round(improvement, digits=1))%")
    end
    
    return (memory_std, memory_opt)
end

# Method overloads for calculate_stats removed — the base calculate_stats now
# accepts Matrix{<:Integer} and Vector{<:Integer} natively, avoiding type-conversion copies.

# Method overloads for haplotypes_to_diploid to handle memory-optimized types

"""
    haplotypes_to_diploid(haplotypes::Matrix{UInt8}) -> Matrix{UInt8}

Convert UInt8 haplotype matrix to diploid genotype matrix.
"""
function haplotypes_to_diploid(haplotypes::Matrix{UInt8})
    return haplotypes_to_diploid_optimized(haplotypes)
end

"""
    haplotypes_to_diploid(haplotypes::SparseGenotypes) -> Union{Matrix{UInt8}, SparseGenotypes}

Convert sparse haplotype matrix to diploid genotype matrix.
"""
function haplotypes_to_diploid(haplotypes::SparseGenotypes)
    return haplotypes_to_diploid_optimized(haplotypes)
end

"""
    haplotypes_to_diploid(haplotypes::BitPackedGenotypes) -> Matrix{UInt8}

Convert bit-packed haplotype matrix to diploid genotype matrix.
"""
function haplotypes_to_diploid(haplotypes::BitPackedGenotypes)
    return haplotypes_to_diploid_optimized(haplotypes)
end

# Additional method overloads for compatibility with standard Julia functions

"""
    cor(bg::BitPackedGenotypes) -> Matrix{Float64}

Compute correlation matrix for BitPackedGenotypes by converting to dense matrix.
"""
function cor(bg::BitPackedGenotypes)
    # Convert to dense matrix for correlation calculation
    dense_matrix = Matrix{Float64}(undef, Int(bg.n_rows), Int(bg.n_cols))
    
    for j in 1:Int(bg.n_cols)
        for i in 1:Int(bg.n_rows)
            dense_matrix[i, j] = bg[i, j] ? 1.0 : 0.0
        end
    end
    
    return cor(dense_matrix)
end

"""
    Matrix{T}(sg::SparseGenotypes) where T

Convert SparseGenotypes to a dense Matrix of type T.
"""
function Matrix{T}(sg::SparseGenotypes) where T
    dense_matrix = zeros(T, Int(sg.n_rows), Int(sg.n_cols))
    
    for j in 1:Int(sg.n_cols)
        start_idx = sg.col_pointers[j]
        end_idx = sg.col_pointers[j + 1] - 1
        
        for idx in start_idx:end_idx
            row = sg.row_indices[idx]
            dense_matrix[row, j] = T(sg.values[idx])
        end
    end
    
    return dense_matrix
end

"""
    Matrix{T}(bg::BitPackedGenotypes) where T

Convert BitPackedGenotypes to a dense Matrix of type T.
"""
function Matrix{T}(bg::BitPackedGenotypes) where T
    dense_matrix = Matrix{T}(undef, Int(bg.n_rows), Int(bg.n_cols))
    
    for j in 1:Int(bg.n_cols)
        for i in 1:Int(bg.n_rows)
            dense_matrix[i, j] = bg[i, j] ? T(1) : T(0)
        end
    end
    
    return dense_matrix
end