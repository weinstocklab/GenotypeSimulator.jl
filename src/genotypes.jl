"""
    simulate_genotypes(params::PopulationParams=HUMAN_PARAMS; 
                      rng::AbstractRNG=Random.GLOBAL_RNG,
                      force_standard::Bool=false) -> (genotypes, positions)

Main simulation function that automatically chooses optimal algorithms based on simulation size.

For large simulations (>500 samples or >1MB sequences), automatically uses memory-optimized algorithms.
Use `force_standard=true` to disable automatic optimization.

# Arguments
- `params::PopulationParams`: Population parameters (default: HUMAN_PARAMS)
- `rng::AbstractRNG`: Random number generator (default: Random.GLOBAL_RNG)
- `force_standard::Bool`: Force use of standard algorithms (default: false)

# Returns
- `genotypes`: Genotype matrix (format chosen automatically based on data characteristics)
- `positions`: Vector of mutation positions
"""
function simulate_genotypes(params::PopulationParams=HUMAN_PARAMS; 
                           rng::AbstractRNG=Random.GLOBAL_RNG,
                           force_standard::Bool=false,
                           demography::Union{Nothing, DemographyModel}=nothing)
    
    # Check if we should use adaptive optimization
    if !force_standard
        # Thresholds for automatic optimization
        large_sample_threshold = 500
        large_sequence_threshold = 1_000_000
        
        if params.sample_size > large_sample_threshold || 
           params.sequence_length > large_sequence_threshold
            
            println("Large simulation detected - using adaptive memory optimization")
            return simulate_genotypes_adaptive(params; rng=rng, demography=demography)
        end
    end
    
    # Standard simulation for smaller problems
    println("Using standard simulation algorithms")
    
    # Build coalescent tree
    root = build_coalescent_tree(params; rng=rng, demography=demography)
    
    # Add mutations
    add_mutations!(root, params.mutation_rate, params.sequence_length; rng=rng)
    
    # Extract genotypes
    genotypes, positions = extract_genotypes(root, params.sample_size, params.sequence_length)
    
    return genotypes, positions
end

"""
    extract_genotypes(root::CoalescentNode, n_samples::Int, sequence_length::Int) -> (Matrix{Int}, Vector{Int})

Extract genotype matrix from the coalescent tree using iterative traversal for better performance.

# Arguments
- `root::CoalescentNode`: Root of the coalescent tree
- `n_samples::Int`: Number of diploid individuals
- `sequence_length::Int`: Length of the sequence

# Returns
- `Matrix{Int}`: Genotype matrix (rows = haplotypes, cols = variant sites)
- `Vector{Int}`: Positions of variant sites
"""
function extract_genotypes(root::CoalescentNode, n_samples::Int, sequence_length::Int)
    # Get all mutation positions using iterative approach
    positions = get_all_mutation_positions_iterative(root)
    n_sites = length(positions)
    
    if n_sites == 0
        return zeros(Int, 2 * n_samples, 0), Int[]
    end
    
    # Create genotype matrix (rows = haplotypes, cols = sites)
    genotypes = zeros(Int, 2 * n_samples, n_sites)
    
    # Create position lookup for O(1) access
    pos_to_idx = Dict{Int, Int}()
    for (idx, pos) in enumerate(positions)
        pos_to_idx[pos] = idx
    end
    
    # Iterative traversal to avoid stack overflow and improve performance
    _extract_genotypes_iterative!(root, genotypes, pos_to_idx, n_samples)
    
    return genotypes, positions
end

"""
    _extract_genotypes_iterative!(root::CoalescentNode, genotypes::Matrix{Int}, 
                                  pos_to_idx::Dict{Int,Int}, n_samples::Int)

Iterative implementation of genotype extraction to avoid recursion overhead.
"""
function _extract_genotypes_iterative!(root::CoalescentNode, genotypes::Matrix{Int}, 
                                      pos_to_idx::Dict{Int,Int}, n_samples::Int)
    # Stack for iterative traversal: (node, inherited_mutations_set)
    stack = Tuple{CoalescentNode, Set{Int}}[]
    push!(stack, (root, Set{Int}()))
    
    while !isempty(stack)
        node, inherited_mutations = pop!(stack)
        
        # Combine inherited mutations with mutations on this branch
        current_mutations = union(inherited_mutations, Set(node.mutations))
        
        if isempty(node.children)  # Leaf node
            hap_id = node.id
            if hap_id <= 2 * n_samples  # Ensure we don't exceed matrix bounds
                for pos in current_mutations
                    site_idx = get(pos_to_idx, pos, 0)
                    if site_idx > 0
                        genotypes[hap_id, site_idx] = 1
                    end
                end
            end
        else
            # Internal node - add children to stack
            for child in node.children
                push!(stack, (child, current_mutations))
            end
        end
    end
end

# ── Interval-aware genotype extraction for ARG ────────────────────────────

"""
    simulate_genotypes_marginal(root, recomb_events, node_intervals, lineage_samples,
                                n_samples, seq_len, mu; rng) -> (Matrix{Int}, Vector{Int})

Correct genotype simulation from an ARG using marginal trees.

For each genomic interval between consecutive recombination breakpoints, the
marginal tree is constant. We trace all sample lineages through the ARG (following
recombination redirects) to build the marginal tree, then place mutations only on
branches that carry a proper subset of samples (below the local MRCA).

This avoids the inflated-tree problem where per-branch mutation placement puts
mutations on branches above the local MRCA at each position.
"""
function simulate_genotypes_marginal(root::CoalescentNode,
                                     recomb_events::Vector{RecombinationEvent},
                                     node_intervals::Dict{Int, Vector{Tuple{Int,Int}}},
                                     lineage_samples::Dict{Int, Vector{Int}},
                                     n_samples::Int, seq_len::Int, mu::Float64;
                                     rng::AbstractRNG=Random.GLOBAL_RNG)
    n_haps = 2 * n_samples

    @inline function interval_contains(ivs::Vector{Tuple{Int,Int}}, pos::Int)
        if length(ivs) == 1
            s, e = ivs[1]
            return s <= pos < e
        end
        lo = 1
        hi = length(ivs)
        while lo <= hi
            mid = (lo + hi) >>> 1
            s, e = ivs[mid]
            if pos < s
                hi = mid - 1
            elseif pos >= e
                lo = mid + 1
            else
                return true
            end
        end
        return false
    end

    @inline function redirected_fragment(rdirs::Vector{Tuple{Int,Int}}, pos::Int)
        if length(rdirs) == 1
            bp, frag_id = rdirs[1]
            return bp <= pos ? frag_id : 0
        end
        lo = 1
        hi = length(rdirs)
        ans = 0
        while lo <= hi
            mid = (lo + hi) >>> 1
            bp, frag_id = rdirs[mid]
            if bp <= pos
                ans = frag_id
                lo = mid + 1
            else
                hi = mid - 1
            end
        end
        return ans
    end

    @inline function interval_length(ivs::Vector{Tuple{Int,Int}}, seq_len::Int)
        total = 0
        for (s, e) in ivs
            lo = max(s, 1)
            hi = min(e, seq_len + 1)
            hi > lo && (total += hi - lo)
        end
        return total
    end

    function sample_interval_position(ivs::Vector{Tuple{Int,Int}}, total_len::Int,
                                      seq_len::Int, rng::AbstractRNG)
        if length(ivs) == 1
            s, e = ivs[1]
            lo = max(s, 1)
            hi = min(e, seq_len + 1)
            return rand(rng, lo:(hi - 1))
        end
        r = rand(rng, 1:total_len)
        for (s, e) in ivs
            lo = max(s, 1)
            hi = min(e, seq_len + 1)
            len = hi - lo
            len <= 0 && continue
            if r <= len
                return lo + r - 1
            end
            r -= len
        end
        return seq_len
    end

    # Build dense node-id lookups. Node IDs are assigned sequentially by the ARG
    # builder, so arrays avoid hash-table work in the interval/haplotype loop.
    nodes_seen = CoalescentNode[]
    max_node_id = 0
    stk = CoalescentNode[root]
    while !isempty(stk)
        n = pop!(stk)
        push!(nodes_seen, n)
        max_node_id = max(max_node_id, n.id)
        for c in n.children
            push!(stk, c)
        end
    end

    node_by_id = Vector{CoalescentNode}(undef, max_node_id)
    parent_id = zeros(Int, max_node_id)
    branch_delta = zeros(Float64, max_node_id)
    node_interval_vec = [Tuple{Int,Int}[] for _ in 1:max_node_id]
    single_interval_start = zeros(Int, max_node_id)
    single_interval_end = zeros(Int, max_node_id)
    recomb_redirects = [Tuple{Int,Int}[] for _ in 1:max_node_id]
    single_redirect_bp = zeros(Int, max_node_id)
    single_redirect_frag = zeros(Int, max_node_id)
    lineage_sample_vec = [Int[] for _ in 1:max_node_id]

    for node in nodes_seen
        node_by_id[node.id] = node
        node.parent === nothing && continue
        parent_id[node.id] = node.parent.id
        branch_delta[node.id] = node.parent.time - node.time
    end

    for (nid, ivs) in node_intervals
        if 1 <= nid <= max_node_id
            if length(ivs) == 1
                single_interval_start[nid], single_interval_end[nid] = ivs[1]
            else
                node_interval_vec[nid] = ivs
            end
        end
    end

    for (nid, samples) in lineage_samples
        1 <= nid <= max_node_id && (lineage_sample_vec[nid] = samples)
    end

    # Build redirect map: lineage_id → [(breakpoint, fragment_id), ...] sorted ascending
    for evt in recomb_events
        if 1 <= evt.left_parent_id <= max_node_id
            push!(recomb_redirects[evt.left_parent_id], (evt.position, evt.right_parent_id))
        end
    end
    for r in recomb_redirects
        sort!(r, by=first)
    end
    for nid in 1:max_node_id
        r = recomb_redirects[nid]
        if length(r) == 1
            single_redirect_bp[nid], single_redirect_frag[nid] = r[1]
            empty!(r)
        end
    end

    all_mutations = Dict{Int, Vector{Int}}()  # position → list of haplotype ids
    full_genome = [(1, seq_len + 1)]
    carriers = Int[]

    function hap_carries_branch(hap::Int, target_id::Int, pos::Int)
        node_id = hap
        safety = 0
        @inbounds while parent_id[node_id] != 0
            safety += 1
            safety > 50000 && return false  # safety valve

            # Check redirects at CURRENT NODE first (handles chained redirects)
            frag_id = 0
            rbp = single_redirect_bp[node_id]
            if rbp != 0
                frag_id = rbp <= pos ? single_redirect_frag[node_id] : 0
            else
                rdirs = recomb_redirects[node_id]
                !isempty(rdirs) && (frag_id = redirected_fragment(rdirs, pos))
            end
            if frag_id != 0
                node_id = frag_id
                continue
            end

            s = single_interval_start[node_id]
            if s != 0
                if !(s <= pos < single_interval_end[node_id])
                    return false
                end
            else
                ivs = node_interval_vec[node_id]
                if !isempty(ivs) && !interval_contains(ivs, pos)
                    return false
                end
            end

            node_id == target_id && return true
            node_id = parent_id[node_id]
        end
        return false
    end

    function collect_branch_carriers!(carriers::Vector{Int}, target_id::Int, pos::Int)
        empty!(carriers)
        @inbounds target_samples = lineage_sample_vec[target_id]
        if length(target_samples) == 1
            push!(carriers, target_samples[1])
        elseif isempty(target_samples)
            for hap in 1:n_haps
                hap_carries_branch(hap, target_id, pos) && push!(carriers, hap)
            end
        else
            for hap in target_samples
                hap_carries_branch(hap, target_id, pos) && push!(carriers, hap)
            end
        end
        return length(carriers)
    end

    # Mutation-driven marginal simulation. Candidate mutations are sampled on
    # each ARG branch over its ancestral material, then thinned to branches that
    # carry a proper subset of samples in the local marginal tree at that site.
    for node in nodes_seen
        nid = node.id
        parent_id[nid] == 0 && continue
        bl = branch_delta[nid]
        bl <= 0 && continue

        s = single_interval_start[nid]
        e = single_interval_end[nid]
        ivs = node_interval_vec[nid]
        if s != 0
            lo = max(s, 1)
            hi = min(e, seq_len + 1)
            material_len = hi - lo
        else
            isempty(ivs) && (ivs = full_genome)
            material_len = interval_length(ivs, seq_len)
        end
        material_len <= 0 && continue

        expected = mu * material_len * bl
        n_mut = rand(rng, Poisson(expected))
        n_mut == 0 && continue

        for _ in 1:n_mut
            pos = s != 0 ? rand(rng, max(s, 1):(min(e, seq_len + 1) - 1)) :
                  sample_interval_position(ivs, material_len, seq_len, rng)
            ns = collect_branch_carriers!(carriers, nid, pos)
            if ns > 0 && ns < n_haps
                muts = get!(all_mutations, pos) do
                    Int[]
                end
                append!(muts, carriers)
            end
        end
    end

    # Build genotype matrix
    positions = sort!(collect(keys(all_mutations)))
    n_sites = length(positions)
    n_sites == 0 && return zeros(Int, n_haps, 0), Int[]

    genotypes = zeros(Int, n_haps, n_sites)
    for (j, pos) in enumerate(positions)
        for hid in all_mutations[pos]
            if 0 < hid <= n_haps
                @inbounds genotypes[hid, j] = 1
            end
        end
    end

    # Filter: keep only polymorphic sites
    poly = Int[]
    for j in 1:n_sites
        col_sum = 0
        @inbounds for i in 1:n_haps
            col_sum += genotypes[i, j]
        end
        if col_sum > 0 && col_sum < n_haps
            push!(poly, j)
        end
    end
    if length(poly) < n_sites
        genotypes = genotypes[:, poly]
        positions = positions[poly]
    end

    return genotypes, positions
end

"""
    extract_genotypes_arg(root, n_samples, seq_len, node_intervals, lineage_samples) -> (Matrix{Int}, Vector{Int})

Extract genotypes from an ARG using top-down traversal matching the top-down
mutation placement. Mutations flow from root to leaves, filtered by node_intervals
at each branch. Fragment leaves route mutations to their original samples via
lineage_samples.
"""
function extract_genotypes_arg(root::CoalescentNode, n_samples::Int, sequence_length::Int,
                               node_intervals::Dict{Int, Vector{Tuple{Int,Int}}},
                               lineage_samples::Dict{Int, Vector{Int}})
    positions = get_all_mutation_positions_iterative(root)
    n_sites = length(positions)
    n_sites == 0 && return zeros(Int, 2 * n_samples, 0), Int[]

    genotypes = zeros(Int, 2 * n_samples, n_sites)
    pos_to_idx = Dict{Int, Int}()
    for (idx, pos) in enumerate(positions)
        pos_to_idx[pos] = idx
    end

    n_haps = 2 * n_samples

    # Stack: (node, inherited_mutations_sorted)
    stack = Tuple{CoalescentNode, Vector{Int}}[]
    push!(stack, (root, Int[]))

    while !isempty(stack)
        node, inherited = pop!(stack)

        # Merge inherited with this node's own mutations (both sorted → O(m))
        own = sort!(copy(node.mutations))
        current = _merge_sorted_unique(inherited, own)

        if isempty(node.children)
            # Leaf node — determine which genotype row(s) to update
            hap_id = node.id
            if hap_id <= n_haps
                for pos in current
                    site_idx = get(pos_to_idx, pos, 0)
                    site_idx > 0 && (genotypes[hap_id, site_idx] = 1)
                end
            else
                # Fragment → route to all samples it represents
                samples = get(lineage_samples, hap_id, Int[])
                for sid in samples
                    if 0 < sid <= n_haps
                        for pos in current
                            site_idx = get(pos_to_idx, pos, 0)
                            site_idx > 0 && (genotypes[sid, site_idx] = 1)
                        end
                    end
                end
            end
        else
            for child in node.children
                child_ivs = get(node_intervals, child.id, nothing)
                if child_ivs === nothing
                    push!(stack, (child, current))
                else
                    push!(stack, (child, _filter_by_intervals(current, child_ivs)))
                end
            end
        end
    end

    # Keep only polymorphic columns, matching the marginal ARG extractor.
    poly = Int[]
    for j in 1:n_sites
        col_sum = 0
        @inbounds for i in 1:2*n_samples
            col_sum += genotypes[i, j]
        end
        if col_sum > 0 && col_sum < 2 * n_samples
            push!(poly, j)
        end
    end
    if length(poly) < n_sites
        genotypes = genotypes[:, poly]
        positions = positions[poly]
    end

    return genotypes, positions
end

"""Merge two sorted vectors into sorted unique result."""
function _merge_sorted_unique(a::Vector{Int}, b::Vector{Int})
    result = Vector{Int}(undef, length(a) + length(b))
    i, j, k = 1, 1, 0
    while i <= length(a) && j <= length(b)
        if a[i] < b[j]
            k += 1; result[k] = a[i]; i += 1
        elseif a[i] > b[j]
            k += 1; result[k] = b[j]; j += 1
        else
            k += 1; result[k] = a[i]; i += 1; j += 1
        end
    end
    while i <= length(a); k += 1; result[k] = a[i]; i += 1; end
    while j <= length(b); k += 1; result[k] = b[j]; j += 1; end
    resize!(result, k)
    return result
end

"""Keep only positions that fall within any of the sorted intervals [s, e)."""
function _filter_by_intervals(positions::Vector{Int}, intervals::Vector{Tuple{Int,Int}})
    isempty(intervals) && return Int[]
    result = Int[]
    iv_idx = 1
    for pos in positions
        while iv_idx <= length(intervals) && intervals[iv_idx][2] <= pos
            iv_idx += 1
        end
        iv_idx > length(intervals) && break
        s, e = intervals[iv_idx]
        s <= pos && pos < e && push!(result, pos)
    end
    return result
end

"""
    get_all_mutation_positions_iterative(root::CoalescentNode) -> Vector{Int}

Get all unique mutation positions using iterative traversal for better performance.
"""
function get_all_mutation_positions_iterative(root::CoalescentNode)
    positions = Set{Int}()
    stack = CoalescentNode[root]
    
    while !isempty(stack)
        node = pop!(stack)
        
        # Add mutations from this node
        union!(positions, node.mutations)
        
        # Add children to stack
        for child in node.children
            push!(stack, child)
        end
    end
    
    return sort(collect(positions))
end

"""
    haplotypes_to_diploid(haplotypes::Matrix{Int}) -> Matrix{Int}

Convert haplotype matrix to diploid genotype matrix with optimized memory access patterns.

# Arguments
- `haplotypes::Matrix{Int}`: Haplotype matrix (rows = haplotypes, cols = sites)

# Returns
- `Matrix{Int}`: Diploid genotype matrix (rows = individuals, cols = sites)
  Values are 0, 1, or 2 representing the number of alternate alleles
"""
function haplotypes_to_diploid(haplotypes::Matrix{Int})
    n_haps, n_sites = size(haplotypes)
    n_individuals = n_haps ÷ 2
    
    diploid_genotypes = zeros(Int, n_individuals, n_sites)
    
    # Optimize memory access pattern by iterating over sites first
    @inbounds for j in 1:n_sites
        for i in 1:n_individuals
            hap1_idx = 2 * i - 1
            hap2_idx = 2 * i
            diploid_genotypes[i, j] = haplotypes[hap1_idx, j] + haplotypes[hap2_idx, j]
        end
    end
    
    return diploid_genotypes
end



"""
    simulate_with_recombination(params::PopulationParams = HUMAN_PARAMS) -> (Matrix{Int}, Vector{Int})

Simulate genotypes using the coalescent with recombination (ARG) with optimized memory management.

# Arguments
- `params::PopulationParams`: Population parameters (defaults to HUMAN_PARAMS)

# Returns
- `Matrix{Int}`: Haplotype genotype matrix
- `Vector{Int}`: Positions of variant sites
"""

"""
Simulate genotypes with recombination using automatically selected algorithms.

For large simulations, automatically uses memory-optimized ARG construction.
Use `force_standard=true` to disable automatic optimization.
"""
function simulate_with_recombination(params::PopulationParams=HUMAN_PARAMS; 
                                   rng::AbstractRNG=Random.GLOBAL_RNG,
                                   force_standard::Bool=false,
                                   demography::Union{Nothing, DemographyModel}=nothing)
    
    # Check if we should use adaptive optimization
    if !force_standard
        large_sample_threshold = 500
        large_sequence_threshold = 1_000_000
        
        if params.sample_size > large_sample_threshold || 
           params.sequence_length > large_sequence_threshold
            
            println("Large ARG simulation detected - using adaptive memory optimization")
            return simulate_with_recombination_adaptive(params; rng=rng, demography=demography)
        end
    end
    
    # Standard ARG simulation
    println("Using standard ARG simulation algorithms")
    println("Simulating genotypes with recombination for $(params.sample_size) individuals...")
    println("Sequence length: $(params.sequence_length) bp")
    println("Mutation rate: $(params.mutation_rate)")
    println("Recombination rate: $(params.recombination_rate)")
    println("Effective population size: $(params.ne)")
    
    # Create recombination map
    println("Creating recombination map...")
    recomb_map = uniform_recombination_map(params.sequence_length, params.recombination_rate)
    
    # Build ARG (Ancestral Recombination Graph)
    println("Building ARG with recombination...")
    root, recomb_events, node_intervals, lineage_samples =
        build_arg_tree(params, recomb_map; rng=rng, demography=demography)
    
    # Simulate genotypes using marginal trees
    println("Placing mutations on marginal trees...")
    genotypes, positions = simulate_genotypes_marginal(
        root, recomb_events, node_intervals, lineage_samples,
        params.sample_size, params.sequence_length, params.mutation_rate; rng=rng)
    
    return genotypes, positions
end

# Legacy function for backward compatibility
function simulate_with_recombination_legacy(params::PopulationParams = HUMAN_PARAMS)
    println("Simulating genotypes with recombination for $(params.sample_size) individuals...")
    println("Sequence length: $(params.sequence_length) bp")
    println("Mutation rate: $(params.mutation_rate)")
    println("Recombination rate: $(params.recombination_rate)")
    println("Effective population size: $(params.ne)")
    
    # Create recombination map
    println("Creating recombination map...")
    recomb_map = uniform_recombination_map(params.sequence_length, params.recombination_rate)
    
    # Build ARG (Ancestral Recombination Graph)
    println("Building ARG with recombination...")
    root, recomb_events = build_arg_tree(params.sample_size, params.ne, 
                                        params.sequence_length, recomb_map)
    
    println("Generated $(length(recomb_events)) recombination events")
    
    # Get local trees for different genomic regions
    println("Extracting local trees...")
    local_trees = get_local_trees(root, recomb_events, params.sequence_length)
    println("Generated $(length(local_trees)) local trees")
    
    # Optimized approach: collect all mutations first, then build matrix once
    println("Adding mutations...")
    all_mutation_data = Tuple{Int, Int, Int}[]  # (position, haplotype_id, tree_index)
    
    for (tree_idx, (start_pos, end_pos, local_root)) in enumerate(local_trees)
        # Create temporary parameters for this region
        region_length = end_pos - start_pos
        if region_length <= 0
            continue
        end
        
        region_params = PopulationParams(
            params.ne,
            params.mutation_rate,
            params.recombination_rate,
            region_length,
            params.sample_size
        )
        
        # Add mutations to this local tree
        add_mutations!(local_root, region_params)
        
        # Extract mutation information without building full matrix
        _collect_mutations_from_tree!(local_root, start_pos - 1, all_mutation_data, params.sample_size)
    end
    
    # Sort mutations by position
    sort!(all_mutation_data, by=x -> x[1])
    
    # Build final genotype matrix efficiently
    if isempty(all_mutation_data)
        return zeros(Int, 2 * params.sample_size, 0), Int[]
    end
    
    println("Building final genotype matrix...")
    all_positions = unique([mut[1] for mut in all_mutation_data])
    n_sites = length(all_positions)
    all_genotypes = zeros(Int, 2 * params.sample_size, n_sites)
    
    # Create position lookup
    pos_to_idx = Dict{Int, Int}()
    for (idx, pos) in enumerate(all_positions)
        pos_to_idx[pos] = idx
    end
    
    # Fill genotype matrix
    for (pos, hap_id, tree_idx) in all_mutation_data
        site_idx = pos_to_idx[pos]
        if hap_id <= 2 * params.sample_size
            all_genotypes[hap_id, site_idx] = 1
        end
    end
    
    println("Simulation complete!")
    println("Generated $(length(all_positions)) variant sites")
    println("Genotype matrix: $(size(all_genotypes))")
    
    return all_genotypes, all_positions
end

"""
    _collect_mutations_from_tree!(root::CoalescentNode, position_offset::Int, 
                                 mutation_data::Vector{Tuple{Int,Int,Int}}, n_samples::Int)

Collect mutation information from a tree without building the full genotype matrix.
"""
function _collect_mutations_from_tree!(root::CoalescentNode, position_offset::Int,
                                      mutation_data::Vector{Tuple{Int,Int,Int}}, n_samples::Int)
    # Stack for iterative traversal: (node, inherited_mutations_set)
    stack = Tuple{CoalescentNode, Set{Int}}[]
    push!(stack, (root, Set{Int}()))
    
    while !isempty(stack)
        node, inherited_mutations = pop!(stack)
        
        # Combine inherited mutations with mutations on this branch
        current_mutations = union(inherited_mutations, Set(node.mutations))
        
        if isempty(node.children)  # Leaf node
            hap_id = node.id
            if hap_id <= 2 * n_samples
                # Record all mutations for this haplotype
                for pos in current_mutations
                    adjusted_pos = pos + position_offset
                    push!(mutation_data, (adjusted_pos, hap_id, 0))
                end
            end
        else
            # Internal node - add children to stack
            for child in node.children
                push!(stack, (child, current_mutations))
            end
        end
    end
end
