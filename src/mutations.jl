"""
    add_mutations!(root::CoalescentNode, params::PopulationParams)

Add mutations to the coalescent tree according to the infinite sites model.

# Arguments
- `root::CoalescentNode`: Root of the coalescent tree
- `params::PopulationParams`: Population parameters including mutation rate and sequence length

# Details
Mutations are added to branches according to a Poisson process with rate
μ * L * t, where μ is the mutation rate, L is the sequence length, and t is the branch length.
"""
function add_mutations!(root::CoalescentNode, params::PopulationParams)
    # Use default RNG and create a reusable scratch buffer
    rng = Random.default_rng()
    scratch = Vector{Int}()
    _add_mutations_with_scratch!(root, params, rng, scratch)
    return nothing
end

"""
    add_mutations!(root::CoalescentNode, params::PopulationParams, rng::AbstractRNG)

Add mutations with explicit RNG for thread safety and reproducibility.
"""
function add_mutations!(root::CoalescentNode, params::PopulationParams, rng::AbstractRNG)
    scratch = Vector{Int}()
    _add_mutations_with_scratch!(root, params, rng, scratch)
    return nothing
end

"""
    add_mutations!(root::CoalescentNode, mu::Float64, sequence_length::Int; rng::AbstractRNG=Random.GLOBAL_RNG)

Add mutations with explicit mutation rate and sequence length parameters.
"""
function add_mutations!(root::CoalescentNode, mu::Float64, sequence_length::Int; rng::AbstractRNG=Random.GLOBAL_RNG)
    # Create temporary PopulationParams with the provided parameters
    temp_params = PopulationParams(10000, mu, 0.0, sequence_length, 100)  # Use dummy values for other fields
    scratch = Vector{Int}()
    _add_mutations_with_scratch!(root, temp_params, rng, scratch)
end

"""
    _add_mutations_with_scratch!(root::CoalescentNode, params::PopulationParams,
                                rng::AbstractRNG, scratch::Vector{Int})

Internal helper for mutation placement with caller-provided scratch buffer.
"""
function _add_mutations_with_scratch!(root::CoalescentNode, params::PopulationParams,
                                      rng::AbstractRNG, scratch::Vector{Int})
    _traverse_and_mutate!(root, params, rng, scratch)
    return nothing
end

"""
    _add_mutations_with_scratch!(root::CoalescentNode, mu::Float64, sequence_length::Int,
                                rng::AbstractRNG, scratch::Vector{Int})

Internal helper for mutation placement with explicit parameters and scratch buffer.
"""
function _add_mutations_with_scratch!(root::CoalescentNode, mu::Float64, sequence_length::Int,
                                      rng::AbstractRNG, scratch::Vector{Int})
    temp_params = PopulationParams(10000, mu, 0.0, sequence_length, 100)
    _add_mutations_with_scratch!(root, temp_params, rng, scratch)
    return nothing
end

"""
    _traverse_and_mutate!(node::CoalescentNode, params::PopulationParams, rng::AbstractRNG, scratch::Vector{Int})

Internal function to traverse tree and add mutations. Uses explicit RNG and reusable scratch buffer
for type stability and performance.
"""
function _traverse_and_mutate!(node::CoalescentNode, params::PopulationParams, rng::AbstractRNG, scratch::Vector{Int})
    if node.parent !== nothing
        # Branch length in generations (type-stable)
        branch_length::Float64 = node.parent.time - node.time
        
        # Expected number of mutations on this branch (type-stable)
        expected_mutations::Float64 = params.mutation_rate * params.sequence_length * branch_length
        
        # Sample actual number of mutations from Poisson distribution
        n_mutations::Int = rand(rng, Poisson(expected_mutations))
        
        # Add mutation positions (infinite sites model - no recurrent mutations)
        if n_mutations > 0
            # Clamp n_mutations to sequence_length to avoid sampling more positions than available
            n_to_sample::Int = min(n_mutations, params.sequence_length)
            
            # Pre-allocate capacity to avoid repeated growth
            mutations_vec = node.mutations
            sizehint!(mutations_vec, length(mutations_vec) + n_to_sample)
            
            if n_to_sample == params.sequence_length
                # If we need all positions, append them directly without creating intermediate vector
                for pos in 1:params.sequence_length
                    push!(mutations_vec, pos)
                end
            else
                # Use reusable scratch buffer for sampling
                resize!(scratch, n_to_sample)
                sample!(rng, 1:params.sequence_length, scratch; replace=false, ordered=true)
                append!(mutations_vec, scratch)
            end
        end
    end
    
    # Recurse to children with same RNG and scratch buffer
    for child in node.children
        _traverse_and_mutate!(child, params, rng, scratch)
    end
    return nothing
end

"""
    count_mutations(root::CoalescentNode) -> Int

Count the total number of mutations in the tree.
"""
function count_mutations(root::CoalescentNode)
    total = length(root.mutations)
    for child in root.children
        total += count_mutations(child)
    end
    return total
end

"""
    get_all_mutation_positions(root::CoalescentNode) -> Vector{Int}

Get all unique mutation positions in the tree, sorted.
"""
function get_all_mutation_positions(root::CoalescentNode)
    positions = Set{Int}()
    
    function collect_positions(node::CoalescentNode)
        union!(positions, node.mutations)
        for child in node.children
            collect_positions(child)
        end
    end
    
    collect_positions(root)
    return sort(collect(positions))
end

# ── Interval-aware mutation placement for ARG ─────────────────────────────

"""
    add_mutations_arg!(root, mu, seq_length, node_intervals; rng)

Place mutations only within each branch's ancestral material.
`node_intervals[id]` maps a node id to the genomic intervals it was ancestral
for at the time it coalesced (became a child).  Nodes missing from the dict
(e.g. the root, or nodes created by the standard coalescent) default to the
full genome.

Uses top-down traversal with effective interval tracking: mutations are only
placed within positions reachable from the root through interval-filtered paths.
This avoids placing mutations on "dead branches" that exist in the ARG tree
but aren't part of the correct marginal tree at those positions.
"""
function add_mutations_arg!(root::CoalescentNode, mu::Float64, sequence_length::Int,
                            node_intervals::Dict{Int, Vector{Tuple{Int,Int}}};
                            rng::AbstractRNG=Random.GLOBAL_RNG)
    # Top-down traversal: track effective intervals (intersection of all ancestor intervals)
    full_genome = [(1, sequence_length + 1)]
    # Stack: (node, effective_intervals_from_above)
    stack = Tuple{CoalescentNode, Vector{Tuple{Int,Int}}}[]
    push!(stack, (root, full_genome))

    while !isempty(stack)
        node, parent_effective = pop!(stack)

        if node.parent !== nothing
            branch_length::Float64 = node.parent.time - node.time

            # Effective intervals for THIS branch = intersection of parent's effective
            # intervals and this node's own node_intervals
            own_ivs = get(node_intervals, node.id, nothing)
            if own_ivs === nothing || isempty(own_ivs)
                effective = parent_effective
            else
                effective = _intersect_intervals_mut(parent_effective, own_ivs)
            end

            if !isempty(effective)
                material_length = sum(e - s for (s, e) in effective)
                expected_mutations::Float64 = mu * material_length * branch_length
                n_mutations::Int = rand(rng, Poisson(expected_mutations))

                if n_mutations > 0
                    n_to_sample = min(n_mutations, material_length)
                    mutations_vec = node.mutations
                    sizehint!(mutations_vec, length(mutations_vec) + n_to_sample)

                    if length(effective) == 1
                        s, e = effective[1]
                        buf = Vector{Int}(undef, n_to_sample)
                        sample!(rng, s:(e - 1), buf; replace=false, ordered=true)
                        append!(mutations_vec, buf)
                    else
                        lengths = [e - s for (s, e) in effective]
                        cum = cumsum(lengths)
                        total_len = cum[end]
                        result = Set{Int}()
                        attempts = 0
                        while length(result) < n_to_sample && attempts < 10 * n_to_sample
                            r = rand(rng, 1:total_len)
                            lo, hi = 1, length(cum)
                            while lo < hi
                                mid = (lo + hi) >>> 1
                                cum[mid] < r ? (lo = mid + 1) : (hi = mid)
                            end
                            s, e = effective[lo]
                            offset = r - (lo > 1 ? cum[lo - 1] : 0)
                            push!(result, s + offset - 1)
                            attempts += 1
                        end
                        append!(mutations_vec, sort!(collect(result)))
                    end
                end
            end

            # Pass effective intervals to children
            for child in node.children
                push!(stack, (child, isempty(effective) ? effective : effective))
            end
        else
            # Root node: pass parent_effective to children
            for child in node.children
                push!(stack, (child, parent_effective))
            end
        end
    end
    return nothing
end

"""Intersect two sorted interval lists [s, e) — used during mutation placement."""
function _intersect_intervals_mut(a::Vector{Tuple{Int,Int}}, b::Vector{Tuple{Int,Int}})
    result = Tuple{Int,Int}[]
    i, j = 1, 1
    while i <= length(a) && j <= length(b)
        s1, e1 = a[i]
        s2, e2 = b[j]
        lo = max(s1, s2)
        hi = min(e1, e2)
        if lo < hi
            push!(result, (lo, hi))
        end
        if e1 <= e2
            i += 1
        else
            j += 1
        end
    end
    return result
end