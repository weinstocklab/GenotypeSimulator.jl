"""
    simulate_coalescent_times(n::Int, ne::Int) -> Vector{Float64}

Simulate coalescent times for a sample of n haplotypes under the standard coalescent model.

# Arguments
- `n::Int`: Number of haplotypes (lineages)
- `ne::Int`: Effective population size

# Returns
- `Vector{Float64}`: Vector of coalescent times (in generations, going backwards in time)
"""
function effective_population_size(demography::Nothing, t::Float64, default_ne::Real)
    return Float64(default_ne)
end

function effective_population_size(demography::DemographyModel, t::Float64, default_ne::Real)
    @inbounds for epoch in demography.epochs
        if epoch.t_start <= t < epoch.t_end
            return epoch.ne
        end
    end
    return demography.epochs[end].ne
end

function next_demography_boundary(demography::Nothing, t::Float64)
    return Inf
end

function next_demography_boundary(demography::DemographyModel, t::Float64)
    @inbounds for epoch in demography.epochs
        if epoch.t_start <= t < epoch.t_end
            return epoch.t_end
        end
    end
    return Inf
end

function sample_piecewise_wait(base_rate::Float64, coal_scale::Float64,
                               current_time::Float64, ne::Real,
                               demography::Union{Nothing, DemographyModel};
                               rng::AbstractRNG=Random.GLOBAL_RNG)
    target_hazard = randexp(rng)
    time = current_time

    while true
        ne_now = effective_population_size(demography, time, ne)
        total_rate = base_rate + coal_scale / ne_now
        total_rate > 0 || error("Total event rate must be positive")
        boundary = next_demography_boundary(demography, time)
        if !isfinite(boundary)
            return target_hazard / total_rate
        end

        interval_len = boundary - time
        if interval_len <= 0
            time = nextfloat(time)
            continue
        end

        interval_hazard = total_rate * interval_len
        if target_hazard < interval_hazard
            return target_hazard / total_rate
        end

        target_hazard -= interval_hazard
        time = boundary
    end
end

function simulate_coalescent_times(n::Int, ne::Int; rng::AbstractRNG=Random.GLOBAL_RNG,
                                   demography::Union{Nothing, DemographyModel}=nothing)
    times = Float64[]
    current_n = n
    current_time = 0.0

    while current_n > 1
        # Rate of coalescence for k lineages: k(k-1)/(4Ne(t))
        coal_scale = current_n * (current_n - 1) / 4.0
        dt = sample_piecewise_wait(0.0, coal_scale, current_time, ne, demography; rng=rng)
        current_time += dt
        push!(times, current_time)

        current_n -= 1
    end

    return times
end

"""
    build_coalescent_tree(sample_size::Int, ne::Int) -> CoalescentNode

Build a coalescent tree for the given sample size and effective population size.

# Arguments
- `sample_size::Int`: Number of diploid individuals (2 * sample_size haplotypes)
- `ne::Int`: Effective population size

# Returns
- `CoalescentNode`: Root node of the coalescent tree
"""
function build_coalescent_tree(sample_size::Int, ne::Int; rng::AbstractRNG=Random.GLOBAL_RNG,
                               demography::Union{Nothing, DemographyModel}=nothing)
    n_haplotypes = 2 * sample_size
    times = simulate_coalescent_times(n_haplotypes, ne; rng=rng, demography=demography)

    # Create leaf nodes (present-day samples)
    nodes = [CoalescentNode(i) for i in 1:n_haplotypes]

    # Build internal nodes by coalescence events
    node_id = n_haplotypes + 1
    active_nodes = copy(nodes)

    for (i, time) in enumerate(times)
        # Randomly choose two nodes to coalesce
        if length(active_nodes) < 2
            break
        end

        idx1, idx2 = sample(rng, 1:length(active_nodes), 2, replace=false)
        child1 = active_nodes[idx1]
        child2 = active_nodes[idx2]

        # Create parent node
        parent = CoalescentNode(node_id, nothing, [child1, child2], time, Int[])
        child1.parent = parent
        child2.parent = parent

        # Update active nodes list - remove higher index first to avoid shifting
        indices_to_remove = sort([idx1, idx2], rev=true)
        for idx in indices_to_remove
            deleteat!(active_nodes, idx)
        end
        push!(active_nodes, parent)

        node_id += 1
    end

    return active_nodes[1]  # Return root node
end

"""
    build_coalescent_tree(params::PopulationParams; rng::AbstractRNG=Random.GLOBAL_RNG) -> CoalescentNode

Build a coalescent tree from population parameters.
"""
function build_coalescent_tree(params::PopulationParams; rng::AbstractRNG=Random.GLOBAL_RNG,
                               demography::Union{Nothing, DemographyModel}=nothing)
    return build_coalescent_tree(params.sample_size, params.ne; rng=rng, demography=demography)
end

"""
    tree_height(root::CoalescentNode) -> Float64

Calculate the height (time to most recent common ancestor) of a coalescent tree.
"""
function tree_height(root::CoalescentNode)
    return root.time
end

"""
    count_nodes(root::CoalescentNode) -> Int

Count the total number of nodes in the tree.
"""
function count_nodes(root::CoalescentNode)
    count = 1
    for child in root.children
        count += count_nodes(child)
    end
    return count
end
