"""
Selection modeling for coalescent simulation

This module implements various forms of natural selection including:
- Directional selection (selective sweeps)
- Balancing selection
- Background selection
- Frequency-dependent selection

Following approaches from Kaplan et al. (1989), Barton (1998), and Etheridge et al. (2006).
"""

"""
    SelectionModel

Abstract type for different selection models.
"""
abstract type SelectionModel end

"""
    NeutralSelection

No selection - standard neutral coalescent.
"""
struct NeutralSelection <: SelectionModel end

"""
    DirectionalSelection

Directional selection model for selective sweeps.

# Fields
- `selection_coefficient::Float64`: Selection coefficient s (s > 0 for beneficial mutations)
- `dominance::Float64`: Dominance coefficient h (0 ≤ h ≤ 1)
- `selected_position::Int`: Genomic position under selection
- `selection_start_time::Float64`: Time when selection begins (generations ago)
- `initial_frequency::Float64`: Initial frequency of beneficial allele
"""
struct DirectionalSelection <: SelectionModel
    selection_coefficient::Float64
    dominance::Float64
    selected_position::Int
    selection_start_time::Float64
    initial_frequency::Float64

    function DirectionalSelection(s, h, pos, start_time, init_freq=1 / (2 * 10000))
        @assert s > 0 "Selection coefficient must be positive for directional selection"
        @assert 0 <= h <= 1 "Dominance coefficient must be between 0 and 1"
        @assert pos > 0 "Selected position must be positive"
        @assert start_time >= 0 "Selection start time must be non-negative"
        @assert 0 < init_freq < 1 "Initial frequency must be between 0 and 1"
        new(s, h, pos, start_time, init_freq)
    end
end

"""
    BalancingSelection

Balancing selection maintaining multiple alleles.

# Fields
- `selection_coefficients::Vector{Float64}`: Selection coefficients for each allele
- `selected_position::Int`: Genomic position under selection
- `equilibrium_frequencies::Vector{Float64}`: Equilibrium allele frequencies
"""
struct BalancingSelection <: SelectionModel
    selection_coefficients::Vector{Float64}
    selected_position::Int
    equilibrium_frequencies::Vector{Float64}

    function BalancingSelection(s_vec, pos, eq_freqs)
        @assert length(s_vec) == length(eq_freqs) "Selection coefficients and frequencies must have same length"
        @assert all(s -> s >= 0, s_vec) "Selection coefficients must be non-negative"
        @assert abs(sum(eq_freqs) - 1.0) < 1e-10 "Equilibrium frequencies must sum to 1"
        @assert all(f -> 0 < f < 1, eq_freqs) "Equilibrium frequencies must be between 0 and 1"
        new(s_vec, pos, eq_freqs)
    end
end

"""
    BackgroundSelection

Background selection against deleterious mutations.

# Fields
- `deleterious_rate::Float64`: Rate of deleterious mutations per bp per generation
- `selection_coefficient::Float64`: Selection coefficient against deleterious alleles (s < 0)
- `dominance::Float64`: Dominance coefficient for deleterious mutations
"""
struct BackgroundSelection <: SelectionModel
    deleterious_rate::Float64
    selection_coefficient::Float64
    dominance::Float64

    function BackgroundSelection(del_rate, s, h=0.0)
        @assert del_rate >= 0 "Deleterious mutation rate must be non-negative"
        @assert s < 0 "Selection coefficient must be negative for deleterious mutations"
        @assert 0 <= h <= 1 "Dominance coefficient must be between 0 and 1"
        new(del_rate, s, h)
    end
end

"""
    SelectionParameters

Container for selection parameters in population simulation.
"""
struct SelectionParameters
    model::SelectionModel
    effective_population_size::Int

    SelectionParameters(model::SelectionModel, ne::Int) = new(model, ne)
end

"""
    allele_frequency_trajectory(selection::DirectionalSelection, ne::Int, 
                               time_points::Vector{Float64}) -> Vector{Float64}

Calculate the deterministic allele frequency trajectory under directional selection.
Uses the logistic growth model for beneficial alleles.
"""
function allele_frequency_trajectory(selection::DirectionalSelection, ne::Int,
    time_points::Vector{Float64})
    s = selection.selection_coefficient
    h = selection.dominance
    p0 = selection.initial_frequency

    frequencies = Float64[]

    for t in time_points
        if t < selection.selection_start_time
            # Before selection starts
            push!(frequencies, p0)
        else
            # During selection phase
            τ = t - selection.selection_start_time

            # Effective selection coefficient
            s_eff = s * (h + (1 - 2 * h) * p0)

            # Logistic growth solution
            if abs(s_eff) < 1e-10
                # Neutral case
                p_t = p0
            else
                exp_term = exp(s_eff * τ)
                p_t = (p0 * exp_term) / (1 - p0 + p0 * exp_term)
            end

            push!(frequencies, min(max(p_t, 1e-10), 1.0 - 1e-10))
        end
    end

    return frequencies
end

"""
    effective_population_size_with_selection(selection::SelectionModel, ne_neutral::Int, 
                                            position::Int) -> Float64

Calculate the effective population size at a given position under selection.
"""
function effective_population_size_with_selection(selection::SelectionModel, ne_neutral::Int,
    position::Int)
    if isa(selection, NeutralSelection)
        return Float64(ne_neutral)
    elseif isa(selection, DirectionalSelection)
        # Reduction in Ne near selected site (Kaplan et al., 1989)
        distance = abs(position - selection.selected_position)
        s = selection.selection_coefficient

        # Approximate reduction based on linkage
        reduction_factor = 1.0 / (1.0 + s * exp(-distance / 1000.0))
        return ne_neutral * reduction_factor
    elseif isa(selection, BackgroundSelection)
        # Background selection reduces effective population size
        U = selection.deleterious_rate * 1000  # Assume 1kb region
        s = abs(selection.selection_coefficient)
        h = selection.dominance

        # Charlesworth (1994) approximation
        reduction = exp(-U * (h + (1 - h) * s) / s)
        return ne_neutral * reduction
    else
        return Float64(ne_neutral)
    end
end

"""
    selection_coefficient_at_position(selection::SelectionModel, position::Int) -> Float64

Get the selection coefficient at a specific genomic position.
"""
function selection_coefficient_at_position(selection::SelectionModel, position::Int)
    if isa(selection, NeutralSelection)
        return 0.0
    elseif isa(selection, DirectionalSelection)
        return position == selection.selected_position ? selection.selection_coefficient : 0.0
    elseif isa(selection, BalancingSelection)
        return position == selection.selected_position ? selection.selection_coefficients[1] : 0.0
    elseif isa(selection, BackgroundSelection)
        return selection.selection_coefficient
    else
        return 0.0
    end
end

"""
    modify_coalescent_rates!(active_lineages::Vector{CoalescentNode}, 
                            selection_params::SelectionParameters,
                            current_time::Float64)

Modify coalescent rates based on selection effects.
"""
function modify_coalescent_rates!(active_lineages::Vector{CoalescentNode},
    selection_params::SelectionParameters,
    current_time::Float64)
    selection = selection_params.model
    ne_base = selection_params.effective_population_size

    if isa(selection, NeutralSelection)
        return Float64(ne_base)
    end

    # For directional selection, modify effective population size based on 
    # allele frequency trajectory
    if isa(selection, DirectionalSelection)
        freq_traj = allele_frequency_trajectory(selection, ne_base, [current_time])
        p = freq_traj[1]

        # Effective population size varies with allele frequency
        # During sweep: Ne_eff ≈ Ne / (4 * p * (1-p) * s)
        if p > 1e-6 && p < 1.0 - 1e-6
            s = selection.selection_coefficient
            ne_eff = ne_base / (1.0 + 4.0 * p * (1 - p) * s)
            return max(ne_eff, ne_base * 0.01)  # Minimum 1% of original Ne
        end
    end

    return Float64(ne_base)
end

"""
    add_selected_mutations!(root::CoalescentNode, params::PopulationParams,
                           selection_params::SelectionParameters)

Add mutations with selection effects to the coalescent tree.
"""
function add_selected_mutations!(root::CoalescentNode, params::PopulationParams,
    selection_params::SelectionParameters)
    selection = selection_params.model

    # First add neutral mutations as usual
    add_mutations!(root, params)

    # Then add selected mutations if applicable
    if isa(selection, DirectionalSelection)
        add_beneficial_mutation!(root, selection, params)
    elseif isa(selection, BackgroundSelection)
        add_deleterious_mutations!(root, selection, params)
    end
end

"""
    add_beneficial_mutation!(root::CoalescentNode, selection::DirectionalSelection,
                            params::PopulationParams)

Add a beneficial mutation at the selected position.
"""
function add_beneficial_mutation!(root::CoalescentNode, selection::DirectionalSelection,
    params::PopulationParams)
    # Find the appropriate time and lineage for the beneficial mutation
    # This is a simplified implementation - in reality, we'd need to track
    # the frequency trajectory more carefully

    selected_pos = selection.selected_position

    # Add the beneficial mutation to a random lineage at the selection start time
    function find_lineage_at_time(node::CoalescentNode, target_time::Float64)
        if node.time <= target_time
            if isempty(node.children)
                # Leaf node - this lineage exists at target time
                return node
            else
                # Internal node - check children
                for child in node.children
                    result = find_lineage_at_time(child, target_time)
                    if result !== nothing
                        return result
                    end
                end
            end
        end
        return nothing
    end

    # Find a lineage that exists at the selection start time
    target_lineage = find_lineage_at_time(root, selection.selection_start_time)

    if target_lineage !== nothing
        # Add the beneficial mutation
        if selected_pos ∉ target_lineage.mutations
            push!(target_lineage.mutations, selected_pos)
            sort!(target_lineage.mutations)
        end
    end
end

"""
    add_deleterious_mutations!(root::CoalescentNode, selection::BackgroundSelection,
                              params::PopulationParams)

Add deleterious mutations throughout the sequence.
"""
function add_deleterious_mutations!(root::CoalescentNode, selection::BackgroundSelection,
    params::PopulationParams)
    del_rate = selection.deleterious_rate

    function traverse_and_add_deleterious(node::CoalescentNode)
        if node.parent !== nothing
            branch_length = node.parent.time - node.time
            expected_del_mutations = del_rate * params.sequence_length * branch_length
            n_del_mutations = rand(Poisson(expected_del_mutations))

            if n_del_mutations > 0
                # Add deleterious mutations (these would be filtered out in reality)
                # For simulation purposes, we track them but they don't contribute to output
                del_positions = sample(1:params.sequence_length,
                    min(n_del_mutations, params.sequence_length),
                    replace=false, ordered=true)
                # Store in a separate field or mark them somehow
                # For now, we'll just reduce the effective mutation rate
            end
        end

        for child in node.children
            traverse_and_add_deleterious(child)
        end
    end

    traverse_and_add_deleterious(root)
end

"""
    simulate_with_selection(params::PopulationParams, 
                           selection_params::SelectionParameters) -> (Matrix{Int}, Vector{Int})

Main function to simulate genotypes with selection effects.
"""
function simulate_with_selection(params::PopulationParams,
    selection_params::SelectionParameters)
    selection = selection_params.model

    println("Simulating genotypes with selection for $(params.sample_size) individuals...")
    println("Selection model: $(typeof(selection))")
    println("Sequence length: $(params.sequence_length) bp")
    println("Mutation rate: $(params.mutation_rate)")
    println("Effective population size: $(params.ne)")

    if isa(selection, DirectionalSelection)
        println("Selection coefficient: $(selection.selection_coefficient)")
        println("Selected position: $(selection.selected_position)")
        println("Selection start time: $(selection.selection_start_time) generations ago")
    elseif isa(selection, BackgroundSelection)
        println("Deleterious mutation rate: $(selection.deleterious_rate)")
        println("Selection against deleterious: $(selection.selection_coefficient)")
    end

    # Build coalescent tree with modified rates
    println("Building coalescent tree with selection...")
    root = build_coalescent_tree_with_selection(params.sample_size, params.ne, selection_params)

    # Add mutations with selection effects
    println("Adding mutations with selection...")
    add_selected_mutations!(root, params, selection_params)

    # Extract genotypes
    println("Extracting genotypes...")
    genotypes, positions = extract_genotypes(root, params.sample_size, params.sequence_length)

    println("Simulation complete!")
    println("Generated $(length(positions)) variant sites")
    println("Genotype matrix: $(size(genotypes))")

    return genotypes, positions
end

"""
    build_coalescent_tree_with_selection(sample_size::Int, ne::Int, 
                                        selection_params::SelectionParameters) -> CoalescentNode

Build a coalescent tree incorporating selection effects on coalescent rates.
"""
function build_coalescent_tree_with_selection(sample_size::Int, ne::Int,
    selection_params::SelectionParameters)
    n_haplotypes = 2 * sample_size

    # Create leaf nodes
    nodes = [CoalescentNode(i) for i in 1:n_haplotypes]
    active_nodes = copy(nodes)

    node_id = n_haplotypes + 1
    current_time = 0.0

    while length(active_nodes) > 1
        k = length(active_nodes)

        # Get effective population size considering selection
        ne_eff = modify_coalescent_rates!(active_nodes, selection_params, current_time)

        # Standard coalescent rate with modified Ne
        rate = k * (k - 1) / (4 * ne_eff)

        # Sample waiting time
        dt = rand(Exponential(1 / rate))
        current_time += dt

        # Choose two lineages to coalesce
        if length(active_nodes) >= 2
            indices_to_remove = sort([rand(1:length(active_nodes)),
                    rand(1:length(active_nodes))], rev=true)
            # Ensure we don't pick the same lineage twice
            while indices_to_remove[1] == indices_to_remove[2]
                indices_to_remove[2] = rand(1:length(active_nodes))
            end
            indices_to_remove = sort(unique(indices_to_remove), rev=true)

            child1 = active_nodes[indices_to_remove[1]]
            child2 = active_nodes[indices_to_remove[2]]

            # Create parent node
            parent = CoalescentNode(node_id, nothing, [child1, child2], current_time, Int[])
            child1.parent = parent
            child2.parent = parent

            # Update active nodes
            for idx in indices_to_remove
                deleteat!(active_nodes, idx)
            end
            push!(active_nodes, parent)

            node_id += 1
        else
            break
        end
    end

    return active_nodes[1]
end