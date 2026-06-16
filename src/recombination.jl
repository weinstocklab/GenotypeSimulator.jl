"""
Recombination modeling for coalescent simulation

This module implements recombination breakpoints and the ancestral recombination graph (ARG)
following approaches from tskit/msprime.
"""

"""
    TimedEvent

Lightweight struct for storing events during ARG construction.
Avoids tuple allocation overhead.
"""
struct TimedEvent
    event_type::Symbol
    time::Float64
    index::Int
end

"""
    RecombinationBreakpoint

Represents a recombination breakpoint in the sequence.
"""
struct RecombinationBreakpoint
    position::Int           # Position in sequence
    time::Float64          # Time when recombination occurred
    left_parent::Int       # Node ID of left parent
    right_parent::Int      # Node ID of right parent
    child::Int             # Node ID of child (recombinant)
end

"""
    RecombinationMap

Represents a recombination rate map across the sequence.
"""
struct RecombinationMap
    positions::Vector{Int}      # Breakpoint positions
    rates::Vector{Float64}      # Recombination rates between positions
    cumulative_rates::Vector{Float64}  # Cumulative recombination rates

    function RecombinationMap(positions::Vector{Int}, rates::Vector{Float64})
        @assert length(positions) == length(rates) + 1 "Positions should be one longer than rates"
        @assert all(diff(positions) .> 0) "Positions must be strictly increasing"
        @assert all(rates .>= 0) "Rates must be non-negative"

        # Calculate cumulative rates for efficient sampling
        cumulative = zeros(Float64, length(rates))
        if length(rates) > 0
            cumulative[1] = rates[1] * (positions[2] - positions[1])
            for i in 2:length(rates)
                cumulative[i] = cumulative[i-1] + rates[i] * (positions[i+1] - positions[i])
            end
        end

        new(positions, rates, cumulative)
    end
end

"""
    uniform_recombination_map(sequence_length::Int, rate::Float64) -> RecombinationMap

Create a uniform recombination map with constant rate across the sequence.
"""
function uniform_recombination_map(sequence_length::Int, rate::Float64)
    positions = [1, sequence_length + 1]
    rates = [rate]
    return RecombinationMap(positions, rates)
end

"""
    periodic_hotspot_recombination_map(sequence_length::Int, rate::Float64;
                                       period_bp::Int=100_000,
                                       hotspot_fraction::Float64=0.05,
                                       hotspot_multiplier::Float64=20.0,
                                       coldspot_multiplier::Float64=0.02) -> RecombinationMap

Create a piecewise-constant recombination map with periodic hotspot/coldspot structure.
The map is rescaled so the genome-wide average recombination rate equals `rate`.
"""
function periodic_hotspot_recombination_map(sequence_length::Int, rate::Float64;
                                            period_bp::Int=100_000,
                                            hotspot_fraction::Float64=0.05,
                                            hotspot_multiplier::Float64=20.0,
                                            coldspot_multiplier::Float64=0.02)
    @assert sequence_length > 0 "sequence_length must be positive"
    @assert period_bp > 1 "period_bp must be > 1"
    @assert 0 < hotspot_fraction < 1 "hotspot_fraction must be in (0,1)"
    @assert hotspot_multiplier >= 0 "hotspot_multiplier must be >= 0"
    @assert coldspot_multiplier >= 0 "coldspot_multiplier must be >= 0"

    avg_multiplier = hotspot_fraction * hotspot_multiplier +
                     (1 - hotspot_fraction) * coldspot_multiplier
    @assert avg_multiplier > 0 "Average multiplier must be positive"

    hotspot_rate = rate * hotspot_multiplier / avg_multiplier
    coldspot_rate = rate * coldspot_multiplier / avg_multiplier
    hotspot_len = max(1, round(Int, hotspot_fraction * period_bp))

    positions = Int[1]
    rates = Float64[]
    start_pos = 1

    while start_pos <= sequence_length
        stop_pos = min(sequence_length + 1, start_pos + period_bp)
        hotspot_stop = min(stop_pos, start_pos + hotspot_len)

        if hotspot_stop > start_pos
            push!(rates, hotspot_rate)
            push!(positions, hotspot_stop)
        end
        if stop_pos > hotspot_stop
            push!(rates, coldspot_rate)
            push!(positions, stop_pos)
        end
        start_pos = stop_pos
    end

    positions[end] == sequence_length + 1 || push!(positions, sequence_length + 1)
    return RecombinationMap(positions, rates)
end

"""
    sample_recombination_position(recomb_map::RecombinationMap) -> Int

Sample a recombination position from the recombination map.
"""
function sample_recombination_position(recomb_map::RecombinationMap)
    return sample_recombination_position(recomb_map, Random.GLOBAL_RNG)
end

"""
    sample_recombination_position(recomb_map::RecombinationMap, rng::AbstractRNG) -> Int

Sample a recombination position from the recombination map using the provided RNG.
"""
function sample_recombination_position(recomb_map::RecombinationMap, rng::AbstractRNG)
    if length(recomb_map.rates) == 0
        return 1
    end

    total_rate = recomb_map.cumulative_rates[end]
    if total_rate == 0
        return rand(rng, 1:recomb_map.positions[end]-1)
    end

    # Sample from cumulative distribution
    u = rand(rng) * total_rate

    # Find the interval
    interval = searchsortedfirst(recomb_map.cumulative_rates, u)
    interval = min(interval, length(recomb_map.rates))

    # Sample position within the interval
    if interval == 1
        cumulative_prev = 0.0
    else
        cumulative_prev = recomb_map.cumulative_rates[interval-1]
    end

    rate = recomb_map.rates[interval]
    if rate == 0
        # Uniform within interval
        pos = rand(rng, recomb_map.positions[interval]:(recomb_map.positions[interval+1]-1))
    else
        # Proportional to rate
        remaining = u - cumulative_prev
        interval_length = recomb_map.positions[interval+1] - recomb_map.positions[interval]
        relative_pos = remaining / (rate * interval_length)
        pos = recomb_map.positions[interval] + floor(Int, relative_pos * interval_length)
    end

    return clamp(pos, 1, recomb_map.positions[end] - 1)
end

"""
    RecombinationEvent

Represents a recombination event during simulation.
"""
struct RecombinationEvent
    time::Float64
    position::Int
    lineage_id::Int
    left_parent_id::Int
    right_parent_id::Int
end

"""
    Lineage

Represents a lineage in the ARG simulation with genomic intervals.
Now includes optimized data structures for faster operations.
"""
mutable struct Lineage
    id::Int
    intervals::Vector{Tuple{Int,Int}}  # (start, end) intervals this lineage covers
    min_start::Int
    max_end::Int
    node::Union{CoalescentNode,Nothing}

    Lineage(id::Int, start::Int, stop::Int) = new(id, [(start, stop)], start, stop, nothing)
end

"""
    update_lineage_bounds!(lineage::Lineage)

Update cached interval bounds for a lineage after interval changes.
"""
function update_lineage_bounds!(lineage::Lineage)
    min_start = typemax(Int)
    max_end = 0
    for (s, e) in lineage.intervals
        if s < min_start
            min_start = s
        end
        if e > max_end
            max_end = e
        end
    end
    lineage.min_start = min_start
    lineage.max_end = max_end
    return nothing
end

"""
    IntervalTreapNode

Treap node storing interval bounds for fast overlap queries.
"""
mutable struct IntervalTreapNode
    start::Int
    stop::Int
    id::Int
    priority::UInt64
    max_end::Int
    min_start::Int
    left::Union{IntervalTreapNode, Nothing}
    right::Union{IntervalTreapNode, Nothing}

    function IntervalTreapNode(start::Int, stop::Int, id::Int)
        prio = treap_priority(start, id)
        new(start, stop, id, prio, stop, start, nothing, nothing)
    end
end

"""
    treap_priority(start::Int, id::Int) -> UInt64

Deterministic priority for the treap to keep results reproducible.
"""
function treap_priority(start::Int, id::Int)
    x = UInt64(start) ⊻ (UInt64(id) * 0x9e3779b97f4a7c15)
    x ⊻= x << 13
    x ⊻= x >> 17
    x ⊻= x << 5
    return x
end

"""
    update_treap_node!(node::IntervalTreapNode)

Refresh subtree bounds after mutations.
"""
function update_treap_node!(node::IntervalTreapNode)
    node.max_end = node.stop
    node.min_start = node.start
    if node.left !== nothing
        node.max_end = max(node.max_end, node.left.max_end)
        node.min_start = min(node.min_start, node.left.min_start)
    end
    if node.right !== nothing
        node.max_end = max(node.max_end, node.right.max_end)
        node.min_start = min(node.min_start, node.right.min_start)
    end
    return nothing
end

function treap_rotate_left(node::IntervalTreapNode)
    right = node.right::IntervalTreapNode
    node.right = right.left
    right.left = node
    update_treap_node!(node)
    update_treap_node!(right)
    return right
end

function treap_rotate_right(node::IntervalTreapNode)
    left = node.left::IntervalTreapNode
    node.left = left.right
    left.right = node
    update_treap_node!(node)
    update_treap_node!(left)
    return left
end

function treap_insert(node::Union{IntervalTreapNode, Nothing}, start::Int, stop::Int, id::Int)
    if node === nothing
        return IntervalTreapNode(start, stop, id)
    end

    if start < node.start || (start == node.start && id < node.id)
        node.left = treap_insert(node.left, start, stop, id)
        if node.left !== nothing && node.left.priority < node.priority
            node = treap_rotate_right(node)
        end
    else
        node.right = treap_insert(node.right, start, stop, id)
        if node.right !== nothing && node.right.priority < node.priority
            node = treap_rotate_left(node)
        end
    end

    update_treap_node!(node)
    return node
end

function treap_merge(left::Union{IntervalTreapNode, Nothing}, right::Union{IntervalTreapNode, Nothing})
    if left === nothing
        return right
    elseif right === nothing
        return left
    end

    if left.priority < right.priority
        left.right = treap_merge(left.right, right)
        update_treap_node!(left)
        return left
    else
        right.left = treap_merge(left, right.left)
        update_treap_node!(right)
        return right
    end
end

function treap_delete(node::Union{IntervalTreapNode, Nothing}, start::Int, id::Int)
    if node === nothing
        return nothing
    end

    if start == node.start && id == node.id
        return treap_merge(node.left, node.right)
    elseif start < node.start || (start == node.start && id < node.id)
        node.left = treap_delete(node.left, start, id)
    else
        node.right = treap_delete(node.right, start, id)
    end

    update_treap_node!(node)
    return node
end

function treap_search_overlap(node::Union{IntervalTreapNode, Nothing}, start::Int, stop::Int, exclude_id::Int)
    if node === nothing
        return nothing
    end

    if node.left !== nothing && node.left.max_end >= start
        candidate = treap_search_overlap(node.left, start, stop, exclude_id)
        if candidate !== nothing
            return candidate
        end
    end

    if node.id != exclude_id && node.start <= stop && node.stop >= start
        return node.id
    end

    if node.right !== nothing && node.start <= stop
        return treap_search_overlap(node.right, start, stop, exclude_id)
    end

    return nothing
end

"""
    LineageManager

Optimized container for managing lineages with O(1) lookup by ID and 
efficient overlap detection.
"""
mutable struct LineageManager
    lineages::Dict{Int, Lineage}     # O(1) lookup by ID
    active_ids::Vector{Int}          # Active lineage IDs
    active_index::Dict{Int, Int}     # Active ID -> index lookup
    interval_root::Union{IntervalTreapNode, Nothing}  # Interval treap root
    
    function LineageManager(n_haplotypes::Int, sequence_length::Int)
        lineages = Dict{Int, Lineage}()
        active_ids = Int[]
        active_index = Dict{Int, Int}()
        interval_root = nothing
        
        for i in 1:n_haplotypes
            lineage = Lineage(i, 1, sequence_length + 1)
            lineages[i] = lineage
            push!(active_ids, i)
            active_index[i] = length(active_ids)
            interval_root = treap_insert(interval_root, lineage.min_start, lineage.max_end, lineage.id)
        end
        
        new(lineages, active_ids, active_index, interval_root)
    end
end

"""
    add_lineage!(manager::LineageManager, lineage::Lineage)

Add a new lineage to the manager.
"""
function add_lineage!(manager::LineageManager, lineage::Lineage)
    manager.lineages[lineage.id] = lineage
    push!(manager.active_ids, lineage.id)
    manager.active_index[lineage.id] = length(manager.active_ids)
    manager.interval_root = treap_insert(manager.interval_root, lineage.min_start, lineage.max_end, lineage.id)
end

"""
    remove_lineages!(manager::LineageManager, ids::Vector{Int})

Remove lineages by their IDs.
"""
function remove_lineages!(manager::LineageManager, ids::Vector{Int})
    for id in ids
        if haskey(manager.lineages, id)
            lineage = manager.lineages[id]
            manager.interval_root = treap_delete(manager.interval_root, lineage.min_start, id)
        end
        if haskey(manager.active_index, id)
            idx = manager.active_index[id]
            last_id = manager.active_ids[end]
            manager.active_ids[idx] = last_id
            pop!(manager.active_ids)
            delete!(manager.active_index, id)
            if last_id != id
                manager.active_index[last_id] = idx
            end
        end
        delete!(manager.lineages, id)
    end
end

"""
    update_lineage_interval!(manager::LineageManager, lineage::Lineage, old_start::Int)

Update interval tree entry after a lineage bounds change.
"""
function update_lineage_interval!(manager::LineageManager, lineage::Lineage, old_start::Int)
    manager.interval_root = treap_delete(manager.interval_root, old_start, lineage.id)
    manager.interval_root = treap_insert(manager.interval_root, lineage.min_start, lineage.max_end, lineage.id)
    return nothing
end

"""
    get_lineage(manager::LineageManager, id::Int) -> Union{Lineage, Nothing}

Get a lineage by ID with O(1) lookup.
Note: For type-stable code, prefer using haskey() check followed by direct dictionary access.
"""
function get_lineage(manager::LineageManager, id::Int)
    return get(manager.lineages, id, nothing)
end

"""
    has_lineage(manager::LineageManager, id::Int) -> Bool

Type-stable check for lineage existence.
"""
function has_lineage(manager::LineageManager, id::Int)
    return haskey(manager.lineages, id)
end

"""
    count_active(manager::LineageManager) -> Int

Get the number of active lineages.
"""
count_active(manager::LineageManager) = length(manager.active_ids)

"""
    overlaps_fast(lineage1::Lineage, lineage2::Lineage) -> Bool

Optimized overlap detection using sorted interval checking.
"""
function overlaps_fast(lineage1::Lineage, lineage2::Lineage)
    # For small numbers of intervals, the simple nested loop is efficient
    # and cache-friendly
    for (s1, e1) in lineage1.intervals
        for (s2, e2) in lineage2.intervals
            if s1 < e2 && s2 < e1  # Intervals overlap
                return true
            end
        end
    end
    return false
end

"""
    find_overlapping_pair(manager::LineageManager) -> Union{Tuple{Int,Int}, Nothing}

Find the first pair of overlapping lineages efficiently.
Returns (id1, id2) or nothing if no overlapping pair exists.
"""
function find_overlapping_pair(manager::LineageManager)
    active = manager.active_ids
    n = length(active)
    
    if n < 2
        return nothing
    end

    for id in active
        lineage = manager.lineages[id]
        candidate_id = treap_search_overlap(manager.interval_root, lineage.min_start, lineage.max_end, id)
        if candidate_id !== nothing
            if overlaps_fast(lineage, manager.lineages[candidate_id])
                return (id, candidate_id)
            end
        end
    end

    return nothing
end

"""
    split_lineage!(lineage::Lineage, position::Int) -> Lineage

Split a lineage at the given position, returning the right part.
"""
function split_lineage!(lineage::Lineage, position::Int)
    new_intervals = Tuple{Int,Int}[]
    right_intervals = Tuple{Int,Int}[]

    for (start, stop) in lineage.intervals
        if position <= start
            # Entire interval goes to right
            push!(right_intervals, (start, stop))
        elseif position >= stop
            # Entire interval stays in left
            push!(new_intervals, (start, stop))
        else
            # Split the interval
            push!(new_intervals, (start, position))
            push!(right_intervals, (position, stop))
        end
    end

    lineage.intervals = new_intervals
    update_lineage_bounds!(lineage)
    right_lineage = Lineage(lineage.id + 1000000, 1, 1)  # Temporary ID
    right_lineage.intervals = right_intervals
    update_lineage_bounds!(right_lineage)

    return right_lineage
end

"""
    simulate_recombination_events(n_lineages::Int, sequence_length::Int, 
                                 recomb_map::RecombinationMap, 
                                 coalescent_times::Vector{Float64}) -> Vector{RecombinationEvent}

Simulate recombination events along with coalescent events.
"""
function simulate_recombination_events(n_lineages::Int, sequence_length::Int,
    recomb_map::RecombinationMap,
    coalescent_times::Vector{Float64},
    rng::AbstractRNG=Random.GLOBAL_RNG)
    events = RecombinationEvent[]
    current_time = 0.0
    lineage_count = n_lineages

    # Calculate total recombination rate
    total_recomb_rate = isempty(recomb_map.cumulative_rates) ? 0.0 : recomb_map.cumulative_rates[end]
    if total_recomb_rate == 0.0
        return events
    end

    coal_event_idx = 1
    n_coal = length(coalescent_times)
    next_coal_time = n_coal > 0 ? coalescent_times[1] : Inf

    while lineage_count > 1 && coal_event_idx <= n_coal
        # Rate of recombination events
        recomb_rate = lineage_count * total_recomb_rate

        # Time to next recombination
        next_recomb_time = current_time + randexp(rng) / recomb_rate

        if next_recomb_time < next_coal_time
            # Recombination happens first
            current_time = next_recomb_time

            # Choose lineage to recombine
            lineage_id = rand(rng, 1:lineage_count)

            # Choose recombination position
            position = sample_recombination_position(recomb_map, rng)

            # Create recombination event
            push!(events, RecombinationEvent(
                current_time,
                position,
                lineage_id,
                lineage_id,  # Will be updated when we build the tree
                lineage_count + 1  # New lineage ID
            ))
            lineage_count += 1  # Recombination increases lineage count
        else
            # Coalescence happens first
            current_time = next_coal_time
            lineage_count -= 1  # Coalescence decreases lineage count

            coal_event_idx += 1
            next_coal_time = coal_event_idx <= n_coal ? coalescent_times[coal_event_idx] : Inf
        end
    end

    return events
end

"""
    build_arg_tree(sample_size::Int, ne::Int, sequence_length::Int, 
                   recomb_map::RecombinationMap) -> (CoalescentNode, Vector{RecombinationEvent})

Build an Ancestral Recombination Graph (ARG) using Hudson's algorithm with
thinning-based coalescence.

Coalescence pairs are drawn uniformly from all active lineages at rate
k*(k-1)/(4Ne).  Pairs that share no ancestral material (no overlapping
genomic bins) are discarded — this is the standard Poisson-thinning
trick that preserves correct per-position coalescence rates.

The genome is divided into ~128 bins. Per-lineage `BitSet` membership
enables O(1) overlap tests via `isdisjoint`.
"""

# Merge two sorted Int vectors into a new sorted UNIQUE vector
function _merge_sorted_int_vecs(a::Vector{Int}, b::Vector{Int})
    result = Vector{Int}()
    sizehint!(result, max(length(a), length(b)))
    i, j = 1, 1
    @inbounds while i <= length(a) && j <= length(b)
        if a[i] < b[j]
            (isempty(result) || result[end] != a[i]) && push!(result, a[i])
            i += 1
        elseif a[i] > b[j]
            (isempty(result) || result[end] != b[j]) && push!(result, b[j])
            j += 1
        else  # equal — take one copy
            (isempty(result) || result[end] != a[i]) && push!(result, a[i])
            i += 1; j += 1
        end
    end
    @inbounds while i <= length(a)
        (isempty(result) || result[end] != a[i]) && push!(result, a[i])
        i += 1
    end
    @inbounds while j <= length(b)
        (isempty(result) || result[end] != b[j]) && push!(result, b[j])
        j += 1
    end
    return result
end

function build_arg_tree(sample_size::Int, ne::Int, sequence_length::Int,
    recomb_map::RecombinationMap; rng::AbstractRNG=Random.GLOBAL_RNG,
    demography::Union{Nothing, DemographyModel}=nothing)
    n_haplotypes = 2 * sample_size
    total_recomb_rate = isempty(recomb_map.cumulative_rates) ? 0.0 : recomb_map.cumulative_rates[end]

    # ── Genomic-bin index (128 bins) ───────────────────────────────────────
    target_bins = 128
    bin_size = max(1, cld(sequence_length, target_bins))
    n_bins   = cld(sequence_length, bin_size)

    @inline pos_to_bin(p::Int) = clamp(cld(p, bin_size), 1, n_bins)
    @inline function bins_for_interval(s::Int, e::Int)
        # interval is [s, e)
        return pos_to_bin(s):pos_to_bin(max(s, e - 1))
    end

    # Per-bin lineage lists (plain vectors — fast append/random access)
    bin_lineages = [Int[] for _ in 1:n_bins]

    # Managed coalescable-bin vector for O(1) random pick + O(1) remove
    coal_bins     = Int[]                  # bins with >= 2 lineages
    coal_bin_idx  = zeros(Int, n_bins)     # bin → index in coal_bins (0 = absent)

    @inline function _add_coal_bin(b::Int)
        if coal_bin_idx[b] == 0
            push!(coal_bins, b)
            coal_bin_idx[b] = length(coal_bins)
        end
    end
    @inline function _del_coal_bin(b::Int)
        idx = coal_bin_idx[b]
        if idx > 0
            last_b = coal_bins[end]
            coal_bins[idx] = last_b
            coal_bin_idx[last_b] = idx
            pop!(coal_bins)
            coal_bin_idx[b] = 0
        end
    end

    # ── index / unindex helpers ────────────────────────────────────────────
    # Per-lineage bin membership (BitSet) — enables O(1) overlap tests.
    lineage_bins = Dict{Int, BitSet}()

    function index_lineage!(id::Int, intervals::Vector{Tuple{Int,Int}})
        bs = BitSet()
        for (s, e) in intervals
            for b in bins_for_interval(s, e)
                push!(bs, b)
            end
        end
        lineage_bins[id] = bs
        for b in bs
            push!(bin_lineages[b], id)
            length(bin_lineages[b]) >= 2 && _add_coal_bin(b)
        end
    end

    function unindex_lineage!(id::Int)
        for b in lineage_bins[id]
            bl = bin_lineages[b]
            idx = findfirst(==(id), bl)
            if idx !== nothing
                bl[idx] = bl[end]
                pop!(bl)
            end
            length(bl) < 2 && _del_coal_bin(b)
        end
        delete!(lineage_bins, id)
    end

    # ── Managed active-lineage vector for O(1) random pick ─────────────
    active_ids  = Int[]                   # active lineage IDs
    active_idx  = Dict{Int,Int}()         # id → index in active_ids

    function add_active!(id::Int)
        push!(active_ids, id)
        active_idx[id] = length(active_ids)
    end
    function remove_active!(id::Int)
        idx = active_idx[id]
        last_id = active_ids[end]
        active_ids[idx] = last_id
        active_idx[last_id] = idx
        pop!(active_ids)
        delete!(active_idx, id)
    end

    # ── Lineage data ───────────────────────────────────────────────────────
    lineage_intervals = Dict{Int, Vector{Tuple{Int,Int}}}()
    lineage_nodes     = Dict{Int, CoalescentNode}()

    nodes = [CoalescentNode(i) for i in 1:n_haplotypes]
    for i in 1:n_haplotypes
        ivs = [(1, sequence_length + 1)]
        lineage_intervals[i] = ivs
        lineage_nodes[i]     = nodes[i]
        index_lineage!(i, ivs)
        add_active!(i)
    end

    recomb_events = RecombinationEvent[]
    node_intervals = Dict{Int, Vector{Tuple{Int,Int}}}()  # ancestral material per node (for mutation placement)
    lineage_samples = Dict{Int, Vector{Int}}()   # lineage_id → sorted list of original sample haplotype ids
    for i in 1:n_haplotypes
        lineage_samples[i] = [i]
    end
    node_id       = n_haplotypes + 1
    current_time  = 0.0

    max_iterations = min(200 * n_haplotypes * max(1, round(Int, 4 * ne * total_recomb_rate + 1)),
                         100_000_000)
    iter = 0

    while !isempty(coal_bins)
        iter += 1
        iter > max_iterations && break

        k = length(active_ids)
        k < 2 && break
        coal_scale = k * (k - 1) / 4.0
        recomb_rate = k * total_recomb_rate
        dt = sample_piecewise_wait(recomb_rate, coal_scale, current_time, ne, demography; rng=rng)
        current_time += dt
        coal_rate = coal_scale / effective_population_size(demography, current_time, ne)
        total_rate = coal_rate + recomb_rate
        total_rate <= 0.0 && break

        if rand(rng) < coal_rate / total_rate
            # ── Coalescence (thinning: uniform pair, skip if disjoint) ──
            i1 = rand(rng, 1:k)
            i2 = rand(rng, 1:(k - 1))
            i2 >= i1 && (i2 += 1)
            id1, id2 = active_ids[i1], active_ids[i2]

            # Poisson thinning: only coalesce if they share ancestral material
            if isdisjoint(lineage_bins[id1], lineage_bins[id2])
                continue
            end

            child1_node = lineage_nodes[id1]
            child2_node = lineage_nodes[id2]

            parent_node = CoalescentNode(node_id, nothing,
                [child1_node, child2_node], current_time, Int[])
            child1_node.parent = parent_node
            child2_node.parent = parent_node

            # Unindex + deactivate old lineages
            unindex_lineage!(id1)
            unindex_lineage!(id2)
            remove_active!(id1)
            remove_active!(id2)

            # Merge intervals (pre-allocate, concat, then merge overlaps)
            ivs1 = lineage_intervals[id1]
            ivs2 = lineage_intervals[id2]
            new_ivs = Vector{Tuple{Int,Int}}(undef, length(ivs1) + length(ivs2))
            copyto!(new_ivs, 1, ivs1, 1, length(ivs1))
            copyto!(new_ivs, length(ivs1) + 1, ivs2, 1, length(ivs2))

            # Merge overlapping/adjacent intervals to prevent unbounded growth
            sort!(new_ivs, by=first)
            n_merged = 1
            for j in 2:length(new_ivs)
                s, e = new_ivs[j]
                ms, me = new_ivs[n_merged]
                if s <= me  # overlapping or adjacent
                    new_ivs[n_merged] = (ms, max(me, e))
                else
                    n_merged += 1
                    new_ivs[n_merged] = (s, e)
                end
            end
            resize!(new_ivs, n_merged)

            # Save ancestral intervals before deletion (needed for mutation placement)
            node_intervals[id1] = lineage_intervals[id1]
            node_intervals[id2] = lineage_intervals[id2]

            delete!(lineage_intervals, id1)
            delete!(lineage_intervals, id2)
            delete!(lineage_nodes, id1)
            delete!(lineage_nodes, id2)

            # Merge sample sets: parent represents all samples from both children
            # NOTE: do NOT delete lineage_samples for id1/id2 — they're needed during
            # extraction to route mutations from fragment leaves back to original samples.
            s1 = get(lineage_samples, id1, Int[])
            s2 = get(lineage_samples, id2, Int[])
            merged_samples = _merge_sorted_int_vecs(s1, s2)
            lineage_samples[node_id] = merged_samples

            lineage_intervals[node_id] = new_ivs
            lineage_nodes[node_id]     = parent_node
            index_lineage!(node_id, new_ivs)
            add_active!(node_id)
            node_id += 1
        else
            # ── Recombination ──────────────────────────────────────────────
            lineage_id = active_ids[rand(rng, 1:length(active_ids))]
            position = sample_recombination_position(recomb_map, rng)

            intervals = lineage_intervals[lineage_id]
            left  = Tuple{Int,Int}[]
            right = Tuple{Int,Int}[]
            for (s, e) in intervals
                if position <= s
                    push!(right, (s, e))
                elseif position >= e
                    push!(left, (s, e))
                else
                    push!(left, (s, position))
                    push!(right, (position, e))
                end
            end

            if !isempty(left) && !isempty(right)
                unindex_lineage!(lineage_id)

                lineage_intervals[lineage_id] = left
                index_lineage!(lineage_id, left)

                right_id   = node_id
                node_id   += 1
                right_node = CoalescentNode(right_id, nothing, CoalescentNode[], lineage_nodes[lineage_id].time, Int[])
                lineage_intervals[right_id] = right
                lineage_nodes[right_id]     = right_node
                lineage_samples[right_id]   = lineage_samples[lineage_id]  # shared reference is fine (read-only)
                index_lineage!(right_id, right)
                add_active!(right_id)

                push!(recomb_events, RecombinationEvent(
                    current_time, position, lineage_id,
                    lineage_id, right_id))
            end
        end
    end

    # Walk up from first leaf to find the MRCA
    root = nodes[1]
    # Save remaining lineages' intervals to node_intervals for extraction filtering
    for id in active_ids
        if haskey(lineage_intervals, id)
            node_intervals[id] = lineage_intervals[id]
        end
    end

    # Connect all remaining lineages under a virtual root.
    # When the simulation ends with k > 1 (non-overlapping lineages), each covers
    # different genomic positions with independent subtrees. A virtual root ensures
    # all subtrees are reachable during genotype extraction.
    remaining_roots = CoalescentNode[]
    for id in active_ids
        push!(remaining_roots, lineage_nodes[id])
    end
    if length(remaining_roots) == 1
        root = remaining_roots[1]
    elseif length(remaining_roots) > 1
        root = CoalescentNode(node_id, nothing, remaining_roots, current_time + 1.0, Int[])
        for child in remaining_roots
            child.parent = root
        end
    else
        # Fallback: walk up from first leaf
        root = nodes[1]
        while root.parent !== nothing
            root = root.parent
        end
    end

    return root, recomb_events, node_intervals, lineage_samples
end

"""
    get_local_trees(root::CoalescentNode, recomb_events::Vector{RecombinationEvent}, 
                   sequence_length::Int) -> Vector{Tuple{Int, Int, CoalescentNode}}

Extract local trees from the ARG for different genomic intervals.
Returns vector of (start_pos, end_pos, local_root) tuples.
"""
function get_local_trees(root::CoalescentNode, recomb_events::Vector{RecombinationEvent},
    sequence_length::Int)
    # Get breakpoints from recombination events
    breakpoints = [1]
    for event in recomb_events
        push!(breakpoints, event.position)
    end
    push!(breakpoints, sequence_length + 1)

    sort!(unique!(breakpoints))

    # For now, return the same tree for all intervals
    # In a full implementation, we would trace through the ARG
    # to get the correct local tree for each interval
    local_trees = Tuple{Int,Int,CoalescentNode}[]

    for i in 1:length(breakpoints)-1
        start_pos = breakpoints[i]
        end_pos = breakpoints[i+1]
        push!(local_trees, (start_pos, end_pos, root))
    end

    return local_trees
end

"""
    build_arg_tree(params::PopulationParams, recomb_map::RecombinationMap; rng::AbstractRNG=Random.GLOBAL_RNG) -> (CoalescentNode, Vector{RecombinationEvent})

Convenience method for build_arg_tree that accepts PopulationParams.
"""
function build_arg_tree(params::PopulationParams, recomb_map::RecombinationMap; rng::AbstractRNG=Random.GLOBAL_RNG,
                        demography::Union{Nothing, DemographyModel}=nothing)
    return build_arg_tree(params.sample_size, params.ne, params.sequence_length, recomb_map;
                          rng=rng, demography=demography)
end
