#!/usr/bin/env julia
push!(LOAD_PATH, "src")
using GenotypeSimulator
using Random

# Use a typed struct instead of Dict{Symbol, Any} for juliac compatibility
mutable struct SimParams
    seed::Int
    has_seed::Bool
    ne::Int
    mu::Float64
    rho::Float64
    seq_length::Int
    samples::Int
    output::String
    format::String
    recombination::Bool
    selection::String
    selection_coeff::Float64
    selected_pos::Int
    has_selected_pos::Bool
    dominance::Float64
    stats::Bool
    quiet::Bool
    help::Bool
end

function default_params()::SimParams
    SimParams(
        0, false,           # seed, has_seed
        10_000,             # ne
        1.25e-8,            # mu
        1e-8,               # rho
        1_000_000,          # seq_length
        100,                # samples
        "simulation",       # output
        "csv",              # format
        false,              # recombination
        "neutral",          # selection
        0.01,               # selection_coeff
        0, false,           # selected_pos, has_selected_pos
        0.5,                # dominance
        false,              # stats
        false,              # quiet
        false               # help
    )
end

function print_usage()
    print(Core.stdout, """
GenotypeSimulator.jl - Coalescent-based genotype simulation

Usage:
    genotype_simulator [OPTIONS]

Options:
    --help, -h              Show this help message
    --seed SEED             Set random seed (default: random)
    --ne NE                 Effective population size (default: 10000)
    --mu MU                 Mutation rate per bp per generation (default: 1.25e-8)
    --rho RHO               Recombination rate per bp per generation (default: 1e-8)
    --length LENGTH         Sequence length in bp (default: 1000000)
    --samples N             Number of diploid individuals (default: 100)
    --output PREFIX         Output file prefix (default: "simulation")
    --format FORMAT         Output format: csv, vcf, plink (default: csv)
    --recombination         Enable recombination (ARG simulation)
    --selection TYPE        Selection type: neutral, directional, balancing, background
    --selection-coeff S     Selection coefficient (default: 0.01)
    --selected-pos POS      Position under selection (default: middle of sequence)
    --dominance H           Dominance coefficient (default: 0.5)
    --stats                 Print detailed statistics
    --quiet                 Suppress progress messages

Examples:
    genotype_simulator --samples 50 --length 500000 --output small_sim
    genotype_simulator --ne 5000 --mu 2e-8 --format vcf --stats
    genotype_simulator --seed 42 --quiet --format plink
""")
end

function parse_args(args::Vector{String})::SimParams
    params = default_params()

    i = 1
    while i <= length(args)
        arg = args[i]

        if arg == "--help" || arg == "-h"
            params.help = true
        elseif arg == "--seed"
            i += 1
            params.seed = parse(Int, args[i])
            params.has_seed = true
        elseif arg == "--ne"
            i += 1
            params.ne = parse(Int, args[i])
        elseif arg == "--mu"
            i += 1
            params.mu = parse(Float64, args[i])
        elseif arg == "--rho"
            i += 1
            params.rho = parse(Float64, args[i])
        elseif arg == "--length"
            i += 1
            params.seq_length = parse(Int, args[i])
        elseif arg == "--samples"
            i += 1
            params.samples = parse(Int, args[i])
        elseif arg == "--output"
            i += 1
            params.output = args[i]
        elseif arg == "--format"
            i += 1
            fmt = lowercase(args[i])
            if fmt == "csv" || fmt == "vcf" || fmt == "plink"
                params.format = fmt
            else
                print(Core.stderr, "Invalid format: ")
                print(Core.stderr, fmt)
                print(Core.stderr, ". Use csv, vcf, or plink.\n")
                return params
            end
        elseif arg == "--recombination"
            params.recombination = true
        elseif arg == "--selection"
            i += 1
            sel = lowercase(args[i])
            if sel == "neutral" || sel == "directional" || sel == "balancing" || sel == "background"
                params.selection = sel
            else
                print(Core.stderr, "Invalid selection type: ")
                print(Core.stderr, sel)
                print(Core.stderr, ". Use neutral, directional, balancing, or background.\n")
                return params
            end
        elseif arg == "--selection-coeff"
            i += 1
            params.selection_coeff = parse(Float64, args[i])
        elseif arg == "--selected-pos"
            i += 1
            params.selected_pos = parse(Int, args[i])
            params.has_selected_pos = true
        elseif arg == "--dominance"
            i += 1
            params.dominance = parse(Float64, args[i])
        elseif arg == "--stats"
            params.stats = true
        elseif arg == "--quiet"
            params.quiet = true
        else
            print(Core.stderr, "Unknown argument: ")
            print(Core.stderr, arg)
            print(Core.stderr, "\n")
        end

        i += 1
    end

    return params
end

function run_simulation(params::SimParams)::Cint
    if params.help
        print_usage()
        return Cint(0)
    end

    # Set random seed if provided
    if params.has_seed
        Random.seed!(params.seed)
        if !params.quiet
            print(Core.stdout, "Random seed set to: ")
            print(Core.stdout, params.seed)
            print(Core.stdout, "\n")
        end
    end

    # Create population parameters
    pop_params = PopulationParams(
        params.ne,
        params.mu,
        params.rho,
        params.seq_length,
        params.samples
    )

    if !params.quiet
        print(Core.stdout, "GenotypeSimulator.jl - Starting simulation\n")
        print(Core.stdout, "Parameters:\n")
        print(Core.stdout, "  Effective population size: ")
        print(Core.stdout, pop_params.ne)
        print(Core.stdout, "\n  Mutation rate: ")
        print(Core.stdout, pop_params.mutation_rate)
        print(Core.stdout, "\n  Recombination rate: ")
        print(Core.stdout, pop_params.recombination_rate)
        print(Core.stdout, "\n  Sequence length: ")
        print(Core.stdout, pop_params.sequence_length)
        print(Core.stdout, " bp\n  Sample size: ")
        print(Core.stdout, pop_params.sample_size)
        print(Core.stdout, " individuals\n  Output format: ")
        print(Core.stdout, params.format)
        print(Core.stdout, "\n\n")
    end

    # Determine selected position
    sel_pos::Int = params.has_selected_pos ? params.selected_pos : pop_params.sequence_length ÷ 2

    # Create selection model (avoiding closures/Box)
    local selection_model::SelectionModel
    if params.selection == "directional"
        selection_model = DirectionalSelection(params.selection_coeff, params.dominance, sel_pos, 1000.0)
    elseif params.selection == "balancing"
        selection_model = BalancingSelection(Float64[params.selection_coeff, -params.selection_coeff], sel_pos, Float64[0.5, 0.5])
    elseif params.selection == "background"
        selection_model = BackgroundSelection(1e-9, -params.selection_coeff, params.dominance)
    else
        selection_model = NeutralSelection()
    end

    selection_params = SelectionParameters(selection_model, pop_params.ne)

    # Run simulation
    local positions::Vector{Int}

    if params.selection != "neutral"
        genotypes_raw, positions = simulate_with_selection(pop_params, selection_params)
    elseif params.recombination
        result = simulate_with_recombination(pop_params)
        genotypes_raw = result[1]
        positions = Vector{Int}(result[2])
    else
        genotypes_raw, positions = simulate_genotypes(pop_params)
    end

    # Convert to smallest dense integer matrix type to save memory.
    # For genotype data (values 0/1), UInt8 uses 1 byte vs Int64's 8 bytes.
    local genotypes_hap::Matrix{UInt8}
    if isa(genotypes_raw, Matrix{UInt8})
        genotypes_hap = genotypes_raw
    elseif isa(genotypes_raw, BitPackedGenotypes)
        genotypes_hap = Matrix{UInt8}(genotypes_raw)
    elseif isa(genotypes_raw, SparseGenotypes)
        genotypes_hap = Matrix{UInt8}(genotypes_raw)
    else
        genotypes_hap = Matrix{UInt8}(genotypes_raw)
    end

    # Free the raw result if it was a different object
    genotypes_raw = nothing
    GC.gc(false)

    # Convert to diploid (values 0/1/2 still fit in UInt8)
    genotypes_diploid = haplotypes_to_diploid(genotypes_hap)

    # Free haploid matrix — diploid is all we need downstream
    genotypes_hap = Matrix{UInt8}(undef, 0, 0)
    GC.gc(false)

    # Calculate statistics
    stats = calculate_stats(genotypes_diploid, positions)

    if params.stats || !params.quiet
        print_stats(stats)
    end

    # Save results
    if params.format == "plink"
        save_genotypes(genotypes_diploid, positions, params.output, format="plink")
    else
        output_file = params.output * "." * params.format
        save_genotypes(genotypes_diploid, positions, output_file, format=params.format)
    end

    if !params.quiet
        print(Core.stdout, "\nSimulation complete!\n")
        if params.format == "plink"
            print(Core.stdout, "Results saved to ")
            print(Core.stdout, params.output)
            print(Core.stdout, ".ped and ")
            print(Core.stdout, params.output)
            print(Core.stdout, ".map\n")
        else
            print(Core.stdout, "Results saved to ")
            print(Core.stdout, params.output)
            print(Core.stdout, ".")
            print(Core.stdout, params.format)
            print(Core.stdout, "\n")
        end
    end

    return Cint(0)
end

function (@main)(args::Vector{String})::Cint
    try
        params = parse_args(args)
        return run_simulation(params)
    catch e
        if isa(e, ArgumentError)
            print(Core.stderr, "Error: ArgumentError\n")
            print(Core.stderr, "Use --help for usage information.\n")
            return Cint(1)
        elseif isa(e, BoundsError)
            print(Core.stderr, "Error: BoundsError\n")
            print(Core.stderr, "Use --help for usage information.\n")
            return Cint(1)
        else
            print(Core.stderr, "Unexpected error occurred: ")
            print(Core.stderr, typeof(e))
            print(Core.stderr, "\n")
            showerror(Core.stderr, e)
            print(Core.stderr, "\n")
            for (exc, bt) in current_exceptions()
                showerror(Core.stderr, exc, bt)
                print(Core.stderr, "\n")
            end
            return Cint(1)
        end
    end
end
