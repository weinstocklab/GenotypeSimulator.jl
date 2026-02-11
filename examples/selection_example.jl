#!/usr/bin/env julia

"""
Selection example for GenotypeSimulator.jl

This script demonstrates how to:
1. Simulate genotypes under different selection models
2. Compare neutral vs. selected scenarios
3. Analyze selection signatures in genetic data
"""

using Pkg
Pkg.activate(".")

using GenotypeSimulator
using Random

# Set random seed for reproducibility
Random.seed!(42)

println("GenotypeSimulator.jl - Selection Example")
println("="^50)

# Define parameters for comparison
params = PopulationParams(
    10_000,     # Effective population size
    1.25e-8,    # Mutation rate per bp per generation
    1e-8,       # Recombination rate per bp per generation
    50_000,     # Sequence length (50kb for faster demo)
    50          # Sample size (50 individuals)
)

println("Parameters:")
println("  Ne: $(params.ne)")
println("  μ: $(params.mutation_rate)")
println("  ρ: $(params.recombination_rate)")
println("  Length: $(params.sequence_length) bp")
println("  Samples: $(params.sample_size) individuals")
println()

# 1. Neutral simulation (baseline)
println("1. Neutral Simulation (Baseline)")
println("-" * 35)

neutral_selection = NeutralSelection()
neutral_params = SelectionParameters(neutral_selection, params.ne)

genotypes_neutral, positions_neutral = simulate_with_selection(params, neutral_params)
genotypes_diploid_neutral = haplotypes_to_diploid(genotypes_neutral)
stats_neutral = calculate_stats(genotypes_diploid_neutral, positions_neutral)

println("Neutral results:")
println("  Variant sites: $(stats_neutral.n_sites)")
println("  Nucleotide diversity (π): $(round(stats_neutral.nucleotide_diversity, digits=6))")
println("  Tajima's D: $(round(tajimas_d(genotypes_diploid_neutral), digits=4))")
println()

# 2. Directional selection (selective sweep)
println("2. Directional Selection (Selective Sweep)")
println("-" * 45)

selected_position = params.sequence_length ÷ 2  # Middle of sequence
directional_selection = DirectionalSelection(
    0.05,                    # Selection coefficient (5% advantage)
    0.5,                     # Dominance coefficient (additive)
    selected_position,       # Position under selection
    2000.0,                  # Selection started 2000 generations ago
    1/(2*params.ne)         # Initial frequency (single copy)
)

directional_params = SelectionParameters(directional_selection, params.ne)

genotypes_directional, positions_directional = simulate_with_selection(params, directional_params)
genotypes_diploid_directional = haplotypes_to_diploid(genotypes_directional)
stats_directional = calculate_stats(genotypes_diploid_directional, positions_directional)

println("Directional selection results:")
println("  Variant sites: $(stats_directional.n_sites)")
println("  Nucleotide diversity (π): $(round(stats_directional.nucleotide_diversity, digits=6))")
println("  Tajima's D: $(round(tajimas_d(genotypes_diploid_directional), digits=4))")
println("  Selected position: $selected_position")
println()

# 3. Background selection
println("3. Background Selection")
println("-" * 25)

background_selection = BackgroundSelection(
    1e-8,      # Deleterious mutation rate
    -0.01,     # Selection coefficient against deleterious mutations
    0.0        # Recessive deleterious mutations
)

background_params = SelectionParameters(background_selection, params.ne)

genotypes_background, positions_background = simulate_with_selection(params, background_params)
genotypes_diploid_background = haplotypes_to_diploid(genotypes_background)
stats_background = calculate_stats(genotypes_diploid_background, positions_background)

println("Background selection results:")
println("  Variant sites: $(stats_background.n_sites)")
println("  Nucleotide diversity (π): $(round(stats_background.nucleotide_diversity, digits=6))")
println("  Tajima's D: $(round(tajimas_d(genotypes_diploid_background), digits=4))")
println()

# 4. Balancing selection
println("4. Balancing Selection")
println("-" * 23)

balancing_selection = BalancingSelection(
    [0.02, 0.02],           # Selection coefficients for both alleles
    selected_position,       # Same position as directional selection
    [0.4, 0.6]              # Equilibrium frequencies
)

balancing_params = SelectionParameters(balancing_selection, params.ne)

genotypes_balancing, positions_balancing = simulate_with_selection(params, balancing_params)
genotypes_diploid_balancing = haplotypes_to_diploid(genotypes_balancing)
stats_balancing = calculate_stats(genotypes_diploid_balancing, positions_balancing)

println("Balancing selection results:")
println("  Variant sites: $(stats_balancing.n_sites)")
println("  Nucleotide diversity (π): $(round(stats_balancing.nucleotide_diversity, digits=6))")
println("  Tajima's D: $(round(tajimas_d(genotypes_diploid_balancing), digits=4))")
println()

# 5. Comparison table
println("5. Selection Effects Comparison")
println("-" * 32)
println("Selection Type\t\tSites\tπ\t\tTajima's D")
println("-" * 55)
println("Neutral\t\t\t$(stats_neutral.n_sites)\t$(round(stats_neutral.nucleotide_diversity, digits=6))\t$(round(tajimas_d(genotypes_diploid_neutral), digits=4))")
println("Directional\t\t$(stats_directional.n_sites)\t$(round(stats_directional.nucleotide_diversity, digits=6))\t$(round(tajimas_d(genotypes_diploid_directional), digits=4))")
println("Background\t\t$(stats_background.n_sites)\t$(round(stats_background.nucleotide_diversity, digits=6))\t$(round(tajimas_d(genotypes_diploid_background), digits=4))")
println("Balancing\t\t$(stats_balancing.n_sites)\t$(round(stats_balancing.nucleotide_diversity, digits=6))\t$(round(tajimas_d(genotypes_diploid_balancing), digits=4))")
println()

# 6. Selection signature analysis
println("6. Selection Signature Analysis")
println("-" * 33)

function analyze_selection_signature(genotypes::Matrix{Int}, positions::Vector{Int}, 
                                   selected_pos::Int, window_size::Int=5000)
    if length(positions) == 0
        return Float64[], Int[]
    end
    
    # Find variants near the selected position
    distances = abs.(positions .- selected_pos)
    nearby_mask = distances .<= window_size
    
    if sum(nearby_mask) == 0
        return Float64[], Int[]
    end
    
    nearby_positions = positions[nearby_mask]
    nearby_genotypes = genotypes[:, nearby_mask]
    
    # Calculate diversity in windows around selected site
    diversities = Float64[]
    window_positions = Int[]
    
    for pos in nearby_positions
        # Calculate diversity for this position
        allele_freq = mean(nearby_genotypes[:, findfirst(==(pos), nearby_positions)]) / 2
        diversity = 2 * allele_freq * (1 - allele_freq)
        push!(diversities, diversity)
        push!(window_positions, pos)
    end
    
    return diversities, window_positions
end

# Analyze diversity around selected site for directional selection
if stats_directional.n_sites > 0
    div_directional, pos_directional = analyze_selection_signature(
        genotypes_diploid_directional, positions_directional, selected_position)
    
    if length(div_directional) > 0
        println("Directional selection signature:")
        println("  Mean diversity near selected site: $(round(mean(div_directional), digits=6))")
        println("  Minimum diversity: $(round(minimum(div_directional), digits=6))")
    end
end

# Analyze diversity for balancing selection
if stats_balancing.n_sites > 0
    div_balancing, pos_balancing = analyze_selection_signature(
        genotypes_diploid_balancing, positions_balancing, selected_position)
    
    if length(div_balancing) > 0
        println("Balancing selection signature:")
        println("  Mean diversity near selected site: $(round(mean(div_balancing), digits=6))")
        println("  Maximum diversity: $(round(maximum(div_balancing), digits=6))")
    end
end

println()

# 7. Save results
println("7. Saving Results")
println("-" * 18)

save_genotypes(genotypes_diploid_neutral, positions_neutral, 
               "neutral_selection.csv", format="csv")
save_genotypes(genotypes_diploid_directional, positions_directional, 
               "directional_selection.csv", format="csv")
save_genotypes(genotypes_diploid_background, positions_background, 
               "background_selection.csv", format="csv")
save_genotypes(genotypes_diploid_balancing, positions_balancing, 
               "balancing_selection.csv", format="csv")

println("Results saved:")
println("  neutral_selection.csv - Neutral evolution")
println("  directional_selection.csv - Selective sweep")
println("  background_selection.csv - Background selection")
println("  balancing_selection.csv - Balancing selection")

println("\nSelection example complete!")
println("\nKey observations:")
println("- Directional selection reduces diversity (negative Tajima's D)")
println("- Background selection reduces overall diversity")
println("- Balancing selection maintains high diversity (positive Tajima's D)")
println("- Selection signatures are detectable in local diversity patterns")