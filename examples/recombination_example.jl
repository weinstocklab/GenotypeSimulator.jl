#!/usr/bin/env julia

"""
Recombination example for GenotypeSimulator.jl

This script demonstrates how to:
1. Simulate genotypes with recombination using ARG
2. Compare results with and without recombination
3. Analyze linkage disequilibrium patterns
"""

using Pkg
Pkg.activate(".")

using GenotypeSimulator
using Random

# Set random seed for reproducibility
Random.seed!(42)

println("GenotypeSimulator.jl - Recombination Example")
println("="^50)

# Define parameters for comparison
params = PopulationParams(
    10_000,     # Effective population size
    1.25e-8,    # Mutation rate per bp per generation
    1e-8,       # Recombination rate per bp per generation
    100_000,    # Sequence length (100kb for faster demo)
    50          # Sample size (50 individuals)
)

println("Parameters:")
println("  Ne: $(params.ne)")
println("  μ: $(params.mutation_rate)")
println("  ρ: $(params.recombination_rate)")
println("  Length: $(params.sequence_length) bp")
println("  Samples: $(params.sample_size) individuals")
println()

# Simulation without recombination
println("1. Simulation WITHOUT recombination (standard coalescent)")
println("-" * 55)

genotypes_no_recomb, positions_no_recomb = simulate_genotypes(params)
genotypes_diploid_no_recomb = haplotypes_to_diploid(genotypes_no_recomb)
stats_no_recomb = calculate_stats(genotypes_diploid_no_recomb, positions_no_recomb)

println("Results without recombination:")
println("  Variant sites: $(stats_no_recomb.n_sites)")
println("  Nucleotide diversity (π): $(round(stats_no_recomb.nucleotide_diversity, digits=6))")
println("  Tajima's D: $(round(tajimas_d(genotypes_diploid_no_recomb), digits=4))")
println()

# Simulation with recombination
println("2. Simulation WITH recombination (ARG)")
println("-" * 40)

genotypes_recomb, positions_recomb = simulate_with_recombination(params)
genotypes_diploid_recomb = haplotypes_to_diploid(genotypes_recomb)
stats_recomb = calculate_stats(genotypes_diploid_recomb, positions_recomb)

println("Results with recombination:")
println("  Variant sites: $(stats_recomb.n_sites)")
println("  Nucleotide diversity (π): $(round(stats_recomb.nucleotide_diversity, digits=6))")
println("  Tajima's D: $(round(tajimas_d(genotypes_diploid_recomb), digits=4))")
println()

# Compare results
println("3. Comparison")
println("-" * 15)
println("Metric\t\t\tNo Recomb\tWith Recomb\tRatio")
println("-" * 55)
println("Variant sites\t\t$(stats_no_recomb.n_sites)\t\t$(stats_recomb.n_sites)\t\t$(round(stats_recomb.n_sites/stats_no_recomb.n_sites, digits=2))")
println("Nucleotide diversity\t$(round(stats_no_recomb.nucleotide_diversity, digits=6))\t$(round(stats_recomb.nucleotide_diversity, digits=6))\t$(round(stats_recomb.nucleotide_diversity/stats_no_recomb.nucleotide_diversity, digits=2))")
println()

# Simple linkage disequilibrium analysis
function calculate_simple_ld(genotypes::Matrix{Int}, positions::Vector{Int}, max_distance::Int=10000)
    n_sites = length(positions)
    if n_sites < 2
        return Float64[], Int[]
    end
    
    r2_values = Float64[]
    distances = Int[]
    
    for i in 1:min(n_sites-1, 100)  # Limit to first 100 sites for speed
        for j in i+1:min(i+20, n_sites)  # Check up to 20 sites ahead
            dist = positions[j] - positions[i]
            if dist > max_distance
                break
            end
            
            # Calculate r² between sites i and j
            x = genotypes[:, i]
            y = genotypes[:, j]
            
            # Simple correlation coefficient squared
            if var(x) > 0 && var(y) > 0
                r = cor(x, y)
                push!(r2_values, r^2)
                push!(distances, dist)
            end
        end
    end
    
    return r2_values, distances
end

println("4. Linkage Disequilibrium Analysis")
println("-" * 35)

if stats_no_recomb.n_sites >= 2
    r2_no_recomb, dist_no_recomb = calculate_simple_ld(genotypes_diploid_no_recomb, positions_no_recomb)
    if length(r2_no_recomb) > 0
        println("Without recombination:")
        println("  Mean r²: $(round(mean(r2_no_recomb), digits=4))")
        println("  Max r²: $(round(maximum(r2_no_recomb), digits=4))")
    end
end

if stats_recomb.n_sites >= 2
    r2_recomb, dist_recomb = calculate_simple_ld(genotypes_diploid_recomb, positions_recomb)
    if length(r2_recomb) > 0
        println("With recombination:")
        println("  Mean r²: $(round(mean(r2_recomb), digits=4))")
        println("  Max r²: $(round(maximum(r2_recomb), digits=4))")
    end
end

println()

# Save results
println("5. Saving Results")
println("-" * 18)

save_genotypes(genotypes_diploid_no_recomb, positions_no_recomb, 
               "no_recombination.csv", format="csv")
save_genotypes(genotypes_diploid_recomb, positions_recomb, 
               "with_recombination.csv", format="csv")

println("Results saved:")
println("  no_recombination.csv - Standard coalescent simulation")
println("  with_recombination.csv - ARG simulation with recombination")

println("\nRecombination example complete!")
println("\nKey observations:")
println("- Recombination typically increases genetic diversity")
println("- ARG simulation produces more realistic linkage patterns")
println("- Recombination breaks up long-range linkage disequilibrium")