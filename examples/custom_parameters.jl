#!/usr/bin/env julia

"""
Example with custom population parameters

This script demonstrates how to:
1. Define custom population parameters
2. Run multiple simulations with different parameters
3. Compare results across different scenarios
"""

using Pkg
Pkg.activate(".")

using GenotypeSimulator
using Random

# Set random seed for reproducibility
Random.seed!(123)

println("GenotypeSimulator.jl - Custom Parameters Example")
println("="^50)

# Define different population scenarios
scenarios = [
    ("Small population", PopulationParams(1_000, 1.25e-8, 1e-8, 100_000, 50)),
    ("Large population", PopulationParams(50_000, 1.25e-8, 1e-8, 100_000, 50)),
    ("High mutation rate", PopulationParams(10_000, 5e-8, 1e-8, 100_000, 50)),
    ("Low mutation rate", PopulationParams(10_000, 2.5e-9, 1e-8, 100_000, 50))
]

results = []

for (name, params) in scenarios
    println("\n" * "="^30)
    println("Scenario: $name")
    println("="^30)
    
    # Run simulation
    genotypes_hap, positions = simulate_genotypes(params)
    genotypes_diploid = haplotypes_to_diploid(genotypes_hap)
    
    # Calculate statistics
    stats = calculate_stats(genotypes_diploid, positions)
    tajima_d = tajimas_d(genotypes_diploid)
    
    # Store results
    push!(results, (name, stats, tajima_d))
    
    # Print summary
    println("Variant sites: $(stats.n_sites)")
    println("Nucleotide diversity (π): $(round(stats.nucleotide_diversity, digits=6))")
    println("Tajima's D: $(round(tajima_d, digits=4))")
end

# Compare results
println("\n" * "="^50)
println("COMPARISON ACROSS SCENARIOS")
println("="^50)
println("Scenario\t\t\tSites\tπ\t\tTajima's D")
println("-"^50)

for (name, stats, tajima_d) in results
    println("$(rpad(name, 20))\t$(stats.n_sites)\t$(round(stats.nucleotide_diversity, digits=6))\t$(round(tajima_d, digits=4))")
end

println("\nKey observations:")
println("- Smaller populations tend to have lower diversity")
println("- Higher mutation rates increase the number of variants")
println("- Tajima's D reflects demographic history and selection")