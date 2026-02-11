#!/usr/bin/env julia

"""
Basic example of using GenotypeSimulator.jl

This script demonstrates how to:
1. Run a basic simulation with default human parameters
2. Calculate population genetics statistics
3. Save results to different file formats
"""

using Pkg
Pkg.activate(".")

using GenotypeSimulator
using Random

# Set random seed for reproducibility
Random.seed!(42)

println("GenotypeSimulator.jl - Basic Example")
println("="^40)

# Run simulation with default human parameters
println("\n1. Running simulation with default parameters...")
genotypes_hap, positions = simulate_genotypes()

# Convert to diploid genotypes
println("\n2. Converting to diploid genotypes...")
genotypes_diploid = haplotypes_to_diploid(genotypes_hap)

# Calculate and display statistics
println("\n3. Calculating population genetics statistics...")
stats = calculate_stats(genotypes_diploid, positions)
print_stats(stats)

# Calculate additional statistics
println("\n4. Additional statistics...")
sfs = allele_frequency_spectrum(genotypes_diploid)
println("Site frequency spectrum (first 10 bins): $(sfs[1:min(10, length(sfs))])")

tajima_d = tajimas_d(genotypes_diploid)
println("Tajima's D: $(round(tajima_d, digits=4))")

# Save results in different formats
println("\n5. Saving results...")
save_genotypes(genotypes_diploid, positions, "simulation_results.csv", format="csv")
save_genotypes(genotypes_diploid, positions, "simulation_results.vcf", format="vcf")
save_genotypes(genotypes_diploid, positions, "simulation_results", format="plink")

println("\nSimulation complete! Check the output files.")
println("- simulation_results.csv (CSV format)")
println("- simulation_results.vcf (VCF format)")  
println("- simulation_results.ped/.map (PLINK format)")