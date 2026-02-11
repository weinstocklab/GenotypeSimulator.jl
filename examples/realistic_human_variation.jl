#!/usr/bin/env julia

"""
Realistic Human Genetic Variation Example

This script demonstrates how to simulate realistic human genetic variation with:
- Common and low-frequency SNPs
- Weak directional selection on low-frequency variants
- Realistic recombination patterns
- Large sample size (1,000 individuals)
- 1Mb genomic region

This example is suitable for:
- Population genetics research
- GWAS simulation studies
- Method development and validation
- Educational purposes
"""

using Pkg
Pkg.activate(".")

using GenotypeSimulator
using Random

# Set seed for reproducibility
Random.seed!(123)

println("="^70)
println("Realistic Human Genetic Variation Simulation")
println("="^70)

# Realistic human population parameters for 1Mb region
params = PopulationParams(
    10_000,     # Effective population size (human)
    1.25e-8,    # Human mutation rate per bp per generation
    1e-8,       # Human recombination rate per bp per generation
    1_000_000,  # 1Mb sequence length
    1000        # Large sample size (1,000 individuals)
)

println("Simulation Parameters:")
println("  Effective population size: $(params.ne)")
println("  Mutation rate: $(params.mutation_rate) per bp per generation")
println("  Recombination rate: $(params.recombination_rate) per bp per generation")
println("  Sequence length: $(params.sequence_length/1e6) Mb")
println("  Sample size: $(params.sample_size) individuals ($(2*params.sample_size) haplotypes)")
println()

# Step 1: Simulate neutral background variation with recombination
println("1. NEUTRAL BACKGROUND VARIATION")
println("-" * 35)
println("Simulating neutral evolution with recombination...")

neutral_genotypes, neutral_positions = simulate_with_recombination(params)
neutral_diploid = haplotypes_to_diploid(neutral_genotypes)
neutral_stats = calculate_stats(neutral_diploid, neutral_positions)

println("\nNeutral simulation results:")
println("  Total variant sites: $(neutral_stats.n_sites)")
println("  Common variants (MAF > 0.05): $(neutral_stats.n_common_variants)")
println("  Rare variants (MAF ≤ 0.05): $(neutral_stats.n_sites - neutral_stats.n_common_variants)")
println("  Nucleotide diversity (π): $(round(neutral_stats.nucleotide_diversity, digits=6))")
println("  Watterson's theta (θw): $(round(neutral_stats.theta_w, digits=6))")
println("  Tajima's D: $(round(tajimas_d(neutral_diploid), digits=4))")

# Calculate allele frequency spectrum for neutral
neutral_sfs = allele_frequency_spectrum(neutral_diploid)
neutral_singletons = length(neutral_positions) > 0 ? neutral_sfs[2] : 0
neutral_doubletons = length(neutral_positions) > 0 ? neutral_sfs[3] : 0

println("  Singletons: $neutral_singletons")
println("  Doubletons: $neutral_doubletons")

# Step 2: Add low-frequency variants under weak directional selection
println("\n2. SELECTION ON LOW-FREQUENCY VARIANTS")
println("-" * 40)

# Define multiple selected sites with realistic weak selection
selected_sites = [200_000, 400_000, 600_000, 800_000]  # Evenly spaced across 1Mb
selection_coefficients = [0.005, 0.008, 0.003, 0.006]  # Weak selection (0.3-0.8%)
dominance_coefficients = [0.5, 0.3, 0.7, 0.5]          # Mix of dominance patterns

println("Simulating selection at $(length(selected_sites)) sites:")
for (i, (pos, s, h)) in enumerate(zip(selected_sites, selection_coefficients, dominance_coefficients))
    dominance_type = h == 0.5 ? "additive" : h < 0.5 ? "recessive" : "dominant"
    println("  Site $(i): position $(pos), s = $(s) ($(dominance_type))")
end

# For demonstration, we'll simulate one selected site in detail
# In practice, you might want to simulate multiple sites simultaneously
selected_pos = selected_sites[1]
selected_s = selection_coefficients[1]
selected_h = dominance_coefficients[1]

println("\nDetailed simulation for site at position $selected_pos:")
println("  Selection coefficient: $selected_s")
println("  Dominance coefficient: $selected_h")

# Create directional selection model
selection_model = DirectionalSelection(
    selected_s,      # Selection coefficient (0.5% advantage)
    selected_h,      # Dominance coefficient
    selected_pos,    # Selected position
    5000.0,          # Selection started 5000 generations ago
    1/(2*params.ne)  # Started from single mutation
)

selection_params = SelectionParameters(selection_model, params.ne)

# Simulate with both selection and recombination
println("Running simulation with selection and recombination...")
selected_genotypes, selected_positions = simulate_with_selection(params, selection_params)
selected_diploid = haplotypes_to_diploid(selected_genotypes)
selected_stats = calculate_stats(selected_diploid, selected_positions)

println("\nSelection simulation results:")
println("  Total variant sites: $(selected_stats.n_sites)")
println("  Common variants (MAF > 0.05): $(selected_stats.n_common_variants)")
println("  Rare variants (MAF ≤ 0.05): $(selected_stats.n_sites - selected_stats.n_common_variants)")
println("  Nucleotide diversity (π): $(round(selected_stats.nucleotide_diversity, digits=6))")
println("  Watterson's theta (θw): $(round(selected_stats.theta_w, digits=6))")
println("  Tajima's D: $(round(tajimas_d(selected_diploid), digits=4))")

# Step 3: Detailed allele frequency spectrum analysis
println("\n3. ALLELE FREQUENCY SPECTRUM ANALYSIS")
println("-" * 38)

if length(selected_positions) > 0
    selected_sfs = allele_frequency_spectrum(selected_diploid)
    
    # Calculate frequency categories
    n_singletons = selected_sfs[2]  # Sites with exactly 1 derived allele
    n_doubletons = selected_sfs[3]  # Sites with exactly 2 derived alleles
    n_rare = sum(selected_sfs[2:21])  # Sites with 1-20 derived alleles (rare: MAF < 0.01)
    n_low_freq = sum(selected_sfs[22:101])  # Sites with 21-100 derived alleles (low freq: 0.01 ≤ MAF < 0.05)
    n_common = sum(selected_sfs[102:end])  # Sites with >100 derived alleles (common: MAF ≥ 0.05)
    
    println("Allele frequency categories:")
    println("  Singletons (n=1): $n_singletons")
    println("  Doubletons (n=2): $n_doubletons")
    println("  Rare variants (MAF < 0.01): $n_rare")
    println("  Low-frequency variants (0.01 ≤ MAF < 0.05): $n_low_freq")
    println("  Common variants (MAF ≥ 0.05): $n_common")
    
    # Calculate ratios
    if n_doubletons > 0
        singleton_doubleton_ratio = n_singletons / n_doubletons
        println("  Singleton/Doubleton ratio: $(round(singleton_doubleton_ratio, digits=2))")
    end
else
    println("No variants generated in this simulation run.")
end

# Step 4: Selection signature analysis
println("\n4. SELECTION SIGNATURE ANALYSIS")
println("-" * 32)

if length(selected_positions) > 0
    # Analyze diversity patterns around the selected site
    window_sizes = [10_000, 25_000, 50_000, 100_000]  # Multiple window sizes
    
    println("Diversity analysis around selected site ($selected_pos):")
    
    for window_size in window_sizes
        # Find variants within window
        distances = abs.(selected_positions .- selected_pos)
        in_window = distances .<= window_size
        n_variants_in_window = sum(in_window)
        
        if n_variants_in_window > 0
            window_positions = selected_positions[in_window]
            window_genotypes = selected_diploid[:, in_window]
            
            # Calculate diversity metrics
            diversities = [2 * mean(window_genotypes[:, j])/2 * (1 - mean(window_genotypes[:, j])/2) 
                          for j in 1:size(window_genotypes, 2)]
            
            mean_diversity = mean(diversities)
            min_diversity = minimum(diversities)
            max_diversity = maximum(diversities)
            
            println("  $(window_size/1000)kb window:")
            println("    Variants: $n_variants_in_window")
            println("    Mean diversity: $(round(mean_diversity, digits=6))")
            println("    Min diversity: $(round(min_diversity, digits=6))")
            println("    Max diversity: $(round(max_diversity, digits=6))")
        else
            println("  $(window_size/1000)kb window: No variants")
        end
    end
    
    # Look for the selected variant itself
    closest_idx = argmin(abs.(selected_positions .- selected_pos))
    distance_to_closest = abs(selected_positions[closest_idx] - selected_pos)
    
    if distance_to_closest <= 5000  # Within 5kb of selected site
        selected_allele_freq = mean(selected_diploid[:, closest_idx]) / 2
        println("\n  Selected allele frequency: $(round(selected_allele_freq, digits=4))")
        println("  Distance to closest variant: $(distance_to_closest) bp")
        
        # Check if this is likely the beneficial allele
        if selected_allele_freq < 0.1  # Low frequency as expected
            println("  ✓ Consistent with beneficial allele under weak selection")
        end
    else
        println("\n  No variant found close to selected position")
        println("  (Closest variant is $(distance_to_closest) bp away)")
    end
else
    println("No variants to analyze.")
end

# Step 5: Comparison with neutral expectation
println("\n5. COMPARISON WITH NEUTRAL EXPECTATION")
println("-" * 40)

println("Metric\t\t\t\tNeutral\t\tWith Selection\tChange")
println("-" * 70)

if length(neutral_positions) > 0 && length(selected_positions) > 0
    println("Total variants\t\t\t$(neutral_stats.n_sites)\t\t$(selected_stats.n_sites)\t\t$(selected_stats.n_sites - neutral_stats.n_sites)")
    println("Common variants\t\t\t$(neutral_stats.n_common_variants)\t\t$(selected_stats.n_common_variants)\t\t$(selected_stats.n_common_variants - neutral_stats.n_common_variants)")
    
    pi_change = ((selected_stats.nucleotide_diversity - neutral_stats.nucleotide_diversity) / neutral_stats.nucleotide_diversity) * 100
    println("Nucleotide diversity (π)\t$(round(neutral_stats.nucleotide_diversity, digits=6))\t$(round(selected_stats.nucleotide_diversity, digits=6))\t$(round(pi_change, digits=1))%")
    
    theta_change = ((selected_stats.theta_w - neutral_stats.theta_w) / neutral_stats.theta_w) * 100
    println("Watterson's theta (θw)\t\t$(round(neutral_stats.theta_w, digits=6))\t$(round(selected_stats.theta_w, digits=6))\t$(round(theta_change, digits=1))%")
    
    tajima_neutral = tajimas_d(neutral_diploid)
    tajima_selected = tajimas_d(selected_diploid)
    tajima_change = tajima_selected - tajima_neutral
    println("Tajima's D\t\t\t$(round(tajima_neutral, digits=4))\t\t$(round(tajima_selected, digits=4))\t\t$(round(tajima_change, digits=4))")
else
    println("Insufficient data for comparison")
end

# Step 6: Save results
println("\n6. SAVING RESULTS")
println("-" * 16)

# Save the selection simulation results
save_genotypes(selected_diploid, selected_positions, "realistic_human_1Mb.vcf", format="vcf")
save_genotypes(selected_diploid, selected_positions, "realistic_human_1Mb.csv", format="csv")
save_genotypes(selected_diploid, selected_positions, "realistic_human_1Mb", format="plink")

# Also save neutral for comparison
save_genotypes(neutral_diploid, neutral_positions, "neutral_human_1Mb.vcf", format="vcf")

println("Results saved:")
println("  realistic_human_1Mb.vcf - Selection simulation (VCF format)")
println("  realistic_human_1Mb.csv - Selection simulation (CSV format)")
println("  realistic_human_1Mb.ped/.map - Selection simulation (PLINK format)")
println("  neutral_human_1Mb.vcf - Neutral comparison (VCF format)")

# Step 7: Summary and interpretation
println("\n7. SUMMARY AND INTERPRETATION")
println("-" * 32)

println("This simulation demonstrates:")
println("✓ Realistic human genetic variation patterns")
println("✓ Effects of weak selection on allele frequency spectra")
println("✓ Recombination breaking up linkage disequilibrium")
println("✓ Large sample sizes revealing rare variant patterns")
println("✓ Selection signatures detectable in local diversity")

println("\nKey findings:")
if length(selected_positions) > 0
    tajima_d_val = tajimas_d(selected_diploid)
    if tajima_d_val < -0.5
        println("• Strong negative Tajima's D suggests recent positive selection")
    elseif tajima_d_val < 0
        println("• Negative Tajima's D suggests mild positive selection or population expansion")
    else
        println("• Tajima's D near zero suggests neutral evolution or balancing selection")
    end
    
    rare_fraction = (selected_stats.n_sites - selected_stats.n_common_variants) / selected_stats.n_sites
    println("• $(round(rare_fraction * 100, digits=1))% of variants are rare (MAF ≤ 0.05)")
    
    if selected_stats.n_sites > neutral_stats.n_sites
        println("• Selection increased total number of segregating sites")
    else
        println("• Selection reduced total number of segregating sites")
    end
else
    println("• No variants generated - try increasing mutation rate or sample size")
end

println("\nThis dataset is suitable for:")
println("• GWAS power calculations")
println("• Population genetics method validation")
println("• Demographic inference studies")
println("• Selection detection algorithm testing")

println("\n" * "="^70)
println("Simulation completed successfully!")
println("="^70)