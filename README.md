# GenotypeSimulator.jl

A Julia package for simulating realistic human genotypes using coalescent theory, inspired by the tskit/msprime ecosystem.

## Features

- **Coalescent simulation**: Generate realistic genealogical trees
- **Recombination modeling**: Ancestral Recombination Graph (ARG) simulation
- **Selection modeling**: Directional, balancing, and background selection
- **Mutation modeling**: Add mutations using the infinite sites model
- **Population genetics statistics**: Calculate π, θw, Tajima's D, and more
- **Multiple output formats**: CSV, VCF, and PLINK formats
- **Efficient implementation**: Optimized for performance with Julia
- **Extensible design**: Modular structure for easy customization

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/GenotypeSimulator.jl")
```

Or clone and develop locally:

```bash
git clone https://github.com/yourusername/GenotypeSimulator.jl
cd GenotypeSimulator.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Quick Start

```julia
using GenotypeSimulator
using Random

# Set seed for reproducibility
Random.seed!(42)

# Run simulation with default human parameters (1Mb, 100 individuals)
genotypes_hap, positions = simulate_genotypes()

# Convert to diploid genotypes
genotypes_diploid = haplotypes_to_diploid(genotypes_hap)

# Calculate population genetics statistics
stats = calculate_stats(genotypes_diploid, positions)
print_stats(stats)

# Save results
save_genotypes(genotypes_diploid, positions, "results.vcf", format="vcf")
```

## Custom Parameters

```julia
# Define custom population parameters
custom_params = PopulationParams(
    5000,      # Effective population size
    2e-8,      # Mutation rate per bp per generation
    1e-8,      # Recombination rate per bp per generation
    500_000,   # Sequence length (500kb)
    50         # Sample size (50 individuals)
)

# Run simulation
genotypes_hap, positions = simulate_genotypes(custom_params)

# Or with recombination (ARG simulation)
genotypes_hap, positions = simulate_with_recombination(custom_params)

# Or with selection
selection_model = DirectionalSelection(0.01, 0.5, 25000, 1000.0)  # 1% advantage, additive
selection_params = SelectionParameters(selection_model, custom_params.ne)
genotypes_hap, positions = simulate_with_selection(custom_params, selection_params)
```

## Realistic Human Genetic Variation Example

Here's a comprehensive example simulating realistic human genetic variation with both common and low-frequency variants under selection:

```julia
using GenotypeSimulator
using Random

# Set seed for reproducibility
Random.seed!(123)

# Realistic human population parameters for 1Mb region
params = PopulationParams(
    10_000,     # Effective population size (human)
    1.25e-8,    # Human mutation rate per bp per generation
    1e-8,       # Human recombination rate per bp per generation
    1_000,  # 1kb sequence length
    10        # Large sample size (1,000 individuals)
)

println("Simulating realistic human genetic variation...")
println("Sample: $(params.sample_size) individuals")
println("Region: $(params.sequence_length/1e6) Mb")

# Step 1: Simulate neutral background variation with recombination
println("\n1. Simulating neutral background variation...")
neutral_genotypes, neutral_positions = simulate_with_recombination(params)
neutral_diploid = haplotypes_to_diploid(neutral_genotypes)
neutral_stats = calculate_stats(neutral_diploid, neutral_positions)

println("Neutral simulation results:")
println("  Total variants: $(neutral_stats.n_sites)")
println("  Common variants (MAF > 0.05): $(neutral_stats.n_common_variants)")
println("  Nucleotide diversity (π): $(round(neutral_stats.nucleotide_diversity, digits=6))")

# Step 2: Add low-frequency variants under weak directional selection
println("\n2. Adding low-frequency variants under selection...")

# Simulate multiple selected sites with weak selection
selected_sites = [200_000, 400_000, 600_000, 800_000]  # 4 selected positions
selection_coefficients = [0.005, 0.008, 0.003, 0.006]  # Weak selection (0.3-0.8%)

all_selected_genotypes = []
all_selected_positions = []

for (i, (pos, s)) in enumerate(zip(selected_sites, selection_coefficients))
    println("  Simulating selection at position $pos (s = $s)...")
    
    # Create directional selection model
    selection_model = DirectionalSelection(
        s,           # Selection coefficient
        0.5,         # Additive (h = 0.5)
        pos,         # Selected position
        5000.0,      # Selection started 5000 generations ago
        1/(2*params.ne)  # Started from single mutation
    )
    
    selection_params = SelectionParameters(selection_model, params.ne)
    
    # Simulate with both selection and recombination
    # Note: This is a simplified approach - in reality, we'd need to model
    # multiple selected sites simultaneously
    selected_genotypes, selected_positions = simulate_with_selection(params, selection_params)
    
    push!(all_selected_genotypes, selected_genotypes)
    push!(all_selected_positions, selected_positions)
end

# Step 3: Analyze the combined variation
println("\n3. Analyzing genetic variation patterns...")

# Combine results (simplified - taking the first selection simulation as example)
combined_genotypes = all_selected_genotypes[1]
combined_positions = all_selected_positions[1]
combined_diploid = haplotypes_to_diploid(combined_genotypes)
combined_stats = calculate_stats(combined_diploid, combined_positions)

println("Combined simulation results:")
println("  Total variants: $(combined_stats.n_sites)")
println("  Common variants (MAF > 0.05): $(combined_stats.n_common_variants)")
println("  Nucleotide diversity (π): $(round(combined_stats.nucleotivity, digits=6))")
println("  Tajima's D: $(round(tajimas_d(combined_diploid), digits=4))")

# Step 4: Analyze allele frequency spectrum
println("\n4. Allele frequency spectrum analysis...")
sfs = allele_frequency_spectrum(combined_diploid)
n_singletons = sfs[2]  # Sites with exactly 1 derived allele
n_doubletons = sfs[3]  # Sites with exactly 2 derived alleles
total_rare = sum(sfs[2:11])  # Sites with 1-10 derived alleles (rare variants)

println("Frequency spectrum:")
println("  Singletons: $n_singletons")
println("  Doubletons: $n_doubletons") 
println("  Rare variants (1-10 copies): $total_rare")
println("  Common variants (>10 copies): $(combined_stats.n_sites - total_rare)")

# Step 5: Selection signature analysis
println("\n5. Selection signature analysis...")
for (i, pos) in enumerate(selected_sites[1:1])  # Analyze first selected site
    window_size = 50_000  # 50kb window
    
    # Find variants near selected site
    distances = abs.(combined_positions .- pos)
    nearby_variants = distances .<= window_size
    
    if sum(nearby_variants) > 0
        nearby_positions = combined_positions[nearby_variants]
        nearby_genotypes = combined_diploid[:, nearby_variants]
        
        # Calculate diversity in windows
        diversities = [2 * mean(nearby_genotypes[:, j])/2 * (1 - mean(nearby_genotypes[:, j])/2) 
                      for j in 1:size(nearby_genotypes, 2)]
        
        mean_diversity = mean(diversities)
        min_diversity = minimum(diversities)
        
        println("  Selected site $pos:")
        println("    Variants in 50kb window: $(sum(nearby_variants))")
        println("    Mean diversity in window: $(round(mean_diversity, digits=6))")
        println("    Minimum diversity: $(round(min_diversity, digits=6))")
        
        # Look for the selected variant itself
        closest_idx = argmin(distances)
        if distances[closest_idx] <= 1000  # Within 1kb
            selected_freq = mean(combined_diploid[:, closest_idx]) / 2
            println("    Selected allele frequency: $(round(selected_freq, digits=4))")
        end
    end
end

# Step 6: Save results in multiple formats
println("\n6. Saving results...")
save_genotypes(combined_diploid, combined_positions, "human_variation.vcf", format="vcf")
save_genotypes(combined_diploid, combined_positions, "human_variation.csv", format="csv")
save_genotypes(combined_diploid, combined_positions, "human_variation", format="plink")

println("\nSimulation complete! Files saved:")
println("  human_variation.vcf - VCF format for analysis")
println("  human_variation.csv - CSV format for data processing")
println("  human_variation.ped/.map - PLINK format for association studies")

# Step 7: Summary statistics comparison
println("\n7. Summary comparison:")
println("Metric\t\t\tNeutral\t\tWith Selection")
println("-" * 50)
println("Total variants\t\t$(neutral_stats.n_sites)\t\t$(combined_stats.n_sites)")
println("Common variants\t\t$(neutral_stats.n_common_variants)\t\t$(combined_stats.n_common_variants)")
println("Nucleotide diversity\t$(round(neutral_stats.nucleotide_diversity, digits=6))\t$(round(combined_stats.nucleotide_diversity, digits=6))")
println("Tajima's D\t\t$(round(tajimas_d(neutral_diploid), digits=4))\t\t$(round(tajimas_d(combined_diploid), digits=4))")

println("\nKey observations:")
println("- Selection reduces local diversity around selected sites")
println("- Low-frequency beneficial variants create negative Tajima's D")
println("- Recombination breaks up linkage disequilibrium")
println("- Large sample sizes reveal rare variant patterns")
```

This example demonstrates:
- **Realistic parameters**: Human-specific mutation and recombination rates
- **Large sample size**: 1,000 individuals to capture rare variants
- **Multiple selected sites**: Weak selection on low-frequency variants
- **Recombination effects**: ARG simulation with realistic crossover rates
- **Comprehensive analysis**: Frequency spectra, selection signatures, and diversity patterns
- **Multiple output formats**: Ready for downstream analysis tools

## Development Workflow

### Using Revise.jl for Development

```julia
using Pkg
Pkg.add("Revise")
using Revise

# Activate the project environment
Pkg.activate(".")

# Load the package with hot reloading
using GenotypeSimulator

# Now you can modify source files and changes will be automatically loaded
```

### Running Tests

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

### Examples

Check out the `examples/` directory for more detailed usage:

- `examples/basic_simulation.jl`: Basic usage example
- `examples/custom_parameters.jl`: Comparing different population scenarios
- `examples/recombination_example.jl`: Recombination and ARG simulation
- `examples/selection_example.jl`: Selection models and signatures
- `examples/realistic_human_variation.jl`: Large-scale human variation with selection

## Compilation to Binary

### Install Julia 1.12 with juliaup

```bash
# Install juliaup
curl -fsSL https://install.julialang.org | sh

# Install Julia 1.12
juliaup add 1.12
juliaup default 1.12
```

### Compile with juliac

Using the Makefile (recommended):

```bash
# Install dependencies and compile
make compile

# Or run individual steps
make deps
make compile

# Install system-wide (Unix/Linux/macOS)
make install
```

Using the Julia script directly:

```bash
julia --project=. scripts/compile.jl
```

This creates a compiled binary `juliac/genotype_simulator` with:
- Ahead-of-time compiled Julia code
- Fast startup times
- Reduced memory footprint
- Full command-line interface

### Command-line Usage

```bash
# Basic simulation
./juliac/genotype_simulator --samples 50 --output my_sim

# Custom parameters
./juliac/genotype_simulator --ne 5000 --mu 2e-8 --format vcf --stats

# With recombination (ARG simulation)
./juliac/genotype_simulator --recombination --rho 1e-8 --format vcf

# With directional selection (selective sweep)
./juliac/genotype_simulator --selection directional --selection-coeff 0.05 --format vcf

# See all options
./juliac/genotype_simulator --help
```

## API Reference

### Core Types

- `PopulationParams`: Parameters for population genetic simulation
- `CoalescentNode`: Node in the coalescent tree
- `PopulationStats`: Population genetics statistics

### Main Functions

- `simulate_genotypes(params)`: Run complete simulation
- `haplotypes_to_diploid(haplotypes)`: Convert haplotype to diploid matrix
- `calculate_stats(genotypes, positions)`: Calculate population genetics statistics
- `save_genotypes(genotypes, positions, filename)`: Save results to file

### Statistics Functions

- `allele_frequency_spectrum(genotypes)`: Calculate site frequency spectrum
- `tajimas_d(genotypes)`: Calculate Tajima's D statistic
- `print_stats(stats)`: Pretty-print statistics

## Comparison with tskit/msprime

| Feature | tskit/msprime | GenotypeSimulator.jl |
|---------|---------------|---------------------|
| Coalescent simulation | ✓ Full featured | ✓ Simplified |
| Mutation models | ✓ Multiple models | ✓ Infinite sites |
| Recombination | ✓ Variable rates | ✓ ARG simulation |
| Population structure | ✓ Complex demography | ○ Single population |
| Selection models | ○ Limited | ✓ Multiple types |
| Performance | ✓ Highly optimized | ✓ Good performance |
| Tree sequences | ✓ Efficient storage | ○ Basic tree |

## Methods

### Theoretical Background

GenotypeSimulator.jl implements coalescent theory-based simulation algorithms to generate realistic patterns of genetic variation. The coalescent process (Kingman, 1982) models the genealogical relationships among sampled sequences by tracing lineages backward in time until they reach common ancestors.

### Standard Coalescent Simulation

The standard coalescent assumes no recombination and models the genealogy as a binary tree. Under the Wright-Fisher model with effective population size Ne, the waiting time until the next coalescence event among k lineages follows an exponential distribution with rate λk = k(k-1)/(4Ne).

**Algorithm 1: Standard Coalescent Tree Construction**
```
Input: sample_size n, effective_population_size Ne
Output: coalescent_tree T

1. Initialize n leaf nodes representing present-day samples
2. Set active_lineages ← {1, 2, ..., n}, current_time ← 0
3. While |active_lineages| > 1:
   4. k ← |active_lineages|
   5. rate ← k(k-1)/(4Ne)
   6. Δt ~ Exponential(1/rate)
   7. current_time ← current_time + Δt
   8. Randomly select two lineages i, j from active_lineages
   9. Create parent node p with time current_time
   10. Set children[p] ← {i, j}, parent[i] ← p, parent[j] ← p
   11. active_lineages ← active_lineages \ {i, j} ∪ {p}
12. Return root node from active_lineages
```

### Ancestral Recombination Graph (ARG) Simulation

The ARG extends the coalescent to include recombination events (Hudson, 1983; Griffiths & Marjoram, 1997). Each lineage is associated with genomic intervals, and recombination events split lineages at breakpoints, creating different local genealogies across the sequence.

**Key Implementation Detail**: To produce theoretically correct segregating site counts independent of recombination rate, mutations must be placed on **marginal trees** (the local genealogy at each position) rather than on the full ARG tree. This prevents spurious "fixed" mutations from branches above the local MRCA at each position. GenotypeSimulator.jl implements this by tracing sample lineages at each genomic interval and placing mutations only on branches carrying a proper subset of samples.

**Algorithm 2: ARG Construction with Sequential Markov Coalescent**
```
Input: sample_size n, Ne, sequence_length L, recombination_rate ρ
Output: ancestral_recombination_graph ARG

1. Initialize n lineages, each covering interval [1, L+1)
2. Create recombination_map with rate ρ
3. Set current_time ← 0
4. While multiple lineages exist:
   5. k ← number of active lineages
   6. total_length ← sum of all lineage interval lengths
   7. λ_coal ← k(k-1)/(4Ne)
   8. λ_recomb ← k × ρ × total_length
   9. Δt ~ Exponential(λ_coal + λ_recomb)
   10. current_time ← current_time + Δt
   11. If Uniform(0,1) < λ_coal/(λ_coal + λ_recomb):
       12. // Coalescence event
       13. Select two overlapping lineages i, j
       14. Create common ancestor for overlapping regions
       15. Update lineage intervals and genealogy
   16. Else:
       17. // Recombination event  
       18. Select lineage i and breakpoint position x
       19. Split lineage i at position x into left and right parts
       20. Add new lineage for right part
21. Return ARG with local trees
```

### Recombination Breakpoint Sampling

Recombination breakpoints are sampled from variable rate maps following the approach of McVean et al. (2004). The implementation uses inverse transform sampling on the cumulative recombination rate function.

**Algorithm 3: Efficient Breakpoint Sampling**
```
Input: recombination_map M = {(s₁,e₁,r₁), (s₂,e₂,r₂), ...}
Output: breakpoint_position x

1. Compute cumulative rates: R₀ = 0
2. For i = 1 to |M|:
   3. Rᵢ ← Rᵢ₋₁ + rᵢ × (eᵢ - sᵢ)
4. Sample u ~ Uniform(0, R|M|)
5. Find interval j where Rⱼ₋₁ < u ≤ Rⱼ using binary search
6. If rⱼ > 0:
   7. relative_position ← (u - Rⱼ₋₁)/(rⱼ × (eⱼ - sⱼ))
   8. x ← sⱼ + relative_position × (eⱼ - sⱼ)
9. Else: x ~ Uniform(sⱼ, eⱼ)
10. Return ⌊x⌋
```

### Mutation Process: Marginal Tree Approach

Mutations are placed on genealogies using the infinite sites model (Kimura, 1969), but with a key improvement for accurate ARG simulation: mutations are placed on **marginal trees** rather than the full inflated ARG tree. This avoids spurious fixed sites that would arise from placing mutations on branches above local MRCAs.

**Algorithm 4: Mutation Placement on Marginal Trees**
```
Input: ARG, mutation_rate μ, sequence_length L, recombination_events
Output: genotypes matrix

1. breakpoints ← {1} ∪ {recombination positions} ∪ {L+1}
2. For each interval [bp_i, bp_{i+1}) between consecutive breakpoints:
   3.   test_pos ← bp_i  // Representative position for this interval
   4.   For each sample haplotype h:
   5.     Trace lineage(h) backward at test_pos through ARG
   6.     Collect all nodes visited, following recombination redirects
   7.   Build marginal_tree with these nodes
8.   For each branch b in marginal_tree carrying proper subset of samples:
   9.     branch_length ← time[parent[b]] - time[b]
   10.    expected_mutations ← μ × (bp_{i+1} - bp_i) × branch_length  
   11.    num_mutations ~ Poisson(expected_mutations)
   12.    Place mutations uniformly in interval and record carriers
13. Combine mutations across all intervals to build complete genotype matrix
```

### Marginal Tree Tracing from ARG

The ARG represents the complete genealogical history, but the marginal tree at each genomic position is determined by tracing sample lineages backward while respecting recombination redirects. This is more efficient than full tree extraction and naturally handles the complex topology of the ARG.

**Algorithm 5: Marginal Tree Tracing**
```
Input: ARG, test_position p, sample lineages
Output: marginal_tree edges and branch lengths

1. For each sample haplotype h:
   2.   node ← leaf[h]
   3.   While node.parent exists:
   4.     Check if node's intervals cover p
   5.     If not covered, stop (lineage coalesced at position p)
   6.     Check if node was recombined at p
   7.     If yes (redirect to fragment at p), follow redirect
   8.     Record branch (node, parent) as part of marginal tree
   9.     Move to parent
10. Collect unique branches into marginal_tree
11. Return branch lengths and sample sets
```

This approach automatically constructs the marginal tree without explicit storage, avoiding the memory overhead of full tree extraction while correctly handling recombination topology.

### Population Genetics Statistics

The package computes standard population genetics summary statistics following established formulations:

**Nucleotide Diversity (π)**:
π = (1/L) × Σᵢ₌₁ˢ 2pᵢ(1-pᵢ)

where S is the number of segregating sites, L is sequence length, and pᵢ is the derived allele frequency at site i.

**Watterson's Estimator (θw)**:
θw = S/(L × aₙ), where aₙ = Σⱼ₌₁ⁿ⁻¹ 1/j

**Tajima's D Statistic**:
D = (π - θw)/√Var(π - θw)

The variance calculation follows Tajima (1989) with correction terms for finite sample size.

**Site Frequency Spectrum (SFS)**:
The SFS counts the number of sites with i derived alleles: SFS[i] = |{sites with exactly i derived alleles}|

### Computational Complexity Analysis

| Algorithm | Time Complexity | Space Complexity | Notes |
|-----------|----------------|------------------|-------|
| Standard Coalescent | O(n log n) | O(n) | Dominated by tree construction |
| ARG Simulation | O(n log n + R log n) | O(n + R) | R = number of recombination events |
| Mutation Placement | O(M) | O(M) | M = total number of mutations |
| Genotype Extraction | O(n × S) | O(n × S) | S = number of segregating sites |
| Statistics Calculation | O(n × S) | O(S) | Vectorized operations in Julia |

**Expected Values**:
- E[R] ≈ ρ × L × log(n) for the number of recombination events
- E[M] ≈ μ × L × 2n for the total number of mutations under the standard coalescent
- E[S] ≈ μ × L × Σⱼ₌₁ⁿ⁻¹ 1/j for segregating sites (Watterson's formula), **independent of ρ** when using marginal tree mutation placement

### Selection Models

The package implements several forms of natural selection that modify the standard neutral coalescent process.

**Directional Selection (Selective Sweeps)**:
Under directional selection, a beneficial mutation with selection coefficient s increases in frequency according to:

dp/dt = sp(1-p)[h + (1-2h)p]

where p is the allele frequency, h is the dominance coefficient, and t is time in generations.

**Algorithm 6: Directional Selection Simulation**
```
Input: selection_coefficient s, dominance h, selected_position, start_time
Output: modified_coalescent_tree

1. Calculate allele frequency trajectory: p(t) using logistic growth
2. Modify effective population size: Ne_eff(t) = Ne / (1 + 4sp(1-p))
3. Build coalescent tree with time-varying Ne_eff
4. Add beneficial mutation at selected position and start_time
5. Propagate mutation through genealogy
```

**Background Selection**:
Deleterious mutations reduce the effective population size according to:
Ne_eff = Ne × exp(-U(h + (1-h)s)/s)

where U is the deleterious mutation rate per genome.

**Balancing Selection**:
Maintains multiple alleles at intermediate frequencies through frequency-dependent selection or overdominance.

### Numerical Stability and Implementation Details

**Random Number Generation**: The implementation uses Julia's high-quality Mersenne Twister PRNG with proper seeding for reproducibility.

**Floating Point Precision**: All time calculations use 64-bit floating point arithmetic. For very large effective population sizes (Ne > 10⁶), care is taken to avoid numerical underflow in exponential rate calculations.

**Memory Management**: The tree structures use mutable references to minimize memory allocation during simulation. Large simulations benefit from Julia's garbage collector optimizations.

**Vectorization**: Statistical calculations leverage Julia's efficient array operations and SIMD instructions where applicable.

## Default Human Parameters

The package includes realistic default parameters for human populations based on empirical estimates:

- Effective population size (Ne): 10,000 (Tenesa et al., 2007)
- Mutation rate: 1.25×10⁻⁸ per bp per generation (Nachman & Crowell, 2000)
- Recombination rate: 1×10⁻⁸ per bp per generation (Kong et al., 2002)
- Sequence length: 1,000,000 bp (1Mb)
- Sample size: 100 individuals (200 haplotypes)

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## References

### Theoretical Foundations

- Kingman, J. F. C. (1982). The coalescent. *Stochastic Processes and their Applications*, 13(3), 235-248. [Foundational coalescent theory]

- Hudson, R. R. (1983). Properties of a neutral allele model with intragenic recombination. *Theoretical Population Biology*, 23(2), 183-201. [Early ARG theory]

- Griffiths, R. C., & Marjoram, P. (1997). An ancestral recombination graph. *Progress in Population Genetics and Human Evolution*, 257, 257-270. [ARG formalization]

- Kimura, M. (1969). The number of heterozygous nucleotide sites maintained in a finite population due to steady flux of mutations. *Genetics*, 61(4), 893-903. [Infinite sites model]

### Computational Methods

- Hudson, R. R. (2002). Generating samples under a Wright–Fisher neutral model of genetic variation. *Bioinformatics*, 18(2), 337-338. [ms simulator - classical implementation]

- Kelleher, J., Etheridge, A. M., & McVean, G. (2016). Efficient coalescent simulation and genealogical analysis for large sample sizes. *PLoS Computational Biology*, 12(5), e1004842. [msprime - modern efficient algorithms]

- Kelleher, J., Thornton, K. R., Ashander, J., & Ralph, P. L. (2018). Efficient pedigree recording for fast population genetics simulation. *PLoS Computational Biology*, 14(11), e1006581. [Tree sequence data structure]

- Haller, B. C., & Messer, P. W. (2019). SLiM 3: forward genetic simulations beyond the Wright–Fisher model. *Molecular Biology and Evolution*, 36(3), 632-637. [Forward simulation comparison]

### Population Genetics Statistics

- Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics*, 123(3), 585-595. [Tajima's D statistic]

- Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. *Theoretical Population Biology*, 7(2), 256-276. [Watterson's theta]

- Fu, Y. X., & Li, W. H. (1993). Statistical tests of neutrality of mutations. *Genetics*, 133(3), 693-709. [Additional neutrality tests]

### Human Genetic Parameters

- Nachman, M. W., & Crowell, S. L. (2000). Estimate of the mutation rate per nucleotide in humans. *Genetics*, 156(1), 297-304. [Human mutation rate estimates]

- Kong, A., Gudbjartsson, D. F., Sainz, J., Jonsdottir, G. M., Gudjonsson, S. A., Richardsson, B., ... & Stefansson, K. (2002). A high-resolution recombination map of the human genome. *Nature Genetics*, 31(3), 241-247. [Human recombination rates]

- McVean, G. A., Myers, S. R., Hunt, S., Deloukas, P., Bentley, D. R., & Donnelly, P. (2004). The fine-scale structure of recombination rate variation in the human genome. *Science*, 304(5670), 581-584. [Fine-scale recombination mapping]

- Tenesa, A., Navarro, P., Hayes, B. J., Duffy, D. L., Clarke, G. M., Goddard, M. E., & Visscher, P. M. (2007). Recent human effective population size estimated from linkage disequilibrium. *Genome Research*, 17(4), 520-526. [Human effective population size]

- Schrider, D. R., & Kern, A. D. (2016). S/HIC: robust identification of soft and hard sweeps using machine learning. *PLoS Genetics*, 12(3), e1005928. [Modern population genetics applications]

### Selection Theory

- Kaplan, N. L., Hudson, R. R., & Langley, C. H. (1989). The "hitchhiking effect" revisited. *Genetics*, 123(4), 887-899. [Selective sweeps and hitchhiking]

- Barton, N. H. (1998). The effect of hitch-hiking on neutral genealogies. *Genetical Research*, 72(2), 123-133. [Selection effects on genealogies]

- Charlesworth, B. (1994). The effect of background selection against deleterious mutations on weakly selected, linked variants. *Genetical Research*, 63(3), 213-227. [Background selection theory]

- Etheridge, A., Pfaffelhuber, P., & Wakolbinger, A. (2006). An approximate sampling formula under genetic hitchhiking. *The Annals of Applied Probability*, 16(2), 685-729. [Mathematical foundations of selection in coalescent]

- Gillespie, J. H. (2000). Genetic drift in an infinite population: the pseudohitchhiking model. *Genetics*, 155(2), 909-919. [Pseudohitchhiking and selection]

- Hermisson, J., & Pennings, P. S. (2005). Soft sweeps: molecular population genetics of adaptation from standing genetic variation. *Genetics*, 169(4), 2335-2352. [Soft vs. hard sweeps]

### Software and Implementation

- Bezanson, J., Edelman, A., Karpinski, S., & Shah, V. B. (2017). Julia: A fresh approach to numerical computing. *SIAM Review*, 59(1), 65-98. [Julia language design]

- Ralph, P., Thornton, K., & Kelleher, J. (2020). Efficiently summarizing relationships in large samples: a general duality between statistics of genealogies and genomes. *Genetics*, 215(3), 779-797. [Tree sequence statistics]

- Baumdicker, F., Bisschop, G., Goldstein, D., Gower, G., Ragsdale, A. P., Tsambos, G., ... & Kelleher, J. (2022). Efficient ancestry and mutation simulation with msprime 1.0. *Genetics*, 220(3), iyab229. [Latest msprime developments]

## Citation

If you use GenotypeSimulator.jl in your research, please cite:

```
GenotypeSimulator.jl: A Julia package for coalescent-based genotype simulation
[Your Name], [Year]
```