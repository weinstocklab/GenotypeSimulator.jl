"""
    SpectralStats

Spectral properties of the genotype matrix (eigenvalue decomposition of the
genetic relatedness matrix, GRM).
"""
struct SpectralStats
    eigenvalues::Vector{Float64}      # All eigenvalues in descending order
    variance_explained::Vector{Float64} # Fraction of variance per component
    n_eff_components::Float64          # Effective number of components (participation ratio)
    eigenvalue_variance::Float64       # Variance of the eigenvalue distribution
    eigenvalue_entropy::Float64        # Shannon entropy of normalised eigenvalues
    tracy_widom_stat::Float64          # Tracy-Widom-like statistic for top eigenvalue
    top1_variance::Float64             # Variance explained by PC1
    top10_variance::Float64            # Cumulative variance explained by top 10 PCs
    marchenko_pastur_ratio::Float64    # Ratio of top eigenvalue to MP upper edge
end

"""
    PopulationStats

Structure to hold population genetics statistics.
"""
struct PopulationStats
    n_sites::Int
    n_common_variants::Int
    mean_allele_freq::Float64
    nucleotide_diversity::Float64
    variant_density::Float64
    segregating_sites::Int
    theta_w::Float64  # Watterson's theta
    spectral::Union{SpectralStats, Nothing}
end

"""
    calculate_stats(genotypes::Matrix{<:Integer}, positions::Vector{<:Integer}) -> PopulationStats

Calculate basic population genetics statistics from genotype data.
Accepts any integer element type (Int, UInt8, etc.) to avoid type-conversion copies.

# Arguments
- `genotypes::Matrix{<:Integer}`: Diploid genotype matrix (individuals × sites)
- `positions::Vector{<:Integer}`: Positions of variant sites

# Returns
- `PopulationStats`: Structure containing various population genetics statistics
"""
function calculate_stats(genotypes::Matrix{<:Integer}, positions::Vector{<:Integer})
    n_individuals, n_sites = size(genotypes)
    
    if n_sites == 0
        return PopulationStats(0, 0, 0.0, 0.0, 0.0, 0, 0.0, nothing)
    end
    
    # Calculate allele frequencies without column copies
    allele_freqs = Vector{Float64}(undef, n_sites)
    for j in 1:n_sites
        s = 0.0
        @inbounds for i in 1:n_individuals
            s += genotypes[i, j]
        end
        allele_freqs[j] = s / (2 * n_individuals)
    end
    
    # Filter for common variants (MAF > 0.05)
    n_common = 0
    for j in 1:n_sites
        maf = min(allele_freqs[j], 1 - allele_freqs[j])
        if maf > 0.05
            n_common += 1
        end
    end
    
    # Nucleotide diversity (π) - average pairwise differences
    pi = 0.0
    for i in 1:n_sites
        p = allele_freqs[i]
        pi += 2 * p * (1 - p)
    end
    
    # Normalize by sequence length
    sequence_length = Int(positions[end]) - Int(positions[1]) + 1
    pi_per_site = pi / sequence_length
    
    # Watterson's theta (based on number of segregating sites)
    n_haplotypes = 2 * n_individuals
    harmonic_number = sum(1/i for i in 1:(n_haplotypes-1))
    theta_w = n_sites / harmonic_number / sequence_length
    
    # Variant density
    variant_density = n_sites / sequence_length
    
    # Spectral analysis of the genotype matrix
    spectral = calculate_spectral_stats(genotypes)

    return PopulationStats(
        n_sites,
        n_common,
        mean(allele_freqs),
        pi_per_site,
        variant_density,
        n_sites,
        theta_w,
        spectral
    )
end

"""
    calculate_spectral_stats(genotypes::Matrix{<:Integer}) -> Union{SpectralStats, Nothing}

Compute spectral statistics from the genetic relatedness matrix (GRM).

Memory-efficient implementation: builds the n×n GRM by streaming over columns
of the genotype matrix one at a time, so peak memory is O(n²) + O(n) instead
of O(n × p).  Uses Float64 for the GRM accumulation.

# Measures returned
- **eigenvalue_variance**: Var(λ) — higher values indicate more dispersed spectrum
- **eigenvalue_entropy**: Shannon entropy of normalised eigenvalues
- **n_eff_components**: Participation ratio — effective dimensionality
- **tracy_widom_stat**: significance of top eigenvalue beyond MP bulk
- **marchenko_pastur_ratio**: λ_max / MP upper edge
- **top1_variance / top10_variance**: Fraction of total variance from top PCs
"""
function calculate_spectral_stats(genotypes::Matrix{<:Integer})::Union{SpectralStats, Nothing}
    n_individuals, n_sites = size(genotypes)

    # Need at least 2 individuals and 2 sites for meaningful spectral analysis
    if n_individuals < 2 || n_sites < 2
        return nothing
    end

    # First pass: compute per-column mean and std to identify polymorphic sites
    # Only O(p) memory for the vectors
    col_means = Vector{Float64}(undef, n_sites)
    col_stds  = Vector{Float64}(undef, n_sites)
    for j in 1:n_sites
        s = 0.0
        @inbounds for i in 1:n_individuals
            s += genotypes[i, j]
        end
        m = s / n_individuals
        col_means[j] = m
        ss = 0.0
        @inbounds for i in 1:n_individuals
            d = genotypes[i, j] - m
            ss += d * d
        end
        col_stds[j] = sqrt(ss / n_individuals)
    end

    # Identify polymorphic columns
    poly_indices = Int[]
    for j in 1:n_sites
        if col_stds[j] > 1e-12
            push!(poly_indices, j)
        end
    end
    n_poly = length(poly_indices)
    if n_poly < 2
        return nothing
    end

    # Build GRM = (1/p) Σ_j z_j z_j^T  by streaming over columns
    # Peak: n×n Float64 GRM + one n-element scratch vector
    grm = zeros(Float64, n_individuals, n_individuals)
    z = Vector{Float64}(undef, n_individuals)

    for j in poly_indices
        m = col_means[j]
        sd = col_stds[j]
        inv_sd = 1.0 / sd
        @inbounds for i in 1:n_individuals
            z[i] = (genotypes[i, j] - m) * inv_sd
        end
        # Rank-1 update: grm += z * z'
        # Use BLAS ger! for efficiency: A := alpha * x * y' + A
        BLAS.ger!(1.0, z, z, grm)
    end
    # Normalise
    rmul!(grm, 1.0 / n_poly)

    # Eigendecomposition (symmetric → real eigenvalues)
    evals_raw = eigvals!(Symmetric(grm))  # in-place, destroys grm
    # Sort descending
    sort!(evals_raw, rev=true)
    # Clamp tiny negatives from numerical noise
    @inbounds for i in eachindex(evals_raw)
        evals_raw[i] = max(evals_raw[i], 0.0)
    end
    evals = evals_raw

    total_var = sum(evals)
    if total_var < 1e-15
        return nothing
    end

    # Normalised spectrum (sums to 1)
    normed = evals ./ total_var

    # Variance of eigenvalues
    eigenvalue_variance = var(evals)

    # Shannon entropy of normalised eigenvalues
    eigenvalue_entropy = 0.0
    for p in normed
        if p > 0
            eigenvalue_entropy -= p * log(p)
        end
    end

    # Participation ratio = (Σλ)² / Σλ² — effective number of components
    sum_sq = sum(x -> x*x, evals)
    n_eff = sum_sq > 0 ? total_var^2 / sum_sq : 1.0

    # Variance explained
    variance_explained = normed
    top1 = normed[1]
    top10 = sum(normed[1:min(10, length(normed))])

    # Marchenko–Pastur edge for the bulk: (1 + √(n/p))²
    gamma = n_individuals / n_poly
    mp_upper = (1.0 + sqrt(gamma))^2
    mp_ratio = evals[1] / mp_upper

    # Tracy-Widom-like centering/scaling
    mu_n = (sqrt(n_individuals - 1) + sqrt(n_poly))^2 / n_poly
    sigma_n = (sqrt(n_individuals - 1) + sqrt(n_poly)) / n_poly *
              (1.0 / sqrt(n_individuals - 1) + 1.0 / sqrt(n_poly))^(1/3)
    tw_stat = sigma_n > 0 ? (evals[1] - mu_n) / sigma_n : 0.0

    return SpectralStats(
        evals,
        variance_explained,
        n_eff,
        eigenvalue_variance,
        eigenvalue_entropy,
        tw_stat,
        top1,
        top10,
        mp_ratio
    )
end

"""
    print_stats(stats::PopulationStats)

Print population genetics statistics in a formatted way.
"""
function print_stats(stats::PopulationStats)
    println("\n" * "="^50)
    println("Population Genetics Statistics")
    println("="^50)
    println("Total variant sites: $(stats.n_sites)")
    println("Common variants (MAF > 0.05): $(stats.n_common_variants)")
    println("Mean allele frequency: $(round(stats.mean_allele_freq, digits=4))")
    println("Nucleotide diversity (π): $(round(stats.nucleotide_diversity, digits=6))")
    println("Watterson's theta (θw): $(round(stats.theta_w, digits=6))")
    println("Variant density: $(round(stats.variant_density * 1000, digits=2)) per kb")

    if stats.spectral !== nothing
        s = stats.spectral
        println("\nSpectral Statistics (GRM eigendecomposition)")
        println("-"^50)
        println("Eigenvalue variance: $(round(s.eigenvalue_variance, digits=6))")
        println("Eigenvalue entropy: $(round(s.eigenvalue_entropy, digits=4))")
        println("Effective dimensionality (participation ratio): $(round(s.n_eff_components, digits=2))")
        println("Tracy-Widom statistic (top eigenvalue): $(round(s.tracy_widom_stat, digits=4))")
        println("Variance explained by PC1: $(round(s.top1_variance * 100, digits=2))%")
        top10_pct = round(s.top10_variance * 100, digits=2)
        println("Cumulative variance by top 10 PCs: $(top10_pct)%")
        println("Top eigenvalue / Marchenko-Pastur edge: $(round(s.marchenko_pastur_ratio, digits=4))")
        n_show = min(5, length(s.eigenvalues))
        println("Top $n_show eigenvalues: $(round.(s.eigenvalues[1:n_show], digits=4))")
    end

    println("="^50)
end

"""
    allele_frequency_spectrum(genotypes::Matrix{Int}) -> Vector{Int}

Calculate the site frequency spectrum (SFS) from genotype data.

# Arguments
- `genotypes::Matrix{Int}`: Diploid genotype matrix

# Returns
- `Vector{Int}`: Site frequency spectrum (count of sites with i derived alleles)
"""
function allele_frequency_spectrum(genotypes::Matrix{Int})
    n_individuals, n_sites = size(genotypes)
    n_haplotypes = 2 * n_individuals
    
    sfs = zeros(Int, n_haplotypes + 1)
    
    for j in 1:n_sites
        # Count derived alleles at this site
        derived_count = sum(genotypes[:, j])
        sfs[derived_count + 1] += 1
    end
    
    return sfs
end

"""
    tajimas_d(genotypes::Matrix{Int}) -> Float64

Calculate Tajima's D statistic.

# Arguments
- `genotypes::Matrix{Int}`: Diploid genotype matrix

# Returns
- `Float64`: Tajima's D value
"""
function tajimas_d(genotypes::Matrix{Int})
    n_individuals, n_sites = size(genotypes)
    n_haplotypes = 2 * n_individuals
    
    if n_sites == 0
        return 0.0
    end
    
    # Calculate π (nucleotide diversity)
    pi = 0.0
    for j in 1:n_sites
        allele_count = sum(genotypes[:, j])
        p = allele_count / n_haplotypes
        pi += 2 * p * (1 - p)
    end
    pi /= n_sites
    
    # Calculate Watterson's theta
    harmonic_number = sum(1/i for i in 1:(n_haplotypes-1))
    theta_w = n_sites / harmonic_number / n_sites  # Per site
    
    # Calculate variance components for Tajima's D
    a1 = harmonic_number
    a2 = sum(1/i^2 for i in 1:(n_haplotypes-1))
    
    b1 = (n_haplotypes + 1) / (3 * (n_haplotypes - 1))
    b2 = 2 * (n_haplotypes^2 + n_haplotypes + 3) / (9 * n_haplotypes * (n_haplotypes - 1))
    
    c1 = b1 - 1/a1
    c2 = b2 - (n_haplotypes + 2)/(a1 * n_haplotypes) + a2/a1^2
    
    e1 = c1 / a1
    e2 = c2 / (a1^2 + a2)
    
    # Variance of Tajima's D
    var_d = e1 * n_sites + e2 * n_sites * (n_sites - 1)
    
    if var_d <= 0
        return 0.0
    end
    
    # Tajima's D
    d = (pi - theta_w) / sqrt(var_d)
    
    return d
end