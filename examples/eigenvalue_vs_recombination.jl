#!/usr/bin/env julia
"""
Eigenvalue spectrum vs. recombination rate analysis
====================================================

Tests whether lower recombination rates increase LD, producing larger
eigenvalue variance in the genetic relatedness matrix (GRM).

Uses the GenotypeSimulator ARG (Ancestral Recombination Graph) to
simulate genotypes with different recombination rates.  For ρ=0 we
fall back to the standard coalescent (no recombination).
"""

using Pkg
Pkg.activate(".")

using GenotypeSimulator
using Random
using Statistics
using Printf
using LinearAlgebra

# ── Parameters ─────────────────────────────────────────────────────────────
const NE            = 10_000
const MUTATION_RATE = 1.25e-8
const SEQ_LENGTH    = 100_000      # 100 kb — keeps ARG tractable at high ρ
const SAMPLE_SIZE   = 80           # 80 diploid individuals (160 haplotypes)
const N_REPLICATES  = 8            # replicates per condition

# Recombination rates to sweep (per bp per generation)
const RECOMB_RATES = [0.0, 1e-10, 5e-10, 1e-9, 5e-9, 1e-8, 5e-8]

# ── Run one replicate ─────────────────────────────────────────────────────
function run_one(ρ::Float64; seed::Int)
    rng = MersenneTwister(seed)
    params = PopulationParams(NE, MUTATION_RATE, ρ, SEQ_LENGTH, SAMPLE_SIZE)

    # Use the module's ARG simulation for ρ>0, standard coalescent for ρ=0
    haplotypes, positions = redirect_stdout(devnull) do
        if ρ > 0
            simulate_with_recombination(params; rng=rng, force_standard=true)
        else
            simulate_genotypes(params; rng=rng, force_standard=true)
        end
    end

    diploid = haplotypes_to_diploid(haplotypes)
    n_sites = size(diploid, 2)

    if n_sites < 3
        return nothing
    end

    stats = calculate_stats(diploid, positions)
    return stats
end

# ── Main sweep ─────────────────────────────────────────────────────────────
function main()
    println("=" ^ 80)
    println("  Eigenvalue spectrum vs. recombination rate")
    println("  (using GenotypeSimulator ARG simulation)")
    println("=" ^ 80)
    println()
    @printf("  Ne = %d  |  μ = %.2e  |  L = %d bp  |  n = %d  |  reps = %d\n",
            NE, MUTATION_RATE, SEQ_LENGTH, SAMPLE_SIZE, N_REPLICATES)
    println()

    # Storage
    ResultTuple = NamedTuple{
        (:eigenvalue_variance, :top1_var, :top10_var, :n_eff, :entropy, :mp_ratio, :n_sites),
        Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int}}
    results = Dict{Float64, Vector{ResultTuple}}()

    for ρ in RECOMB_RATES
        results[ρ] = ResultTuple[]
        @printf("ρ = %.1e :  ", ρ)
        flush(stdout)
        for rep in 1:N_REPLICATES
            seed = Int(hash((ρ, rep)) % (2^30))
            stats = run_one(ρ; seed=seed)
            if stats !== nothing && stats.spectral !== nothing
                s = stats.spectral
                push!(results[ρ], (
                    eigenvalue_variance = s.eigenvalue_variance,
                    top1_var            = s.top1_variance,
                    top10_var           = s.top10_variance,
                    n_eff               = s.n_eff_components,
                    entropy             = s.eigenvalue_entropy,
                    mp_ratio            = s.marchenko_pastur_ratio,
                    n_sites             = stats.n_sites,
                ))
                print(".")
            else
                print("x")
            end
            flush(stdout)
        end
        println()
    end

    # ── Summary table ──────────────────────────────────────────────────────
    println()
    println("─" ^ 110)
    @printf("%-12s  %7s  %16s  %10s  %10s  %10s  %10s  %10s\n",
            "ρ", "#sites", "EigVar(±SE)", "Top1%", "Top10%", "N_eff", "Entropy", "MP ratio")
    println("─" ^ 110)

    for ρ in RECOMB_RATES
        recs = results[ρ]
        if isempty(recs)
            @printf("%-12s  %7s\n", @sprintf("%.1e", ρ), "N/A")
            continue
        end
        mn_sites = round(Int, mean(r.n_sites for r in recs))
        ev    = mean(r.eigenvalue_variance for r in recs)
        ev_se = length(recs) > 1 ? std([r.eigenvalue_variance for r in recs]) / sqrt(length(recs)) : 0.0
        t1  = mean(r.top1_var  for r in recs) * 100
        t10 = mean(r.top10_var for r in recs) * 100
        ne  = mean(r.n_eff     for r in recs)
        ent = mean(r.entropy   for r in recs)
        mp  = mean(r.mp_ratio  for r in recs)
        @printf("%-12s  %7d  %8.4f ± %.4f  %9.2f%%  %9.2f%%  %10.2f  %10.4f  %10.4f\n",
                @sprintf("%.1e", ρ), mn_sites, ev, ev_se, t1, t10, ne, ent, mp)
    end
    println("─" ^ 110)

    # ── Trend analysis ─────────────────────────────────────────────────────
    println()
    println("Trend analysis (Spearman rank-correlation with ρ):")

    ρ_vals  = Float64[]
    ev_vals = Float64[]
    t1_vals = Float64[]
    ne_vals = Float64[]
    ent_vals = Float64[]

    for ρ in RECOMB_RATES
        recs = results[ρ]
        if !isempty(recs)
            push!(ρ_vals, ρ)
            push!(ev_vals,  mean(r.eigenvalue_variance for r in recs))
            push!(t1_vals,  mean(r.top1_var            for r in recs))
            push!(ne_vals,  mean(r.n_eff               for r in recs))
            push!(ent_vals, mean(r.entropy              for r in recs))
        end
    end

    if length(ρ_vals) >= 3
        function rank_corr(x, y)
            n = length(x)
            rx = sortperm(sortperm(x)) .|> Float64
            ry = sortperm(sortperm(y)) .|> Float64
            mx, my = mean(rx), mean(ry)
            num = sum((rx .- mx) .* (ry .- my))
            den = sqrt(sum((rx .- mx).^2) * sum((ry .- my).^2))
            return den > 0 ? num / den : 0.0
        end

        @printf("  EigVar  vs ρ:  r_s = %+.3f  (expect negative  — low ρ → high EigVar)\n",
                rank_corr(ρ_vals, ev_vals))
        @printf("  Top1%%   vs ρ:  r_s = %+.3f  (expect negative  — low ρ → dominant PC1)\n",
                rank_corr(ρ_vals, t1_vals))
        @printf("  N_eff   vs ρ:  r_s = %+.3f  (expect positive  — high ρ → flatter spectrum)\n",
                rank_corr(ρ_vals, ne_vals))
        @printf("  Entropy vs ρ:  r_s = %+.3f  (expect positive  — high ρ → uniform eigenvalues)\n",
                rank_corr(ρ_vals, ent_vals))
    end

    # ── Interpretation ─────────────────────────────────────────────────────
    println()
    println("=" ^ 80)
    println("Interpretation")
    println("=" ^ 80)
    println()
    println("  Lower ρ → fewer independent genealogical segments → stronger LD:")
    println("    • GRM eigenvalue variance is HIGHER (concentrated spectrum)")
    println("    • PC1 captures more variance")
    println("    • Effective dimensionality (N_eff) is LOWER")
    println("    • Eigenvalue entropy is LOWER")
    println()
    println("  Higher ρ → many independent segments → genealogies decorrelate:")
    println("    • GRM approaches Wishart null → Marchenko–Pastur spectrum")
    println("    • Eigenvalue variance DECREASES toward the MP prediction")
    println("    • N_eff increases toward n (uniform spectrum)")
    println()
    println("  The Marchenko–Pastur (MP) ratio reports λ_max relative to the")
    println("  theoretical upper edge of the null spectrum. Values ≫ 1 indicate")
    println("  significant covariance structure (LD) beyond sampling noise.")

    # ── Per-replicate detail ───────────────────────────────────────────────
    println()
    println("Per-replicate eigenvalue variances:")
    for ρ in RECOMB_RATES
        recs = results[ρ]
        vals = [r.eigenvalue_variance for r in recs]
        if !isempty(vals)
            @printf("  ρ=%.1e : mean=%.4f  std=%.4f  [%s]\n",
                    ρ, mean(vals), length(vals) > 1 ? std(vals) : 0.0,
                    join([@sprintf("%.4f", v) for v in vals], ", "))
        end
    end
end

main()
