using Test
using GenotypeSimulator
using Random

@testset "GenotypeSimulator.jl Tests" begin

    @testset "Types and Parameters" begin
        # Test PopulationParams constructor
        params = PopulationParams(1000, 1e-8, 1e-8, 10000, 10)
        @test params.ne == 1000
        @test params.mutation_rate == 1e-8
        @test params.sequence_length == 10000
        @test params.sample_size == 10

        # Test parameter validation
        @test_throws AssertionError PopulationParams(-1000, 1e-8, 1e-8, 10000, 10)  # negative Ne
        @test_throws AssertionError PopulationParams(1000, -1e-8, 1e-8, 10000, 10)  # negative mutation rate
        @test_throws AssertionError PopulationParams(1000, 1e-8, 1e-8, -10000, 10)  # negative sequence length
        @test_throws AssertionError PopulationParams(1000, 1e-8, 1e-8, 10000, -10)  # negative sample size

        # Test CoalescentNode
        node = CoalescentNode(1)
        @test node.id == 1
        @test node.parent === nothing
        @test isempty(node.children)
        @test node.time == 0.0
        @test isempty(node.mutations)
    end

    @testset "Coalescent Simulation" begin
        Random.seed!(42)

        # Test coalescent times
        times = simulate_coalescent_times(4, 1000)
        @test length(times) == 3  # n-1 coalescent events
        @test all(times[i] < times[i+1] for i in 1:length(times)-1)  # Times should be increasing

        # Test tree building
        root = build_coalescent_tree(2, 1000)  # 2 individuals = 4 haplotypes
        @test count_nodes(root) == 7  # 4 leaves + 3 internal nodes
        @test tree_height(root) > 0
    end

    @testset "Mutations" begin
        Random.seed!(42)

        # Build a simple tree
        root = build_coalescent_tree(2, 1000)
        params = PopulationParams(1000, 1e-6, 1e-8, 1000, 2)  # High mutation rate for testing

        # Add mutations
        add_mutations!(root, params)

        # Check that mutations were added
        total_mutations = count_mutations(root)
        @test total_mutations >= 0

        # Get mutation positions
        positions = get_all_mutation_positions(root)
        @test length(positions) <= total_mutations
        @test all(1 <= pos <= params.sequence_length for pos in positions)
    end

    @testset "Genotype Extraction" begin
        Random.seed!(42)

        # Simple test with known tree structure
        params = PopulationParams(1000, 1e-6, 1e-8, 1000, 2)
        root = build_coalescent_tree(params.sample_size, params.ne)
        add_mutations!(root, params)

        genotypes_hap, positions = extract_genotypes(root, params.sample_size, params.sequence_length)

        @test size(genotypes_hap, 1) == 2 * params.sample_size  # Correct number of haplotypes
        @test size(genotypes_hap, 2) == length(positions)  # Correct number of sites
        @test all(g in [0, 1] for g in genotypes_hap)  # Binary genotypes

        # Test diploid conversion
        genotypes_diploid = haplotypes_to_diploid(genotypes_hap)
        @test size(genotypes_diploid, 1) == params.sample_size  # Correct number of individuals
        @test size(genotypes_diploid, 2) == length(positions)  # Same number of sites
        @test all(g in [0, 1, 2] for g in genotypes_diploid)  # Diploid genotypes
    end

    @testset "Statistics" begin
        # Create test genotype data
        genotypes = [0 1 2; 1 0 1; 2 2 0]  # 3 individuals, 3 sites
        positions = [100, 200, 300]

        stats = calculate_stats(genotypes, positions)

        @test stats.n_sites == 3
        @test stats.mean_allele_freq ≈ 0.5  # Average allele frequency
        @test stats.nucleotide_diversity > 0
        @test stats.variant_density > 0

        # Test site frequency spectrum
        sfs = allele_frequency_spectrum(genotypes)
        @test length(sfs) == 7  # 2 * n_individuals + 1
        @test sum(sfs) == 3  # Total number of sites

        # Test Tajima's D (should not crash)
        d = tajimas_d(genotypes)
        @test isa(d, Float64)
    end

    @testset "Full Simulation" begin
        Random.seed!(42)

        # Test with small parameters for speed
        params = PopulationParams(100, 1e-7, 1e-8, 1000, 5)
        genotypes_hap, positions = simulate_genotypes(params)

        @test size(genotypes_hap, 1) == 2 * params.sample_size
        @test length(positions) == size(genotypes_hap, 2)
        @test all(g in [0, 1] for g in genotypes_hap)

        # Convert to diploid and test
        genotypes_diploid = haplotypes_to_diploid(genotypes_hap)
        @test size(genotypes_diploid, 1) == params.sample_size
        @test all(g in [0, 1, 2] for g in genotypes_diploid)

        # Calculate statistics (should not crash)
        stats = calculate_stats(genotypes_diploid, positions)
        @test stats.n_sites >= 0
    end

    @testset "Recombination" begin
        Random.seed!(42)

        # Test recombination map
        recomb_map = uniform_recombination_map(1000, 1e-8)
        @test length(recomb_map.positions) == 2
        @test recomb_map.positions[1] == 1
        @test recomb_map.positions[2] == 1001
        @test length(recomb_map.rates) == 1
        @test recomb_map.rates[1] == 1e-8

        # Test recombination position sampling
        pos = sample_recombination_position(recomb_map)
        @test 1 <= pos <= 1000

        # Test ARG simulation with small parameters
        params = PopulationParams(50, 1e-7, 1e-7, 500, 3)
        genotypes_hap, positions = simulate_with_recombination(params)

        @test size(genotypes_hap, 1) == 2 * params.sample_size
        @test length(positions) == size(genotypes_hap, 2)
        @test all(g in [0, 1] for g in genotypes_hap)

        # Convert to diploid
        genotypes_diploid = haplotypes_to_diploid(genotypes_hap)
        @test size(genotypes_diploid, 1) == params.sample_size
        @test all(g in [0, 1, 2] for g in genotypes_diploid)

        # Calculate statistics
        stats = calculate_stats(genotypes_diploid, positions)
        @test stats.n_sites >= 0
    end

    @testset "ARG Summary Statistics" begin
        # Validate that segregating sites roughly match theory and are stable across rho.
        ne = 1000
        mu = 2e-6
        seq_len = 1000
        n_samples = 5
        n_haps = 2 * n_samples
        harmonic = sum(1.0 / i for i in 1:(n_haps - 1))
        expected_s = 4 * ne * mu * seq_len * harmonic

        function avg_segregating_sites(rho::Float64)
            params = PopulationParams(ne, mu, rho, seq_len, n_samples)
            n_reps = 8
            total_sites = 0
            for rep in 1:n_reps
                rng = MersenneTwister(1000 + rep + hash(rho))
                _, positions = simulate_with_recombination(params; rng=rng)
                total_sites += length(positions)
            end
            return total_sites / n_reps
        end

        avg_rho0 = avg_segregating_sites(0.0)
        avg_rho1 = avg_segregating_sites(1e-7)

        @test avg_rho0 > 0.4 * expected_s
        @test avg_rho0 < 1.6 * expected_s
        @test avg_rho1 > 0.4 * expected_s
        @test avg_rho1 < 1.6 * expected_s
        @test abs(avg_rho0 - avg_rho1) < 0.6 * expected_s
    end

    @testset "Selection" begin
        Random.seed!(42)

        # Test selection models
        neutral = NeutralSelection()
        @test isa(neutral, NeutralSelection)

        directional = DirectionalSelection(0.01, 0.5, 500, 1000.0)
        @test directional.selection_coefficient == 0.01
        @test directional.dominance == 0.5
        @test directional.selected_position == 500

        background = BackgroundSelection(1e-8, -0.01, 0.0)
        @test background.deleterious_rate == 1e-8
        @test background.selection_coefficient == -0.01

        balancing = BalancingSelection([0.01, 0.01], 500, [0.5, 0.5])
        @test length(balancing.selection_coefficients) == 2
        @test balancing.selected_position == 500

        # Test selection parameters
        selection_params = SelectionParameters(neutral, 1000)
        @test selection_params.effective_population_size == 1000

        # Test allele frequency trajectory
        traj = allele_frequency_trajectory(directional, 1000, [0.0, 500.0, 1500.0])
        @test length(traj) == 3
        @test all(0 <= f <= 1 for f in traj)

        # Test effective population size calculation
        ne_eff = effective_population_size_with_selection(neutral, 1000, 500)
        @test ne_eff == 1000.0

        ne_eff_dir = effective_population_size_with_selection(directional, 1000, 500)
        @test ne_eff_dir <= 1000.0  # Should be reduced near selected site

        # Test selection simulation with small parameters
        params = PopulationParams(50, 1e-7, 1e-8, 500, 3)
        neutral_params = SelectionParameters(neutral, params.ne)

        genotypes_hap, positions = simulate_with_selection(params, neutral_params)

        @test size(genotypes_hap, 1) == 2 * params.sample_size
        @test length(positions) == size(genotypes_hap, 2)
        @test all(g in [0, 1] for g in genotypes_hap)

        # Convert to diploid
        genotypes_diploid = haplotypes_to_diploid(genotypes_hap)
        @test size(genotypes_diploid, 1) == params.sample_size
        @test all(g in [0, 1, 2] for g in genotypes_diploid)

        # Calculate statistics
        stats = calculate_stats(genotypes_diploid, positions)
        @test stats.n_sites >= 0
    end

    @testset "Edge Cases" begin
        # Test with no mutations
        genotypes = zeros(Int, 10, 0)
        positions = Int[]
        stats = calculate_stats(genotypes, positions)

        @test stats.n_sites == 0
        @test stats.nucleotide_diversity == 0.0

        # Test with single individual
        params = PopulationParams(100, 1e-8, 1e-8, 1000, 1)
        genotypes_hap, positions = simulate_genotypes(params)
        @test size(genotypes_hap, 1) == 2
    end
end