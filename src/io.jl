"""
    save_genotypes(genotypes::Matrix{Int}, positions::Vector{Int}, filename::String; format::String="csv")

Save genotype data to file in various formats.

# Arguments
- `genotypes::Matrix{Int}`: Genotype matrix (individuals × sites)
- `positions::Vector{Int}`: Positions of variant sites
- `filename::String`: Output filename
- `format::String`: Output format ("csv", "vcf", or "plink")
"""
function save_genotypes(genotypes::Matrix{<:Integer}, positions::Vector{<:Integer}, filename::String; format::String="csv")
    if format == "csv"
        save_csv(genotypes, positions, filename)
    elseif format == "vcf"
        save_vcf(genotypes, positions, filename)
    elseif format == "plink"
        save_plink(genotypes, positions, filename)
    else
        error("Unsupported format: $format. Use 'csv', 'vcf', or 'plink'")
    end
end

"""
    save_csv(genotypes::Matrix{<:Integer}, positions::Vector{<:Integer}, filename::String)

Save genotypes in wide-format CSV, streaming directly to disk (no DataFrame).
Columns: Position_1, Position_2, ..., one row per individual.
"""
function save_csv(genotypes::Matrix{<:Integer}, positions::Vector{<:Integer}, filename::String)
    n_individuals, n_sites = size(genotypes)

    open(filename, "w") do io
        # Header row: position labels
        for j in 1:n_sites
            j > 1 && print(io, ',')
            print(io, positions[j])
        end
        println(io)

        # One row per individual
        for i in 1:n_individuals
            for j in 1:n_sites
                j > 1 && print(io, ',')
                print(io, Int(genotypes[i, j]))
            end
            println(io)
        end
    end

    println("Genotypes saved to $filename (CSV format)")
end

"""
    save_vcf(genotypes::Matrix{Int}, positions::Vector{Int}, filename::String)

Save genotypes in VCF format.
"""
function save_vcf(genotypes::Matrix{<:Integer}, positions::Vector{<:Integer}, filename::String)
    n_individuals, n_sites = size(genotypes)
    
    open(filename, "w") do io
        # VCF header
        println(io, "##fileformat=VCFv4.2")
        println(io, "##source=GenotypeSimulator.jl")
        println(io, "##contig=<ID=1,length=1000000>")
        println(io, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        
        # Column headers
        print(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for i in 1:n_individuals
            print(io, "\tIND$i")
        end
        println(io)
        
        # Variant records
        for j in 1:n_sites
            print(io, "1\t$(positions[j])\t.\tA\tT\t.\tPASS\t.\tGT")
            
            for i in 1:n_individuals
                gt = genotypes[i, j]
                if gt == 0
                    print(io, "\t0/0")
                elseif gt == 1
                    print(io, "\t0/1")
                else
                    print(io, "\t1/1")
                end
            end
            println(io)
        end
    end
    
    println("Genotypes saved to $filename (VCF format)")
end

"""
    save_plink(genotypes::Matrix{Int}, positions::Vector{Int}, filename::String)

Save genotypes in PLINK format (.ped and .map files).
"""
function save_plink(genotypes::Matrix{<:Integer}, positions::Vector{<:Integer}, filename::String)
    n_individuals, n_sites = size(genotypes)
    
    # Remove extension if provided
    base_name = replace(filename, r"\.[^.]*$" => "")
    
    # Save .map file (marker information)
    open("$(base_name).map", "w") do io
        for j in 1:n_sites
            println(io, "1\tSNP$j\t0\t$(positions[j])")
        end
    end
    
    # Save .ped file (pedigree and genotype data)
    open("$(base_name).ped", "w") do io
        for i in 1:n_individuals
            print(io, "FAM$i\tIND$i\t0\t0\t0\t-9")  # Family, Individual, Father, Mother, Sex, Phenotype
            
            for j in 1:n_sites
                gt = genotypes[i, j]
                if gt == 0
                    print(io, "\tA A")
                elseif gt == 1
                    print(io, "\tA T")
                else
                    print(io, "\tT T")
                end
            end
            println(io)
        end
    end
    
    println("Genotypes saved to $(base_name).ped and $(base_name).map (PLINK format)")
end

"""
    load_genotypes(filename::String; format::String="csv") -> (Matrix{Int}, Vector{Int})

Load genotype data from file.

# Arguments
- `filename::String`: Input filename
- `format::String`: Input format ("csv" only for now)

# Returns
- `Matrix{Int}`: Genotype matrix
- `Vector{Int}`: Positions vector
"""
function load_genotypes(filename::String; format::String="csv")
    if format == "csv"
        return load_csv(filename)
    else
        error("Unsupported format: $format. Only 'csv' is currently supported for loading.")
    end
end

"""
    load_csv(filename::String) -> (Matrix{Int}, Vector{Int})

Load genotypes from CSV format.
"""
function load_csv(filename::String)
    df = CSV.read(filename, DataFrame)
    
    # Get unique individuals and positions
    individuals = sort(unique(df.individual))
    positions = sort(unique(df.position))
    
    n_individuals = length(individuals)
    n_sites = length(positions)
    
    # Create genotype matrix
    genotypes = zeros(Int, n_individuals, n_sites)
    
    for row in eachrow(df)
        i_idx = findfirst(==(row.individual), individuals)
        j_idx = findfirst(==(row.position), positions)
        genotypes[i_idx, j_idx] = row.genotype
    end
    
    return genotypes, positions
end