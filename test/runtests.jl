using SeqFold
using Test
using LinearAlgebra
using Aqua

function _randna(n)
    String(rand("AGTC", n))
end

function _randconds()
    MeltingConditions(
        100*rand(), 100*rand(), 100*rand(), 5*rand(), 5*rand(),5*rand(), 5*rand()
    )
end

function _gc(seq::String)
    seq = uppercase(seq)
    if !all(in("AGTC"), seq)
        throw(ArgumentError("Invalid DNA sequence"))
    end
    count(in("GC"), seq) / length(seq)
end

function verify_tm_cache(seq)
    seq2 = SeqFold.complement(seq)
    n = length(seq)
    
    for conds in [_randconds() for _ in 1:10]
        cache = tm_cache(seq, seq2, conditions=conds)
        for _ in 1:10
            i = rand(1:n-2)
            j = min(i+rand(1:n÷2), n)
            direct_tm = tm(SubString(seq, i:j), SubString(seq2, i:j), conditions=conds)
            cache_tm = cache[i, j]
            @test isapprox(direct_tm, cache_tm, atol=0.1)
        end
    end
end

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(SeqFold)
end

@testset "MeltingConditions" begin
    @testset "Throws" begin
        @test_throws DomainError MeltingConditions(-1, 0, 0, 0, 0, 0, 0)
        @test_throws ArgumentError MeltingConditions(0, 0, 0, 0, 0, 0, 0)
    end
    @testset "Defaults" begin
        @test SeqFold.PCR_CONDITIONS == MeltingConditions(250.0, 0.0, 0.0, 50.0, 2.0, 1.5, 0.2) == MeltingConditions(:pcr)
        @test SeqFold.STD_CONDITIONS == MeltingConditions(25.0, 25.0, 50.0, 0.0, 0.0, 0.0, 0.0) == MeltingConditions(:std)
    end
end

@testset "Complement" begin

    seq1 = "AAGGTTCC"
    seq1c = "TTCCAAGG"
    @test SeqFold.complement(seq1) == seq1c

    for _ in 1:100
        n = rand(5:55)
        seq = _randna(n)
        @test length(SeqFold.complement(seq)) == n
        @test seq == SeqFold.complement(SeqFold.complement(seq))
    end

    bad_nucl_compl = SeqFold.complement("some bad nucleotides")
    @test issubset(Set(bad_nucl_compl), Set("AGTC."))
end

@testset "tm" begin
    @testset "Throws" begin
        @test_throws MethodError tm()
        @test_throws ArgumentError tm("BAD NUCLEOTIDES")
        @test_throws ArgumentError tm("A")
        @test_throws DimensionMismatch tm("ACT", "ACGGT")
    end

    @testset "Reference TMs" begin
        common_conds = MeltingConditions(1000.0, 1000.0, 0.0, 0.0, 2.0, 0.0, 0.0)
        abs_tolerance = 2.0
        Mgs = [0.5, 1.5, 3.0]
        tm_test_cases = Dict(
            # Values from the reference paper (Owczarzy et al., 2008)
            "TTGTAGTCAT" => [27.2, 30.9, 32.3],
            "ATCGTCTGGA" => [35.1, 38.9, 40.6],
            "GATGCGCTCG" => [44.0, 47.7, 48.8],
            "TTCTACCTATGTGAT" => [43.3, 46.1, 47.8],
            "TGATTCTACCTATGTGATTT" => [51.8, 54.6, 55.8],
            "GTTCTATACTCTTGAAGTTGATTAC" => [57.2, 59.7, 60.8],
            "CTTAAGATATGAGAACTTCAACTAATGTGT" => [61.0, 63.1, 64.1],

            # Values, calculated for several random sequences, using BioPython:
            # Bio.SeqUtils.MeltingTemp.Tm_NN( sequence, saltcorr=7, Na=0, K=0, Tris=2, Mg=Mg, dNTPs=0, dnac1=1000, dnac2=1000)
            "GAATGATTGC" => [29.6, 33.6, 35.7],
            "GATCAGGTCA" => [33.0, 37.0, 39.0],
            "AGCTGGTCGC" => [45.6, 49.6, 51.5],
            "GTCAAATAGG" => [26.7, 30.7, 32.7],
            "GGAGATCTTT" => [28.2, 32.2, 34.2],
            "CTGCCAACTG" => [38.2, 42.1, 44.1],
            "CAAATGAGACAGCCG" => [53.2, 56.1, 57.5],
            "CCGAACTCACAGTAT" => [50.1, 53.0, 54.5],
            "GCGGTTCTGGAGCTA" => [57.1, 59.9, 61.3],
            "TGGTGGGTACAGGGA" => [58.0, 60.8, 62.3],
            "GGACGGCGCTCACAT" => [61.6, 64.4, 65.8],
            "CGGCTTAGGAAGGCG" => [58.9, 61.7, 63.1],
            "TGTGTCTTCGATAAGATCAC" => [55.9, 58.4, 59.6],
            "TTTGAGGATCCACCAGTCGG" => [63.6, 65.9, 67.1],
            "GTTGGCCACCCACGGCGCGG" => [76.3, 78.2, 79.2],
            "TTGTTCACAGTGCGTTACAA" => [59.4, 61.9, 63.2],
            "ATGCCTCTTCTCGCCGTTGA" => [65.9, 68.2, 69.4],
            "GTGAATTGTTTGGAGTGGAT" => [57.3, 59.8, 61.1],
            "ATAGGCTCAGATAAGTCCGGGCTAT" => [66.0, 68.0, 69.0],
            "TTCAGAAGCATTATACGCCAAGCCA" => [65.9, 68.0, 69.1],
            "GGCCCTCAGATGGAGTCAGCATTGC" => [71.1, 72.9, 73.8],
            "TGGACTAACGTCACGTGGTTTCTGG" => [67.9, 69.8, 70.8],
            "ATCGATTGTACGGAATATCTGGACG" => [63.3, 65.4, 66.4],
            "GCCCTTCGTGGTAACCCCCCAATCT" => [72.1, 74.0, 74.9],
            "CTGATAGCTAAGACGTCTACCTAAGTCCGA" => [67.0, 68.8, 69.7],
            "GACGTCTACTCAAATGATCGAATGCTCGTT" => [67.3, 69.1, 70.0],
            "GTCTGCTAACTCGGTGTACGCTTCGTTAAA" => [68.8, 70.6, 71.5],
            "GGTGCTTCTTCTTGCTCTCCGAACACAATC" => [70.0, 71.7, 72.5],
            "AATGGAGTAAAGCTCTATAGCTTCTAAGCT" => [64.2, 66.0, 67.0],
            "GCTTCCACCTAGGTGAGCCCCTGTAAGTAA" => [71.6, 73.3, 74.1],
            "GGGAGATGGAGGGTCGGAAAGTCCGAATTACTCGGAGGAA" => [76.9, 78.2, 78.9],
            "TCTCCACAATGCTTAGCTCGATAGCGCAAGAGATTCGCATGATTACGACG" => [76.9, 78.1, 78.7],
            "TGTCTTTGCTAATGGCCCATTCAATAACCAGGCGCACCCAATATTCGATTCCTAGCAACA" => [78.7, 79.8, 80.4],
        )

        for (seq, expected_tms) in tm_test_cases
            seq_complement = SeqFold.complement(seq)
            for (Mg, expected_tm) in zip(Mgs, expected_tms)
                calc_tm = tm(seq; conditions=common_conds, Mg=Mg)
                calc_tm_comp = tm(seq; conditions=common_conds, Mg=Mg)
                @test calc_tm == tm(lowercase(seq); conditions=common_conds, Mg=Mg) == tm(uppercasefirst(seq); conditions=common_conds, Mg=Mg)
                @test calc_tm == calc_tm_comp
                @test isapprox(calc_tm, expected_tm; atol=abs_tolerance)
            end
        end
    end
end

@testset "gc_cache" begin
    @testset "Basic functionality" begin
        for _ in 1:100
            n = rand(5:55)
            seq = _randna(n)
            M = gc_cache(seq)
            @test size(M) == (n, n)
            @test count(isinf, M) == (n^2 - n) ÷ 2

            isgc = Float64[c in "GC" for c in seq]
            @test isgc == diag(M)

            for _ in 1:n
                i = rand(1:n-1)
                j = min(n, i + rand(1:n÷2))

                @test M[i,j] == _gc(seq[i:j])
            end
        end
    end
end

@testset "tm_cache" begin
    @testset "Throws" begin
        @test_throws MethodError tm_cache()
        @test_throws ArgumentError tm_cache("BAD NUCLEOTIDES")
        @test_throws ArgumentError tm_cache("A")
        @test_throws DimensionMismatch tm_cache("ACT", "ACGGT")
    end

    @testset "Verification" begin
        for seq in [_randna(rand(5:50)) for _ in 1:10]
            verify_tm_cache(seq)
        end
    end
end

# @testset "fold tests" begin
    
# end
