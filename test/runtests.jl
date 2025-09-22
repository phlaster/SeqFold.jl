using SeqFold
using Test
using LinearAlgebra
using Aqua
using Logging
using JET

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
        cache = SeqFold.tm_cache(seq, seq2, conditions=conds)
        @test count(isnan, cache) == (n^2 - n) ÷ 2 + n
        for _ in 1:10
            i = rand(1:n-2)
            j = min(i+rand(1:n÷2), n)
            direct_tm = tm(seq[i:j], seq2[i:j], conditions=conds)
            cache_tm = cache[i, j]
            @test isapprox(direct_tm, cache_tm, atol=0.1)
        end
    end
end

function verify_dg_cache(seq)
    n = length(seq)
    for temp in rand(-20:1e-5:120, 10)
        cache = SeqFold.dg_cache(seq, temp=temp)
        @test count(isnan, cache) == (n^2 - n) ÷ 2 + n
        last_dg = dg(seq, temp=temp)
        last_cache = cache[1, end]
        @test isapprox(last_dg, last_cache, atol=0.1)
    end
end

@testset "Aqua.jl" begin
    Aqua.test_all(SeqFold)
end

@testset "JET.jl" begin
    JET.test_package(SeqFold; target_defined_modules = true)
end

@testset "tm.jl" begin
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
                M = SeqFold.gc_cache(seq)
                @test size(M) == (n, n)
                @test count(isnan, M) == (n^2 - n) ÷ 2

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
            @test_throws MethodError SeqFold.tm_cache()
            @test_throws ArgumentError SeqFold.tm_cache("BAD NUCLEOTIDES")
            @test_throws ArgumentError SeqFold.tm_cache("A")
            @test_throws DimensionMismatch SeqFold.tm_cache("ACT", "ACGGT")
        end

        @testset "Verification" begin
            for seq in [_randna(rand(5:50)) for _ in 1:10]
                verify_tm_cache(seq)
            end
        end
    end
end

@testset "fold.jl" begin
    @testset "dg" begin
        @test_throws ArgumentError dg("bad sequence")
        @test_throws ArgumentError dg("ATGCATGACGATUU")
        @test dg("ATGGATTTAGATAGAT") isa Float64

        @test_throws ArgumentError fold("A")
        @test_throws ArgumentError fold("TTTUUU")
        @test_throws ArgumentError fold("ATGCATGACGATUU")
        @test_throws ArgumentError dg("")
        @test_throws ArgumentError dg("ATCG", temp=-100.0)
        @test_throws ArgumentError dg("ATCG", temp=1000.0)
    end

    @testset "dg_cache" begin
        @testset "Throws" begin
            @test_throws MethodError SeqFold.dg_cache()
            @test_throws ArgumentError SeqFold.dg_cache("BAD NUCLEOTIDES")
            @test_throws ArgumentError SeqFold.dg_cache("A")
        end

        @testset "Verification" begin
            for seq in [_randna(rand(10:50)) for _ in 1:10]
                verify_dg_cache(seq)
            end
        end
    end


    @testset "fold throws" begin
        @test_throws ArgumentError fold("bad sequence")
        @test_throws ArgumentError fold("A")
        @test_throws ArgumentError fold("TTTUUU")
        @test_throws ArgumentError fold("ATGCATGACGATUU")
        @test_throws ArgumentError fold("")
        @test_throws ArgumentError fold("ATCG", temp=-100.0)
        @test_throws ArgumentError fold("ATCG", temp=1000.0)
    end

    @testset "_cache throws" begin
        @test_throws ArgumentError SeqFold._cache("bad sequence", 37)
        @test_throws ArgumentError SeqFold._cache("A", 37)
        @test_throws ArgumentError SeqFold._cache("TTTUUU", 37)
        @test_throws ArgumentError SeqFold._cache("ATGCATGACGATUU", 37)
        @test_throws ArgumentError SeqFold._cache("", 37)
        @test_throws ArgumentError SeqFold._cache("ATCG", -100.0)
        @test_throws ArgumentError SeqFold._cache("ATCG", 1000.0)
    end

    @testset "fold_dna" begin
        # unafold's estimates for free energy estimates of DNA oligos
        unafold_dgs = Dict(
            "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC" => -10.94,  # three branched structure
            "GGGAGGTCGCTCCAGCTGGGAGGAGCGTTGGGGGTATATACCCCCAACACCGGTACTGATCCGGTGACCTCCC" => -23.4,  # four branched structure
            "CGCAGGGAUACCCGCG" => -3.8,
            "TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT" => -6.85,
            "GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA" => -15.50,
            "TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA" => -18.10,
            "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA" => -3.65,
        )

        for (seq, ufold) in unafold_dgs
            d = dg(seq)

            # accepting a 60% difference
            delta = abs(0.6 * min(d, ufold))
            @test isapprox(d, ufold, atol=delta)
        end
    end

    @testset "fold_rna" begin
        # unafold's estimates for free energy estimates of RNA oligos
        # most tests available at https://github.com/jaswindersingh2/SPOT-RNA/blob/master/sample_inputs/batch_seq.fasta  
        unafold_dgs = Dict(
            "ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA" => -9.5,
            "AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC" => -10.1,
            "UUGGAGUACACAACCUGUACACUCUUUC" => -4.3,
            "AGGGAAAAUCCC" => -3.3,
            "GCUUACGAGCAAGUUAAGCAAC" => -4.6,
            "UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA" => -32.8,
            "GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG" => -20.7,
            "GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA" => -31.4,
            "CAGCGCGGCGGGCGGGAGUCCGGCGCGCCCUCCAUCCCCGGCGGCGUCGGCAAGGAGUAG" => -18.26,
        )

        for (seq, ufold) in unafold_dgs
            d = dg(seq)

            # accepting a 30% difference
            delta = abs(0.3 * min(d, ufold))
            @test isapprox(d, ufold, atol=delta)
        end
    end

    @testset "dot_bracket" begin
        seq = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
        structs = fold(seq)

        @test dot_bracket(seq, structs) == "((((((((.((((......))))..((((.......)))).))))))))"

        seq = "ACGCTCACCGTGCCCAGTGAGCGA"
        structs = fold(seq)
        @test length(seq) == length(dot_bracket(seq, structs))
    end

    @testset "_multibranch" begin
        seq = "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"  # three branch

        structs = fold(seq)
        @test any(s -> occursin("BIFURCATION", s.desc) && (8, 42) in s.ij, structs)

        seq = "CAGCGCGGCGGGCGGGAGUCCGGCGCGCCCUCCAUCCCCGGCGGCGUCGGCAAGGAGUAG"

        structs = fold(seq)
        @test any(s -> occursin("BIFURCATION", s.desc) && (3, 57) in s.ij, structs)
    end

    @testset "_pair" begin
        seq = "ATGGAATAGTG"
        @test SeqFold._pair(seq, 1, 2, 10, 11) == "AT/TG"
    end

    @testset "_stack" begin
        seq = "GCUCAGCUGGGAGAGC"
        temp = 37+273.15

        @test isapprox(SeqFold._stack(seq, 2, 3, 15, 14, temp, SeqFold.RNA_ENERGIES), -2.1, atol=0.1)
    end

    @testset "_bulge" begin
        # mock bulge of CAT on one side and AG on other
        # from pg 429 of SantaLucia, 2004
        seq = "ACCCCCATCCTTCCTTGAGTCAAGGGGCTCAA"

        pair_dg = SeqFold._bulge(seq, 6, 8, 19, 18, 37+273.15, SeqFold.DNA_ENERGIES)
        @test isapprox(pair_dg, 3.22, atol=0.4)
    end

    @testset "_hairpin" begin
        # hairpin = "CCTTGG"
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 12
        j = 17
        temp = 37+273.15
        hairpin_dg = SeqFold._hairpin(seq, i, j, temp, SeqFold.DNA_ENERGIES)
        # this differs from Unafold
        @test isapprox(hairpin_dg, 4.3, atol=1.0)

        # from page 428 of SantaLucia, 2004
        # hairpin = "CGCAAG"
        seq = "ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 4
        j = 9
        hairpin_dg = SeqFold._hairpin(seq, i, j, temp, SeqFold.DNA_ENERGIES)
        @test isapprox(hairpin_dg, 0.67, atol=0.1)

        seq = "CUUUGCACG"
        i = 1
        j = 9
        hairpin_dg = SeqFold._hairpin(seq, i, j, temp, SeqFold.RNA_ENERGIES)
        @test isapprox(hairpin_dg, 4.5, atol=0.2)
    end

    @testset "_internal_loop" begin
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 7
        j = 22
        temp = 37+273.15
        dg_val = SeqFold._internal_loop(seq, i, i + 4, j, j - 4, temp, SeqFold.DNA_ENERGIES)
        @test isapprox(dg_val, 3.5, atol=0.1)
    end

    @testset "_w!" begin
        test_cases = [
            ("GCUCAGCUGGGAGAGC", -3.8),
            ("CCUGCUUUGCACGCAGG", -6.4),
            ("GCGGUUCGAUCCCGC", -4.2)
        ]

        for (seq, expected_e) in test_cases
            i = 1
            j = length(seq)
            temp = 37 + 273.15
            n = length(seq)
            v_cache = [fill(SeqFold.STRUCT_DEFAULT, n) for _ in 1:n]
            w_cache = [fill(SeqFold.STRUCT_DEFAULT, n) for _ in 1:n]
            str = SeqFold._w!(seq, i, j, temp, v_cache, w_cache, SeqFold.RNA_ENERGIES)
            @test isapprox(str.e, expected_e, atol=0.2)
        end
    end
end

@testset "utils.jl" begin
    @testset "complement" begin
        @test isempty(SeqFold.complement(""))

        @test_throws ErrorException SeqFold.complement("Σωμϵ υνικωδε ħℯρ€")

        seq1 = "AaGgtTCCc"
        seq1c = "TTCCAAGGG"
        @test SeqFold.complement(seq1) == seq1c

        for _ in 1:20
            n = rand(1:55)
            seq = _randna(n)
            @test length(SeqFold.complement(seq)) == n
            @test seq == SeqFold.complement(SeqFold.complement(seq))
        end

        @test SeqFold.complement.(collect("AGCTNagctn")) == collect("TCGANTCGAN")

        bad_nucl_compl = SeqFold.complement("some bad nucleotides")
        @test issubset(bad_nucl_compl, "AGTCN")
    end

    @testset "revcomp" begin
        @test isempty(SeqFold.revcomp(""))
        @test_throws ErrorException SeqFold.revcomp("Σωμϵ υνικωδε ħℯρ€")

        seq1 = "AaGgtTCCc"
        seq1rc = "GGGAACCTT"
        @test SeqFold.revcomp(seq1) == seq1rc

        for _ in 1:20
            n = rand(1:55)
            seq = _randna(n)
            @test length(SeqFold.revcomp(seq)) == n
            @test seq == SeqFold.revcomp(SeqFold.revcomp(seq))
        end

        @test SeqFold.revcomp.(collect("AGCTNagctn")) == collect("TCGANTCGAN")

        bad_nucl_compl = SeqFold.revcomp("some bad nucleotides")
        @test issubset(bad_nucl_compl, "AGTCN")
    end

    @testset "gc_content" begin 
        @test isnan(SeqFold.gc_content(""))

        for _ in 1:20
            n = rand(1:55)
            seq = _randna(n)
            gc = SeqFold.gc_content(seq)
            @test isapprox(gc, _gc(seq), atol=1e-5)
        end
    end
end
