using Pkg
Pkg.activate(;temp=false)
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PKG_OFFLINE"] = true

Pkg.add(url="https://github.com/phlaster/SeqFold.jl")
Pkg.add(["PythonCall", "CSV", "DataFrames", "Plots", "StatsBase"])

using SeqFold
using PythonCall
using Random
using CSV
using DataFrames
using Statistics
using Printf
using Plots
using StatsBase

seqfold_tm = PythonCall.pyimport("seqfold"=>"tm")
seqfold_fold = PythonCall.pyimport("seqfold"=>"fold")
bio_tm = PythonCall.pyimport("Bio.SeqUtils.MeltingTemp" => "Tm_NN")

generate_random_dna(n) = String(rand("AGTC", n))


function benchmark_function(func, sequence, repetitions::Int)
    fastest_time = Inf
    for _ in 1:repetitions
        start_time = time_ns()
        try
            func(sequence)
            elapsed = (time_ns() - start_time) / 1e9
            fastest_time = min(fastest_time, elapsed)
        catch
            return Inf
        end
    end
    return fastest_time
end

function run_benchmark(batch_size=100)
    lengths = [5, 10, 20, 40, 80]
    n_sequences = [2^12, 2^11, 2^10, 2^9, 2^8]  
    results = DataFrame(
        length = Int[],
        function_name = String[],
        fastest_time_us = Float64[],
        avg_time_us = Float64[],
        median_time_us = Float64[],
        p90_time_us = Float64[]
    )
    Random.seed!(12345)
    
    julia_tm_func = tm
    python_seqfold_tm_func = seqfold_tm
    biopython_tm_func = bio_tm
    
    julia_fold_func = fold
    python_seqfold_fold_func = seqfold_fold

    for (i, len) in enumerate(lengths)
        total_sequences = n_sequences[i]
        println("\n- Benchmarking length $len with $total_sequences sequences...")
        
        all_julia_tm_times = Float64[]
        all_python_seqfold_tm_times = Float64[]
        all_biopython_tm_times = Float64[]
        
        all_julia_fold_times = Float64[]
        all_python_seqfold_fold_times = Float64[]
        
        num_batches = ceil(Int, total_sequences / batch_size)
        
        for batch_idx in 1:num_batches
            batch_start = (batch_idx - 1) * batch_size + 1
            batch_end = min(batch_start + batch_size - 1, total_sequences)
            batch_size_actual = batch_end - batch_start + 1
            println("  → Processing batch $batch_idx/$num_batches ($batch_start:$batch_end)...")
            
            sequences = [generate_random_dna(len) for _ in 1:batch_size_actual]
            
            println("    - Benchmarking Tm functions...")
            println("      • SeqFold.jl tm...")
            julia_tm_times = Vector{Float64}(undef, batch_size_actual)
            for j in 1:batch_size_actual
                julia_tm_times[j] = benchmark_function(julia_tm_func, sequences[j], 3)
            end
            append!(all_julia_tm_times, julia_tm_times)
            
            println("      • seqfold (Python) tm...")
            python_seqfold_tm_times = Vector{Float64}(undef, batch_size_actual)
            for j in 1:batch_size_actual
                python_seqfold_tm_times[j] = benchmark_function(python_seqfold_tm_func, sequences[j], 3)
            end
            append!(all_python_seqfold_tm_times, python_seqfold_tm_times)
            println("      • Biopython tm...")
            biopython_tm_times = Vector{Float64}(undef, batch_size_actual)
            for j in 1:batch_size_actual
                biopython_tm_times[j] = benchmark_function(biopython_tm_func, sequences[j], 3)
            end
            append!(all_biopython_tm_times, biopython_tm_times)
            
            println("    - Benchmarking fold functions...")
            
            println("      • SeqFold.jl fold...")
            julia_fold_times = Vector{Float64}(undef, batch_size_actual)
            for j in 1:batch_size_actual
                julia_fold_times[j] = benchmark_function(julia_fold_func, sequences[j], 3)
            end
            append!(all_julia_fold_times, julia_fold_times)
            
            println("      • seqfold (Python) fold...")
            python_seqfold_fold_times = Vector{Float64}(undef, batch_size_actual)
            for j in 1:batch_size_actual
                python_seqfold_fold_times[j] = benchmark_function(python_seqfold_fold_func, sequences[j], 3)
            end
            append!(all_python_seqfold_fold_times, python_seqfold_fold_times)
            
            println("\n      Batch statistics for length $len:")
            println("      TM Functions:")
            println(@sprintf("        SeqFold.jl: min=%.2f μs, avg=%.2f μs", 
                minimum(julia_tm_times) * 1e6, 
                mean(julia_tm_times) * 1e6))
            println(@sprintf("        seqfold (Python): min=%.2f μs, avg=%.2f μs", 
                minimum(python_seqfold_tm_times) * 1e6, 
                mean(python_seqfold_tm_times) * 1e6))
            println(@sprintf("        Biopython: min=%.2f μs, avg=%.2f μs", 
                minimum(biopython_tm_times) * 1e6, 
                mean(biopython_tm_times) * 1e6))
            
            println("      Fold Functions:")
            println(@sprintf("        SeqFold.jl: min=%.2f μs, avg=%.2f μs", 
                minimum(julia_fold_times) * 1e6, 
                mean(julia_fold_times) * 1e6))
            println(@sprintf("        seqfold (Python): min=%.2f μs, avg=%.2f μs", 
                minimum(python_seqfold_fold_times) * 1e6, 
                mean(python_seqfold_fold_times) * 1e6))
            println("      " * "-"^50)
        end
        
        all_julia_tm_times_us = all_julia_tm_times * 1e6
        all_python_seqfold_tm_times_us = all_python_seqfold_tm_times * 1e6
        all_biopython_tm_times_us = all_biopython_tm_times * 1e6
        
        all_julia_fold_times_us = all_julia_fold_times * 1e6
        all_python_seqfold_fold_times_us = all_python_seqfold_fold_times * 1e6
        
        push!(results, (
            len,
            "SeqFold.jl (tm)",
            minimum(all_julia_tm_times_us),
            mean(all_julia_tm_times_us),
            median(all_julia_tm_times_us),
            quantile(all_julia_tm_times_us, 0.9)
        ))
        
        push!(results, (
            len,
            "seqfold (Python) (tm)",
            minimum(all_python_seqfold_tm_times_us),
            mean(all_python_seqfold_tm_times_us),
            median(all_python_seqfold_tm_times_us),
            quantile(all_python_seqfold_tm_times_us, 0.9)
        ))
        
        push!(results, (
            len,
            "Biopython (tm)",
            minimum(all_biopython_tm_times_us),
            mean(all_biopython_tm_times_us),
            median(all_biopython_tm_times_us),
            quantile(all_biopython_tm_times_us, 0.9)
        ))
        
        push!(results, (
            len,
            "SeqFold.jl (fold)",
            minimum(all_julia_fold_times_us),
            mean(all_julia_fold_times_us),
            median(all_julia_fold_times_us),
            quantile(all_julia_fold_times_us, 0.9)
        ))
        
        push!(results, (
            len,
            "seqfold (Python) (fold)",
            minimum(all_python_seqfold_fold_times_us),
            mean(all_python_seqfold_fold_times_us),
            median(all_python_seqfold_fold_times_us),
            quantile(all_python_seqfold_fold_times_us, 0.9)
        ))
        
        CSV.write("assets/benchmark.csv", results)
        
        println("\n- Summary for length $len:")
        println("  TM Functions:")
        println(@sprintf("    SeqFold.jl: min=%.2f μs, avg=%.2f μs", 
            results[end-4, :fastest_time_us], 
            results[end-4, :avg_time_us]))
        println(@sprintf("    seqfold (Python): min=%.2f μs, avg=%.2f μs", 
            results[end-3, :fastest_time_us], 
            results[end-3, :avg_time_us]))
        println(@sprintf("    Biopython: min=%.2f μs, avg=%.2f μs", 
            results[end-2, :fastest_time_us], 
            results[end-2, :avg_time_us]))
        
        println("  Fold Functions:")
        println(@sprintf("    SeqFold.jl: min=%.2f μs, avg=%.2f μs", 
            results[end-1, :fastest_time_us], 
            results[end-1, :avg_time_us]))
        println(@sprintf("    seqfold (Python): min=%.2f μs, avg=%.2f μs", 
            results[end, :fastest_time_us], 
            results[end, :avg_time_us]))
        println("-" ^ 70)
    end
    
    println("\n Benchmarking complete. Results saved to assets/benchmark.csv")
    return results
end

function plot_benchmark()
    results = CSV.read("assets/benchmark.csv", DataFrame)

    results.fastest_time_us ./= 1e6
    results.avg_time_us ./= 1e6
    results.median_time_us ./= 1e6
    results.p90_time_us ./= 1e6

    tm_speedup_20 = results[results.function_name .== "seqfold (Python) (tm)" .&& results.length .== 20, :avg_time_us][1] / 
                    results[results.function_name .== "SeqFold.jl (tm)" .&& results.length .== 20, :avg_time_us][1]

    fold_speedup_20 = results[results.function_name .== "seqfold (Python) (fold)" .&& results.length .== 20, :avg_time_us][1] / 
                    results[results.function_name .== "SeqFold.jl (fold)" .&& results.length .== 20, :avg_time_us][1]

    p = plot(
        dpi=300,
        size=(800, 600),
        framestyle=:box,
        grid=true,
        minorgrid=true,
        gridalpha=0.3,
        gridlinewidth=1.5,
        gridcolor=:gray80,
        foreground_color_axis=:black,
        background_color_inside=:white,
        legend=:topleft,
        legendfontsize=10,
        title="SeqFold.jl vs seqfold Performance Comparison",
        titlefontsize=14,
        xlabel="Sequence Length (nt)",
        ylabel="Time (s)",
        xscale=:log10,
        yscale=:log10,
        xticks=([5, 10, 20, 40, 80], ["5", "10", "20", "40", "80"]),
        yticks=10.0 .^ collect(-6:0),
        ylims=(1e-6, 1.5),
    )

    biopython_tm = results[results.function_name .== "Biopython (tm)", :]
    plot!(
        biopython_tm.length,
        biopython_tm.avg_time_us,
        fillalpha=0.2,
        label="Biopython: Tm_NN",
        color="#359024FF",
        markershape=:utriangle,
        markersize=6,
        linewidth=2,
        markerstrokewidth=0.5
    )

    python_tm = results[results.function_name .== "seqfold (Python) (tm)", :]
    plot!(
        python_tm.length,
        python_tm.avg_time_us,
        fillalpha=0.2,
        label="Py: tm",
        color="#4063D8FF",
        markershape=:diamond,
        markersize=6,
        linewidth=2,
        markerstrokewidth=0.5
    )

    julia_tm = results[results.function_name .== "SeqFold.jl (tm)", :]
    plot!(
        julia_tm.length,
        julia_tm.avg_time_us,
        fillalpha=0.2,
        label="Jl: tm, $(round(tm_speedup_20, digits=1))× than seqfold",
        color="#9558B2FF",
        markershape=:circle,
        markersize=6,
        linewidth=2,
        markerstrokewidth=0.5
    )

    python_fold = results[results.function_name .== "seqfold (Python) (fold)", :]
    plot!(
        python_fold.length,
        python_fold.avg_time_us,
        fillalpha=0.2,
        label="Py: fold",
        color="#4063D8FF",
        markershape=:diamond,
        markersize=6,
        linewidth=2,
        markerstrokewidth=0.5,
        linestyle=:dash
    )

    julia_fold = results[results.function_name .== "SeqFold.jl (fold)", :]
    plot!(
        julia_fold.length,
        julia_fold.avg_time_us,
        fillalpha=0.2,
        label="Jl: fold, $(round(fold_speedup_20, digits=1))× than seqfold",
        color="#9558B2FF",
        markershape=:circle,
        markersize=6,
        linewidth=2,
        markerstrokewidth=0.5,
        linestyle=:dash
    )

    savefig(p, "assets/benchmark.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark()
    plot_benchmark()
end