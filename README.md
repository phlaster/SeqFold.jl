# <div align="center"> <img src="docs/src/assets/logo.png" alt="SeqFold.jl" width="500"></img></div>
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://phlaster.github.io/SeqFold.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://phlaster.github.io/SeqFold.jl/dev/)
[![Build Status](https://github.com/phlaster/SeqFold.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/phlaster/SeqFold.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)


## Introduction

`SeqFold.jl` is a high-performance pure Julia reimplementation of [`seqfold`](https://github.com/Lattice-Automation/seqfold) Python library for predicting nucleic acid secondary structures and calculating melting temperatures, which is, in turn, an implementation of the `Zuker, 1981` dynamic programming algorithm, the basis for [UNAFold/mfold](https://www.unafold.org/), with energy functions from `SantaLucia, 2004` (DNA) and `Turner, 2009` (RNA).

## Motivation

Secondary structure prediction is essential for:
- Designing PCR primers with minimal secondary structure;
- Creating oligos for genome editing techniques like MAGE;
- Tuning ribosome binding site (RBS) expression rates;
- Analyzing potential off-target binding in CRISPR applications.

`SeqFold.jl` provides a minimalist open-source alternative to proprietary solutions [with on-demand access](https://vfold.missouri.edu/software.html). The package has [well-documented API](https://phlaster.github.io/SeqFold.jl/stable/) and lightweight architecture with no 3rd party dependencies.

Apart from sequence folding utilities several functions for $T_m$ calculation are provided with fine-grained ionic conditions control.

[Julia programming language](https://julialang.org/) is well-known for its speed and this package outperforms both the [original Python package](https://github.com/Lattice-Automation/seqfold) and Biopython `Tm_NN` reference, with little to no difference in results to the latter:
![seqfold vs SeqFold.jl](docs/src/assets/benchmark.png)



## Basic Usage

### Installation

First option (inside Julia REPL):
```julia
julia> ]
pkg> add SeqFold
```
Second option (from Julia scripts):
```julia
using Pkg
Pkg.add("SeqFold")
```
either way you should now be able to call the package:
```julia
using SeqFold
```

### Melting Temperature Calculation

```julia
julia> seq = "GGGAGGTCAGCAAACCTGAACCTGTTGAGATGTTGACGTCAGGAAACCCT";

julia> tm(seq) # melting at PCR conditions, same as tm(seq, conditions=:pcr)
80.4

julia> MeltingConditions(:pcr) # here they are as a preset inherited from seqfold library
MeltingConditions (PCR preset)
  • NEB PCR buffer conditions for Taq DNA Polymerase
  • seq1 concentration: 250.0 nM (typical primer concentration)
  • seq2 concentration: 0.0 nM (asymmetric PCR)
  • Mg²⁺: 1.5 mM (optimal for Taq DNA Polymerase per NEB guidelines)
  • K⁺: 50.0 mM
  • Tris: 2.0 mM
  • dNTPs: 0.2 mM

julia> tm(seq, Mg=5) # altering default conditions
81.4

julia> MeltingConditions(:std) # second preset, inherited from seqfold library
MeltingConditions (standard preset)
  • Standard hybridization buffer
  • seq1 concentration: 25.0 nM
  • seq2 concentration: 25.0 nM
  • Na⁺: 50.0 mM

julia> tm(seq, conditions=:std, Na=150) # altering chosen preset
79.8

julia> custom_conds = MeltingConditions( # Specifying custom conditions
       seq1_conc=30,seq2_conc=20,Na=5,K=5,Tris=15,Mg=0,dNTPs=0)
MeltingConditions (custom)
  • seq1 concentration: 30.0 nM
  • seq2 concentration: 20.0 nM
  • Na⁺: 5.0 mM
  • K⁺: 5.0 mM
  • Tris: 15.0 mM

julia> tm(seq, conditions=custom_conds, Tris=50) # flexible adjustments
68.4
```

### Secondary Structure Prediction

```julia
# Get the list of possible secondary structures for a sequence
julia> folded = fold(seq) # calculated at 37°C by default, Julia uses 1-based indexing
12-element Vector{SeqFold.Structure}:
    3   40   -1.3  STACK:GA/CT    
    4   39   -1.4  STACK:AGG/TGC  
    6   37   -1.5  STACK:GT/CA    
    7   36   -1.3  STACK:TC/AG    
    8   35   -1.5  STACK:CA/GT    
    9   34   -2.0  STACK:AGC/TTG  
   11   32   -1.5  STACK:CA/GT    
   12   31    2.8  INTERIOR_LOOP:4/2
   16   29   -1.3  STACK:CT/GA    
   17   28   -0.1  STACK:TGA/AGT  
   19   26   -1.0  STACK:AA/TT    
   20   25    3.5  HAIRPIN:AC/TG 

julia> dot_bracket(seq, folded) # get the dot-bracket notation
"..((.((((.((...((.((....)).)).)).)))).)).........."

julia> folded_hot = fold(seq, temp=80)
4-element Vector{SeqFold.Structure}:
    7   19   -0.4  STACK:TC/AG    
    8   18   -0.5  STACK:CA/GT    
    9   17   -0.4  STACK:AG/TC    
   10   16    3.6  HAIRPIN:GC/CC 

julia> dot_bracket(seq, folded_hot) # at a higher temperature number of secondary structures is decreased
"......((((.....))))..............................."

julia> dg(folded) # get the total free energy of the predicted secondary structures
-6.6

julia> dg(folded_hot)
2.3

julia> dg(seq, temp=4) # this method calls the computationally intensive `fold` internally
-16.2
```
## Calling from Python scripts
Install [`juliacall`](pypi.org/project/juliacall) package through you favourite package manager (pip, Conda, uv):
```bash
$ pip install juliacall
```
In Python REPL:
```py
Python 3.13.7 | packaged by conda-forge | (main, Sep  3 2025, 14:30:35) [GCC 14.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from juliacall import Main as jl # this will download Julia binaries
>>> jl.seval("""using Pkg; Pkg.add(url="https://github.com/phlaster/SeqFold.jl")""")
>>> jl.seval("using SeqFold") # this exports tm, fold, dg, ...
>>> jl.fold("ACACGCAGCGGGTGTGGCTGTGGTGC")
Julia:
4-element Vector{SeqFold.Structure}:
    1   13   -1.9  STACK_DE:AC/TG 
    2   12   -2.3  STACK:CAC/GGG  
    4   10   -2.2  STACK:CG/GC    
    5    9    3.5  HAIRPIN:GC/CG  
>>> jl.dg("ACACGCAGCGGGTGTGGCTGTGGTGC")
-2.9
```
Remember, that calling a function in Julia for the first time takes more time due to JIT compiliation.

## Citations

Papers, which helped in developing the original `seqfold` package:
* **Nussinov, 1980:**
  >Nussinov, Ruth, and Ann B. Jacobson. "Fast algorithm for predicting the secondary structure of single-stranded RNA." Proceedings of the National Academy of Sciences 77.11 (1980): 6309-6313.
* **Zuker, 1981:**
  >Zuker, Michael, and Patrick Stiegler. "Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information." Nucleic acids research 9.1 (1981): 133-148.
* **Jaeger, 1989:**
  >Jaeger, John A., Douglas H. Turner, and Michael Zuker. "Improved predictions of secondary structures for RNA." Proceedings of the National Academy of Sciences 86.20 (1989): 7706-7710.
* **SantaLucia, 2004:**
  >SantaLucia Jr, John, and Donald Hicks. "The thermodynamics of DNA structural motifs." Annu. Rev. Biophys. Biomol. Struct. 33 (2004): 415-440.
* **Turner, 2009:**
  >Turner, Douglas H., and David H. Mathews. "NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure." Nucleic acids research 38.suppl_1 (2009): D280-D282.
* **Ward, 2017:**
  >Ward, M., Datta, A., Wise, M., & Mathews, D. H. (2017). Advanced multi-loop algorithms for RNA secondary structure prediction reveal that the simplest model is best. Nucleic acids research, 45(14), 8541-8550.
