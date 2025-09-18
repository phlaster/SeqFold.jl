# Melting Temperature Calculation

[`SeqFold.jl`](https://github.com/phlaster/SeqFold.jl) provides functions for calculating DNA melting temperatures using nearest-neighbor thermodynamics.

## Showcase

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

## Public API

```@docs
tm
MeltingConditions
SeqFold.tm_cache
```